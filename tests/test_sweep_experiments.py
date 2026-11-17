"""Numerical, format, real bcftools/SF2, and optional simulation integration tests."""

import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from pysrc.sweep_experiments import benchmark, evo, inputs, simulation
from pysrc.sweep_experiments.__main__ import main


class SequenceReference(inputs.Reference):
    chrom = "7"
    accession = "NC_000007.14"

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

    def fetch(self, start, end):
        return self.sequence[start:end]


class ExperimentTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def vcf(self):
        return self.write("sites.vcf", """##fileformat=VCFv4.2
##contig=<ID=7,length=10000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestor">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb
7\t100\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|0\t0|1
7\t200\t.\tA\tG\t.\tPASS\tAA=G\tGT\t0|1\t1|1
7\t300\t.\tA\tG\t.\tPASS\tAA=A\tGT\t./.\t0|1
7\t400\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|1\t0|0
7\t400\t.\tA\tG,T\t.\tPASS\tAA=A\tGT\t0|1\t0|0
7\t500\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0
7\t600\t.\tA\tG\t.\tPASS\tAA=A\tGT\t1|1\t1|1
""")

    @unittest.skipUnless(shutil.which("bcftools"), "bcftools not installed")
    def test_real_bcftools_mask_counts_duplicate_positions_and_ancestry(self):
        mask = self.write("mask.bed", "7\t99\t450\n")
        sites = list(inputs.variants(self.vcf(), True, mask))
        self.assertEqual([(s.position, s.x, s.n) for s in sites],
                         [(100, 1, 4), (200, 1, 4), (300, 1, 2)])
        unfolded = list(inputs.variants(self.root / "sites.vcf", False, mask))
        self.assertEqual([s.x for s in unfolded], [1, 1, 1])
        with self.assertRaisesRegex(ValueError, "INFO/AA"):
            list(inputs.variants(self.root / "sites.vcf", False))

    def test_interval_genetic_distances_and_bed_boundaries(self):
        sites = [inputs.Variant("7", p, "A", "G", 1, 4, True) for p in (100, 200, 500)]
        map_file = self.write("map.tsv", "position\tcM\n1\t0\n1001\t0.001\n")
        rec = self.root / "rec.tsv"
        inputs.write_recombination(rec, sites, map_file, 1e-8)
        with rec.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for observed, expected in zip(rows, (0, .0001, .0003)):
            self.assertAlmostEqual(float(observed["rate"]), expected)
        bmap = self.write("b.bedgraph", "7\t99\t200\t0.4\n7\t200\t500\t0.8\n")
        output = self.root / "b.tsv"
        inputs.write_bvalues(output, sites, bmap)
        self.assertEqual(output.read_text().splitlines()[1:], ["100\t0.4", "200\t0.4", "500\t0.8"])
        with self.assertRaisesRegex(ValueError, "does not span"):
            inputs.interpolate(1100, [1, 1000], [0, 1])
        with self.assertRaisesRegex(ValueError, "does not cover"):
            inputs.write_bvalues(output, [inputs.Variant("7", 501, "A", "G", 1, 4, True)], bmap)

    def test_window_offsets_and_reverse_complement_probabilities(self):
        reference = SequenceReference("ACGT" * 10)
        calls = []

        class Predictor:
            def log_probs(self, sequence, target):
                calls.append((sequence, target))
                return {b: math.log(p) for b, p in zip("ACGT", (.1, .2, .3, .4))}

        output = self.root / "evo.tsv"
        count = evo.score_sites([inputs.Variant("7", 17, "A", "G", 1, 4, True)],
                                reference, Predictor(), output, 8, [2, 5])
        self.assertEqual(count, 4)
        self.assertEqual([target for _, target in calls], [2, 2, 5, 5])
        with output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertAlmostEqual(float(rows[0]["delta_logp"]), math.log(3))
        self.assertAlmostEqual(float(rows[1]["delta_logp"]), math.log(.2 / .4))
        self.assertAlmostEqual(float(rows[1]["p_A"]), .4)
        self.assertEqual(rows[1]["strand"], "-")
        for window in evo.windows(17, 40, 8, [2, 5]):
            genomic_target = (window.start + window.target if window.strand == "+"
                              else window.end - 1 - window.target)
            self.assertEqual(genomic_target, 16)
        self.assertEqual(evo.windows(1, 40, 8, [2]), [])
        with self.assertRaisesRegex(ValueError, "offset"):
            evo.windows(17, 40, 8, [0])

    def test_model_predicts_at_previous_token_not_target_token(self):
        class Scalar:
            def item(self):
                return -.5

        class Distribution:
            def float(self):
                return self

            def log_softmax(self, dim):
                return self

            def __getitem__(self, key):
                return Scalar()

        seen = []

        class Logits:
            def __getitem__(self, key):
                seen.append(key)
                return Distribution()

        class Tokenizer:
            def tokenize(self, sequence):
                return list(range(len(sequence)))

        class Model:
            tokenizer = Tokenizer()

            def __call__(self, inputs):
                return (Logits(), None), None

        class Torch:
            long = "long"

            def tensor(self, values, **kwargs):
                return values

            def inference_mode(self):
                from contextlib import nullcontext
                return nullcontext()

        predictor = evo.EvoPredictor.__new__(evo.EvoPredictor)
        predictor.model, predictor.torch = Model(), Torch()
        self.assertEqual(predictor.log_probs("ACGT", 2)["A"], -.5)
        self.assertEqual(seen, [(0, 1)])

    def test_scoring_accepts_mismatching_ref_and_preserves_simulated_alleles(self):
        reference = SequenceReference("ACGT" * 10)
        calls = []

        def predict(sequence, target):
            calls.append((sequence, target))
            return {b: math.log(p) for b, p in zip("ACGT", (.1, .2, .3, .4))}

        output = self.root / "mismatching-ref.tsv"
        count = evo.score_sites(
            [inputs.Variant("chr7", 17, "C", "G", 1, 4, True)], reference,
            SimpleNamespace(log_probs=predict), output, 8, [2])
        self.assertEqual(count, 2)
        self.assertEqual([sequence[target] for sequence, target in calls], ["A", "T"])
        with output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([(row["ref"], row["alt"], row["status"]) for row in rows],
                         [("C", "G", "ok"), ("C", "G", "ok")])
        self.assertAlmostEqual(float(rows[0]["delta_logp"]), math.log(.3 / .2))
        self.assertAlmostEqual(float(rows[1]["delta_logp"]), math.log(.2 / .3))

    def test_reference_check_keeps_chromosome_and_position_validation(self):
        reference = SequenceReference("ACGT")
        for chrom, position in (("7", 1), ("chr7", 4)):
            reference.check(inputs.Variant(chrom, position, "C", "G", 1, 4, True))
        with self.assertRaisesRegex(ValueError, "contig"):
            reference.check(inputs.Variant("9", 1, "A", "G", 1, 4, True))
        for position in (0, 5):
            with self.subTest(position=position):
                with self.assertRaisesRegex(ValueError, "position outside"):
                    reference.check(inputs.Variant("7", position, "A", "G", 1, 4, True))

    def test_bad_map_values_and_external_frequency_mask_fail_explicitly(self):
        bad = self.write("bad.tsv", "position\tcM\n1\t0\n2\tnan\n")
        with self.assertRaisesRegex(ValueError, "finite"):
            inputs.cumulative_map(bad)
        frequency = self.write("external.freq", "position\tx\tn\tfolded\n100\t4\t4\t0\n")
        sites = inputs.read_frequency(frequency, "7", False)
        mask = self.write("mask.bed", "7\t0\t99\n")
        with self.assertRaisesRegex(ValueError, "noncallable"):
            inputs.check_mask(sites, mask)
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            inputs.read_frequency(frequency, "7", True)

    def test_rocm_model_selection_and_explicit_opt_in(self):
        class FakeOOM(RuntimeError):
            pass

        calls = []

        class FakeEvo:
            def __init__(self, name):
                calls.append(name)
                self.model = SimpleNamespace(eval=lambda: None)

        torch = SimpleNamespace(
            version=SimpleNamespace(hip="test"), OutOfMemoryError=FakeOOM,
            cuda=SimpleNamespace(is_available=lambda: True, device_count=lambda: 1,
                                 mem_get_info=lambda i: (24 * 1024**3, 24 * 1024**3)))
        with patch.dict(sys.modules, {"torch": torch, "evo2": SimpleNamespace(Evo2=FakeEvo)}):
            with self.assertRaisesRegex(RuntimeError, "experimental"):
                evo.EvoPredictor("auto", False, 8)
            with patch.object(evo.EvoPredictor, "log_probs", return_value={"A": -.5}):
                predictor = evo.EvoPredictor("auto", True, 8)
                self.assertEqual(predictor.selected, "evo2_7b")
                self.assertEqual(calls, ["evo2_7b"])
                self.assertEqual(predictor.attempts[0]["model"], "evo2_40b")

    def test_only_oom_causes_model_downgrade(self):
        class FakeOOM(RuntimeError):
            pass

        calls = []

        class FakeEvo:
            def __init__(self, name):
                calls.append(name)
                if name == "evo2_40b":
                    raise FakeOOM("test memory pressure")
                self.model = SimpleNamespace(eval=lambda: None)

        torch = SimpleNamespace(
            version=SimpleNamespace(hip=None), OutOfMemoryError=FakeOOM,
            cuda=SimpleNamespace(is_available=lambda: True, device_count=lambda: 1,
                                 get_device_capability=lambda i: (9, 0),
                                 mem_get_info=lambda i: (200 * 1024**3, 200 * 1024**3),
                                 empty_cache=lambda: None))
        with patch.dict(sys.modules, {"torch": torch, "evo2": SimpleNamespace(Evo2=FakeEvo)}), \
                patch.object(evo.EvoPredictor, "log_probs", return_value={"A": -.5}):
            predictor = evo.EvoPredictor("auto", False, 8)
            self.assertEqual(calls, ["evo2_40b", "evo2_7b"])
            self.assertEqual(predictor.selected, "evo2_7b")
        def broken_backend(name):
            raise ImportError("missing kernel")

        with patch.dict(sys.modules, {"torch": torch, "evo2": SimpleNamespace(Evo2=broken_backend)}):
            with self.assertRaisesRegex(ImportError, "missing kernel"):
                evo.EvoPredictor("auto", False, 8)

    @unittest.skipUnless(importlib.util.find_spec("biocommons"), "optional SeqRepo not installed")
    def test_real_seqrepo_lookup_uses_assembly_and_exact_chromosome_intervals(self):
        from biocommons.seqrepo import SeqRepo
        import pysam

        directory = self.root / "seqrepo"
        # Use pysam's real BGZF compressor to build the fixture without a bgzip CLI.
        def compress(command):
            pysam.tabix_compress(command[-1], command[-1] + ".gz", force=True)

        with patch("biocommons.seqrepo.fastadir.fabgz._find_bgzip", return_value="bgzip"), \
                patch("biocommons.seqrepo.fastadir.fabgz.subprocess.check_call", side_effect=compress):
            repo = SeqRepo(str(directory), writeable=True)
            repo.store("ACGT" * 100, [
                {"namespace": "NCBI", "alias": "NC_000007.14"},
            ])
            repo.store("TGCA" * 100, [
                {"namespace": "GRCh37", "alias": "7"},
                {"namespace": "NCBI", "alias": "NC_000007.13"},
                {"namespace": "GRCh38.p14", "alias": "9"},
            ])
            repo.commit()
        reference = inputs.Reference(directory, "testing")
        self.assertEqual((reference.chrom, reference.accession, reference.length),
                         ("7", "NC_000007.14", 159345973))
        self.assertEqual(reference.fetch(1, 5), "CGTA")
        self.assertEqual(reference.fetch(0, 1), "A")
        self.assertEqual(reference.fetch(399, 400), "T")
        with patch.object(reference.repo, "fetch", wraps=reference.repo.fetch) as fetch:
            self.assertEqual(reference.fetch(3, 4), "T")
            fetch.assert_called_once_with("NC_000007.14", namespace="NCBI", start=3, end=4)
            fetch.reset_mock()
            reference.check(inputs.Variant("chr7", 4, "C", "G", 1, 4, True))
            for start, end in ((-1, 1), (0, 0), (2, 1), (0, reference.length + 1)):
                with self.subTest(start=start, end=end):
                    with self.assertRaisesRegex(ValueError, "outside"):
                        reference.fetch(start, end)
            fetch.assert_not_called()
        with self.assertRaisesRegex(ValueError, "Truncated"):
            reference.fetch(399, 401)
        with self.assertRaisesRegex(KeyError, "NC_000009.12"):
            inputs.Reference(directory, "final")
        reference.check(inputs.Variant("chr7", 1, "A", "G", 1, 4, True))
        reference.check(inputs.Variant("7", 1, "C", "G", 1, 4, True))
        with self.assertRaisesRegex(ValueError, "contig"):
            reference.check(inputs.Variant("9", 1, "A", "G", 1, 4, True))

    def test_estimate_b_uses_ratio_of_means_and_matched_replicates(self):
        path = self.write("pi.tsv", "scenario\treplicate\tchrom\tstart\tend\tpi\n"
                          "neutral\t1\t7\t0\t100\t0.01\nneutral\t2\t7\t0\t100\t0.03\n"
                          "bgs\t1\t7\t0\t100\t0.005\nbgs\t2\t7\t0\t100\t0.006\n")
        output = self.root / "b.tsv"
        simulation.estimate_b(path, output, bootstrap=100)
        with output.open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertAlmostEqual(float(row["B"]), .275)

    def test_held_out_calibration_and_conservative_ties(self):
        lines = ["case\tscenario\tsplit\tconfiguration\tscore\tpeak_position\t"
                 "sweep_position\tlocalization_error_bp"]
        lines.extend(f"c{i}\tneutral\tcalibration\tphysical\t{i}\t1\t\t" for i in range(19))
        lines.extend(["e1\tbgs\tevaluation\tphysical\t20\t1\t\t",
                      "e2\tsweep\tevaluation\tphysical\t18\t1\t2\t1"])
        summary = self.write("summary.tsv", "\n".join(lines) + "\n")
        output = self.root / "evaluation.json"
        benchmark.evaluate(summary, output, .05)
        report = json.loads(output.read_text())
        self.assertEqual([r["detected"] for r in report["results"]], [True, False])
        self.assertEqual(report["results"][0]["empirical_p"], .05)

    @unittest.skipUnless(importlib.util.find_spec("msprime"), "optional msprime not installed")
    def test_neutral_simulation_reproducibility_and_reference_alignment(self):
        reference = SequenceReference("ACGT" * 5000)
        manifests = [simulation.simulate(self.root / name, reference, 100, 5100,
                                        2, 4, 1000, 1e-6, 1e-8, 41)
                     for name in ("first", "second")]
        first = json.loads(manifests[0].read_text())
        self.assertEqual([case["split"] for case in first["cases"]], ["calibration", "evaluation"])
        for case in first["cases"]:
            a = manifests[0].parent / case["vcf"]
            b = manifests[1].parent / case["vcf"]
            self.assertEqual(a.read_text(), b.read_text())
            sites = list(inputs.variants(a, False))
            self.assertTrue(sites)
            for site in sites:
                reference.check(site)
                self.assertEqual(reference.fetch(site.position - 1, site.position), site.ref)

    @unittest.skipUnless(shutil.which("bcftools"), "bcftools not installed")
    def test_real_sweepfinder_manifest_runner(self):
        executable = os.environ.get("SWEEPFINDER2_TEST_BINARY", shutil.which("SweepFinder2"))
        if not executable:
            self.skipTest("set SWEEPFINDER2_TEST_BINARY to enable real SF2 integration")
        self.vcf()
        self.write("mask.bed", "7\t0\t10000\n")
        self.write("map.tsv", "position\tcM\n1\t0\n10000\t0.01\n")
        manifest = self.write("manifest.json", json.dumps({"cases": [{
            "id": "test", "vcf": "sites.vcf", "mask": "mask.bed",
            "recombination_map": "map.tsv", "dataset": "testing",
            "scenario": "unknown", "split": "exploratory"}]}))
        output = self.root / "run"
        main(["run", "--manifest", str(manifest), "--output", str(output),
              "--configurations", "physical,constant,map", "--grid", "10",
              "--sweepfinder", executable])
        self.assertTrue((output / "test" / "background.spect").is_file())
        self.assertEqual(json.loads((output / "run.json").read_text())["status"], "complete")
        with (output / "summary.tsv").open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 3)
        sentinel = output / "test" / "keep.txt"
        sentinel.write_text("unrelated file")
        (output / "summary.tsv").write_text("stale summary")
        with self.assertLogs(benchmark.LOG, level="WARNING") as messages:
            main(["run", "--manifest", str(manifest), "--output", str(output),
                  "--configurations", "physical,constant,map", "--grid", "10",
                  "--sweepfinder", executable])
        self.assertTrue(any("Reusing output directory" in message for message in messages.output))
        self.assertEqual(sentinel.read_text(), "unrelated file")
        self.assertEqual(json.loads((output / "run.json").read_text())["status"], "complete")
        with (output / "summary.tsv").open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 3)

    @unittest.skipUnless(shutil.which("bcftools"), "bcftools not installed")
    def test_real_sweepfinder_bmap_mode(self):
        executable = os.environ.get("SWEEPFINDER2_TEST_BINARY", shutil.which("SweepFinder2"))
        if not executable:
            self.skipTest("set SWEEPFINDER2_TEST_BINARY to enable real SF2 integration")
        self.vcf()
        self.write("mask.bed", "7\t0\t450\n")
        self.write("map.tsv", "position\tcM\n1\t0\n10000\t0.01\n")
        self.write("bmap.bedgraph", "7\t0\t10000\t0.6\n")
        frequency = self.write("external.freq", "position\tx\tn\tfolded\n"
                               "100\t1\t4\t0\n150\t4\t4\t0\n200\t1\t4\t0\n"
                               "250\t4\t4\t0\n300\t1\t2\t0\n350\t4\t4\t0\n")
        spectrum = self.root / "external.spect"
        benchmark.run_command([executable, "-f", str(frequency), str(spectrum)],
                              self.root / "spectrum.log")
        manifest = self.write("manifest.json", json.dumps({"cases": [{
            "id": "test", "vcf": "sites.vcf", "mask": "mask.bed",
            "recombination_map": "map.tsv", "bmap": "bmap.bedgraph",
            "sf2_frequency": frequency.name, "background_spectrum": spectrum.name,
            "n1": 10000, "n2": 10000, "divergence_time": 100000,
            "dataset": "testing", "scenario": "unknown", "split": "exploratory"}]}))
        output = self.root / "run-b"
        main(["run", "--manifest", str(manifest), "--output", str(output),
              "--configurations", "bmap", "--grid", "10", "--unfolded",
              "--sweepfinder", executable])
        self.assertTrue((output / "test" / "bmap.sf2.out").is_file())

    @unittest.skipUnless(importlib.util.find_spec("msprime") and shutil.which("bcftools"),
                         "requires msprime and bcftools")
    def test_simulation_to_scanning_to_held_out_evaluation(self):
        executable = os.environ.get("SWEEPFINDER2_TEST_BINARY", shutil.which("SweepFinder2"))
        if not executable:
            self.skipTest("set SWEEPFINDER2_TEST_BINARY to enable real SF2 integration")
        with self.assertLogs("msprime", level="INFO"):
            manifest = simulation.simulate(self.root / "sim", SequenceReference("ACGT" * 2500),
                                           100, 5100, 40, 4, 1000, 1e-6, 1e-8, 41)
        output = self.root / "scan"
        with self.assertLogs(benchmark.LOG, level="INFO"):
            main(["run", "--manifest", str(manifest), "--output", str(output),
                  "--configurations", "physical,map", "--grid", "5",
                  "--sweepfinder", executable])
        evaluation = self.root / "evaluation.json"
        main(["evaluate", "--summary", str(output / "summary.tsv"),
              "--output", str(evaluation)])
        report = json.loads(evaluation.read_text())
        self.assertEqual(len(report["results"]), 40)
        self.assertEqual({row["scenario"] for row in report["results"]}, {"neutral"})


if __name__ == "__main__":
    unittest.main()
