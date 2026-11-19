"""Numerical, format, real bcftools/SF2, and optional simulation integration tests."""

import csv
from contextlib import redirect_stderr
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from pysrc.sweep_experiments import benchmark, ntv3, inputs, simulation
from pysrc.sweep_experiments.__main__ import main, parser


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

        output = self.root / "ntv3.tsv"
        count = ntv3.score_sites([inputs.Variant("7", 17, "A", "G", 1, 4, True)],
                                reference, Predictor(), output, 8, [2, 5])
        self.assertEqual(count, 4)
        self.assertEqual([target for _, target in calls], [2, 2, 5, 5])
        with output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertAlmostEqual(float(rows[0]["delta_logp"]), math.log(3))
        self.assertAlmostEqual(float(rows[1]["delta_logp"]), math.log(.2 / .4))
        self.assertAlmostEqual(float(rows[1]["p_A"]), .4)
        self.assertEqual(rows[1]["strand"], "-")
        for window in ntv3.windows(17, 40, 8, [2, 5]):
            genomic_target = (window.start + window.target if window.strand == "+"
                              else window.end - 1 - window.target)
            self.assertEqual(genomic_target, 16)
        self.assertEqual(ntv3.windows(1, 40, 8, [2]), [])
        with self.assertRaisesRegex(ValueError, "offset"):
            ntv3.windows(17, 40, 8, [-1])
        for offsets in ([], [2, 2], [8]):
            with self.subTest(offsets=offsets):
                with self.assertRaises(ValueError):
                    ntv3.windows(17, 40, 8, offsets)

    def test_masked_offsets_include_both_window_endpoints(self):
        windows = ntv3.windows(17, 40, 8, [0, 7])
        self.assertEqual([w.target for w in windows], [0, 0, 7, 7])
        for window in windows:
            genomic_target = (window.start + window.target if window.strand == "+"
                              else window.end - 1 - window.target)
            self.assertEqual(genomic_target, 16)

    def test_scoring_progress_counts_variants_including_skipped_contexts(self):
        reference = SequenceReference("ACGT" * 7 + "N" * 12)
        sites = [inputs.Variant("7", position, "A", "G", 1, 4, True)
                 for position in (17, 1, 33)]
        predictor = SimpleNamespace(log_probs=Mock(
            return_value={base: math.log(.25) for base in "ACGT"}))
        output = self.root / "ntv3.tsv"
        progress = io.StringIO()

        with redirect_stderr(progress):
            count = ntv3.score_sites(sites, reference, predictor, output, 8, [2])

        self.assertIn("Scoring SNPs", progress.getvalue())
        self.assertIn("3/3", progress.getvalue())
        self.assertIn("variant", progress.getvalue())
        self.assertEqual(count, 2)
        self.assertEqual(predictor.log_probs.call_count, 2)
        with output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["status"] for row in rows],
                         ["ok", "ok", "contig_edge", "ambiguous_context", "ambiguous_context"])

    def test_scoring_accepts_mismatching_ref_and_preserves_simulated_alleles(self):
        reference = SequenceReference("ACGT" * 10)
        calls = []

        def predict(sequence, target):
            calls.append((sequence, target))
            return {b: math.log(p) for b, p in zip("ACGT", (.1, .2, .3, .4))}

        output = self.root / "mismatching-ref.tsv"
        count = ntv3.score_sites(
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

    def test_ntv3_cli_and_early_context_validation(self):
        required = ["run", "--manifest", str(self.root / "missing.json"),
                    "--output", str(self.root / "output")]
        args = parser().parse_args(required)
        self.assertFalse(args.ntv3)
        self.assertFalse(args.hyenadna)
        self.assertIsNone(args.model_revision)
        with patch("sys.stderr"):
            for legacy in (["--evo"], ["--allow-rocm"], ["--model", "evo2_7b"],
                           ["--evo-stride", "1"]):
                with self.subTest(legacy=legacy), self.assertRaises(SystemExit):
                    parser().parse_args(required + legacy)
        manifest = self.write("manifest.json", json.dumps({"cases": [{
            "id": "test", "vcf": "sites.vcf", "split": "exploratory"}]}))
        args.manifest = manifest
        args.ntv3 = True
        with self.assertRaisesRegex(ValueError, "--seqrepo"):
            benchmark.run(args)
        args.seqrepo = self.root
        for context in (0, 127, 129, 1_048_704):
            args.context = context
            with self.subTest(context=context), self.assertRaisesRegex(ValueError, "multiple of 128"):
                benchmark.run(args)
        self.assertFalse(args.output.exists())

    def test_ntv3_runner_records_new_outputs_provenance_and_stride(self):
        self._check_glm_runner("ntv3")

    def test_hyenadna_runner_records_distinct_outputs_and_exports(self):
        self._check_glm_runner("hyenadna")

    def test_hyenadna_cli_rejects_conflicting_flags_and_invalid_context(self):
        manifest = self.write("manifest.json", json.dumps({"cases": [{
            "id": "test", "vcf": "sites.vcf", "split": "exploratory"}]}))
        required = ["run", "--manifest", str(manifest), "--output", str(self.root / "run"),
                    "--hyenadna"]
        with patch("sys.stderr"), self.assertRaises(SystemExit):
            parser().parse_args(required + ["--ntv3"])
        with self.assertRaisesRegex(ValueError, "--seqrepo"):
            main(required)
        for flags, message in ((["--context", "1000001"], "context must be"),
                               (["--offsets", "0,64"], "offsets must satisfy"),
                               (["--offsets", "64,64"], "unique")):
            with self.subTest(flags=flags), self.assertRaisesRegex(ValueError, message):
                main(required + ["--seqrepo", str(self.root), *flags])
        self.assertFalse((self.root / "run").exists())

    def _check_glm_runner(self, backend):
        from pysrc.sweep_experiments import hyenadna

        manifest = self.write("manifest.json", json.dumps({"cases": [{
            "id": name, "vcf": "sites.vcf", "mask": "mask.bed", "dataset": "testing",
            "scenario": "unknown", "split": "exploratory"} for name in ("first", "second")]}))
        self.write("sites.vcf", "fixture")
        self.write("mask.bed", "7\t0\t2000\n")
        executable = self.write("sf2", "fixture")
        sites = [inputs.Variant("7", p, "A", "G", 1, 4, True) for p in (257, 513, 769)]
        reference = SequenceReference("ACGT" * 512)
        module = ntv3 if backend == "ntv3" else hyenadna
        metadata = {"model": module.MODEL, "revision": module.REVISION}
        predictor_name = "NTv3Predictor" if backend == "ntv3" else "HyenaDNAPredictor"
        offsets = [0, 64, 127] if backend == "ntv3" else [1, 64, 127]
        stride = "--ntv3-stride" if backend == "ntv3" else "--glm-stride"
        calls = []

        class Predictor:
            def log_probs(self, sequence, target):
                calls.append((sequence, target))
                return {b: math.log(p) for b, p in zip("ACGT", (.1, .2, .3, .4))}

            def metadata(self):
                return metadata

        def run_command(command, log):
            Path(command[-1]).write_text("location\tLR\n500\t2\n")

        with patch.object(benchmark.shutil, "which", return_value=str(executable)), \
                patch.object(benchmark.subprocess, "run",
                             return_value=SimpleNamespace(stdout="bcftools test\n")), \
                patch.object(benchmark, "variants", return_value=sites), \
                patch.object(benchmark, "Reference", return_value=reference), \
                patch.object(benchmark, "run_command", side_effect=run_command), \
                patch.object(module, predictor_name, return_value=Predictor()) as load:
            main(["run", "--manifest", str(manifest), "--output", str(self.root / "run"),
                  "--configurations", "physical", f"--{backend}", "--seqrepo", str(self.root),
                  "--context", "128", "--offsets", ",".join(map(str, offsets)), stride, "2"])
        if backend == "ntv3":
            load.assert_called_once_with(ntv3.MODEL, ntv3.REVISION, 128)
        else:
            load.assert_called_once_with(hyenadna.REVISION, 128, offsets)
        self.assertEqual(len(calls), 24)
        output = self.root / "run"
        raw = output / "first" / f"{backend}.tsv"
        with raw.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["position"] for row in rows}, {"257", "769"})
        self.assertEqual({row["status"] for row in rows}, {"ok"})
        self.assertFalse((output / "first" / "evo.tsv").exists())
        provenance = json.loads((output / "run.json").read_text())
        self.assertEqual(provenance["status"], "complete")
        self.assertEqual(provenance[backend], metadata)
        self.assertEqual(provenance["arguments"]["model"], module.MODEL)
        self.assertEqual(provenance["arguments"]["model_revision"], module.REVISION)
        regions = output / "regions" / "first"
        self.assertTrue((regions / "physical.candidates.bed").is_file())
        self.assertEqual((regions / "physical.sweeps.bed").read_text(), "")
        self.assertEqual(len((regions / f"{backend}.sites.bed").read_text().splitlines()), 2)
        for line in (regions / "physical.candidates.bed").read_text().splitlines():
            self.assertEqual(len(line.split("\t")), 3)
        export_metadata = json.loads((output / "regions" / "metadata.json").read_text())
        self.assertEqual(export_metadata["status"], "complete")
        self.assertEqual(export_metadata["cases"]["first"]["raw_glm"]["sha256"],
                         benchmark.file_metadata(raw)["sha256"])
        with (output / "summary.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["configuration"] for row in rows],
                         ["physical", f"{backend}_mean_abs_delta"] * 2)
        self.assertAlmostEqual(float(rows[1]["score"]), (math.log(3) + math.log(2)) / 2)

    def test_runner_automatically_exports_calibrated_ghist_sweep_regions(self):
        cases = []
        for i in range(20):
            name = f"case-{i}"
            self.write(f"{name}.vcf", name)
            cases.append({"id": name, "vcf": f"{name}.vcf", "mask": "mask.bed",
                          "dataset": "testing", "scenario": "neutral" if i < 19 else "unknown",
                          "split": "calibration" if i < 19 else "exploratory"})
        manifest = self.write("manifest.json", json.dumps({"cases": cases}))
        self.write("mask.bed", "7\t0\t1000\n")
        executable = self.write("sf2", "fixture")
        sites = [inputs.Variant("7", p, "A", "G", 1, 4, True) for p in (100, 500, 900)]

        def run_command(command, log):
            path = Path(command[-1])
            index = int(path.parent.name.split("-")[-1])
            path.write_text(f"location\tLR\n100\t1\n500\t{index + 1}\n900\t1\n")

        with patch.object(benchmark.shutil, "which", return_value=str(executable)), \
                patch.object(benchmark.subprocess, "run",
                             return_value=SimpleNamespace(stdout="bcftools test\n")), \
                patch.object(benchmark, "variants", return_value=sites), \
                patch.object(benchmark, "run_command", side_effect=run_command):
            main(["run", "--manifest", str(manifest), "--output", str(self.root / "run"),
                  "--configurations", "physical"])
        output = self.root / "run"
        self.assertEqual((output / "regions" / "case-19" / "physical.sweeps.bed").read_text(),
                         "7\t300\t700\n")
        for i in range(19):
            self.assertEqual((output / "regions" / f"case-{i}" / "physical.sweeps.bed").read_text(), "")
        with (output / "regions" / "case-19" / "physical.regions.tsv").open() as handle:
            row = next(row for row in csv.DictReader(handle, delimiter="\t")
                       if row["kind"] == "prediction")
        self.assertEqual(float(row["empirical_p"]), .05)
        self.assertEqual(row["n_null"], "19")
        self.assertEqual(row["glm_status"], "disabled")
        self.assertEqual(json.loads((output / "run.json").read_text())["status"], "complete")

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
        for configuration in ("physical", "constant", "map"):
            directory = output / "regions" / "test"
            candidates = directory / f"{configuration}.candidates.bed"
            self.assertTrue(candidates.read_text().strip())
            for line in candidates.read_text().splitlines():
                chrom, start, end = line.split("\t")
                self.assertEqual(chrom, "7")
                self.assertTrue(99 <= int(start) < int(end) <= 500)
            self.assertEqual((directory / f"{configuration}.sweeps.bed").read_text(), "")
            with (directory / f"{configuration}.sf2.tsv").open() as handle:
                evidence = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(evidence)
            self.assertTrue(all(row["status"] == "uncalibrated" and row["empirical_p"] == "."
                                for row in evidence))
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
