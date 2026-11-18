"""GHIST BED3 format, grid-region decisions, and held-out SF2 calibration."""

import csv
import json
from pathlib import Path
import tempfile
import unittest

from pysrc.sweep_experiments import benchmark, regions
from pysrc.sweep_experiments.calibration import empirical_p
from pysrc.sweep_experiments.__main__ import main


class RegionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.run = self.root / "run"
        self.case = self.run / "target"
        self.case.mkdir(parents=True)
        self.mask = self.root / "mask.bed"
        self.mask.write_text("7\t0\t100\n")
        self.scan = self.case / "physical.sf2.out"
        self.scan.write_text("location\tLR\n1\t1\n25\t30\n50\t25\n75\t1\n100\t40\n")
        self.provenance = {
            "status": "complete",
            "arguments": {"configurations": "physical", "ntv3": False,
                          "grid": 5, "unfolded": False, "constant_rate": 1e-8},
            "sweepfinder_binary": {"sha256": "binary"},
            "inputs": {
                "target:vcf": {"sha256": "target"},
                "target:mask": benchmark.file_metadata(self.mask)},
            "cases": {"target": {"chrom": "7", "start0": 0, "end0": 100,
                                 "split": "exploratory", "scenario": "unknown",
                                 "mask": str(self.mask)}}}
        self.summary = self.run / "summary.tsv"
        self.write_summary(self.summary, [("target", "unknown", "exploratory", "physical", 40)])
        self.save()

    def save(self):
        (self.run / "run.json").write_text(json.dumps(self.provenance))

    def write_summary(self, path, rows):
        with path.open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["case", "scenario", "split", "configuration", "score"])
            writer.writerows(rows)

    def null_source(self, values=range(1, 20)):
        directory = self.root / "null"
        directory.mkdir(exist_ok=True)
        summary = directory / "summary.tsv"
        self.write_summary(summary, [(f"null-{i}", "neutral", "calibration", "physical", v)
                                     for i, v in enumerate(values)])
        source = {
            "status": "complete", "arguments": dict(self.provenance["arguments"]),
            "sweepfinder_binary": {"sha256": "binary"},
            "inputs": {f"null-{i}:vcf": {"sha256": f"null-{i}"} for i, _ in enumerate(values)}}
        (directory / "run.json").write_text(json.dumps(source))
        return summary

    def exported(self, name):
        return self.run / "regions" / "target" / name

    def bed(self, name):
        rows = [line.split("\t") for line in self.exported(name).read_text().splitlines()]
        for row in rows:
            self.assertEqual(len(row), 3)
            self.assertEqual(row[0], "7")
            self.assertTrue(0 <= int(row[1]) < int(row[2]) <= 100)
        return [(int(row[1]), int(row[2])) for row in rows]

    def test_pvalues_use_region_maxima_conservative_ties_and_never_zero(self):
        null = list(range(1, 20))
        self.assertEqual(empirical_p(20, null), .05)
        self.assertEqual(empirical_p(19, null), .1)
        self.assertEqual(empirical_p(0, null), 1)
        for bad in (float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                empirical_p(bad, null)
        with self.assertRaises(ValueError):
            empirical_p(1, [])

    def test_calibrated_grid_merging_and_strict_ghist_bed3(self):
        result = regions.export_regions(self.run, null_summary=self.null_source())
        self.assertEqual(self.bed("physical.sweeps.bed"), [(13, 62), (87, 100)])
        self.assertEqual(self.bed("physical.candidates.bed"), [(13, 37), (87, 100)])
        evidence = regions.read_table(self.exported("physical.sf2.tsv"))
        self.assertEqual([float(row["empirical_p"]) for row in evidence], [1, .05, .05, 1, .05])
        self.assertEqual([int(row["significant"]) for row in evidence], [0, 1, 1, 0, 1])
        decisions = regions.read_table(self.exported("physical.regions.tsv"))
        predictions = [row for row in decisions if row["kind"] == "prediction"]
        self.assertEqual([float(row["peak_clr"]) for row in predictions], [30, 40])
        self.assertTrue(all(row["status"] == "calibrated" for row in predictions))
        self.assertTrue(all(row["glm_status"] == "disabled" for row in predictions))
        metadata = json.loads((self.run / "regions" / "metadata.json").read_text())
        self.assertEqual(metadata["status"], "complete")
        self.assertEqual(metadata["cases"]["target"]["prediction_counts"]["physical"], 2)
        self.assertEqual(result["alpha"], .05)

    def test_no_null_and_insufficient_nulls_do_not_manufacture_predictions(self):
        for source, status in ((None, "uncalibrated"),
                               (self.null_source([1, 2]), "insufficient_nulls")):
            with self.subTest(status=status), self.assertLogs(regions.LOG, level="WARNING"):
                regions.export_regions(self.run, null_summary=source)
            self.assertEqual(self.bed("physical.sweeps.bed"), [])
            self.assertTrue(self.bed("physical.candidates.bed"))
            rows = regions.read_table(self.exported("physical.sf2.tsv"))
            self.assertTrue(all(row["status"] == status for row in rows))
            if source is None:
                self.assertTrue(all(row["empirical_p"] == "." for row in rows))
            else:
                self.assertAlmostEqual(float(rows[1]["empirical_p"]), 1 / 3)

    def test_mask_gaps_split_regions_and_chr_prefix_is_removed(self):
        self.mask.write_text("chr7\t0\t30\nchr7\t40\t100\n")
        self.provenance["cases"]["target"]["chrom"] = "chr7"
        self.provenance["inputs"]["target:mask"] = benchmark.file_metadata(self.mask)
        self.save()
        regions.export_regions(self.run, null_summary=self.null_source())
        self.assertEqual(self.bed("physical.sweeps.bed"), [(13, 30), (40, 62), (87, 100)])

    def test_raw_glm_is_preserved_and_annotation_averages_per_snp(self):
        self.provenance["arguments"]["ntv3"] = True
        self.save()
        raw = self.case / "ntv3.tsv"
        raw.write_text("chrom\tposition\tstrand\toffset\tdelta_logp\tstatus\n"
                       "7\t20\t+\t0\t-1\tok\n7\t20\t-\t0\t3\tok\n"
                       "7\t30\t+\t0\t6\tok\n7\t50\t+\t0\t\tambiguous_context\n"
                       "7\t100\t+\t0\t2\tok\n7\t1\t\t\t\tcontig_edge\n")
        original = raw.read_bytes()
        regions.export_regions(self.run, null_summary=self.null_source())
        self.assertEqual(raw.read_bytes(), original)
        self.assertEqual(self.bed("ntv3.sites.bed"), [(19, 20), (29, 30), (99, 100)])
        rows = regions.read_table(self.exported("physical.regions.tsv"))
        first = next(row for row in rows if row["kind"] == "prediction")
        self.assertEqual(first["glm_snps"], "2")
        self.assertEqual(first["glm_rows"], "3")
        self.assertEqual(float(first["glm_mean_abs_delta"]), 4)
        self.assertEqual(float(first["glm_max_abs_delta"]), 6)
        metadata = json.loads((self.run / "regions" / "metadata.json").read_text())
        self.assertEqual(metadata["cases"]["target"]["raw_glm"], benchmark.file_metadata(raw))

    def test_glm_half_open_boundaries_and_empty_sites(self):
        raw = self.case / "ntv3.tsv"
        raw.write_text("chrom\tposition\tdelta_logp\tstatus\n"
                       "7\t1\t1\tok\n7\t10\t2\tok\n7\t11\t3\tok\n")
        glm = regions.GLMEvidence(raw, "7", self.case / "sites.bed")
        self.assertEqual(glm.annotate(0, 10)["glm_snps"], 2)
        self.assertEqual(glm.annotate(10, 11)["glm_mean_abs_delta"], 3)
        self.assertEqual(glm.annotate(20, 30)["glm_status"], "no_scored_sites")
        raw.write_text("chrom\tposition\tdelta_logp\tstatus\n7\t1\tnan\tok\n")
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            regions.GLMEvidence(raw, "7", self.case / "bad.bed")

    def test_calibration_cases_are_not_called_or_self_scored(self):
        self.provenance["cases"]["target"]["split"] = "calibration"
        self.save()
        regions.export_regions(self.run, null_summary=self.null_source())
        self.assertEqual(self.bed("physical.sweeps.bed"), [])
        rows = regions.read_table(self.exported("physical.sf2.tsv"))
        self.assertTrue(all(row["status"] == "calibration_only" and row["empirical_p"] == "."
                            for row in rows))

    def test_current_run_calibration_and_configuration_separation(self):
        rows = [("target", "unknown", "exploratory", "physical", 40)]
        for i in range(19):
            rows.extend([(f"null-{i}", "bgs", "calibration", "physical", i),
                         (f"null-{i}", "bgs", "calibration", "constant", 1000)])
            self.provenance["inputs"][f"null-{i}:vcf"] = {"sha256": f"null-{i}"}
        self.write_summary(self.summary, rows)
        self.save()
        regions.export_regions(self.run)
        self.assertEqual(self.bed("physical.sweeps.bed"), [(13, 62), (87, 100)])
        self.assertEqual(regions.null_maxima(regions.read_table(self.summary))["constant"], [1000] * 19)

    def test_null_summary_rejects_bad_truth_duplicates_and_nonfinite_scores(self):
        base = {"case": "null", "scenario": "neutral", "split": "calibration",
                "configuration": "physical", "score": "1"}
        for rows in ([base, base], [{**base, "scenario": "sweep"}],
                     [{**base, "score": "nan"}], [{**base, "score": "-1"}]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                regions.null_maxima(rows)
        self.assertEqual(regions.null_maxima([{**base, "split": "evaluation", "score": "999"}]), {})

    def test_external_controls_must_be_independent_and_match_scan_settings(self):
        summary = self.null_source()
        path = summary.parent / "run.json"
        source = json.loads(path.read_text())
        for key in ("grid", "constant_rate", "unfolded"):
            bad = json.loads(json.dumps(source))
            bad["arguments"][key] = "mismatch"
            path.write_text(json.dumps(bad))
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "mismatched"):
                regions.export_regions(self.run, null_summary=summary)
        for change in ("reuse_target", "duplicate_controls", "unfinished"):
            bad = json.loads(json.dumps(source))
            if change == "reuse_target":
                bad["inputs"]["null-0:vcf"]["sha256"] = "target"
            elif change == "duplicate_controls":
                bad["inputs"]["null-0:vcf"]["sha256"] = "null-1"
            else:
                bad["status"] = "running"
            path.write_text(json.dumps(bad))
            with self.subTest(change=change), self.assertRaises(ValueError):
                regions.export_regions(self.run, null_summary=summary)

    def test_bad_scan_and_changed_mask_fail_explicitly(self):
        for contents in ("location\tLR\n", "location\tLR\n1\tnan\n",
                         "location\tLR\n0\t1\n", "location\tLR\n10\t1\n9\t2\n",
                         "location\tLR\n10\t-1\n", "location\tLR\n101\t1\n"):
            self.scan.write_text(contents)
            with self.subTest(contents=contents), self.assertRaises(ValueError):
                regions.scan_cells(self.scan, 0, 100)
        self.mask.write_text("7\t0\t50\n")
        with self.assertRaisesRegex(ValueError, "mask changed"):
            regions.export_regions(self.run)

    def test_plateaus_fractional_positions_and_single_grid_point(self):
        self.scan.write_text("location\tLR\n1\t1\n1.1\t2\n1.2\t2\n100\t1\n")
        points = regions.scan_cells(self.scan, 0, 100)
        self.assertEqual(regions.peak_indices(points), [1])
        for point in points:
            self.assertTrue(0 <= point.start < point.end <= 100)
        self.scan.write_text("location\tLR\n50.5\t1\n")
        self.assertEqual(regions.scan_cells(self.scan, 0, 100), [regions.Point(50.5, 1, 0, 100)])

    def test_reexport_cli_replaces_calls_and_updates_provenance(self):
        summary = self.null_source()
        main(["export-regions", "--run-dir", str(self.run), "--null-summary", str(summary)])
        self.assertTrue(self.bed("physical.sweeps.bed"))
        main(["export-regions", "--run-dir", str(self.run), "--null-summary", str(summary),
              "--alpha", "0.01"])
        self.assertEqual(self.bed("physical.sweeps.bed"), [])
        provenance = json.loads((self.run / "run.json").read_text())
        self.assertEqual(provenance["region_exports"]["alpha"], .01)
        for alpha in (0, 1, float("nan")):
            with self.subTest(alpha=alpha), self.assertRaises(ValueError):
                regions.export_regions(self.run, alpha)
        self.provenance["status"] = "running"
        self.save()
        with self.assertRaisesRegex(ValueError, "complete run"):
            regions.export_regions(self.run, require_complete=True)


if __name__ == "__main__":
    unittest.main()
