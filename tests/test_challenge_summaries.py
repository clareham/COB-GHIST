"""Numerical scaling, exploratory-estimator, and combined-report regressions."""

import contextlib
import io
import json
import math
import tempfile
import unittest
from pathlib import Path

from pysrc.plotting import demographic_curves as curves
from pysrc.plotting import genetic_history, relatedness, report, sweep_detection
from pysrc.plotting.common import plt


class SummaryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(temporary.name)

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def metadata(self, challenge="growth", rate=1e-8):
        return self.write(f"genetic_history/{challenge}/msmc2/{challenge}.msmc2.run/assumptions.json",
                          json.dumps({"options": {"mutation_rate": rate}}))

    def psmc(self):
        return self.write("genetic_history/growth/psmc/growth.psmc",
                          "RD\t0\nTR\t.08\t.01\nRS\t0\t0\t99\n//\n"
                          "RD\t1\nTR\t.04\t.01\nRS\t0\t0\t2\nRS\t1\t2\t3\n//\n"
                          "RD\t2\nTR\t.9\t.01\nRS\t0\t0\t99\n")

    def msmc(self, challenge="growth"):
        return self.write(f"genetic_history/{challenge}/msmc2/{challenge}.msmc2.final.txt",
                          "time_index\tleft_time_boundary\tright_time_boundary\t"
                          "lambda_00\tlambda_01\tlambda_11\n"
                          "0\t0\t.00001\t1000\t100\t2000\n"
                          "1\t.00001\tinf\t1000\t500\t2000\n")

    def test_psmc_scaling_uses_final_theta_and_100bp_bins(self):
        path = self.psmc()
        self.metadata()
        curve = curves.scaled_psmc(path, curves.mutation_rate(path))
        self.assertAlmostEqual(curve.left[1], 40000)
        self.assertAlmostEqual(curve.size[0], 20000)
        self.assertAlmostEqual(curve.size[1], 30000)
        self.assertTrue(math.isinf(curve.right[-1]))
        figure = genetic_history.plot_psmc(path)
        self.assertIn("Generations", figure.axes[0].get_xlabel())
        self.assertAlmostEqual(figure.axes[0].lines[0].get_ydata()[1], 30000)

    def test_psmc_bin_metadata_and_missing_theta(self):
        path = self.psmc()
        self.write(str(path.relative_to(self.root)) + ".assumptions.json", '{"window_size_bp": 50}')
        self.assertAlmostEqual(curves.scaled_psmc(path, 1e-8).left[1], 80000)
        bad = self.write("bad.psmc", "RD\t0\nRS\t0\t0\t1\n//\n")
        with self.assertRaisesRegex(ValueError, "TR theta"):
            curves.scaled_psmc(bad, 1e-8)

    def test_msmc_scaling_never_inverts_cross_population_rate(self):
        path = self.msmc()
        scaled = curves.scaled_msmc(path, 1e-8)
        self.assertEqual(len(scaled), 2)
        self.assertAlmostEqual(scaled[0].left[1], 1000)
        self.assertAlmostEqual(scaled[0].size[0], 50000)
        self.assertAlmostEqual(scaled[1].size[0], 25000)
        self.assertTrue(all("lambda_01" not in curve.label for curve in scaled))
        self.metadata()
        self.assertAlmostEqual(genetic_history.plot_msmc(path).axes[0].lines[0].get_xdata()[1], 1000)

    def test_scaling_missing_or_conflicting_metadata_is_explicit(self):
        path = self.psmc()
        self.assertIsNone(curves.mutation_rate(path))
        self.metadata()
        self.write(str(path.relative_to(self.root)) + ".assumptions.json", '{"mutation_rate": 2e-8}')
        with self.assertRaisesRegex(ValueError, "conflicting"):
            curves.mutation_rate(path)
        self.metadata(rate=-1)
        with self.assertRaisesRegex(ValueError, "positive"):
            curves.mutation_rate(path)

    def test_scaling_rejects_nonfinite_rates_and_nonobject_metadata(self):
        path = self.psmc()
        for rate in (0, -1, math.inf, math.nan):
            with self.subTest(rate=rate):
                with self.assertRaisesRegex(ValueError, "positive"):
                    curves.scaled_psmc(path, rate)
                with self.assertRaisesRegex(ValueError, "positive"):
                    curves.scaled_msmc(self.msmc(), rate)
        metadata = self.metadata()
        metadata.write_text("[]")
        with self.assertRaisesRegex(ValueError, "JSON object"):
            curves.mutation_rate(path)
        metadata.write_text('{"options": []}')
        with self.assertRaisesRegex(ValueError, "JSON object"):
            curves.mutation_rate(path)

    def test_known_two_epoch_growth_and_decline_targets(self):
        for growth, sizes, expected in ((True, [400] * 4 + [100] * 4, 4),
                                        (False, [100] * 4 + [400] * 4, .25)):
            with self.subTest(growth=growth):
                curve = curves.Curve(list(range(8)), list(range(1, 9)), sizes, "test", Path("test"))
                estimates = curves.change_estimates(curve, growth, "test")
                self.assertEqual(estimates[0].value, 4)
                self.assertAlmostEqual(estimates[1].value, expected)
                self.assertEqual(estimates[0].status, "exploratory")

    def test_change_fit_rejects_flat_wrong_direction_and_infinite_terminal_bin(self):
        curve = curves.Curve(list(range(8)), list(range(1, 8)) + [math.inf],
                             [100] * 8, "test", Path("test"))
        self.assertTrue(all(item.value is None for item in curves.change_estimates(curve, True, "test")))
        curve.size = [100] * 4 + [400] * 4
        self.assertTrue(all(item.value is None for item in curves.change_estimates(curve, True, "test")))
        curve.right = [math.inf] * 8
        self.assertIn("Insufficient", curves.change_estimates(curve, False, "test")[0].method)

    def test_change_fit_enforces_twenty_percent_directional_threshold(self):
        for growth, recent, available in ((True, 119, False), (True, 120, True),
                                           (False, 81, False), (False, 80, True)):
            with self.subTest(growth=growth, recent=recent):
                curve = curves.Curve(list(range(8)), list(range(1, 9)),
                                     [recent] * 4 + [100] * 4, "test", Path("test"))
                estimate = curves.change_estimates(curve, growth, "test")[0]
                self.assertEqual(estimate.value is not None, available)

    def test_rccr_half_crossing_is_sustained_and_bracketed(self):
        rows = [{"left_generations": str(time), "relative_cross_coalescence": str(value)}
                for time, value in ((0, .1), (10, .6), (20, .4), (30, .7), (40, .9))]
        estimate = curves.rccr_split(rows, "test")
        self.assertEqual(estimate.value, 30)
        self.assertEqual(estimate.interval, (20, 30))
        self.assertIn("not a confidence", estimate.method)
        for values in ((.6, .7, .8), (.1, .2, .4), (.1, .2, .6)):
            rows = [{"left_generations": str(i), "relative_cross_coalescence": str(value)}
                    for i, value in enumerate(values)]
            self.assertIsNone(curves.rccr_split(rows, "test").value)

    def test_genetic_summary_separates_population_and_cross_curves(self):
        path = self.msmc("secondary_contact")
        self.metadata("secondary_contact")
        self.write("genetic_history/secondary_contact/msmc2/secondary_contact.msmc2.run/rccr.tsv",
                   "left_generations\trelative_cross_coalescence\n0\t.1\n10\t.7\n20\t.8\n")
        summary = genetic_history.compose_challenge(path.parents[1])
        self.assertEqual(len(summary.figure.axes), 2)
        self.assertEqual(len(summary.figure.axes[0].lines), 2)
        self.assertEqual(summary.figure.axes[0].get_xlim()[0], 0)
        self.assertEqual(summary.figure.axes[1].get_xlim()[0], 0)
        estimates = {item.target: item for item in summary.estimates}
        self.assertAlmostEqual(estimates["mainland_population_size"].value, 50000)
        self.assertAlmostEqual(estimates["island_population_size"].value, 25000)
        self.assertIsNone(estimates["migration_rate"].value)
        self.assertIsNone(estimates["generations_since_migration"].value)
        self.assertEqual(estimates["generations_since_split"].value, 10)

    def test_relatedness_threshold_boundaries(self):
        for value, label in ((.354, "duplicate/MZ-like"), (.177, "first degree"),
                             (.0884, "second degree"), (.0442, "third degree"),
                             (-.01, "more distant/unrelated")):
            self.assertEqual(relatedness.relationship_degree(value), label)

    def test_relatedness_combined_key_and_named_pair_assignments(self):
        self.write("relatedness/relatedness3/admixture/test.2.Q", ".25 .75\n.8 .2\n")
        path = self.write("relatedness/relatedness3/king_plink/test.kin0",
                          "#IID1\tIID2\tKINSHIP\nA\tB\t.25\nA\tC\t-.1\n")
        self.write("relatedness/relatedness3/ngsrelate/test.ngsrelate.res",
                   "a\tb\trab\ttheta\tKING\n0\t1\t.5\t.25\t.24\n")
        summary = relatedness.compose_challenge(path.parents[1])
        self.assertEqual(len(summary.figure.axes), 3)
        self.assertEqual(len(summary.assignments), 2)
        self.assertEqual(summary.assignments[0]["individual_1"], "A")
        self.assertEqual(summary.assignments[0]["degree_proxy"], "first degree")
        self.assertEqual(summary.figure.axes[0].get_legend_handles_labels()[1],
                         ["Component 1", "Component 2"])
        self.assertTrue(any("PLINK KINSHIP" == label for label in
                            summary.figure.axes[1].get_legend_handles_labels()[1]))

    def test_sweep_interval_union_peak_and_unidentifiable_targets(self):
        self.write("sweep_detection/onesweep/sweepfinder2/a.sf2.out",
                   "location\tLR\talpha\n100\t2\t1\n200\t5\t2\n")
        path = self.write("sweep_detection/onesweep/flexsweep/run.flexsweep/predictions.tsv",
                          "chrom\tstart\tend\tprob_sweep\n7\t0\t100\t.95\n"
                          "7\t50\t150\t.99\n7\t500\t600\t.92\n")
        bed = self.write("sweep_detection/onesweep/flexsweep/run.flexsweep/sweeps.bed",
                         "7\t0\t100\n7\t50\t150\n7\t500\t600\n")
        self.assertEqual(sweep_detection.merged_intervals(bed), [("7", 0, 150), ("7", 500, 600)])
        summary = sweep_detection.compose_challenge(path.parents[2])
        estimates = {item.target: item for item in summary.estimates}
        self.assertEqual(estimates["maximum_CLR_position"].value, 200)
        self.assertEqual(estimates["candidate_sweep_regions"].value, 2)
        self.assertEqual(estimates["candidate_region_1_position"].value, 100)
        self.assertEqual(estimates["candidate_region_1_position"].interval, (0, 150))
        self.assertIsNone(estimates["sweep_time"].value)
        self.assertIsNone(estimates["selection_coefficient"].value)
        self.assertEqual(len(summary.assignments), 2)
        self.assertEqual(len(summary.figure.axes), 2)

    def test_empty_sweep_bed_means_zero_candidates_not_missing_data(self):
        path = self.write("sweep_detection/a/flexsweep/a.flexsweep/predictions.tsv",
                          "chrom\tstart\tend\tprob_sweep\n7\t0\t100\t.1\n")
        self.write("sweep_detection/a/flexsweep/a.flexsweep/sweeps.bed", "")
        summary = sweep_detection.compose_challenge(path.parents[2])
        count = next(item for item in summary.estimates if item.target == "candidate_sweep_regions")
        self.assertEqual(count.value, 0)
        self.assertEqual(summary.assignments, [])

    def test_report_exports_reproducible_artifacts_and_excludes_them_on_rerun(self):
        self.write("relatedness/a/admixture/a.2.Q", ".25 .75\n")
        output, count, failures = report.build_report(self.root)
        self.assertEqual((count, failures), (2, 0))
        self.assertTrue((self.root / "plot_summaries/relatedness/a.png").read_bytes().startswith(b"\x89PNG"))
        records = json.loads((self.root / "challenge_estimates.json").read_text())
        self.assertIn("relatedness/a", records["challenges"])
        self.assertIn("Challenge overview", output.read_text())
        self.assertIn("Download combined graphic", output.read_text())
        report.build_report(self.root)
        self.assertIn("1 files; 2 plots; 0 errors", output.read_text())

    def test_summary_failure_remains_visible_and_other_challenges_continue(self):
        path = self.psmc()
        self.metadata()
        self.write(str(path.relative_to(self.root)) + ".assumptions.json", '{"mutation_rate": 2e-8}')
        self.write("relatedness/a/admixture/a.2.Q", ".25 .75\n")
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            output, _, failures = report.build_report(self.root)
        self.assertGreater(failures, 0)
        self.assertIn("Challenge overview failed", output.read_text())
        self.assertIn("conflicting", stderr.getvalue())
        records = json.loads((self.root / "challenge_estimates.json").read_text())
        self.assertIn("error", records["challenges"]["genetic_history/growth"])
        self.assertIn("estimates", records["challenges"]["relatedness/a"])


if __name__ == "__main__":
    unittest.main()
