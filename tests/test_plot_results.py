"""Synthetic plotting and HTML integration tests; no inference binaries required."""

import base64
import contextlib
import io
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path

from matplotlib import pyplot as plt

from pysrc import plot_results as plots


class ReportParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images = []
        self.links = []
        self.rows = 0

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "img":
            self.images.append(attributes["src"])
        if tag == "a":
            self.links.append(attributes["href"])
        if tag == "tr":
            self.rows += 1


class PlotTests(unittest.TestCase):
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

    def test_psmc_uses_last_complete_iteration(self):
        path = self.write("curve.psmc", "RD\t0\nRS\t0\t0\t1\n//\n"
                          "RD\t1\nRS\t0\t0\t3\nRS\t1\t2\t4\n//\n"
                          "RD\t2\nRS\t0\t0\t99\n")
        figure = plots.plot_psmc(path)
        self.assertEqual(list(figure.axes[0].lines[0].get_ydata()), [3, 4])
        self.assertEqual(list(figure.axes[0].lines[0].get_xdata()), [0, 2])

    def test_msmc_keeps_cross_rate_and_native_units(self):
        path = self.write("msmc2/run.final.txt",
                          "time_index\tleft_time_boundary\tright_time_boundary\t"
                          "lambda_00\tlambda_01\tlambda_11\n"
                          "0\t0\t1\t4\t2\t8\n1\t1\tinf\t5\t3\t9\n")
        figure = plots.plot_msmc(path)
        self.assertEqual(len(figure.axes[0].lines), 3)
        self.assertEqual(list(figure.axes[0].lines[1].get_ydata()), [2, 3])
        self.assertIn("mutations per site", figure.axes[0].get_xlabel())

    def test_rccr_does_not_clip_values_above_one(self):
        path = self.write("rccr.tsv", "left_generations\trelative_cross_coalescence\n"
                          "0\t0.2\n100\t1.2\n")
        figure = plots.plot_rccr(path)
        self.assertGreater(figure.axes[0].get_ylim()[1], 1.2)

    def test_admixture_stacks_in_input_order(self):
        figure = plots.plot_admixture(self.write("test.2.Q", ".25 .75\n.8 .2\n"))
        ax = figure.axes[0]
        self.assertEqual([patch.get_height() for patch in ax.containers[0]], [.25, .8])
        self.assertEqual([patch.get_y() for patch in ax.containers[1]], [.25, .8])
        with self.assertRaisesRegex(ValueError, "sum to one"):
            plots.plot_admixture(self.write("bad.Q", ".2 .2\n"))

    def test_relatedness_preserves_negative_estimates_and_all_pairs(self):
        path = self.write("sample.kin0", "#IID1\tIID2\tKINSHIP\nA\tB\t-0.1\nA\tC\t0.25\n")
        figure = plots.plot_relatedness(path)
        self.assertEqual(sum(bar.get_height() for bar in figure.axes[0].patches), 2)
        self.assertLess(figure.axes[0].get_xlim()[0], -0.1)
        path = self.write("sample.ngsrelate.res", "a\tb\trab\ttheta\tKING\n0\t1\t.5\t.25\t.24\n")
        self.assertEqual(len(plots.plot_relatedness(path).axes), 3)

    def test_sweep_plots_use_coordinates_and_probabilities(self):
        sf = self.write("test.sf2.out", "location\tLR\talpha\n100\t2\t.1\n200\t5\t.2\n")
        self.assertEqual(list(plots.plot_sweepfinder(sf).axes[0].lines[0].get_xdata()), [100, 200])
        flex = self.write("flexsweep/predictions.tsv",
                          "chrom\tstart\tend\tprob_sweep\n1\t0\t100\t.9\n2\t0\t200\t.2\n")
        figure = plots.plot_flexsweep(flex)
        self.assertEqual(len(figure.axes[0].lines), 2)
        self.assertEqual(list(figure.axes[0].lines[0].get_xdata()), [50])

    def test_bed_empty_and_overlapping_intervals(self):
        figure = plots.plot_bed(self.write("empty.bed", ""))
        self.assertEqual(len(figure.axes[0].patches), 0)
        figure = plots.plot_bed(self.write("overlap.bed", "1\t0\t10\n1\t5\t15\n"))
        self.assertEqual(figure.axes[0].patches[0].get_height(), 2)
        self.assertEqual(figure.axes[1].patches[0].get_height(), 20)

    def test_input_and_iteration_diagnostics(self):
        path = self.write("input.multihetsep.txt", "1\t10\t10\tAC\n1\t15\t5\tCA\n")
        figure = plots.plot_multihetsep(path)
        self.assertEqual(figure.axes[0].patches[0].get_height(), 2)
        self.assertEqual(figure.axes[1].patches[0].get_height(), 15)
        path = self.write("msmc2/test.loop.txt", ".01\t-100\t0,inf\t2\n.02\t-90\t0,inf\t3\n")
        self.assertEqual(list(plots.plot_msmc_loop(path).axes[0].lines[0].get_ydata()), [-100, -90])
        path = self.write("test.normalise.tsv", "metric\tvalue\nretained\t10\nexcluded\t2\n")
        self.assertEqual(len(plots.plot_normalisation(path).axes[0].patches), 2)

    def test_fsc_folded_and_joint_spectra_and_parameters(self):
        path = self.write("folded.obs", "1 observations\nd0_0 d0_1 d0_2\n100 20 5\n")
        figure = plots.plot_sfs(path)
        self.assertEqual([bar.get_height() for bar in figure.axes[0].patches], [100, 20, 5])
        path = self.write("joint.obs", "1 observations\nd0_0 d0_1\n"
                          "d1_0 100 2\nd1_1 3 4\n")
        figure = plots.plot_sfs(path)
        self.assertEqual(figure.axes[0].images[0].get_array().tolist(), [[100, 2], [3, 4]])
        path = self.write("test.bestlhoods", "NPOP\tTSPLIT\tMaxEstLhood\n1000\t200\t-1234\n")
        self.assertEqual(list(plots.plot_fsc_parameters(path).axes[0].lines[0].get_ydata()),
                         [1000, 200, -1234])

    def test_report_embeds_plots_escapes_metadata_and_inventories_every_file(self):
        self.write("relatedness/a/admixture/test.2.Q", ".25 .75\n")
        self.write("relatedness/a/options.json", '{"note": "<script>alert(1)</script>"}')
        self.write("relatedness/a/model.keras", "not a real model")
        self.write("relatedness/a/quality.svg", '<svg xmlns="http://www.w3.org/2000/svg"></svg>')
        self.write("relatedness/a/a & b.log", "diagnostic")
        output, count, failures = plots.build_report(self.root)
        document = output.read_text()
        parser = ReportParser()
        parser.feed(document)
        self.assertEqual((count, failures), (2, 0))
        self.assertEqual(parser.rows, 6)
        self.assertEqual(len(parser.images), 2)
        for image in parser.images:
            self.assertTrue(image.startswith("data:image/"))
            self.assertTrue(base64.b64decode(image.split(",", 1)[1], validate=True))
        self.assertNotIn("<script>", document)
        self.assertIn("&lt;script&gt;", document)
        self.assertIn("relatedness/a/a%20%26%20b.log", parser.links)
        self.assertIn("Supporting file", document)
        # Reruns must not include the previous report in its own inventory.
        plots.build_report(self.root)
        self.assertNotIn('href="plots.html"', output.read_text())
        self.assertEqual(plt.get_fignums(), [])

    def test_report_records_errors_and_continues(self):
        self.write("test.2.Q", "nan 1\n")
        self.write("broken.json", "{")
        self.write("good.normalise.tsv", "metric\tvalue\nretained\t10\n")
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            output, count, failures = plots.build_report(self.root)
        self.assertEqual((count, failures), (1, 2))
        self.assertIn("ERROR", output.read_text())
        self.assertIn("test.2.Q", stderr.getvalue())

    def test_missing_empty_and_malformed_inputs_are_explicit_errors(self):
        with self.assertRaisesRegex(ValueError, "does not exist"):
            plots.build_report(self.root / "missing")
        with self.assertRaisesRegex(ValueError, "No result files"):
            plots.build_report(self.root)
        with self.assertRaisesRegex(ValueError, "no data"):
            plots.read_table(self.write("empty.tsv", "metric\tvalue\n"))
        with self.assertRaisesRegex(ValueError, "row width"):
            plots.read_table(self.write("ragged.tsv", "a\tb\n1\n"))
        with self.assertRaisesRegex(ValueError, "finite"):
            plots.number("inf")


if __name__ == "__main__":
    unittest.main()
