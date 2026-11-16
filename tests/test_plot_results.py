"""Synthetic plotting and HTML integration tests; no inference binaries required."""

import base64
import contextlib
import io
import subprocess
import sys
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path

from matplotlib import pyplot as plt

from pysrc.plotting import common, genetic_history, relatedness, report, sweep_detection


class ReportParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images = []
        self.links = []
        self.rows = 0
        self.ancestors = []
        self.headings = []
        self.ids = []
        self.heading = None

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag in {"section", "article"}:
            self.ancestors.append(attributes["data-path"])
        if "id" in attributes:
            self.ids.append(attributes["id"])
        if tag in {"h1", "h2", "h3", "h4", "h5", "h6"} or attributes.get("role") == "heading":
            level = int(attributes["aria-level"]) if tag == "div" else int(tag[1])
            self.heading = (level, tuple(self.ancestors), [])
        if tag == "img":
            self.images.append(attributes["src"])
        if tag == "a":
            self.links.append(attributes["href"])
        if tag == "tr":
            self.rows += 1

    def handle_data(self, data):
        if self.heading is not None:
            self.heading[2].append(data)

    def handle_endtag(self, tag):
        if tag in {"section", "article"}:
            self.ancestors.pop()
        if self.heading is not None and tag in {"h1", "h2", "h3", "h4", "h5", "h6", "div"}:
            level, ancestors, text = self.heading
            self.headings.append((level, ancestors, "".join(text)))
            self.heading = None


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

    def test_cli_supports_module_and_direct_script_execution(self):
        repository = Path(__file__).resolve().parents[1]
        results = self.root / "results"
        self.write("results/test.normalise.tsv", "metric\tvalue\nretained\t10\n")
        invocations = [
            (repository, ["-m", "pysrc.plot_results"]),
            (repository, ["pysrc/plot_results.py"]),
            (repository / "pysrc", ["plot_results.py"]),
        ]
        for index, (cwd, arguments) in enumerate(invocations):
            with self.subTest(arguments=arguments, cwd=cwd):
                output = self.root / f"report-{index}.html"
                completed = subprocess.run(
                    [sys.executable, *arguments, "--results-dir", str(results),
                     "--output", str(output)],
                    cwd=cwd, capture_output=True, text=True, timeout=60,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertIn("1 plots, 0 errors", completed.stdout)
                self.assertIn("test.normalise.tsv", output.read_text())

    def test_psmc_uses_last_complete_iteration(self):
        path = self.write("curve.psmc", "RD\t0\nRS\t0\t0\t1\n//\n"
                          "RD\t1\nRS\t0\t0\t3\nRS\t1\t2\t4\n//\n"
                          "RD\t2\nRS\t0\t0\t99\n")
        figure = genetic_history.plot_psmc(path)
        self.assertEqual(list(figure.axes[0].lines[0].get_ydata()), [3, 4])
        self.assertEqual(list(figure.axes[0].lines[0].get_xdata()), [0, 2])

    def test_msmc_keeps_cross_rate_and_native_units(self):
        path = self.write("msmc2/run.final.txt",
                          "time_index\tleft_time_boundary\tright_time_boundary\t"
                          "lambda_00\tlambda_01\tlambda_11\n"
                          "0\t0\t1\t4\t2\t8\n1\t1\tinf\t5\t3\t9\n")
        figure = genetic_history.plot_msmc(path)
        self.assertEqual(len(figure.axes[0].lines), 3)
        self.assertEqual(list(figure.axes[0].lines[1].get_ydata()), [2, 3])
        self.assertIn("mutations per site", figure.axes[0].get_xlabel())

    def test_rccr_does_not_clip_values_above_one(self):
        path = self.write("rccr.tsv", "left_generations\trelative_cross_coalescence\n"
                          "0\t0.2\n100\t1.2\n")
        figure = genetic_history.plot_rccr(path)
        self.assertGreater(figure.axes[0].get_ylim()[1], 1.2)

    def test_admixture_stacks_in_input_order(self):
        figure = relatedness.plot_admixture(self.write("test.2.Q", ".25 .75\n.8 .2\n"))
        ax = figure.axes[0]
        self.assertEqual([patch.get_height() for patch in ax.containers[0]], [.25, .8])
        self.assertEqual([patch.get_y() for patch in ax.containers[1]], [.25, .8])
        with self.assertRaisesRegex(ValueError, "sum to one"):
            relatedness.plot_admixture(self.write("bad.Q", ".2 .2\n"))

    def test_relatedness_preserves_negative_estimates_and_all_pairs(self):
        path = self.write("sample.kin0", "#IID1\tIID2\tKINSHIP\nA\tB\t-0.1\nA\tC\t0.25\n")
        figure = relatedness.plot_relatedness(path)
        self.assertEqual(sum(bar.get_height() for bar in figure.axes[0].patches), 2)
        self.assertLess(figure.axes[0].get_xlim()[0], -0.1)
        path = self.write("sample.ngsrelate.res", "a\tb\trab\ttheta\tKING\n0\t1\t.5\t.25\t.24\n")
        self.assertEqual(len(relatedness.plot_relatedness(path).axes), 3)

    def test_sweep_plots_use_coordinates_and_probabilities(self):
        sf = self.write("test.sf2.out", "location\tLR\talpha\n100\t2\t.1\n200\t5\t.2\n")
        self.assertEqual(list(sweep_detection.plot_sweepfinder(sf).axes[0].lines[0].get_xdata()),
                         [100, 200])
        flex = self.write("flexsweep/predictions.tsv",
                          "chrom\tstart\tend\tprob_sweep\n1\t0\t100\t.9\n2\t0\t200\t.2\n")
        figure = sweep_detection.plot_flexsweep(flex)
        self.assertEqual(len(figure.axes[0].lines), 2)
        self.assertEqual(list(figure.axes[0].lines[0].get_xdata()), [50])

    def test_bed_empty_and_overlapping_intervals(self):
        figure = common.plot_bed(self.write("empty.bed", ""))
        self.assertEqual(len(figure.axes[0].patches), 0)
        figure = common.plot_bed(self.write("overlap.bed", "1\t0\t10\n1\t5\t15\n"))
        self.assertEqual(figure.axes[0].patches[0].get_height(), 2)
        self.assertEqual(figure.axes[1].patches[0].get_height(), 20)

    def test_input_and_iteration_diagnostics(self):
        path = self.write("input.multihetsep.txt", "1\t10\t10\tAC\n1\t15\t5\tCA\n")
        figure = genetic_history.plot_multihetsep(path)
        self.assertEqual(figure.axes[0].patches[0].get_height(), 2)
        self.assertEqual(figure.axes[1].patches[0].get_height(), 15)
        path = self.write("msmc2/test.loop.txt", ".01\t-100\t0,inf\t2\n.02\t-90\t0,inf\t3\n")
        self.assertEqual(list(genetic_history.plot_msmc_loop(path).axes[0].lines[0].get_ydata()),
                         [-100, -90])
        path = self.write("test.normalise.tsv", "metric\tvalue\nretained\t10\nexcluded\t2\n")
        self.assertEqual(len(common.plot_normalisation(path).axes[0].patches), 2)

    def test_fsc_folded_and_joint_spectra_and_parameters(self):
        path = self.write("folded.obs", "1 observations\nd0_0 d0_1 d0_2\n100 20 5\n")
        figure = genetic_history.plot_sfs(path)
        self.assertEqual([bar.get_height() for bar in figure.axes[0].patches], [100, 20, 5])
        path = self.write("joint.obs", "1 observations\nd0_0 d0_1\n"
                          "d1_0 100 2\nd1_1 3 4\n")
        figure = genetic_history.plot_sfs(path)
        self.assertEqual(figure.axes[0].images[0].get_array().tolist(), [[100, 2], [3, 4]])
        path = self.write("test.bestlhoods", "NPOP\tTSPLIT\tMaxEstLhood\n1000\t200\t-1234\n")
        self.assertEqual(list(genetic_history.plot_fsc_parameters(path).axes[0].lines[0].get_ydata()),
                         [1000, 200, -1234])

    def test_report_embeds_plots_escapes_metadata_and_inventories_every_file(self):
        self.write("relatedness/a/admixture/test.2.Q", ".25 .75\n")
        self.write("relatedness/a/options.json", '{"note": "<script>alert(1)</script>"}')
        self.write("relatedness/a/model.keras", "not a real model")
        self.write("relatedness/a/quality.svg", '<svg xmlns="http://www.w3.org/2000/svg"></svg>')
        self.write("relatedness/a/a & b.log", "diagnostic")
        output, count, failures = report.build_report(self.root)
        document = output.read_text()
        parser = ReportParser()
        parser.feed(document)
        self.assertEqual((count, failures), (3, 0))
        self.assertEqual(document.count("Supporting file (linked; no plot adapter)"), 4)
        self.assertEqual(len(parser.images), 3)
        for image in parser.images:
            self.assertTrue(image.startswith("data:image/"))
            self.assertTrue(base64.b64decode(image.split(",", 1)[1], validate=True))
        self.assertNotIn("<script>", document)
        self.assertIn("&lt;script&gt;", document)
        self.assertIn("relatedness/a/a%20%26%20b.log", parser.links)
        self.assertIn("Supporting file", document)
        # Reruns must not include the previous report in its own inventory.
        report.build_report(self.root)
        self.assertNotIn('href="plots.html"', output.read_text())
        self.assertEqual(plt.get_fignums(), [])

    def test_report_records_errors_and_continues(self):
        self.write("test.2.Q", "nan 1\n")
        self.write("broken.json", "{")
        self.write("good.normalise.tsv", "metric\tvalue\nretained\t10\n")
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            output, count, failures = report.build_report(self.root)
        self.assertEqual((count, failures), (1, 2))
        self.assertIn("ERROR", output.read_text())
        self.assertIn("test.2.Q", stderr.getvalue())

    def test_missing_empty_and_malformed_inputs_are_explicit_errors(self):
        with self.assertRaisesRegex(ValueError, "does not exist"):
            report.build_report(self.root / "missing")
        with self.assertRaisesRegex(ValueError, "No result files"):
            report.build_report(self.root)
        with self.assertRaisesRegex(ValueError, "no data"):
            common.read_table(self.write("empty.tsv", "metric\tvalue\n"))
        with self.assertRaisesRegex(ValueError, "row width"):
            common.read_table(self.write("ragged.tsv", "a\tb\n1\n"))
        with self.assertRaisesRegex(ValueError, "finite"):
            common.number("inf")

    def test_headings_follow_every_directory_and_keep_runs_separate(self):
        first = "genetic_history/growth/msmc2/growth.msmc2.run"
        second = "genetic_history/bottleneck/msmc2/bottleneck.msmc2.run"
        self.write(f"{first}/inputs/contig_0000.multihetsep.txt", "1\t10\t10\tAC\n")
        self.write(f"{second}/inputs/contig_0000.multihetsep.txt", "1\t20\t20\tAC\n")
        self.write(f"{first}/within1.log", "diagnostic")
        self.write("execution_report.html", "<html>pipeline report</html>")
        self.write("sweep_detection/one & two/flexsweep/run/training/model.keras", "model")
        output, count, failures = report.build_report(self.root)
        parser = ReportParser()
        parser.feed(output.read_text())
        self.assertEqual((count, failures), (5, 0))
        path = f"{first}/inputs/contig_0000.multihetsep.txt"
        ancestors = tuple("/".join(path.split("/")[:i]) for i in range(1, 7))
        self.assertIn((7, ancestors, "contig_0000.multihetsep.txt"), parser.headings)
        self.assertIn((2, ("genetic_history",), "genetic_history"), parser.headings)
        self.assertIn((3, ("genetic_history", "genetic_history/growth"), "growth"), parser.headings)
        self.assertIn((4, ancestors[:3], "msmc2"), parser.headings)
        self.assertIn((5, ancestors[:4], "growth.msmc2.run"), parser.headings)
        self.assertIn((6, ancestors[:5], "inputs"), parser.headings)
        self.assertIn((2, ("execution_report.html",), "execution_report.html"), parser.headings)
        log_path = f"{first}/within1.log"
        self.assertIn((6, ancestors[:4] + (log_path,), "within1.log"), parser.headings)
        run_files = [parents for _, parents, title in parser.headings
                     if title == "contig_0000.multihetsep.txt"]
        self.assertEqual(len(run_files), 2)
        self.assertNotEqual(run_files[0], run_files[1])
        self.assertEqual(len(parser.ids), len(set(parser.ids)))
        self.assertEqual({link[1:] for link in parser.links if link.startswith("#")}, set(parser.ids))
        self.assertNotIn("<h7>", output.read_text())
        self.assertIn("one &amp; two", output.read_text())

    def test_custom_output_preserves_relative_source_links(self):
        self.write("results/relatedness/a/admixture/test.2.Q", ".25 .75\n")
        output, count, failures = report.build_report(self.root / "results", self.root / "custom.html")
        parser = ReportParser()
        parser.feed(output.read_text())
        self.assertEqual((count, failures), (2, 0))
        self.assertIn("results/relatedness/a/admixture/test.2.Q", parser.links)

    def test_dispatch_uses_challenge_modules_and_shared_qc(self):
        cases = [
            ("genetic_history/growth/psmc/a.psmc", genetic_history.plot_psmc),
            ("genetic_history/growth/msmc2/a.final.txt", genetic_history.plot_msmc),
            ("relatedness/relatedness3/admixture/a.3.Q", relatedness.plot_admixture),
            ("relatedness/relatedness3/king_plink/a.kin0", relatedness.plot_relatedness),
            ("sweep_detection/onesweep/sweepfinder2/a.sf2.out", sweep_detection.plot_sweepfinder),
            ("sweep_detection/onesweep/flexsweep/run/predictions.tsv",
             sweep_detection.plot_flexsweep),
            ("genetic_history/growth/normalise_vcf/a.normalise.tsv", common.plot_normalisation),
            ("sweep_detection/onesweep/flexsweep/run/sweeps.bed", common.plot_bed),
        ]
        for path, expected in cases:
            with self.subTest(path=path):
                self.assertIs(report.plotter_for(Path(path)), expected)
        self.assertIsNone(report.plotter_for(Path("model.keras")))


if __name__ == "__main__":
    unittest.main()
