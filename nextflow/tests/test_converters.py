"""Small synthetic regression cases; requires cyvcf2, not any inference binary."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

from callable_regions import callable_bases_in_range, load_mask, resolve_mask, subtract_mask
from msmc_survey import combine, read_rates
from vcf_sample_to_psmcfa import convert as psmc
from vcf_to_multihetsep import convert as multihetsep
from vcf_to_sweepfinder_freq import convert as sweepfinder


class ConverterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def write(self, name, content):
        path = self.directory / name
        path.write_text(content)
        return str(path)

    def vcf(self, records, contigs=(("1", 300),), samples=("sample_1", "sample_2")):
        header = "##fileformat=VCFv4.2\n"
        header += "".join(f"##contig=<ID={name},length={length}>\n" for name, length in contigs)
        header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        header += '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">\n'
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        header += "\t".join(samples) + "\n"
        return self.write("input.vcf", header + records)

    def test_psmc_retains_invariant_and_trailing_windows(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n")
        output = self.directory / "output.psmcfa"
        psmc(path, output, assume_callable=True)
        self.assertEqual(output.read_text(), ">1\nTKT\n")

    def test_psmc_missing_genotype_masks_window_and_selects_sample(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\n"
                        "1\t250\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t.|.\n")
        output = self.directory / "output.psmcfa"
        psmc(path, output, sample="sample_2", assume_callable=True)
        self.assertEqual(output.read_text(), ">1\nTKN\n")

    def test_psmc_requires_explicit_callability(self):
        with self.assertRaisesRegex(ValueError, "assume-callable"):
            psmc(self.vcf(""), self.directory / "output")

    def test_psmc_masks_stacked_positions(self):
        # GHIST VCFs stack several mutations at one coordinate as separate
        # records; such windows are masked rather than aborting the run.
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n"
                        "1\t250\t.\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\n"
                        "1\t250\t.\tC\tG\t.\tPASS\t.\tGT\t1|0\t0|0\n")
        output = self.directory / "output.psmcfa"
        psmc(path, output, assume_callable=True)
        self.assertEqual(output.read_text(), ">1\nTKN\n")

    def test_unsorted_positions_are_rejected(self):
        path = self.vcf("1\t250\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n"
                        "1\t150\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|0\n")
        with self.assertRaisesRegex(ValueError, "sorted"):
            psmc(path, self.directory / "output.psmcfa", assume_callable=True)

    def test_multihetsep_preserves_unphased_ambiguity(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0|0\n")
        output = self.directory / "output.txt"
        multihetsep(path, output, assume_callable=True)
        self.assertEqual(output.read_text(), "1\t150\t150\tAGAA,GAAA\n")

    def test_multihetsep_splits_contigs(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n"
                        "2\t50\t.\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|1\n",
                        contigs=(("1", 300), ("2", 100)))
        files, ordered = multihetsep(path, self.directory / "split", assume_callable=True, split=True)
        self.assertEqual(len(files), 2)
        self.assertEqual(ordered, ["sample_1", "sample_2"])
        self.assertEqual(Path(files[1]).read_text(), "2\t50\t50\tCCCT\n")

    def test_multihetsep_excludes_stacked_positions_from_distance(self):
        path = self.vcf("1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\n"
                        "1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t0|0\n"
                        "1\t150\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|0\n")
        output = self.directory / "output.txt"
        multihetsep(path, output, assume_callable=True)
        self.assertEqual(output.read_text(), "1\t150\t149\tCTCC\n")

    def test_mask_absent_contig_is_not_callable(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n")
        mask = self.write("mask.bed", "2\t0\t300\n")
        with self.assertRaisesRegex(ValueError, "No callable segregating"):
            multihetsep(path, self.directory / "output", mask_path=mask)

    def test_subtract_mask_removes_excluded_sites(self):
        mask = {"1": ([0, 100], [50, 200])}
        excluded = {"1": ([10, 49, 150], [20, 60, 151]), "2": ([0], [5])}
        self.assertEqual(subtract_mask(mask, excluded),
                         {"1": ([0, 20, 100, 151], [10, 49, 150, 200])})

    def test_resolve_mask_assumed_callable_minus_exclusions(self):
        excluded = self.write("excluded.bed", "1\t149\t150\n")
        self.assertIsNone(resolve_mask({"1": 300}, assume_callable=True))
        self.assertEqual(resolve_mask({"1": 300}, exclude_path=excluded, assume_callable=True),
                         {"1": ([0, 150], [149, 300])})
        empty = self.write("empty.bed", "")
        self.assertIsNone(resolve_mask({"1": 300}, exclude_path=empty, assume_callable=True))

    def test_psmc_excluded_site_masks_window(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n")
        excluded = self.write("excluded.bed", "1\t249\t250\n")
        output = self.directory / "output.psmcfa"
        psmc(path, output, assume_callable=True, exclude_path=excluded)
        self.assertEqual(output.read_text(), ">1\nTKN\n")

    def test_multihetsep_excluded_site_reduces_distance(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n")
        excluded = self.write("excluded.bed", "1\t99\t100\n")
        output = self.directory / "output.txt"
        multihetsep(path, output, assume_callable=True, exclude_path=excluded)
        self.assertEqual(output.read_text(), "1\t150\t149\tAGAA\n")

    def test_overlapping_mask_intervals_are_not_double_counted(self):
        mask = load_mask(self.write("mask.bed", "1\t0\t20\n1\t10\t30\n1\t40\t50\n"))
        self.assertEqual(callable_bases_in_range(mask["1"], 0, 50), 40)

    def test_sweepfinder_folds_without_ancestral_tag(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|1\n")
        output = self.directory / "output.freq"
        sweepfinder(path, output, folded=True)
        self.assertEqual(output.read_text(), "position\tx\tn\tfolded\n150\t1\t4\t1\n")

    def test_sweepfinder_skips_stacked_and_multiallelic_sites(self):
        path = self.vcf("1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|1\n"
                        "1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|1\n"
                        "1\t120\t.\tA\tG,T\t.\tPASS\t.\tGT\t1|2\t0|1\n"
                        "1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|1\n")
        output = self.directory / "output.freq"
        sweepfinder(path, output, folded=True)
        self.assertEqual(output.read_text(), "position\tx\tn\tfolded\n150\t1\t4\t1\n")

    def test_sweepfinder_polarizes_using_ancestral_allele(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\tAA=G\tGT\t1|1\t0|1\n")
        output = self.directory / "output.freq"
        sweepfinder(path, output)
        self.assertTrue(output.read_text().endswith("150\t1\t4\t0\n"))

    def test_sweepfinder_rejects_unpolarized_unfolded_input(self):
        path = self.vcf("1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n")
        with self.assertRaisesRegex(ValueError, "INFO/AA"):
            sweepfinder(path, self.directory / "output.freq")

    def test_rccr_uses_time_weighted_within_rates(self):
        header = "time_index\tleft_time_boundary\tright_time_boundary\tlambda\n"
        cross = self.write("cross.txt", header + "0\t0\t2\t1\n1\t2\tinf\t1\n")
        within = self.write("within.txt", header + "0\t0\t1\t1\n1\t1\t2\t3\n2\t2\tinf\t2\n")
        output, rccr = self.directory / "combined.txt", self.directory / "rccr.tsv"
        combine(cross, within, within, output, rccr, 1e-8)
        self.assertEqual(output.read_text().splitlines()[1], "0\t0.0\t2.0\t2.0\t1.0\t2.0")
        self.assertTrue(rccr.read_text().splitlines()[1].endswith("\t0.5"))

    def test_msmc_rejects_nan_time_boundary(self):
        path = self.write("rates.txt",
                          "time_index\tleft_time_boundary\tright_time_boundary\tlambda\n"
                          "0\t0\tnan\t1\n")
        with self.assertRaisesRegex(ValueError, "Invalid MSMC2 rate interval"):
            read_rates(path)


if __name__ == "__main__":
    unittest.main()
