#!/usr/bin/env python3
"""Prepare a challenge-matched fastsimcoal2 SFS and demographic model."""

import argparse
import gzip
import io
import json
import math
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


def fail(message):
    raise SystemExit(f"fastsimcoal2 survey: {message}")


def read_options(value):
    try:
        return json.loads(value)
    except json.JSONDecodeError as exc:
        fail(f"invalid options JSON: {exc}")


def open_text(path):
    """Open a plain or (b)gzip-compressed text file for reading.

    Decompression is delegated to an external `gzip` process when available:
    it runs in parallel with this parser and is several times faster than
    Python's gzip module on the multi-GB challenge VCFs.
    """
    if not str(path).endswith(".gz"):
        return open(path, encoding="utf-8")
    if shutil.which("gzip") is None:
        return gzip.open(path, "rt", encoding="utf-8")
    process = subprocess.Popen(["gzip", "-cd", str(path)], stdout=subprocess.PIPE)
    return io.TextIOWrapper(process.stdout, encoding="utf-8")


def parse_vcf_header(path, prefixes):
    contigs = {}
    samples = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##contig=<"):
                attrs = dict(
                    item.split("=", 1)
                    for item in line.rstrip()[10:-1].split(",")
                    if "=" in item
                )
                if "ID" in attrs and "length" in attrs:
                    contigs[attrs["ID"]] = int(attrs["length"])
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    if not samples or not contigs:
        fail("VCF must have samples and ##contig ID/length declarations")
    if prefixes:
        if len(prefixes) != 2:
            fail("population_prefixes must be empty or contain exactly two prefixes")
        groups = [[i for i, sample in enumerate(samples) if sample.startswith(prefix)]
                  for prefix in prefixes]
        unassigned = [
            sample for sample in samples
            if not any(sample.startswith(prefix) for prefix in prefixes)
        ]
        if unassigned or not all(groups):
            fail(
                "population_prefixes do not partition the VCF samples; "
                f"unassigned={unassigned[:5]}, group_sizes={[len(g) for g in groups]}"
            )
        if sum(map(len, groups)) != len(samples):
            fail("a VCF sample matches more than one population prefix")
    else:
        groups = [list(range(len(samples)))]
    return contigs, samples, groups


def iter_vcf_records(path):
    """Yield per-position record groups; stacked records share CHROM/POS."""
    group = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            record = line.rstrip("\n").split("\t")
            if group and record[:2] != group[0][:2]:
                yield group
                group = []
            group.append(record)
    if group:
        yield group


def load_exclusions(paths, contigs):
    """Load and merge zero or more exclusion BEDs into {contig: [(start, end), ...]}."""
    intervals = defaultdict(list)
    if isinstance(paths, (str, Path)):
        paths = [paths]
    for path in paths or []:
        _collect_exclusions(path, contigs, intervals)
    merged = defaultdict(list)
    for chrom, spans in intervals.items():
        for start, end in sorted(spans):
            if merged[chrom] and start <= merged[chrom][-1][1]:
                merged[chrom][-1] = (merged[chrom][-1][0], max(end, merged[chrom][-1][1]))
            else:
                merged[chrom].append((start, end))
    return merged


def _collect_exclusions(path, contigs, intervals):
    bed = Path(path)
    if not bed.is_file():
        fail(f"configured exclusion BED does not exist: {bed}")
    for line_number, line in enumerate(bed.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            fail(f"invalid BED row {line_number} in {bed}")
        chrom, start, end = fields[:3]
        if chrom not in contigs:
            fail(f"BED contig {chrom!r} is absent from the VCF header")
        try:
            start, end = int(start), int(end)
        except ValueError:
            fail(f"invalid BED coordinates at {bed}:{line_number}")
        if start < 0 or end <= start:
            fail(f"invalid BED interval at {bed}:{line_number}")
        start, end = min(start, contigs[chrom]), min(end, contigs[chrom])
        if end > start:
            intervals[chrom].append((start, end))


def in_intervals(position0, spans):
    return any(start <= position0 < end for start, end in spans)


def excluded_length(contigs, exclusions):
    return sum(end - start for spans in exclusions.values() for start, end in spans)


def get_alleles(record, samples):
    if len(record) < 10:
        fail("VCF data row has no sample genotypes")
    chrom, pos, _identifier, ref, alt, _qual, filt, _info, fmt, *calls = record
    if chrom not in samples["contigs"]:
        fail(f"record contig {chrom!r} has no ##contig length")
    try:
        position = int(pos)
    except ValueError:
        fail(f"invalid VCF position {pos!r} on {chrom}")
    if position < 1 or position > samples["contigs"][chrom]:
        fail(f"VCF position outside declared contig bounds: {chrom}:{position}")
    # Multiallelic sites cannot enter a biallelic SFS; the caller drops them
    # from the callable length like any other rejected record.
    if "," in alt or len(ref) != 1 or len(alt) != 1 or alt == ".":
        return chrom, position, None, filt
    fmt_fields = fmt.split(":")
    if "GT" not in fmt_fields:
        fail(f"site has no GT field: {chrom}:{position}")
    gt_index = fmt_fields.index("GT")
    alleles = []
    for sample_call in calls:
        values = sample_call.split(":")
        gt = values[gt_index] if gt_index < len(values) else "."
        genotype = re.split(r"[/|]", gt)
        if len(genotype) != 2 or any(a not in ("0", "1") for a in genotype):
            return chrom, position, None, filt
        alleles.extend(int(a) for a in genotype)
    return chrom, position, alleles, filt


def render_fsc_model(challenge, hap_counts, options, mutation_rate):
    pop_count = len(hap_counts)
    if pop_count not in (1, 2):
        fail("only one- and two-population survey models are supported")
    size_bounds = options.get("population_size_bounds", [100, 1_000_000])
    time_bounds = options.get("time_bounds", [10, 100_000])
    migration_bounds = options.get("migration_rate_bounds", [1e-8, 0.1])
    if len(size_bounds) != 2 or len(time_bounds) != 2:
        fail("population_size_bounds and time_bounds must be two-element arrays")
    if size_bounds[0] <= 0 or size_bounds[1] <= size_bounds[0]:
        fail("invalid population_size_bounds")
    if time_bounds[0] <= 0 or time_bounds[1] <= time_bounds[0]:
        fail("invalid time_bounds")

    common = [
        "// Priors are survey search ranges; population sizes are diploid individuals.",
        "[PARAMETERS]",
    ]
    # fsc reads the template as alternating "//comment" + data blocks, so
    # keep exactly one comment line ahead of every block (see the fsc27
    # manual's .tpl examples).
    tpl = ["//Number of population samples (demes)", str(pop_count)]
    if pop_count == 1:
        tpl += ["//Population effective sizes (number of genes)", "HAPSIZE",
                "//Sample sizes", str(hap_counts[0]),
                "//Growth rates: negative growth implies population expansion", "0",
                "//Number of migration matrices", "0",
                "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix"]
        if challenge == "bottleneck":
            fraction_bounds = options.get("fraction_bounds", [0.01, 1.0])
            tpl += ["1 historical event", "TCHANGE 0 0 0 RESIZE 0 0"]
            common += [
                f"1 NDIP unif {size_bounds[0]} {size_bounds[1]} output bounded",
                f"1 FRAC unif {fraction_bounds[0]} {fraction_bounds[1]} output bounded",
                f"1 TCHANGE unif {time_bounds[0]} {time_bounds[1]} output bounded",
                "[COMPLEX PARAMETERS]",
                "0 HAPSIZE = 2*NDIP hide",
                "0 RESIZE = 1/FRAC hide",
            ]
        elif challenge == "growth":
            scale_bounds = options.get("growth_scale_bounds", [1.01, 100.0])
            tpl += ["1 historical event", "TCHANGE 0 0 0 RESIZE 0 0"]
            common += [
                f"1 NDIP unif {size_bounds[0]} {size_bounds[1]} output bounded",
                f"1 FRAC unif {scale_bounds[0]} {scale_bounds[1]} output bounded",
                f"1 TCHANGE unif {time_bounds[0]} {time_bounds[1]} output bounded",
                "[COMPLEX PARAMETERS]",
                "0 HAPSIZE = 2*NDIP hide",
                "0 RESIZE = 1/FRAC hide",
            ]
        else:
            fail(f"no single-population model is defined for {challenge!r}")
    else:
        if len(migration_bounds) != 2 or migration_bounds[0] <= 0 or migration_bounds[1] <= migration_bounds[0]:
            fail("migration_rate_bounds must be two increasing positive values")
        tpl += [
            "//Population effective sizes (number of genes)", "HAPA", "HAPB",
            "//Sample sizes", str(hap_counts[0]), str(hap_counts[1]),
            "//Growth rates: negative growth implies population expansion", "0", "0",
        ]
        if challenge == "secondary_contact":
            ratio_bounds = options.get("contact_age_ratio_bounds", [0.01, 0.99])
            if (len(ratio_bounds) != 2 or ratio_bounds[0] <= 0
                    or ratio_bounds[1] >= 1 or ratio_bounds[1] <= ratio_bounds[0]):
                fail("contact_age_ratio_bounds must be an increasing interval inside (0, 1)")
            tpl += [
                "//Number of migration matrices", "2",
                "//Migration after contact, symmetric in both directions",
                "0 MIG", "MIG 0",
                "//Migration matrix during isolation",
                "0 0", "0 0",
                "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
                "2 historical event",
                "TCONTACT 0 0 0 1 0 1",
                "TSPLIT 1 0 1 1 0 1",
            ]
            common += [
                f"1 NPOPA unif {size_bounds[0]} {size_bounds[1]} output bounded",
                f"1 NPOPB unif {size_bounds[0]} {size_bounds[1]} output bounded",
                f"1 TSPLIT unif {time_bounds[0]} {time_bounds[1]} output bounded",
                f"1 RATIO unif {ratio_bounds[0]} {ratio_bounds[1]} output bounded",
                f"1 MIG unif {migration_bounds[0]} {migration_bounds[1]} output bounded",
                "[COMPLEX PARAMETERS]",
                "0 HAPA = 2*NPOPA hide",
                "0 HAPB = 2*NPOPB hide",
                "0 TCONTACT = TSPLIT*RATIO output",
            ]
        elif challenge == "split_migration":
            tpl += [
                "//Number of migration matrices", "2",
                "//Symmetric migration until the split",
                "0 MIG", "MIG 0",
                "//No migration in the ancestral population",
                "0 0", "0 0",
                "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
                "1 historical event",
                "TSPLIT 1 0 1 1 0 1",
            ]
            common += [
                f"1 NPOPA unif {size_bounds[0]} {size_bounds[1]} output bounded",
                f"1 NPOPB unif {size_bounds[0]} {size_bounds[1]} output bounded",
                f"1 TSPLIT unif {time_bounds[0]} {time_bounds[1]} output bounded",
                f"1 MIG unif {migration_bounds[0]} {migration_bounds[1]} output bounded",
                "[COMPLEX PARAMETERS]",
                "0 HAPA = 2*NPOPA hide",
                "0 HAPB = 2*NPOPB hide",
            ]
        else:
            fail(f"no two-population model is defined for {challenge!r}")
    tpl += [
        "//Number of independent loci [chromosome]", "1 0",
        "//Per chromosome: Number of linkage blocks", "1",
        "//per Block: data type, num loci, rec. rate and mut rate",
        f"FREQ 1 0 {mutation_rate:.12g}",
    ]
    return "\n".join(tpl) + "\n", "\n".join(common) + "\n"


def prepare(args):
    options = read_options(args.options)
    mutation_rate = float(options.get("mutation_rate", 1.25e-8))
    if not math.isfinite(mutation_rate) or mutation_rate <= 0:
        fail("mutation_rate must be a positive finite per-generation rate")
    prefixes = options.get("population_prefixes", [])
    if not isinstance(prefixes, list):
        fail("population_prefixes must be a list")
    contigs, samples, groups = parse_vcf_header(args.vcf, prefixes)
    if len(groups) == 2 and options.get("assume_ref_ancestral") is not True:
        fail("joint DAF requires the configured REF-as-ancestral assumption")
    if len(groups) == 2 and options.get("migration_mode", "symmetric") != "symmetric":
        fail("only symmetric two-population migration is implemented")
    if len(groups) == 2 and options.get("ancestral_size_rule", "population_0") != "population_0":
        fail("the fsc split template currently retains population 0 size as ancestral size")
    exclusion_beds = [path for path in (options.get("exclude_bed"), args.exclude_bed) if path]
    exclusions = load_exclusions(exclusion_beds, contigs)
    if options.get("assume_callable") is not True and not options.get("exclude_bed"):
        fail("set assume_callable=true or configure an exclusion BED")
    callable_length = sum(contigs.values()) - excluded_length(contigs, exclusions)
    if callable_length <= 0:
        fail("no callable sequence remains after applying the configured mask")
    hap_counts = [2 * len(group) for group in groups]
    spectrum = (
        [[0] * (hap_counts[0] // 2 + 1)] if len(groups) == 1
        else [[0] * (hap_counts[0] + 1) for _ in range(hap_counts[1] + 1)]
    )
    callable_variant_records = 0
    rejected_sites = 0
    previous = {}
    metadata = {"contigs": contigs}
    for records in iter_vcf_records(args.vcf):
        chrom, position, alleles, filt = get_alleles(records[0], metadata)
        if position <= previous.get(chrom, 0):
            fail(f"VCF must be sorted within each contig: {chrom}:{position}")
        previous[chrom] = position
        if len(records) > 1:
            # Stacked mutations at one coordinate are not a biallelic site.
            alleles = None
        if in_intervals(position - 1, exclusions.get(chrom, [])):
            continue
        if filt not in ("PASS", ".") or alleles is None:
            callable_length -= 1
            rejected_sites += 1
            continue
        callable_variant_records += 1
        alt_counts = []
        for group in groups:
            alt_counts.append(sum(alleles[2 * i] + alleles[2 * i + 1] for i in group))
        if len(groups) == 1:
            folded = min(alt_counts[0], hap_counts[0] - alt_counts[0])
            spectrum[0][folded] += 1
        else:
            spectrum[alt_counts[1]][alt_counts[0]] += 1
    if callable_length < callable_variant_records:
        fail("callable sequence length is shorter than the retained variant count")
    if len(groups) == 1:
        spectrum[0][0] += callable_length - callable_variant_records
    else:
        spectrum[0][0] += callable_length - callable_variant_records
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / args.challenge
    if len(groups) == 1:
        with open(f"{prefix}_MAFpop0.obs", "w", encoding="utf-8") as handle:
            handle.write("1 observations\n")
            handle.write(" ".join(f"d0_{i}" for i in range(len(spectrum[0]))) + "\n")
            handle.write(" ".join(map(str, spectrum[0])) + "\n")
    else:
        with open(f"{prefix}_jointDAFpop1_0.obs", "w", encoding="utf-8") as handle:
            handle.write("1 observations\n")
            handle.write(" ".join(f"d0_{i}" for i in range(hap_counts[0] + 1)) + "\n")
            for i, row in enumerate(spectrum):
                handle.write(f"d1_{i} " + " ".join(map(str, row)) + "\n")
    tpl, est = render_fsc_model(args.challenge, hap_counts, options, mutation_rate)
    (outdir / f"{args.challenge}.tpl").write_text(tpl, encoding="utf-8")
    (outdir / f"{args.challenge}.est").write_text(est, encoding="utf-8")
    (outdir / "survey_manifest.json").write_text(json.dumps({
        "tool": "fastsimcoal2",
        "challenge": args.challenge,
        "options": options,
        "contig_lengths_bp": contigs,
        "samples": samples,
        "population_prefixes": prefixes,
        "population_sample_counts": [len(group) for group in groups],
        "population_haplotype_counts": hap_counts,
        "sfs": "folded one-dimensional" if len(groups) == 1 else "unfolded joint two-dimensional (REF assumed ancestral)",
        "callable_length_bp": callable_length,
        "retained_biallelic_snps": callable_variant_records,
        "rejected_sites": rejected_sites,
        "excluded_bed": options.get("exclude_bed"),
        "excluded_sites_bed": args.exclude_bed,
        "mutation_rate_per_generation_per_bp": mutation_rate,
        "assumptions": {
            "ploidy": 2,
            "callability": "VCF ##contig lengths minus excluded BED intervals and rejected VCF records",
            "missing_or_filtered_sites": "excluded from callable sequence",
            "multiallelic_or_stacked_sites": "excluded from the SFS and from callable sequence",
            "population_size_units": "diploid individuals; model converts sizes to haploid gene copies",
            "growth_and_bottleneck": "instantaneous size change at T generations before present",
            "two_population_migration": "symmetric backward migration; split event merges population 1 into population 0",
            "joint_sfs_polarization": "VCF REF allele assumed ancestral",
            "ancestral_size": "population 0 size is retained at the split",
        },
    }, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--challenge", required=True,
                        choices=("bottleneck", "growth", "secondary_contact", "split_migration"))
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--options", default="{}")
    parser.add_argument("--exclude-bed", default=None,
                        help="additional BED of positions removed from the VCF (normalise_vcf.sh output)")
    args = parser.parse_args()
    prepare(args)


if __name__ == "__main__":
    main()
