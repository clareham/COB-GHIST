#!/usr/bin/env python3
"""Prepare phased-per-contig VCFs and a challenge-specific mrpast model."""

import argparse
import json
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def fail(message):
    raise SystemExit(f"mrpast survey: {message}")


def read_options(value):
    try:
        return json.loads(value)
    except json.JSONDecodeError as exc:
        fail(f"invalid options JSON: {exc}")


def load_exclusions(path, contigs):
    intervals = defaultdict(list)
    if not path:
        return intervals
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
    merged = defaultdict(list)
    for chrom, spans in intervals.items():
        for start, end in sorted(spans):
            if merged[chrom] and start <= merged[chrom][-1][1]:
                merged[chrom][-1] = (merged[chrom][-1][0], max(end, merged[chrom][-1][1]))
            else:
                merged[chrom].append((start, end))
    return merged


def inside_mask(position0, spans):
    return any(start <= position0 < end for start, end in spans)


def parse_header(vcf_path, prefixes):
    contigs = {}
    header = []
    samples = []
    with open(vcf_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                header.append(line)
                break
            header.append(line)
            if line.startswith("##contig=<"):
                attrs = dict(
                    item.split("=", 1)
                    for item in line.rstrip()[10:-1].split(",")
                    if "=" in item
                )
                if "ID" in attrs and "length" in attrs:
                    contigs[attrs["ID"]] = int(attrs["length"])
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
    return contigs, header, samples, groups


def emit_contig_vcfs(vcf_path, outdir, contigs, header, samples, exclusions):
    vcf_dir = outdir / "vcf"
    vcf_dir.mkdir(parents=True, exist_ok=True)
    contig_names = sorted(contigs)
    handles = []
    outputs = []
    record_counts = defaultdict(int)
    try:
        # Emit a minimal, well-formed header: only GT is carried over, so
        # the source INFO/FORMAT definitions (one of which is unterminated
        # in the GHIST drops and makes htslib swallow following lines) are
        # not reproduced.
        fileformat = next((line for line in header if line.startswith("##fileformat=")),
                          "##fileformat=VCFv4.2\n")
        column_line = next(line for line in header if line.startswith("#CHROM"))
        for index, chrom in enumerate(contig_names):
            output = vcf_dir / f"contig_{index:03d}.vcf"
            outputs.append(output)
            handle = open(output, "w", encoding="utf-8")
            handles.append(handle)
            handle.write(fileformat)
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write(f"##contig=<ID={chrom},length={contigs[chrom]}>\n")
            handle.write(column_line)
        index_by_contig = {chrom: index for index, chrom in enumerate(contig_names)}
        previous = {}
        dropped = defaultdict(int)

        def usable(fields):
            """Return GT-only fields for a complete, biallelic, PASS SNP record, else None."""
            ref, alt = fields[3], fields[4]
            if len(ref) != 1 or len(alt) != 1 or alt == "." or "," in alt:
                return None
            if fields[6] not in ("PASS", "."):
                return None
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                fail(f"record has no GT field: {fields[0]}:{fields[1]}")
            gt_idx = format_fields.index("GT")
            genotypes = []
            for sample_call in fields[9:]:
                values = sample_call.split(":")
                gt = values[gt_idx] if gt_idx < len(values) else "."
                alleles = gt.replace("|", "/").split("/")
                if len(alleles) != 2 or any(a not in ("0", "1") for a in alleles):
                    return None
                genotypes.append(gt)
            return fields[:7] + [".", "GT"] + genotypes

        def emit(group):
            # Stacked records at one coordinate, multiallelic sites, filtered
            # records and incomplete genotypes are dropped: tsinfer needs
            # one clean biallelic site per position.
            fields = group[0]
            chrom = fields[0]
            cleaned = None if len(group) > 1 else usable(fields)
            if cleaned is None:
                dropped[chrom] += 1
                return
            handles[index_by_contig[chrom]].write("\t".join(cleaned) + "\n")
            record_counts[chrom] += 1

        group = []
        with open(vcf_path, encoding="utf-8") as source:
            for line in source:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 + len(samples):
                    fail(f"VCF row has the wrong number of samples at {fields[0]}:{fields[1]}")
                chrom, pos = fields[0], int(fields[1])
                if chrom not in contigs or not 1 <= pos <= contigs[chrom]:
                    fail(f"record outside declared contig bounds: {chrom}:{pos}")
                if group and (chrom, pos) != (group[0][0], int(group[0][1])):
                    emit(group)
                    group = []
                if not group:
                    if pos <= previous.get(chrom, 0):
                        fail(f"VCF must be sorted within each contig: {chrom}:{pos}")
                    previous[chrom] = pos
                if inside_mask(pos - 1, exclusions.get(chrom, [])):
                    continue
                group.append(fields)
            if group:
                emit(group)
    finally:
        for handle in handles:
            handle.close()
    if dropped:
        print(f"[mrpast input] dropped stacked/multiallelic/incomplete sites per contig: "
              f"{dict(dropped)}", file=sys.stderr)
    retained = [
        (path, chrom, index)
        for index, (chrom, path) in enumerate(zip(contig_names, outputs))
        if record_counts[chrom] > 0
    ]
    for path, chrom, _index in (
        (path, chrom, index)
        for index, (chrom, path) in enumerate(zip(contig_names, outputs))
        if record_counts[chrom] == 0
    ):
        path.unlink()
    if not retained:
        fail("no per-contig VCF records remain after applying the mask")
    return retained


def rate_bounds(options):
    size_bounds = options.get("population_size_bounds", [100, 1_000_000])
    if len(size_bounds) != 2 or size_bounds[0] <= 0 or size_bounds[1] <= size_bounds[0]:
        fail("population_size_bounds must contain two increasing positive sizes")
    return 1.0 / (2.0 * size_bounds[1]), 1.0 / (2.0 * size_bounds[0])


def parameter(ground_truth, bounds, index):
    return {
        "ground_truth": float(ground_truth),
        "lb": float(bounds[0]),
        "ub": float(bounds[1]),
        "index": index,
    }


def make_model(challenge, names, options):
    rate_min, rate_max = rate_bounds(options)
    size_mid = math.sqrt(
        options.get("population_size_bounds", [100, 1_000_000])[0]
        * options.get("population_size_bounds", [100, 1_000_000])[1]
    )
    base_rate = 1.0 / (2.0 * size_mid)
    model = {"ploidy": 2, "pop_count": len(names), "pop_names": names}
    coal_entries = []
    coal_parameters = []
    migration_entries = []
    migration_parameters = []
    time_entries = []

    if len(names) == 1:
        coal_entries = [
            {"epoch": 0, "deme": names[0], "rate": {"param": 1}},
            {"epoch": 1, "deme": names[0], "rate": {"param": 2}},
        ]
        coal_parameters = [parameter(base_rate, [rate_min, rate_max], i) for i in (1, 2)]
        bounds = options.get("time_bounds", [10, 100_000])
        time_entries = [parameter(5000, bounds, 1)]
    else:
        current = [
            {"epoch": 0, "deme": names[0], "rate": {"param": 1}},
            {"epoch": 0, "deme": names[1], "rate": {"param": 2}},
        ]
        ancestral = [
            {"deme": names[0], "rate": {"param": 3}},
            {"deme": names[1], "rate": {"param": 3}},
        ]
        migration_bounds = options.get("migration_rate_bounds", [1e-8, 0.1])
        if len(migration_bounds) != 2 or migration_bounds[0] <= 0 or migration_bounds[1] <= migration_bounds[0]:
            fail("migration_rate_bounds must contain two increasing positive rates")
        initial_migration = float(options.get("initial_migration_rate", 0.001))
        if not migration_bounds[0] <= initial_migration <= migration_bounds[1]:
            fail("initial_migration_rate must lie within migration_rate_bounds")
        ancestral_migration = float(options.get("ancestral_migration_rate", 0.1))
        if ancestral_migration <= 0:
            fail("ancestral_migration_rate must be positive")
        coal_parameters = [
            parameter(base_rate, [rate_min, rate_max], 1),
            parameter(base_rate, [rate_min, rate_max], 2),
            parameter(base_rate, [rate_min, rate_max], 3),
        ]
        for source, dest in ((names[0], names[1]), (names[1], names[0])):
            migration_entries.append({
                "epoch": 0, "source": source, "dest": dest, "rate": {"param": 1}
            })
            migration_entries.append({
                "epoch": 1 if challenge == "split_migration" else 2,
                "source": source, "dest": dest, "rate": ancestral_migration
            })
        migration_parameters = [parameter(initial_migration, migration_bounds, 1)]
        if challenge == "secondary_contact":
            current += [
                {"epoch": 1, "deme": names[0], "rate": {"param": 3}},
                {"epoch": 1, "deme": names[1], "rate": {"param": 4}},
            ]
            ancestral = [
                {"deme": names[0], "rate": {"param": 5}},
                {"deme": names[1], "rate": {"param": 5}},
            ]
            coal_entries = current
            coal_parameters += [
                parameter(base_rate, [rate_min, rate_max], 4),
                parameter(base_rate, [rate_min, rate_max], 5),
            ]
            migration_entries = [
                entry for entry in migration_entries if entry["epoch"] != 2
            ]
            migration_entries += [
                {
                    "epoch": 2, "source": names[0], "dest": names[1],
                    "rate": ancestral_migration,
                },
                {
                    "epoch": 2, "source": names[1], "dest": names[0],
                    "rate": ancestral_migration,
                },
            ]
            recent_bounds = options.get("recent_time_bounds", [10, 9000])
            old_bounds = options.get("old_time_bounds", [10000, 100000])
            if (len(recent_bounds) != 2 or len(old_bounds) != 2
                    or recent_bounds[0] <= 0 or recent_bounds[1] <= recent_bounds[0]
                    or old_bounds[0] <= recent_bounds[1] or old_bounds[1] <= old_bounds[0]):
                fail("recent_time_bounds must be strictly more recent than old_time_bounds")
            time_entries = [parameter(1000, recent_bounds, 1), parameter(20000, old_bounds, 2)]
        else:
            time_entries = [parameter(10000, options.get("time_bounds", [10, 100000]), 1)]
        coal_entries = current + [
            {"epoch": 2 if challenge == "secondary_contact" else 1,
             "deme": entry["deme"], "rate": entry["rate"]}
            for entry in ancestral
        ]
    model["coalescence"] = {"entries": coal_entries, "parameters": coal_parameters}
    if migration_entries:
        model["migration"] = {"entries": migration_entries, "parameters": migration_parameters}
    model["epochTimeSplit"] = time_entries
    return model


def write_rate_maps(outdir, contigs, exclusions, inference_contigs, options):
    prefix_dir = outdir / "rate_maps"
    prefix_dir.mkdir(parents=True, exist_ok=True)
    neutral_rate = float(options.get("tree_sampling_rate", 1e-10))
    masked_rate = float(options.get("masked_tree_sampling_rate", 1e-8))
    if neutral_rate <= 0 or masked_rate <= 1e-9:
        fail("tree_sampling_rate must be positive and masked_tree_sampling_rate must exceed 1e-9")
    maps = []
    for _vcf_path, chrom, index in inference_contigs:
        spans = exclusions.get(chrom, [])
        boundaries = {0, contigs[chrom]}
        for start, end in spans:
            boundaries.add(start)
            boundaries.add(end)
        points = sorted(boundaries)
        path = prefix_dir / f"contig_{index:03d}.txt"
        with open(path, "w", encoding="utf-8") as handle:
            first = True
            for start, end in zip(points, points[1:]):
                masked = inside_mask(start, spans)
                rate = masked_rate if masked else neutral_rate
                # mrpast's parser seeds position 0 from the first row; use 1
                # as the first explicit breakpoint to avoid a duplicate zero.
                if first and start == 0:
                    start = 1
                if end <= start:
                    continue
                handle.write(f"{start} {end} {rate:.12g}\n")
                first = False
        maps.append(path)
    return str(prefix_dir / "contig_")


def run(args):
    options = read_options(args.options)
    mut_rate = float(options.get("mutation_rate", 1.25e-8))
    recomb_rate = float(options.get("recombination_rate", 1e-8))
    if not math.isfinite(mut_rate) or mut_rate <= 0:
        fail("mutation_rate must be positive and finite")
    if not math.isfinite(recomb_rate) or recomb_rate <= 0:
        fail("recombination_rate must be positive and finite")
    prefixes = options.get("population_prefixes", [])
    contigs, header, samples, groups = parse_header(args.vcf, prefixes)
    exclusions = load_exclusions(options.get("exclude_bed"), contigs)
    if options.get("assume_ref_ancestral") is not True:
        fail("tsinfer workflow requires the documented REF-as-ancestral assumption")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    inference_contigs = emit_contig_vcfs(args.vcf, outdir, contigs, header, samples, exclusions)
    for vcf_path, _chrom, _index in inference_contigs:
        vcz_path = vcf_path.with_suffix(".vcz")
        subprocess.run(
            ["vcf2zarr", "convert", str(vcf_path), str(vcz_path)],
            check=True,
        )
    names = [prefix.rstrip("_") for prefix in prefixes] if prefixes else ["population_0"]
    model = make_model(args.challenge, names, options)
    model_path = outdir / f"{args.challenge}.yaml"
    model_path.write_text(json.dumps(model, indent=2) + "\n", encoding="utf-8")
    popmap_path = outdir / "population_map.json"
    popmap = {
        "mapping": [group for group in groups],
        "names": names,
    }
    popmap_path.write_text(json.dumps(popmap, indent=2) + "\n", encoding="utf-8")
    arg_prefix = str(outdir / "args" / "arg_")
    Path(arg_prefix).parent.mkdir(parents=True, exist_ok=True)
    jobs = int(options.get("jobs", 1))
    arginfer = [
        "mrpast", "arginfer", "-j", str(jobs), "--mut-rate", str(mut_rate),
        "--recomb-rate", str(recomb_rate), "--tool", "tsinfer",
        str(outdir / "vcf" / "contig_"), arg_prefix, str(popmap_path),
    ]
    subprocess.run(arginfer, check=True)
    arg_files = sorted((outdir / "args").glob("*.trees"))
    if not arg_files or any(path.stat().st_size == 0 for path in arg_files):
        fail("mrpast arginfer did not produce nonempty .trees files")
    # mrpast process writes into --out-dir but does not create it (the
    # documented workflow runs `mkdir -p` first).
    (outdir / "solve").mkdir(parents=True, exist_ok=True)
    process = [
        "mrpast", "process", "-j", str(jobs),
        "--num-times", str(int(options.get("num_times", 50))),
        "--replicates", str(int(options.get("replicates", 1))),
        "--tree-sample-rate", str(int(options.get("tree_sample_rate", 10000))),
        "--bootstrap", "none", "--solve", "--out-dir", str(outdir / "solve"),
    ]
    if options.get("exclude_bed"):
        rate_prefix = write_rate_maps(outdir, contigs, exclusions, inference_contigs, options)
        process += ["--rate-maps", rate_prefix, "--rate-map-threshold", "1e-9"]
    process += [str(model_path), arg_prefix]
    subprocess.run(process, check=True)
    results = sorted((outdir / "solve").glob("*.out.json"))
    if not results:
        fail("mrpast process --solve produced no *.out.json solver result")
    scored = []
    for path in results:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        likelihood = data.get("negLL")
        if isinstance(likelihood, (int, float)) and math.isfinite(likelihood):
            scored.append((likelihood, path))
    if not scored:
        fail("mrpast solver results contain no finite negLL values")
    best = min(scored, key=lambda item: item[0])[1]
    with open(outdir / "best-fit.txt", "w", encoding="utf-8") as handle:
        subprocess.run(["mrpast", "show", str(best)], check=True, stdout=handle)
    (outdir / "survey_manifest.json").write_text(json.dumps({
        "tool": "mrpast",
        "version": "0.5",
        "challenge": args.challenge,
        "options": options,
        "vcf": str(args.vcf),
        "contigs": contigs,
        "samples": samples,
        "population_prefixes": prefixes,
        "population_sample_counts": [len(group) for group in groups],
        "population_map": str(popmap_path),
        "model": model,
        "arg_files": [str(path) for path in arg_files],
        "best_result": str(best),
        "best_negLL": min(scored, key=lambda item: item[0])[0],
        "mutation_rate_per_generation_per_bp": mut_rate,
        "recombination_rate_per_generation_per_bp": recomb_rate,
        "excluded_bed": options.get("exclude_bed"),
        "assumptions": {
            "ploidy": 2,
            "ancestral_state": "VCF REF allele, as used by mrpast's tsinfer integration",
            "per_contig_inference": "multi-contig VCF is split using declared contig lengths before inference",
            "size_change": "instantaneous piecewise-constant coalescence-rate change",
            "ancestral_two_population_model": "shared ancestral coalescence rate plus fixed symmetric migration approximation",
            "mask_sampling": "configured BED excluded from VCF variants and marginal-tree sampling",
        },
    }, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--challenge", required=True,
                        choices=("bottleneck", "growth", "secondary_contact", "split_migration"))
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--options", default="{}")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
