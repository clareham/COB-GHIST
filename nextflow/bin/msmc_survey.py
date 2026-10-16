#!/usr/bin/env python3
"""Run within-population MSMC2 and, for two populations, cross/RCCR analysis."""

import argparse
import csv
import json
import math
import shlex
import shutil
import subprocess
from pathlib import Path

from cyvcf2 import VCF

from vcf_to_multihetsep import convert


def read_rates(path):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty MSMC2 rates: {path}")
    rates = [(float(row["left_time_boundary"]), float(row["right_time_boundary"]),
              float(row["lambda"])) for row in rows]
    for index, (left, right, rate) in enumerate(rates):
        if (not math.isfinite(left) or left < 0 or math.isnan(right) or right <= left or
                not math.isfinite(rate) or rate < 0):
            raise ValueError(f"Invalid MSMC2 rate interval in {path}")
        if index and not math.isclose(left, rates[index - 1][1]):
            raise ValueError(f"Noncontiguous MSMC2 intervals in {path}")
    return rates


def mean_rate(rates, left, right):
    if math.isinf(right):
        return rates[-1][2]
    total = 0.0
    covered = 0.0
    for start, end, rate in rates:
        width = max(0.0, min(end, right) - max(start, left))
        total += width * rate
        covered += width
    if not math.isclose(covered, right - left):
        raise ValueError("Within-population rate intervals do not cover cross-population interval")
    return total / covered


def combine(cross_path, first_path, second_path, output, rccr_output, mutation_rate):
    cross, first, second = (read_rates(path) for path in (cross_path, first_path, second_path))
    with open(output, "w") as handle, open(rccr_output, "w") as rccr:
        handle.write("time_index\tleft_time_boundary\tright_time_boundary\tlambda_00\tlambda_01\tlambda_11\n")
        rccr.write("time_index\tleft_generations\tright_generations\trelative_cross_coalescence\n")
        for index, (left, right, rate) in enumerate(cross):
            within1, within2 = mean_rate(first, left, right), mean_rate(second, left, right)
            denominator = within1 + within2
            if denominator <= 0:
                raise ValueError("Undefined relative cross-coalescence: both within rates are zero")
            handle.write(f"{index}\t{left}\t{right}\t{within1}\t{rate}\t{within2}\n")
            rccr.write(f"{index}\t{left / mutation_rate}\t{right / mutation_rate}\t{2 * rate / denominator}\n")


def run(vcf_path, prefix, options, mask=None, population_map=None):
    threads = options["threads"]
    mutation_rate = options["mutation_rate"]
    if (not isinstance(threads, int) or threads < 1 or
            not math.isfinite(mutation_rate) or mutation_rate <= 0):
        raise ValueError("Positive thread count and mutation rate required")
    binary = next((path for name in ("msmc2", "msmc2_Linux", "msmc2_Darwin")
                   if (path := shutil.which(name))), None)
    if binary is None:
        raise ValueError("No MSMC2 executable found on PATH")
    samples = VCF(vcf_path).samples
    prefixes = options.get("population_prefixes", [])
    if population_map:
        with open(population_map) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            mapping = {row["sample"]: row["pop"] for row in reader}
        if set(samples) - set(mapping):
            raise ValueError("Population map does not cover VCF samples")
    elif prefixes:
        mapping = {}
        for sample in samples:
            matches = [value for value in prefixes if sample.startswith(value)]
            if len(matches) != 1:
                raise ValueError(f"Sample {sample} must match exactly one configured population prefix")
            mapping[sample] = matches[0]
    else:
        mapping = dict.fromkeys(samples, "population")
    populations = list(dict.fromkeys(mapping[sample] for sample in samples))
    if not 1 <= len(populations) <= 2:
        raise ValueError("Survey workflow supports one or two populations")
    count = 2 if len(populations) == 1 else 1
    selected = []
    for population in populations:
        members = [sample for sample in samples if mapping[sample] == population]
        if len(members) < count:
            raise ValueError(f"Need {count} diploid samples for population {population}")
        selected.extend(members[:count])
    work = Path(f"{prefix}.msmc2.run")
    work.mkdir()
    inputs, ordered = convert(vcf_path, work / "inputs", selected, mask,
                              options.get("assume_callable", False), split=True)
    groups = [[hap for i, sample in enumerate(ordered) if mapping[sample] == population
               for hap in (2 * i, 2 * i + 1)] for population in populations]
    extra = shlex.split(options.get("extra_args", ""))
    if any(arg in ("-I", "--pairIndices", "-o", "--outFilePrefix", "-t", "--nrThreads") or
           arg.startswith(("--pairIndices=", "--outFilePrefix=", "--nrThreads=")) for arg in extra):
        raise ValueError("Pair selection, output prefix and threads are managed by msmc_survey")
    runs = {}
    for index, group in enumerate(groups):
        pairs = ",".join(map(str, group))
        runs[f"within{index + 1}"] = (pairs, [])
    if len(groups) == 2:
        pairs = ",".join(f"{a}-{b}" for a in groups[0] for b in groups[1])
        runs["cross"] = (pairs, ["-s"])
    for name, (pairs, flags) in runs.items():
        run_prefix = work / name
        subprocess.run([binary, *extra, "-t", str(threads), "-I", pairs, *flags,
                        "-o", str(run_prefix), *inputs], check=True)
        read_rates(f"{run_prefix}.final.txt")
    final = f"{prefix}.msmc2.final.txt"
    if len(groups) == 2:
        combine(work / "cross.final.txt", work / "within1.final.txt", work / "within2.final.txt",
                final, work / "rccr.tsv", mutation_rate)
    else:
        shutil.copyfile(work / "within1.final.txt", final)
    with open(work / "assumptions.json", "w") as handle:
        json.dump({"options": options, "sample_populations": {sample: mapping[sample] for sample in ordered},
                   "haplotype_groups": groups, "mask": mask,
                   "rccr_combination": "time-weighted within rates on cross-run intervals"}, handle, indent=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf")
    parser.add_argument("prefix")
    parser.add_argument("options_json")
    parser.add_argument("--mask")
    parser.add_argument("--population-map")
    args = parser.parse_args()
    with open(args.options_json) as handle:
        options = json.load(handle)
    run(args.vcf, args.prefix, options, args.mask, args.population_map)


if __name__ == "__main__":
    main()
