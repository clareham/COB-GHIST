"""Neutral, reference-anchored controls and replicate-based B estimates."""

import csv
import json
import math
from pathlib import Path
import random

from .inputs import Reference, cumulative_map, interpolate


def simulate(output: Path, reference: Reference, start: int, end: int,
             replicates: int, samples: int, ne: float, mutation_rate: float,
             rate: float, seed: int, map_path: Path | None = None) -> Path:
    import msprime

    if (replicates < 2 or samples < 2 or not 0 <= start < end - 1 < reference.length
            or not math.isfinite(ne) or ne <= 0 or not math.isfinite(mutation_rate)
            or mutation_rate <= 0 or not math.isfinite(rate) or rate < 0
            or seed < 1 or seed + 2 * replicates >= 2**32):
        raise ValueError("Invalid simulation parameters; need >=2 replicates and samples")
    sequence = reference.fetch(start, end)
    if any(base not in "ACGT" for base in sequence):
        raise ValueError("Choose an unambiguous ACGT reference interval for neutral simulation")
    output.mkdir(parents=True, exist_ok=False)
    mask = output / "accessible.bed"
    mask.write_text(f"{reference.chrom}\t{start}\t{end}\n")
    map_file = output / "recombination.tsv"
    if map_path:
        positions, cm = cumulative_map(map_path)
        boundaries = [start + 1, *(p for p in positions if start + 1 < p < end), end]
        genetic = [interpolate(p, positions, cm) for p in boundaries]
        breakpoints = [0, *(p - start - 1 for p in boundaries[1:-1]), end - start]
        # Rightmost boundary is between the final base and its successor.
        rates = [(genetic[i + 1] - genetic[i]) / 100 /
                 (boundaries[i + 1] - boundaries[i]) for i in range(len(boundaries) - 1)]
        recombination = msprime.RateMap(position=breakpoints, rate=rates)
    else:
        boundaries = [start + 1, end]
        genetic = [0.0, (end - start - 1) * rate * 100]
        recombination = rate
    with map_file.open("w") as handle:
        handle.write("position\tcM\n")
        for position, cm_value in zip(boundaries, genetic):
            handle.write(f"{position}\t{cm_value:.12g}\n")
    cases = []
    for replicate in range(replicates):
        ancestry = msprime.sim_ancestry(
            samples=samples, ploidy=2, population_size=ne, sequence_length=end - start,
            recombination_rate=recombination, random_seed=seed + 2 * replicate)
        trees = msprime.sim_mutations(
            ancestry, rate=mutation_rate, model=msprime.JC69(), discrete_genome=True,
            random_seed=seed + 2 * replicate + 1)
        name = f"neutral-{replicate:04d}"
        path = output / f"{name}.vcf"
        with path.open("w") as handle:
            handle.write('##fileformat=VCFv4.2\n'
                         f'##contig=<ID={reference.chrom},length={reference.length}>\n'
                         '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">\n'
                         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                         '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
            handle.write("\t".join(f"sample_{i}" for i in range(samples)) + "\n")
            for variant in trees.variants():
                local = int(variant.site.position)
                ancestor = variant.site.ancestral_state
                ref = sequence[local]
                # JC69 is symmetric: a per-site nucleotide permutation preserves the
                # simulated process while conditioning its root on the reference.
                shift = "ACGT".index(ref) - "ACGT".index(ancestor)
                labels = ["ACGT"[("ACGT".index(a) + shift) % 4]
                          for a in variant.alleles if a is not None]
                alternatives = [a for a in labels if a != ref]
                if not alternatives:
                    continue
                ordered = [ref, *alternatives]
                remap = [ordered.index(a) for a in labels]
                genotypes = [remap[int(g)] for g in variant.genotypes]
                gt = [f"{genotypes[i]}|{genotypes[i + 1]}"
                      for i in range(0, len(genotypes), 2)]
                handle.write(f"{reference.chrom}\t{start + local + 1}\t.\t{ref}\t"
                             f"{','.join(alternatives)}\t.\tPASS\tAA={ref}\tGT\t"
                             + "\t".join(gt) + "\n")
        cases.append({"id": name, "vcf": path.name, "mask": mask.name,
                      "recombination_map": map_file.name, "dataset": "testing"
                      if reference.chrom == "7" else "final", "scenario": "neutral",
                      "split": "calibration" if replicate < replicates // 2 else "evaluation",
                      "sweep_position": None})
    manifest = output / "manifest.json"
    manifest.write_text(json.dumps({
        "cases": cases,
        "simulation": {"engine": "msprime", "reference": reference.accession,
                       "start0": start, "end0": end, "samples_diploid": samples, "ne": ne,
                       "mutation_rate": mutation_rate, "constant_rate": rate,
                       "source_map": str(map_path) if map_path else None, "seed": seed,
                       "msprime_version": msprime.__version__,
                       "sequence_model": "reference-conditioned JC69; not a context mutation model"}
    }, indent=2) + "\n")
    return manifest


def estimate_b(diversity: Path, output: Path, bootstrap: int = 1000, seed: int = 1) -> None:
    """Matched no-sweep replicates: scenario, replicate, chrom, start, end, pi."""
    if bootstrap < 1:
        raise ValueError("bootstrap must be positive")
    groups: dict[tuple[str, int, int], dict[str, dict[str, float]]] = {}
    with diversity.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            scenario = row["scenario"]
            if scenario not in ("neutral", "bgs"):
                raise ValueError("B estimates accept only neutral and bgs no-sweep controls")
            key = row["chrom"], int(row["start"]), int(row["end"])
            value = float(row["pi"])
            if key[1] < 0 or key[2] <= key[1] or not math.isfinite(value) or value < 0:
                raise ValueError("Invalid diversity interval or value")
            bucket = groups.setdefault(key, {}).setdefault(scenario, {})
            if row["replicate"] in bucket:
                raise ValueError("Duplicate diversity replicate")
            bucket[row["replicate"]] = value
    if not groups:
        raise ValueError("No diversity data")
    rng = random.Random(seed)
    with output.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chrom", "start", "end", "B", "lower95", "upper95", "replicates"])
        for (chrom, start, end), scenarios in sorted(groups.items()):
            neutral, bgs = scenarios.get("neutral", {}), scenarios.get("bgs", {})
            if len(neutral) < 2 or neutral.keys() != bgs.keys():
                raise ValueError("Need >=2 matched neutral/BGS replicate IDs per window")
            ids = sorted(neutral)

            def ratio(sample: list[str]) -> float:
                denominator = sum(neutral[i] for i in sample)
                if denominator <= 0:
                    raise ValueError("Zero neutral diversity in B estimate/bootstrap")
                return sum(bgs[i] for i in sample) / denominator

            values = sorted(ratio(rng.choices(ids, k=len(ids))) for _ in range(bootstrap))
            writer.writerow([chrom, start, end, ratio(ids),
                             values[int(.025 * (bootstrap - 1))],
                             values[int(.975 * (bootstrap - 1))], len(ids)])
