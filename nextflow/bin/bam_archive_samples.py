#!/usr/bin/env python3
"""Select diploid BAMs from a GHIST testing archive using popfile.txt."""

import argparse
import csv
import subprocess
from pathlib import Path


def select(archive, destination, output, maximum=2):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    subprocess.run(["tar", "xzf", str(archive), "-C", str(destination)], check=True)
    popfiles = list(destination.rglob("popfile.txt"))
    if len(popfiles) != 1:
        raise ValueError("Archive must contain exactly one popfile.txt")
    with open(popfiles[0]) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows or any(not row.get("sample") or not row.get("pop") for row in rows):
        raise ValueError("popfile.txt requires sample and pop columns and nonempty rows")
    if len({row["sample"] for row in rows}) != len(rows):
        raise ValueError("Duplicate sample identifiers in popfile.txt")
    populations = list(dict.fromkeys(row["pop"] for row in rows))
    if len(populations) > 2:
        raise ValueError("Survey BAM preparation supports at most two populations")
    selected = []
    for population in populations:
        selected.append(next(row for row in rows if row["pop"] == population))
    selected = selected[:maximum]
    for row in rows:
        if len(selected) >= maximum:
            break
        if row not in selected:
            selected.append(row)
    if len(selected) != maximum:
        raise ValueError(f"Archive needs at least {maximum} individuals")
    bams = {}
    for bam in destination.rglob("*.bam"):
        header = subprocess.run(["samtools", "view", "-H", str(bam)], check=True,
                                capture_output=True, text=True).stdout
        samples = {field[3:] for line in header.splitlines() if line.startswith("@RG\t")
                   for field in line.split("\t") if field.startswith("SM:")}
        if len(samples) != 1:
            raise ValueError(f"{bam}: expected exactly one read-group sample name")
        sample = samples.pop()
        if sample in bams:
            raise ValueError(f"Multiple BAMs for sample {sample}; merge before running")
        bams[sample] = str(bam)
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "pop", "bam"])
        for row in selected:
            if row["sample"] not in bams:
                raise ValueError(f"No BAM read-group SM matches sample {row['sample']}")
            writer.writerow([row["sample"], row["pop"], bams[row["sample"]]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive")
    parser.add_argument("destination")
    parser.add_argument("output")
    parser.add_argument("--max-samples", type=int, choices=[1, 2], default=2)
    args = parser.parse_args()
    select(args.archive, args.destination, args.output, args.max_samples)


if __name__ == "__main__":
    main()
