"""Reference lookup, streamed VCF conversion, and SweepFinder2 map adapters."""

from bisect import bisect_right
from dataclasses import dataclass
from itertools import groupby
import csv
import math
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Iterator


@dataclass(frozen=True)
class Variant:
    chrom: str
    position: int
    ref: str
    alt: str
    x: int
    n: int
    folded: bool


class Reference:
    """Pinned primary chromosome sequences from a locally populated SeqRepo."""

    ACCESSIONS = {"testing": ("7", "NC_000007.14", 159345973),
                  "final": ("9", "NC_000009.12", 138394717)}

    def __init__(self, root: Path, dataset: str):
        from biocommons.seqrepo import SeqRepo

        self.chrom, self.accession, self.length = self.ACCESSIONS[dataset]
        self.repo = SeqRepo(str(root), writeable=False)
        # Versioned accessions avoid ambiguous GRCh37/GRCh38 aliases.
        self.fetch(0, 1)

    def fetch(self, start: int, end: int) -> str:
        if not 0 <= start < end <= self.length:
            raise ValueError(f"Reference interval outside {self.accession}: {start}:{end}")
        sequence = self.repo.fetch(self.accession, start=start, end=end).upper()
        if len(sequence) != end - start:
            raise ValueError(f"Truncated SeqRepo sequence for {self.accession}")
        return sequence

    def check(self, variant: Variant) -> None:
        if variant.chrom.removeprefix("chr") != self.chrom:
            raise ValueError(f"VCF contig {variant.chrom} does not match {self.accession}")
        if self.fetch(variant.position - 1, variant.position) != variant.ref:
            raise ValueError(
                f"REF mismatch at {variant.chrom}:{variant.position} against {self.accession}; "
                "chromosome-like simulations are not necessarily GRCh38 sequence. "
                "Do not attach human contexts to incompatible simulated alleles."
            )


def variants(vcf: Path, folded: bool, mask: Path | None = None,
             bcftools: str = "bcftools") -> Iterator[Variant]:
    """Stream biallelic SNPs, excluding duplicate positions and missing allele copies."""
    with tempfile.TemporaryDirectory(prefix="sweep-vcf-") as directory:
        selected = Path(directory) / "selected.bcf"
        command = [bcftools, "view", "-Ob", "-o", str(selected)]
        if mask is not None:
            command.extend(["-T", str(mask), "--targets-overlap", "0"])
        subprocess.run([*command, str(vcf)], check=True)
        # An absent AA tag is acceptable only in folded mode.
        with tempfile.TemporaryFile(mode="w+t") as errors:
            process = subprocess.Popen(
                [bcftools, "query", "-u", "-f",
                 "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AA[\\t%GT]\\n", str(selected)],
                stdout=subprocess.PIPE, stderr=errors, text=True,
            )
            if process.stdout is None:
                raise RuntimeError("bcftools query has no output pipe")
            try:
                previous: tuple[str, int] | None = None
                contig: str | None = None
                for key, records in groupby(
                        (line.rstrip("\n").split("\t") for line in process.stdout),
                        key=lambda row: (row[0], int(row[1]))):
                    rows = list(records)
                    if contig is None:
                        contig = key[0]
                    if key[0] != contig or (previous is not None and key <= previous):
                        raise ValueError("Provide one sorted contig per experiment case")
                    previous = key
                    if len(rows) != 1:
                        continue
                    chrom, pos, ref, alt, ancestor, *genotypes = rows[0]
                    if ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1:
                        continue
                    alleles = [int(a) for gt in genotypes for a in re.split(r"[/|]", gt)
                               if a != "."]
                    if any(a not in (0, 1) for a in alleles):
                        raise ValueError(f"Non-biallelic genotype at {chrom}:{pos}")
                    n, alternate = len(alleles), sum(alleles)
                    if not n or alternate in (0, n):
                        continue
                    if folded:
                        x = min(alternate, n - alternate)
                    elif ancestor in (ref, alt):
                        x = alternate if ancestor == ref else n - alternate
                    else:
                        raise ValueError(f"Missing/invalid INFO/AA at {chrom}:{pos}")
                    yield Variant(chrom, int(pos), ref, alt, x, n, folded)
                returncode = process.wait()
                if returncode:
                    errors.seek(0)
                    raise RuntimeError(f"bcftools query failed ({returncode}): {errors.read()}")
            finally:
                process.stdout.close()
                if process.poll() is None:
                    process.terminate()
                    process.wait()


def write_frequency(path: Path, sites: list[Variant]) -> None:
    if not sites:
        raise ValueError("No eligible polymorphic SNPs")
    with path.open("w") as handle:
        handle.write("position\tx\tn\tfolded\n")
        for site in sites:
            handle.write(f"{site.position}\t{site.x}\t{site.n}\t{int(site.folded)}\n")


def read_frequency(path: Path, chrom: str, folded: bool) -> list[Variant]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sites = [Variant(chrom, int(row["position"]), "N", "N", int(row["x"]),
                     int(row["n"]), bool(int(row["folded"]))) for row in rows]
    if (not sites or any(site.position < 1 or site.n < 2 or not 0 <= site.x <= site.n
                         or site.folded != folded for site in sites)
            or any(b.position <= a.position for a, b in zip(sites, sites[1:]))
            or any(row["folded"] not in ("0", "1") for row in rows)):
        raise ValueError(f"{path}: invalid or inconsistent per-site frequency input")
    return sites


def check_mask(sites: list[Variant], path: Path) -> None:
    intervals: list[tuple[int, int]] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, *_ = line.split()
            if chrom != sites[0].chrom:
                continue
            left, right = int(start), int(end)
            if left < 0 or right <= left:
                raise ValueError(f"{path}: invalid callable interval")
            intervals.append((left, right))
    merged: list[tuple[int, int]] = []
    for left, right in sorted(intervals):
        if merged and left <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(right, merged[-1][1])
        else:
            merged.append((left, right))
    starts = [left for left, _ in merged]
    for site in sites:
        index = bisect_right(starts, site.position - 1) - 1
        if index < 0 or site.position - 1 >= merged[index][1]:
            raise ValueError(f"Frequency input contains a noncallable position: {site.position}")


def cumulative_map(path: Path) -> tuple[list[int], list[float]]:
    """Read position/cM TSV with a header, using 1-based physical coordinates."""
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    positions = [int(row["position"]) for row in rows]
    values = [float(row["cM"]) for row in rows]
    if (len(rows) < 2 or positions[0] < 1
            or any(b <= a for a, b in zip(positions, positions[1:]))
            or any(not math.isfinite(v) or v < 0 for v in values)
            or any(b < a for a, b in zip(values, values[1:]))):
        raise ValueError(f"{path}: expected increasing positions and nondecreasing finite cM")
    return positions, values


def interpolate(position: int, positions: list[int], values: list[float]) -> float:
    if not positions[0] <= position <= positions[-1]:
        raise ValueError(f"Map does not span position {position}")
    index = min(bisect_right(positions, position) - 1, len(positions) - 2)
    fraction = (position - positions[index]) / (positions[index + 1] - positions[index])
    return values[index] + fraction * (values[index + 1] - values[index])


def write_recombination(path: Path, sites: list[Variant], map_path: Path | None,
                        constant_rate: float) -> None:
    if not math.isfinite(constant_rate) or constant_rate < 0:
        raise ValueError("Recombination rate must be finite and nonnegative")
    coordinates, values = cumulative_map(map_path) if map_path else ([], [])
    previous: float | None = None
    with path.open("w") as handle:
        handle.write("position\trate\n")
        for site in sites:
            cm = (interpolate(site.position, coordinates, values) if map_path
                  else site.position * constant_rate * 100)
            distance = 0.0 if previous is None else cm - previous
            handle.write(f"{site.position}\t{distance:.12g}\n")
            previous = cm


def write_bvalues(path: Path, sites: list[Variant], bedgraph: Path) -> None:
    """Sample headerless chromosome/start/end/B intervals (0-based half-open)."""
    starts: list[int] = []
    ends: list[int] = []
    values: list[float] = []
    for line in bedgraph.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start_text, end_text, value_text = line.split()
        if chrom != sites[0].chrom:
            continue
        start, end, value = int(start_text), int(end_text), float(value_text)
        if (start < 0 or end <= start or (ends and start < ends[-1])
                or not math.isfinite(value) or not 0 < value <= 1):
            raise ValueError("B intervals must be sorted, disjoint, and have 0 < B <= 1")
        starts.append(start)
        ends.append(end)
        values.append(value)
    with path.open("w") as handle:
        handle.write("position\tbvalue\n")
        for site in sites:
            index = bisect_right(starts, site.position - 1) - 1
            if index < 0 or site.position - 1 >= ends[index]:
                raise ValueError(f"B map does not cover {site.chrom}:{site.position}")
            handle.write(f"{site.position}\t{values[index]:.12g}\n")
