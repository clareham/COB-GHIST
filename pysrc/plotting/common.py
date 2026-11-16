"""Shared parsing, plotting setup, and cross-challenge QC plots."""

import csv
import math
from collections import Counter
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

Plotter = Callable[[Path], Figure]


def number(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Expected a finite number, got {value!r}")
    return result


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Missing table header")
        rows = list(reader)
    if not rows:
        raise ValueError("Table contains no data rows")
    if any(None in row or None in row.values() for row in rows):
        raise ValueError("Table row width does not match the header")
    return rows


def new_plot() -> tuple[Figure, Axes]:
    return plt.subplots(figsize=(10, 4.5), layout="constrained")


def plot_normalisation(path: Path) -> Figure:
    rows = read_table(path)
    fig, ax = new_plot()
    ax.barh([row["metric"] for row in rows], [number(row["value"]) for row in rows])
    ax.set(xlabel="Count", title="VCF normalisation")
    return fig


def plot_bed(path: Path) -> Figure:
    """Aggregate interval lengths without claiming overlapping windows are merged."""
    counts: Counter[str] = Counter()
    lengths: Counter[str] = Counter()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, *_ = line.split()
            start_value, end_value = int(start), int(end)
            if start_value < 0 or end_value <= start_value:
                raise ValueError("Invalid BED interval")
            counts[chrom] += 1
            lengths[chrom] += end_value - start_value
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), layout="constrained")
    axes[0].bar(list(counts), list(counts.values()))
    axes[1].bar(list(lengths), list(lengths.values()))
    axes[0].set(xlabel="Contig", ylabel="Intervals")
    axes[1].set(xlabel="Contig", ylabel="Summed interval length (bp; overlaps counted)")
    fig.suptitle("BED intervals" if counts else "BED intervals: no intervals")
    return fig


def plotter_for(path: Path) -> Plotter | None:
    if path.name.endswith(".normalise.tsv"):
        return plot_normalisation
    if path.suffix == ".bed":
        return plot_bed
    return None
