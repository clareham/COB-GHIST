"""Sweep-detection challenge plots: SweepFinder2 and FlexSweep."""

from pathlib import Path

from .common import Plotter, new_plot, number, read_table
from matplotlib.figure import Figure


def plot_sweepfinder(path: Path) -> Figure:
    rows = read_table(path)
    fig, ax = new_plot()
    ax.plot([number(row["location"]) for row in rows], [number(row["LR"]) for row in rows])
    ax.set(xlabel="Position (bp)", ylabel="Composite likelihood ratio",
           title="SweepFinder2 scan (no significance threshold fitted)")
    return fig


def plot_flexsweep(path: Path) -> Figure:
    rows = read_table(path)
    fig, ax = new_plot()
    for chrom in dict.fromkeys(row["chrom"] for row in rows):
        group = [row for row in rows if row["chrom"] == chrom]
        probabilities = [number(row["prob_sweep"]) for row in group]
        if any(value < 0 or value > 1 for value in probabilities):
            raise ValueError("Sweep probabilities must be in [0, 1]")
        ax.plot([(number(row["start"]) + number(row["end"])) / 2 for row in group],
                probabilities, label=chrom)
    ax.set(xlabel="Window midpoint (bp)", ylabel="Predicted sweep probability",
           ylim=(0, 1), title="FlexSweep predictions")
    ax.legend(title="Contig")
    return fig


def plotter_for(path: Path) -> Plotter | None:
    if path.name.endswith(".sf2.out"):
        return plot_sweepfinder
    if path.name == "predictions.tsv" and "flexsweep" in path.parts:
        return plot_flexsweep
    return None
