"""Relatedness challenge plots: ADMIXTURE, KING, and ngsRelate."""

import math
from pathlib import Path

from .common import Plotter, new_plot, number, plt, read_table
from matplotlib.figure import Figure


def plot_admixture(path: Path) -> Figure:
    with path.open(encoding="utf-8") as handle:
        rows = [[number(value) for value in line.split()] for line in handle if line.strip()]
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("Empty or ragged ADMIXTURE Q matrix")
    if any(any(value < 0 or value > 1 for value in row) or
           not math.isclose(sum(row), 1, abs_tol=1e-4) for row in rows):
        raise ValueError("Ancestry proportions must be in [0, 1] and sum to one")
    fig, ax = new_plot()
    bottom = [0.0] * len(rows)
    for k in range(len(rows[0])):
        values = [row[k] for row in rows]
        ax.bar(range(1, len(rows) + 1), values, bottom=bottom, width=1, label=f"Component {k + 1}")
        bottom = [a + b for a, b in zip(bottom, values)]
    ax.set(xlabel="Individual (Q-file row order; sample IDs not published)",
           ylabel="Ancestry proportion", ylim=(0, 1), title="ADMIXTURE")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    return fig


def plot_relatedness(path: Path) -> Figure:
    rows = read_table(path)
    columns = ["KINSHIP"] if "KINSHIP" in rows[0] else ["rab", "theta", "KING"]
    fig, axes = plt.subplots(1, len(columns), figsize=(5 * len(columns), 4.5),
                             layout="constrained", squeeze=False)
    for ax, column in zip(axes[0], columns):
        values = [number(row[column]) for row in rows]
        ax.hist(values, bins=60)
        ax.set(xlabel=column, ylabel="Number of pairs",
               title=f"{column} ({len(values):,} pairs)")
    fig.suptitle("Pairwise relatedness estimates (not relationship classifications)")
    return fig


def plotter_for(path: Path) -> Plotter | None:
    if path.suffix == ".Q":
        return plot_admixture
    if path.name.endswith((".kin0", ".ngsrelate.res")):
        return plot_relatedness
    return None
