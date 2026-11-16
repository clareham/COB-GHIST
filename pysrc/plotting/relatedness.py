"""Relatedness challenge plots: ADMIXTURE, KING, and ngsRelate."""

import math
from collections import Counter
from pathlib import Path

from .common import Plotter, new_plot, number, plt, read_table
from .summary import ChallengeSummary, Estimate, overview_axes
from matplotlib.figure import Figure


def read_ancestry(path: Path) -> list[list[float]]:
    with path.open(encoding="utf-8") as handle:
        rows = [[number(value) for value in line.split()] for line in handle if line.strip()]
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("Empty or ragged ADMIXTURE Q matrix")
    if any(any(value < 0 or value > 1 for value in row) or
           not math.isclose(sum(row), 1, abs_tol=1e-4) for row in rows):
        raise ValueError("Ancestry proportions must be in [0, 1] and sum to one")
    return rows


def draw_ancestry(ax, rows: list[list[float]]) -> None:
    bottom = [0.0] * len(rows)
    for k in range(len(rows[0])):
        values = [row[k] for row in rows]
        ax.bar(range(1, len(rows) + 1), values, bottom=bottom, width=1, label=f"Component {k + 1}")
        bottom = [a + b for a, b in zip(bottom, values)]
    ax.set(xlabel="Individual (Q-file row order; sample IDs not published)",
           ylabel="Ancestry proportion", ylim=(0, 1), title="ADMIXTURE")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))


def plot_admixture(path: Path) -> Figure:
    fig, ax = new_plot()
    draw_ancestry(ax, read_ancestry(path))
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


def relationship_degree(kinship: float) -> str:
    for boundary, label in ((.354, "duplicate/MZ-like"), (.177, "first degree"),
                            (.0884, "second degree"), (.0442, "third degree")):
        if kinship >= boundary:
            return label
    return "more distant/unrelated"


def compose_challenge(challenge: Path) -> ChallengeSummary:
    ancestry = sorted((challenge / "admixture").glob("*.Q"))
    king = sorted((challenge / "king_plink").glob("*.kin0"))
    ngs = sorted((challenge / "ngsrelate").glob("*.ngsrelate.res"))
    figure, axes = overview_axes(len(ancestry) + 2, challenge.name + ": relatedness")
    result = ChallengeSummary(figure)
    for ax, path in zip(axes, ancestry):
        rows = read_ancestry(path)
        draw_ancestry(ax, rows)
        result.estimates.extend([
            Estimate("individuals_in_Q", len(rows), "individuals", "Q-matrix row count",
                     path.relative_to(challenge).as_posix(), "descriptive"),
            Estimate("ancestry_components", len(rows[0]), "configured components",
                     "Q-matrix column count; configured K, not inferred population count",
                     path.relative_to(challenge).as_posix(), "descriptive"),
        ])
    histogram = axes[-2]
    counts: Counter[str] = Counter()
    distributions = []
    for path in king + ngs:
        rows = read_table(path)
        source = path.relative_to(challenge).as_posix()
        columns = ["KINSHIP"] if path in king else ["theta", "KING"]
        for column in columns:
            values = [number(row[column]) for row in rows]
            distributions.append((values, ("PLINK " if path in king else "ngsRelate ") + column))
        if path in king:
            for row in rows:
                kinship = number(row["KINSHIP"])
                degree = relationship_degree(kinship)
                counts[degree] += 1
                result.assignments.append({"individual_1": row["#IID1"], "individual_2": row["IID2"],
                                           "kinship": kinship, "degree_proxy": degree, "source": source})
            for label, count in sorted(Counter(
                    relationship_degree(number(row["KINSHIP"])) for row in rows).items()):
                result.estimates.append(Estimate(
                    label, count, "pairs", "KING kinship thresholds 0.354/0.177/0.0884/0.0442; "
                    "provisional degree bins, not pedigree labels.", source))
    if distributions:
        low = min(min(values) for values, _ in distributions)
        high = max(max(values) for values, _ in distributions)
        if low == high:
            low, high = low - .01, high + .01
        bins = [low + (high - low) * index / 60 for index in range(61)]
        for values, label in distributions:
            histogram.hist(values, bins=bins, histtype="step", label=label)
    histogram.set(xlabel="Kinship coefficient", ylabel="Number of pairs",
                  title="Kinship estimators (rab is deliberately not mixed with kinship)")
    if king or ngs:
        histogram.legend(fontsize=8)
    axes[-1].bar(list(counts), list(counts.values()), label="PLINK KING threshold classification")
    axes[-1].set(ylabel="Number of pairs", title="Exploratory pairwise degree estimates")
    axes[-1].tick_params(axis="x", rotation=15)
    axes[-1].legend()
    result.estimates.append(Estimate(
        "pedigree_relationships", None, "relationships",
        "Degree bins cannot distinguish parent-offspring from siblings or reconstruct a pedigree; "
        "requires IBD/pedigree inference and validation.", "", "unavailable"))
    result.notes.extend([
        "Pairwise degree proxies with published individual IDs are exported in challenge_estimates.json. "
        "No IDs are guessed for Q rows or ngsRelate numeric indices, and pairs are not matched across tools.",
        "KING cutoffs are conventional exploratory thresholds, not species-/cohort-calibrated. "
        "Dog inbreeding and population structure can bias labels. No parent/offspring, sibling, or "
        "ancestry population names are asserted. Negative coefficients are retained.",
        "These outputs have no temporal axis; generation scaling is not applicable.",
    ])
    return result


def plotter_for(path: Path) -> Plotter | None:
    if path.suffix == ".Q":
        return plot_admixture
    if path.name.endswith((".kin0", ".ngsrelate.res")):
        return plot_relatedness
    return None
