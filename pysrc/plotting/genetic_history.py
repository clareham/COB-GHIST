"""Demographic-history challenge plots: PSMC, MSMC2, and fastsimcoal2."""

from collections import Counter
from pathlib import Path

from .common import Plotter, new_plot, number, plt, read_table
from matplotlib.figure import Figure


def plot_psmc(path: Path) -> Figure:
    """Use RS records from the last complete iteration, in native PSMC units."""
    rounds: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] | None = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "RD":
                current = []
            elif fields[0] == "RS" and current is not None:
                current.append((number(fields[2]), number(fields[3])))
            elif fields[0] == "//" and current:
                rounds.append(current)
                current = None
    if not rounds:
        raise ValueError("No complete PSMC iteration containing RS records")
    times, sizes = zip(*rounds[-1])
    if any(time < 0 for time in times) or any(size <= 0 for size in sizes):
        raise ValueError("PSMC times must be nonnegative and relative sizes positive")
    fig, ax = new_plot()
    ax.step(times, sizes, where="post")
    ax.set(xlabel="Time before present (2N0 generations)",
           ylabel="Relative effective population size (N / N0)",
           title=f"PSMC: final complete iteration ({len(rounds)} available)")
    return fig


def plot_msmc(path: Path) -> Figure:
    rows = read_table(path)
    columns = [key for key in rows[0] if key.startswith("lambda")]
    if not columns:
        raise ValueError("Missing MSMC2 coalescence-rate columns")
    left = [number(row["left_time_boundary"]) for row in rows]
    fig, ax = new_plot()
    for column in columns:
        rates = [number(row[column]) for row in rows]
        if any(rate < 0 for rate in rates):
            raise ValueError("Coalescence rates cannot be negative")
        # The last interval may end at infinity; do not invent a finite end.
        ax.step(left, rates, where="post", label=column)
    ax.set(xlabel="Time before present (mutations per site)",
           ylabel="Coalescence rate (per mutation unit)", title="MSMC2 coalescence rates")
    ax.legend()
    return fig


def plot_rccr(path: Path) -> Figure:
    rows = read_table(path)
    fig, ax = new_plot()
    ax.step([number(row["left_generations"]) for row in rows],
            [number(row["relative_cross_coalescence"]) for row in rows], where="post")
    ax.axhline(0.5, color="grey", linestyle="--", label="0.5 reference (not a fitted split time)")
    ax.set(xlabel="Generations before present", ylabel="Relative cross-coalescence rate",
           title="MSMC2 relative cross-coalescence")
    ax.legend()
    return fig


def plot_multihetsep(path: Path) -> Figure:
    counts: Counter[str] = Counter()
    callable_bases: Counter[str] = Counter()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            chrom, position, distance, pattern = line.split()
            if int(position) <= 0 or int(distance) < 0 or not pattern:
                raise ValueError("Invalid multihetsep record")
            counts[chrom] += 1
            callable_bases[chrom] += int(distance)
    if not counts:
        raise ValueError("Empty multihetsep input")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), layout="constrained")
    axes[0].bar(list(counts), list(counts.values()))
    axes[1].bar(list(callable_bases), list(callable_bases.values()))
    axes[0].set(xlabel="Contig", ylabel="Segregating records")
    axes[1].set(xlabel="Contig", ylabel="Callable bases through last record")
    fig.suptitle("MSMC2 input diagnostics")
    return fig


def plot_msmc_loop(path: Path) -> Figure:
    with path.open(encoding="utf-8") as handle:
        likelihoods = [number(line.split()[1]) for line in handle if line.strip()]
    if not likelihoods:
        raise ValueError("Empty MSMC2 iteration log")
    fig, ax = new_plot()
    ax.plot(range(1, len(likelihoods) + 1), likelihoods, marker=".")
    ax.set(xlabel="Iteration", ylabel="Log likelihood", title="MSMC2 convergence")
    return fig


def plot_fsc_parameters(path: Path) -> Figure:
    with path.open(encoding="utf-8") as handle:
        lines = [line.split() for line in handle if line.strip()]
    if len(lines) < 2 or any(len(row) != len(lines[0]) for row in lines[1:]):
        raise ValueError("Invalid fastsimcoal2 best-likelihood table")
    fig, ax = new_plot()
    for index, row in enumerate(lines[1:], 1):
        ax.plot(lines[0], [number(value) for value in row], marker="o", label=f"Fit {index}")
    ax.set(yscale="symlog", ylabel="Fitted value (symlog; parameter-specific units)",
           title="fastsimcoal2 parameters and likelihoods")
    ax.tick_params(axis="x", rotation=45)
    ax.legend()
    return fig


def plot_sfs(path: Path) -> Figure:
    with path.open(encoding="utf-8") as handle:
        lines = [line.split() for line in handle if line.strip()]
    if len(lines) < 3:
        raise ValueError("Missing observed spectrum header or counts")
    joint = lines[2][0].startswith("d1_")
    values = [[number(value) for value in (row[1:] if joint else row)] for row in lines[2:]]
    if not values[0] or any(len(row) != len(values[0]) for row in values):
        raise ValueError("Ragged observed spectrum")
    if any(value < 0 for row in values for value in row):
        raise ValueError("Observed spectrum counts cannot be negative")
    fig, ax = new_plot()
    if not joint:
        ax.bar(range(len(values[0])), values[0])
        ax.set(xlabel="Allele-count bin", ylabel="Observed sites")
    else:
        image = ax.imshow(values, origin="lower", aspect="auto")
        fig.colorbar(image, ax=ax, label="Observed sites")
        ax.set(xlabel="Population 0 derived allele count", ylabel="Population 1 derived allele count")
    ax.set_title("Observed site-frequency spectrum (includes invariant bin)")
    return fig


def plotter_for(path: Path) -> Plotter | None:
    name = path.name
    if name.endswith(".psmc"):
        return plot_psmc
    if name.endswith(".final.txt") and "msmc2" in path.parts:
        return plot_msmc
    if name == "rccr.tsv":
        return plot_rccr
    if name.endswith(".multihetsep.txt"):
        return plot_multihetsep
    if name.endswith(".loop.txt") and "msmc2" in path.parts:
        return plot_msmc_loop
    if path.suffix == ".bestlhoods":
        return plot_fsc_parameters
    if path.suffix == ".obs":
        return plot_sfs
    return None
