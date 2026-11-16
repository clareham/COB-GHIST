"""Demographic-history challenge plots: PSMC, MSMC2, and fastsimcoal2."""

from collections import Counter
from pathlib import Path

from .common import Plotter, new_plot, number, plt, read_table
from .demographic_curves import (
    change_estimates, mutation_rate, rccr_split, read_psmc, recent_size,
    scaled_msmc, scaled_psmc,
)
from .summary import ChallengeSummary, Estimate, overview_axes
from matplotlib.figure import Figure


def plot_psmc(path: Path) -> Figure:
    """Use RS records from the last complete iteration, in native PSMC units."""
    times, sizes, _ = read_psmc(path)
    fig, ax = new_plot()
    rate = mutation_rate(path)
    if rate is not None:
        curve = scaled_psmc(path, rate)
        ax.step(curve.left, curve.size, where="post", label=curve.label)
        ax.set(xlabel="Generations before present", ylabel="Diploid effective population size",
               title=f"PSMC: final complete iteration (mu={rate:g}; recorded/pipeline bin size)")
        ax.legend()
        return fig
    ax.step(times, sizes, where="post")
    ax.set(xlabel="Time before present (2N0 generations)",
           ylabel="Relative effective population size (N / N0)",
           title="PSMC: final complete iteration (mutation-rate metadata unavailable)")
    return fig


def plot_msmc(path: Path) -> Figure:
    rows = read_table(path)
    columns = [key for key in rows[0] if key.startswith("lambda")]
    if not columns:
        raise ValueError("Missing MSMC2 coalescence-rate columns")
    left = [number(row["left_time_boundary"]) for row in rows]
    fig, ax = new_plot()
    rate = mutation_rate(path)
    if rate is not None:
        left = [time / rate for time in left]
    for column in columns:
        rates = [number(row[column]) for row in rows]
        if any(rate < 0 for rate in rates):
            raise ValueError("Coalescence rates cannot be negative")
        # The last interval may end at infinity; do not invent a finite end.
        ax.step(left, rates, where="post", label=column)
    ax.set(xlabel=("Generations before present" if rate is not None
                   else "Time before present (mutations per site; scaling unavailable)"),
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


def compose_challenge(challenge: Path) -> ChallengeSummary:
    """Compare within-population Ne curves and RCCR without pooling estimators."""
    psmc = sorted((challenge / "psmc").glob("*.psmc"))
    msmc = sorted((challenge / "msmc2").glob("*.msmc2.final.txt"))
    rccr = sorted((challenge / "msmc2").glob("*.run/rccr.tsv"))
    figure, axes = overview_axes(2 if rccr else 1, challenge.name + ": demographic history")
    result = ChallengeSummary(figure)
    curves = []
    for path in psmc + msmc:
        rate = mutation_rate(path)
        source = path.relative_to(challenge).as_posix()
        if rate is None:
            result.notes.append(f"{source}: no published mutation rate; omitted from scaled overlay.")
            continue
        result.notes.append(f"{source}: mu={rate:g} per base per generation (published survey assumption).")
        curves.extend([scaled_psmc(path, rate)] if path in psmc else scaled_msmc(path, rate))
    for curve in curves:
        source = curve.source.relative_to(challenge).as_posix()
        times = curve.left
        sizes = curve.size
        if curve.right[-1] != float("inf"):
            times = times + [curve.right[-1]]
            sizes = sizes + [sizes[-1]]
        population_names = {
            "secondary_contact": ("mainland", "island"),
            "split_migration": ("AFR", "EUR"),
            "split_isolation_low_coverage": ("pop1", "pop2"),
        }
        label = curve.label
        if label in {"MSMC2 lambda_00", "MSMC2 lambda_11"} and challenge.name in population_names:
            label = "MSMC2 " + population_names[challenge.name][int(label.endswith("11"))]
        line, = axes[0].step(times, sizes, where="post", label=label)
        if challenge.name.startswith("bottleneck") or challenge.name == "growth":
            estimates = change_estimates(curve, challenge.name == "growth", source)
            result.estimates.extend(estimates)
            for estimate in estimates:
                if estimate.target == "generations" and isinstance(estimate.value, (int, float)):
                    axes[0].axvline(estimate.value, color=line.get_color(), linestyle=":",
                                    label=curve.label + " two-epoch change proxy")
        elif curve.label in {"MSMC2 lambda_00", "MSMC2 lambda_11"}:
            populations = {
                "secondary_contact": ("mainland_population_size", "island_population_size"),
                "split_migration": ("AFR_population_size", "EUR_population_size"),
                "split_isolation_low_coverage": ("pop1_population_size", "pop2_population_size"),
            }
            if challenge.name in populations:
                target = populations[challenge.name][int(curve.label.endswith("11"))]
                result.estimates.append(Estimate(
                    target, recent_size(curve), "diploid individuals",
                    "Median Ne of the first three MSMC2 within-population bins; "
                    "recent-size proxy, not a fitted census population size.", source))
    ax = axes[0]
    ax.set(xlabel="Generations before present", ylabel="Diploid effective population size (Ne)",
           title="Within-population histories (cross-population rates are excluded)",
           xscale="symlog", yscale="log")
    ax.set_xlim(left=0)
    if curves:
        ax.legend(fontsize=8)
    else:
        ax.text(.5, .5, "No scalable demographic curves available", ha="center", transform=ax.transAxes)
    for path in rccr:
        rows = read_table(path)
        source = path.relative_to(challenge).as_posix()
        axes[1].step([number(row["left_generations"]) for row in rows],
                     [number(row["relative_cross_coalescence"]) for row in rows],
                     where="post", label="MSMC2 " + path.parent.name)
        estimate = rccr_split(rows, source)
        result.estimates.append(estimate)
        if isinstance(estimate.value, (float, int)):
            axes[1].axvline(estimate.value, color="grey", linestyle=":",
                            label="Sustained half-crossing proxy")
    if rccr:
        axes[1].axhline(.5, color="black", linestyle="--", label="RCCR = 0.5")
        axes[1].set(xlabel="Generations before present", ylabel="Relative cross-coalescence",
                    title="Split-time proxy; migration confounds this crossing", xscale="symlog")
        axes[1].set_xlim(left=0)
        axes[1].legend(fontsize=8)
    targets = {
        "secondary_contact": {"mainland_population_size": "diploid individuals",
                              "island_population_size": "diploid individuals",
                              "generations_since_split": "generations",
                              "generations_since_migration": "generations",
                              "migration_rate": "per generation"},
        "split_migration": {"AFR_population_size": "diploid individuals",
                            "EUR_population_size": "diploid individuals",
                            "generations_since_split": "generations", "migration_rate": "per generation"},
        "split_isolation_low_coverage": {"pop1_population_size": "diploid individuals",
                                        "pop2_population_size": "diploid individuals",
                                        "generations_since_split": "generations"},
    }
    if challenge.name.startswith("bottleneck") or challenge.name == "growth":
        targets[challenge.name] = {"generations": "generations",
                                  ("post_growth_scaling" if challenge.name == "growth"
                                   else "post_decline_fraction"): "ratio"}
    for target, units in targets.get(challenge.name, {}).items():
        if not any(item.target == target for item in result.estimates):
            result.estimates.append(Estimate(
                target, None, units, "Not identifiable from the available curves; requires "
                "a fitted demographic model or additional data.", "", "unavailable"))
    result.notes.extend([
        "PSMC: N0 = theta / (4 * mu * bin_bp), t = 2*N0*t_RS, Ne = N0*lambda_RS. "
        "Uses TR theta from the final complete iteration and 100-bp pipeline bins "
        "unless window_size_bp is explicitly recorded. Mu is shared with this challenge's MSMC2 run.",
        "MSMC2: generations = boundary / mu; Ne = 1/(2*mu*within_rate). "
        "lambda_00/population 0 and lambda_11/population 1 follow the published population order.",
        "All time axes are generations, not calendar years. Survey mutation rates, "
        "low coverage, and background selection can bias estimates. No confidence intervals "
        "or combined consensus are inferred. Infinite terminal bins are not used in change fits.",
    ])
    changes = [item.value for item in result.estimates
               if item.target == "generations" and isinstance(item.value, (int, float))]
    if changes and max(changes) / min(changes) > 4:
        result.notes.append("WARNING: per-tool change-time proxies differ by more than fourfold. "
                            "This disagreement is evidence of instability, not a basis for averaging.")
    return result


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
