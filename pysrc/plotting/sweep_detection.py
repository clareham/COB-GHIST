"""Sweep-detection challenge plots: SweepFinder2 and FlexSweep."""

from pathlib import Path

from .common import Plotter, new_plot, number, read_table
from .summary import ChallengeSummary, Estimate, overview_axes
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


def merged_intervals(path: Path) -> list[tuple[str, int, int]]:
    intervals = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.split()
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError(f"{path}: invalid sweep BED interval")
            intervals.append((chrom, start, end))
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def compose_challenge(challenge: Path) -> ChallengeSummary:
    scans = sorted((challenge / "sweepfinder2").glob("*.sf2.out"))
    predictions = sorted((challenge / "flexsweep").glob("*.flexsweep/predictions.tsv"))
    figure, axes = overview_axes(2, challenge.name + ": sweep candidates")
    result = ChallengeSummary(figure)
    for path in scans:
        rows = read_table(path)
        positions = [number(row["location"]) for row in rows]
        scores = [number(row["LR"]) for row in rows]
        axes[0].plot(positions, scores, label="SweepFinder2 CLR")
        peak = max(range(len(scores)), key=scores.__getitem__)
        axes[0].axvline(positions[peak], linestyle=":", color="grey", label="Maximum CLR candidate")
        result.estimates.append(Estimate(
            "maximum_CLR_position", positions[peak], "bp (scan coordinate)",
            "Global maximum of the published scan; no null-calibrated significance threshold. "
            "A grid maximum is not a confirmed selected site.",
            path.relative_to(challenge).as_posix()))
    for path in predictions:
        rows = read_table(path)
        for chrom in dict.fromkeys(row["chrom"] for row in rows):
            group = [row for row in rows if row["chrom"] == chrom]
            probabilities = [number(row["prob_sweep"]) for row in group]
            if any(not 0 <= value <= 1 for value in probabilities):
                raise ValueError(f"{path}: invalid sweep probabilities")
            axes[1].plot([(number(row["start"]) + number(row["end"])) / 2 for row in group],
                         probabilities, label=f"FlexSweep contig {chrom}")
        bed = path.parent / "sweeps.bed"
        if not bed.is_file():
            result.notes.append(f"{bed.relative_to(challenge)} missing; no thresholded intervals inferred.")
            continue
        merged = merged_intervals(bed)
        result.estimates.append(Estimate(
            "candidate_sweep_regions", len(merged), "merged regions",
            "Merge overlapping or adjacent published FlexSweep BED windows. "
            "Regions reflect the configured classifier threshold, not confirmed sweep counts.",
            bed.relative_to(challenge).as_posix()))
        for index, (chrom, start, end) in enumerate(merged, 1):
            group = [row for row in rows if row["chrom"] == chrom and
                     int(row["start"]) < end and int(row["end"]) > start]
            if not group:
                raise ValueError(f"{bed}: BED interval has no overlapping prediction windows")
            peak = max(group, key=lambda row: number(row["prob_sweep"]))
            position = (int(peak["start"]) + int(peak["end"])) / 2
            result.estimates.append(Estimate(
                f"candidate_region_{index}_position", position, f"bp on contig {chrom}",
                "Midpoint of highest-probability overlapping window, not a resolved causal variant. "
                "Interval is the merged window extent (0-based half-open), not a confidence interval.",
                path.relative_to(challenge).as_posix(), interval=(float(start), float(end))))
            result.assignments.append({"chrom": chrom, "start": start, "end": end,
                                       "peak_window_midpoint": position,
                                       "peak_probability": number(peak["prob_sweep"]),
                                       "source": bed.relative_to(challenge).as_posix()})
            axes[1].axvspan(start, end, alpha=.12, color="grey",
                            label="Published candidate window union" if index == 1 else None)
    axes[0].set(xlabel="Scan position (bp; contig not recorded by SweepFinder2)",
                ylabel="Composite likelihood ratio", title="SweepFinder2: uncalibrated scan peaks")
    axes[1].set(xlabel="Window midpoint (bp)", ylabel="Sweep probability", ylim=(0, 1),
                title="FlexSweep: classifier probabilities and candidate intervals")
    for ax in axes:
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    if not predictions:
        axes[1].text(.5, .5, "No FlexSweep predictions available", ha="center", transform=axes[1].transAxes)
    for target, units in (("sweep_time", "generations"), ("selection_coefficient", "per generation"),
                          ("hard_or_soft_sweep", "sweep type")):
        result.estimates.append(Estimate(
            target, None, units, "Not identifiable from published CLR/probability outputs; "
            "requires a calibrated parameter/type estimator.", "", "unavailable"))
    result.notes.extend([
        "SweepFinder2 scan coordinates and FlexSweep contig coordinates are kept on separate axes; "
        "no contig identity or cross-tool candidate agreement is assumed.",
        "A high classifier probability or CLR is not a calibrated significance level. "
        "FlexSweep BED windows are merged for reporting, not declared distinct biological events.",
        "SweepFinder2 alpha is not directly reported as a selection coefficient. "
        "No sweep age is read from model-training time ranges or genomic coordinates.",
    ])
    return result


def plotter_for(path: Path) -> Plotter | None:
    if path.name.endswith(".sf2.out"):
        return plot_sweepfinder
    if path.name == "predictions.tsv" and "flexsweep" in path.parts:
        return plot_flexsweep
    return None
