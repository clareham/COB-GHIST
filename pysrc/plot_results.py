"""Plot published survey outputs and assemble an offline HTML report.

Run from the repository root with ``python -m pysrc.plot_results``.
"""

import argparse
import base64
import csv
import html
import io
import json
import math
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Callable
from urllib.parse import quote

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.figure import Figure


DEFAULT_RESULTS = Path(__file__).resolve().parents[1] / "nextflow" / "results"


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


def new_plot() -> tuple[Figure, plt.Axes]:
    return plt.subplots(figsize=(10, 4.5), layout="constrained")


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


def plotter_for(path: Path) -> Callable[[Path], Figure] | None:
    name = path.name
    if name.endswith(".psmc"):
        return plot_psmc
    if name.endswith(".final.txt") and "msmc2" in path.parts:
        return plot_msmc
    if name == "rccr.tsv":
        return plot_rccr
    if path.suffix == ".Q":
        return plot_admixture
    if name.endswith((".kin0", ".ngsrelate.res")):
        return plot_relatedness
    if name.endswith(".sf2.out"):
        return plot_sweepfinder
    if name == "predictions.tsv" and "flexsweep" in path.parts:
        return plot_flexsweep
    if name.endswith(".normalise.tsv"):
        return plot_normalisation
    if path.suffix == ".bed":
        return plot_bed
    if name.endswith(".multihetsep.txt"):
        return plot_multihetsep
    if name.endswith(".loop.txt") and "msmc2" in path.parts:
        return plot_msmc_loop
    if path.suffix == ".bestlhoods":
        return plot_fsc_parameters
    if path.suffix == ".obs":
        return plot_sfs
    return None


def embed_figure(figure: Figure) -> str:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=110)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f'<img alt="Result plot" src="data:image/png;base64,{encoded}">'


def build_report(results: Path, output: Path | None = None) -> tuple[Path, int, int]:
    results = results.resolve()
    if not results.is_dir():
        raise ValueError(f"Results directory does not exist: {results}")
    output = output.resolve() if output is not None else results / "plots.html"
    files = sorted(path for path in results.rglob("*") if path.is_file() and path.resolve() != output)
    if not files:
        raise ValueError(f"No result files found in {results}")
    groups: dict[str, list[str]] = {}
    inventory = []
    failures = 0
    plots = 0
    for path in files:
        relative = path.relative_to(results)
        label = html.escape(relative.as_posix())
        link = quote(Path(os.path.relpath(path, output.parent)).as_posix(), safe="/")
        source = f'<a href="{link}">{label}</a>'
        status = "Supporting file (linked; no plot adapter)"
        content = ""
        plotter = plotter_for(path)
        try:
            if plotter is not None:
                try:
                    figure = plotter(path)
                    content = embed_figure(figure)
                finally:
                    plt.close("all")
                plots += 1
                status = "Plotted"
            elif path.suffix.lower() in {".svg", ".png", ".jpg", ".jpeg"}:
                mime = {".svg": "image/svg+xml", ".png": "image/png",
                        ".jpg": "image/jpeg", ".jpeg": "image/jpeg"}[path.suffix.lower()]
                encoded = base64.b64encode(path.read_bytes()).decode("ascii")
                content = f'<img alt="{label}" src="data:{mime};base64,{encoded}">'
                plots += 1
                status = "Existing plot embedded"
            elif path.suffix == ".json":
                metadata = json.loads(path.read_text(encoding="utf-8"))
                content = "<details><summary>Metadata</summary><pre>" + html.escape(
                    json.dumps(metadata, indent=2)) + "</pre></details>"
                status = "Metadata displayed"
        except (ValueError, KeyError, IndexError, OSError, csv.Error) as exc:
            failures += 1
            status = f"ERROR: {exc}"
            content = f'<p class="error">{html.escape(status)}</p>'
            print(f"{relative}: {status}", file=sys.stderr)
        if content:
            group = "/".join(relative.parts[:2]) if len(relative.parts) > 2 else "Pipeline"
            groups.setdefault(group, []).append(f"<article><h3>{source}</h3>{content}</article>")
        inventory.append(f"<tr><td>{source}</td><td>{path.stat().st_size:,}</td>"
                         f"<td>{html.escape(status)}</td></tr>")
    sections = []
    navigation = []
    for index, (group, articles) in enumerate(groups.items()):
        navigation.append(f'<a href="#group-{index}">{html.escape(group)}</a>')
        sections.append(f'<section id="group-{index}"><h2>{html.escape(group)}</h2>'
                        + "".join(articles) + "</section>")
    document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GHIST survey results</title><style>
body{font-family:system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem}
img{max-width:100%;height:auto}article{border-top:1px solid #ccc;padding:1rem 0}
h3,td{overflow-wrap:anywhere}table{border-collapse:collapse;width:100%}
td,th{border:1px solid #ccc;padding:.4rem;text-align:left}pre{overflow:auto}
nav{display:flex;flex-wrap:wrap;gap:1rem}.error{color:#b00020}
</style></head><body><h1>GHIST survey results</h1>
<p>Exploratory outputs, not scientifically validated estimates. Missing tools are
not imputed. PSMC uses native relative units; MSMC2 uses mutation-scaled time and
rates (cross-population rates are not effective population sizes). RCCR uses the
generations provided by the pipeline. Ancestry and relatedness plots do not
assign population or relationship labels. BED length sums count overlaps.</p>
<p>All plots and metadata are embedded for offline viewing. Source links require
the results tree. Models, binary intermediates, logs and other unsupported formats
are inventoried below, not executed or deserialized.</p>
"""
    document += f"<p>{len(files)} files; {plots} plots; {failures} errors.</p>"
    document += "<nav>" + "".join(navigation) + "</nav>" + "".join(sections)
    document += "<h2>Complete file inventory</h2><table><thead><tr><th>File</th>"
    document += "<th>Bytes</th><th>Status</th></tr></thead><tbody>"
    document += "".join(inventory) + "</tbody></table></body></html>"
    output.write_text(document, encoding="utf-8")
    return output, plots, failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, help="Default: <results-dir>/plots.html")
    args = parser.parse_args()
    try:
        output, plots, failures = build_report(args.results_dir, args.output)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"plot-results: {exc}\n")
    print(f"Wrote {output} ({plots} plots, {failures} errors)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
