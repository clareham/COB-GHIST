"""Discover result files, dispatch challenge plots, and render their directory tree."""

import base64
import csv
import html
import io
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path
from urllib.parse import quote

from . import common, genetic_history, relatedness, sweep_detection
from .common import Plotter, plt
from .summary import ChallengeSummary
from matplotlib.figure import Figure


def plotter_for(path: Path) -> Plotter | None:
    for module in (common, genetic_history, relatedness, sweep_detection):
        plotter = module.plotter_for(path)
        if plotter is not None:
            return plotter
    return None


def embed_figure(figure: Figure) -> str:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=110)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f'<img alt="Result plot" src="data:image/png;base64,{encoded}">'


def heading(level: int, title: str) -> str:
    # Deep run directories still need accurate accessible heading levels.
    if level <= 6:
        return f"<h{level}>{title}</h{level}>"
    return f'<div role="heading" aria-level="{level}">{title}</div>'


def render_tree(entries: dict[Path, tuple[str, str]], directory: Path = Path(),
                level: int = 2, summaries: dict[Path, str] | None = None) -> tuple[str, str]:
    """Return nested result sections and a matching directory navigation tree."""
    content = [summaries[directory]] if summaries and directory in summaries else []
    navigation = []
    children = set()
    for path, (source, body) in entries.items():
        if path.parent == directory:
            content.append(f'<article data-path="{html.escape(path.as_posix(), quote=True)}">'
                           + heading(level, source) + body + "</article>")
        elif directory in path.parents:
            children.add(directory / path.relative_to(directory).parts[0])
    for child in sorted(children):
        anchor = "dir-" + child.as_posix().encode("utf-8").hex()
        title = html.escape(child.name)
        nested_content, nested_navigation = render_tree(entries, child, level + 1, summaries)
        content.append(f'<section id="{anchor}" data-path="'
                       f'{html.escape(child.as_posix(), quote=True)}">'
                       + heading(level, title) + nested_content + "</section>")
        navigation.append(f'<li><a href="#{anchor}">{title}</a>'
                          + (f"<ul>{nested_navigation}</ul>" if nested_navigation else "") + "</li>")
    return "".join(content), "".join(navigation)


def summary_html(summary: ChallengeSummary, image_link: str) -> str:
    rows = []
    for estimate in summary.estimates:
        value = (f"{estimate.value:,.6g}" if isinstance(estimate.value, (float, int))
                 else str(estimate.value) if estimate.value is not None else "Unavailable")
        bracket = (f" [{estimate.interval[0]:,.6g}, {estimate.interval[1]:,.6g}]"
                   if estimate.interval is not None else "")
        rows.append("<tr>" + "".join(f"<td>{html.escape(str(item))}</td>" for item in (
            estimate.target, value + bracket, estimate.units, estimate.status,
            estimate.method, estimate.source)) + "</tr>")
    return (heading(4, "Challenge overview and target estimates")
            + f'<p><a href="{image_link}" download>Download combined graphic (PNG)</a></p>'
            + embed_figure(summary.figure)
            + "<table class=\"estimates\"><thead><tr><th>Target / summary</th><th>Estimate / bracket</th>"
            "<th>Units</th><th>Status</th><th>Method / limitation</th><th>Source</th></tr></thead><tbody>"
            + "".join(rows) + "</tbody></table><ul>"
            + "".join(f"<li>{html.escape(note)}</li>" for note in summary.notes) + "</ul>")


def build_summaries(results: Path, output: Path) -> tuple[dict[Path, str], int, int]:
    summaries = {}
    records = {}
    failures = 0
    for module in (genetic_history, relatedness, sweep_detection):
        category = module.__name__.rsplit(".", 1)[-1]
        category_dir = results / category
        if not category_dir.is_dir():
            continue
        for challenge in sorted(path for path in category_dir.iterdir() if path.is_dir()):
            relative = challenge.relative_to(results)
            try:
                try:
                    summary = module.compose_challenge(challenge)
                    image = results / "plot_summaries" / category / f"{challenge.name}.png"
                    image.parent.mkdir(parents=True, exist_ok=True)
                    summary.figure.savefig(image, dpi=150)
                    link = quote(Path(os.path.relpath(image, output.parent)).as_posix(), safe="/")
                    summaries[relative] = summary_html(summary, link)
                    records[relative.as_posix()] = {
                        "estimates": [asdict(item) for item in summary.estimates],
                        "notes": summary.notes,
                        "assignments": summary.assignments,
                    }
                finally:
                    plt.close("all")
            except (ValueError, KeyError, IndexError, OSError, csv.Error) as exc:
                failures += 1
                message = f"Challenge overview failed: {exc}"
                print(f"{relative}: {message}", file=sys.stderr)
                summaries[relative] = f'<p class="error">{html.escape(message)}</p>'
                records[relative.as_posix()] = {"error": message}
    if records:
        (results / "challenge_estimates.json").write_text(json.dumps({
            "schema_version": 1,
            "warning": "Exploratory curve-derived summaries, not validated challenge submissions. "
                       "Brackets are bin/window extents, not confidence intervals. "
                       "Sources are relative to each challenge directory.",
            "challenges": records,
        }, indent=2, allow_nan=False), encoding="utf-8")
    return summaries, len(records) - failures, failures


def build_report(results: Path, output: Path | None = None) -> tuple[Path, int, int]:
    results = results.resolve()
    if not results.is_dir():
        raise ValueError(f"Results directory does not exist: {results}")
    output = output.resolve() if output is not None else results / "plots.html"
    files = sorted(path for path in results.rglob("*")
                   if path.is_file() and path.resolve() != output
                   and path.relative_to(results).parts[0] not in
                   {"plot_summaries", "challenge_estimates.json"})
    if not files:
        raise ValueError(f"No result files found in {results}")
    entries = {}
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
        file_link = f'<a href="{link}">{html.escape(path.name)}</a>'
        entries[relative] = (file_link, content or f"<p>{html.escape(status)}</p>")
        inventory.append(f"<tr><td>{source}</td><td>{path.stat().st_size:,}</td>"
                         f"<td>{html.escape(status)}</td></tr>")
    summaries, combined_count, summary_failures = build_summaries(results, output)
    plots += combined_count
    failures += summary_failures
    sections, navigation = render_tree(entries, summaries=summaries)
    document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GHIST survey results</title><style>
body{font-family:system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem}
img{max-width:100%;height:auto}article{border-top:1px solid #ccc;padding:1rem 0}
h2,h3,h4,h5,h6,[role=heading],td{overflow-wrap:anywhere}
[role=heading]{font-weight:bold;margin:1rem 0}section{margin-left:1rem}
table{border-collapse:collapse;width:100%}
.estimates{font-size:.8rem}table{display:block;overflow-x:auto}
td,th{border:1px solid #ccc;padding:.4rem;text-align:left}pre{overflow:auto}
nav ul{padding-left:1.25rem}nav>ul{columns:3 18rem}.error{color:#b00020}
</style></head><body><h1>GHIST survey results</h1>
<p>Exploratory outputs, not scientifically validated estimates. Challenge overviews
compose compatible outputs with legends and separate panels for different quantities.
Demographic time is converted to generations using published mutation assumptions
and PSMC theta/bin size; within-population Ne curves are compared, not cross rates.
Files lacking scaling metadata retain explicitly labeled native units.
RCCR uses the generations provided by the pipeline. Target estimates come from the
underlying numbers, not image pixels. Unidentifiable targets remain unavailable;
bin/window brackets are not confidence intervals. BED length sums count overlaps.</p>
<p>All plots and metadata are embedded for offline viewing. Source links require
the results tree. Models, binary intermediates, logs and other unsupported formats
are linked in their original directories and inventoried below, not executed or
deserialized. Headings follow category, challenge, tool, run and nested directories.</p>
"""
    document += f"<p>{len(files)} files; {plots} plots; {failures} errors.</p>"
    if summaries:
        estimates_link = quote(Path(os.path.relpath(
            results / "challenge_estimates.json", output.parent)).as_posix(), safe="/")
        document += f'<p><a href="{estimates_link}">Download target estimates and pair/region assignments (JSON)</a></p>'
    document += '<nav aria-label="Results directories"><ul>' + navigation + "</ul></nav>" + sections
    document += "<h2>Complete file inventory</h2><table><thead><tr><th>File</th>"
    document += "<th>Bytes</th><th>Status</th></tr></thead><tbody>"
    document += "".join(inventory) + "</tbody></table></body></html>"
    output.write_text(document, encoding="utf-8")
    return output, plots, failures
