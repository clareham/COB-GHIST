"""Discover result files, dispatch challenge plots, and render their directory tree."""

import base64
import csv
import html
import io
import json
import os
import sys
from pathlib import Path
from urllib.parse import quote

from . import common, genetic_history, relatedness, sweep_detection
from .common import Plotter, plt
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
                level: int = 2) -> tuple[str, str]:
    """Return nested result sections and a matching directory navigation tree."""
    content = []
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
        nested_content, nested_navigation = render_tree(entries, child, level + 1)
        content.append(f'<section id="{anchor}" data-path="'
                       f'{html.escape(child.as_posix(), quote=True)}">'
                       + heading(level, title) + nested_content + "</section>")
        navigation.append(f'<li><a href="#{anchor}">{title}</a>'
                          + (f"<ul>{nested_navigation}</ul>" if nested_navigation else "") + "</li>")
    return "".join(content), "".join(navigation)


def build_report(results: Path, output: Path | None = None) -> tuple[Path, int, int]:
    results = results.resolve()
    if not results.is_dir():
        raise ValueError(f"Results directory does not exist: {results}")
    output = output.resolve() if output is not None else results / "plots.html"
    files = sorted(path for path in results.rglob("*") if path.is_file() and path.resolve() != output)
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
    sections, navigation = render_tree(entries)
    document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GHIST survey results</title><style>
body{font-family:system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem}
img{max-width:100%;height:auto}article{border-top:1px solid #ccc;padding:1rem 0}
h2,h3,h4,h5,h6,[role=heading],td{overflow-wrap:anywhere}
[role=heading]{font-weight:bold;margin:1rem 0}section{margin-left:1rem}
table{border-collapse:collapse;width:100%}
td,th{border:1px solid #ccc;padding:.4rem;text-align:left}pre{overflow:auto}
nav ul{padding-left:1.25rem}nav>ul{columns:3 18rem}.error{color:#b00020}
</style></head><body><h1>GHIST survey results</h1>
<p>Exploratory outputs, not scientifically validated estimates. Missing tools are
not imputed. PSMC uses native relative units; MSMC2 uses mutation-scaled time and
rates (cross-population rates are not effective population sizes). RCCR uses the
generations provided by the pipeline. Ancestry and relatedness plots do not
assign population or relationship labels. BED length sums count overlaps.</p>
<p>All plots and metadata are embedded for offline viewing. Source links require
the results tree. Models, binary intermediates, logs and other unsupported formats
are linked in their original directories and inventoried below, not executed or
deserialized. Headings follow category, challenge, tool, run and nested directories.</p>
"""
    document += f"<p>{len(files)} files; {plots} plots; {failures} errors.</p>"
    document += '<nav aria-label="Results directories"><ul>' + navigation + "</ul></nav>" + sections
    document += "<h2>Complete file inventory</h2><table><thead><tr><th>File</th>"
    document += "<th>Bytes</th><th>Status</th></tr></thead><tbody>"
    document += "".join(inventory) + "</tbody></table></body></html>"
    output.write_text(document, encoding="utf-8")
    return output, plots, failures
