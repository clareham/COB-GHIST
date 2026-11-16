"""Generate a survey plotting report as a module or standalone script."""

import argparse
from pathlib import Path

if __package__:
    from .plotting.report import build_report
else:
    from plotting.report import build_report


DEFAULT_RESULTS = Path(__file__).resolve().parents[1] / "nextflow" / "results"


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
