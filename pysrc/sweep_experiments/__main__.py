"""Command-line entry point; heavyweight dependencies are loaded only when used."""

import argparse
import logging
from pathlib import Path

from . import benchmark
from .ntv3 import MODEL


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description="Exploratory SweepFinder2 and genomic LM experiments")
    sub = root.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run a JSON experiment manifest")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True,
                     help="Output directory; existing generated files are overwritten")
    run.add_argument("--configurations", default="physical,constant",
                     help="Comma-separated physical,constant,map,bmap")
    run.add_argument("--sweepfinder", default="SweepFinder2")
    run.add_argument("--bcftools", default="bcftools")
    run.add_argument("--grid", type=int, default=1000)
    run.add_argument("--constant-rate", type=float, default=1e-8)
    run.add_argument("--unfolded", action="store_true", help="Require valid INFO/AA")
    glm = run.add_mutually_exclusive_group()
    glm.add_argument("--ntv3", action="store_true", help="Enable bidirectional masked-base scoring")
    glm.add_argument("--hyenadna", action="store_true", help="Enable HyenaDNA large causal next-base scoring")
    run.add_argument("--seqrepo", type=Path)
    run.add_argument("--model", choices=[MODEL], help="Explicit NTv3 checkpoint (legacy option)")
    run.add_argument("--model-revision",
                     help="Override the selected model's pinned Hub revision (weights, tokenizer, and code)")
    run.add_argument("--context", type=int, default=8192,
                     help="Window length in bases (default: 8192; NTv3 requires a multiple of 128)")
    run.add_argument("--offsets", default="1024,4096,7168",
                     help="0-based target offsets in oriented context; HyenaDNA requires offsets >= 1")
    run.add_argument("--glm-stride", "--ntv3-stride", type=int, default=1,
                     help="Score every Nth eligible SNP with the selected model")
    run.add_argument("--window-size", type=int, default=100000,
                     help="Physical bins for exploratory genomic LM aggregation")
    run.add_argument("--alpha", type=float, default=.05,
                     help="Per-configuration, region-maximum empirical p-value cutoff for sweep BEDs")
    run.add_argument("--null-summary", type=Path,
                     help="Completed independent run's summary.tsv containing matched calibration controls")
    simulation = sub.add_parser("simulate", help="Generate reference-anchored neutral msprime controls")
    simulation.add_argument("--seqrepo", type=Path, required=True)
    simulation.add_argument("--dataset", choices=["testing", "final"], default="testing")
    simulation.add_argument("--start", type=int, required=True, help="0-based inclusive")
    simulation.add_argument("--end", type=int, required=True, help="0-based exclusive")
    simulation.add_argument("--replicates", type=int, default=100)
    simulation.add_argument("--samples", type=int, default=25, help="Diploid individuals")
    simulation.add_argument("--ne", type=float, default=10000)
    simulation.add_argument("--mutation-rate", type=float, default=1.2e-8)
    simulation.add_argument("--constant-rate", type=float, default=1e-8)
    simulation.add_argument("--recombination-map", type=Path)
    simulation.add_argument("--seed", type=int, default=1729)
    simulation.add_argument("--output", type=Path, required=True)
    bvalues = sub.add_parser("estimate-b", help="Estimate B and bootstrap intervals from replicate diversity")
    bvalues.add_argument("--diversity", type=Path, required=True)
    bvalues.add_argument("--bootstrap", type=int, default=1000)
    bvalues.add_argument("--seed", type=int, default=1729)
    bvalues.add_argument("--output", type=Path, required=True)
    evaluation = sub.add_parser("evaluate", help="Calibrate maxima using held-out no-sweep replicates")
    evaluation.add_argument("--summary", type=Path, required=True)
    evaluation.add_argument("--alpha", type=float, default=.05)
    evaluation.add_argument("--output", type=Path, required=True)
    regions = sub.add_parser("export-regions", help="Re-export BED evidence and sweep calls from a completed run")
    regions.add_argument("--run-dir", type=Path, required=True)
    regions.add_argument("--alpha", type=float, default=.05)
    regions.add_argument("--null-summary", type=Path)
    return root


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.command == "run":
        benchmark.run(args)
    elif args.command == "simulate":
        from .inputs import Reference
        from .simulation import simulate

        simulate(args.output, Reference(args.seqrepo, args.dataset), args.start, args.end,
                 args.replicates, args.samples, args.ne, args.mutation_rate,
                 args.constant_rate, args.seed, args.recombination_map)
    elif args.command == "estimate-b":
        from .simulation import estimate_b

        estimate_b(args.diversity, args.output, args.bootstrap, args.seed)
    elif args.command == "evaluate":
        benchmark.evaluate(args.summary, args.output, args.alpha)
    else:
        from .regions import export_regions

        export_regions(args.run_dir, args.alpha, args.null_summary, require_complete=True)


if __name__ == "__main__":
    main()
