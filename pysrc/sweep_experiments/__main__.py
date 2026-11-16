"""Command-line entry point; heavyweight dependencies are loaded only when used."""

import argparse
import logging
from pathlib import Path

from . import benchmark


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description="Exploratory SweepFinder2/Evo 2 experiments")
    sub = root.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run a JSON experiment manifest")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True, help="New directory; never overwritten")
    run.add_argument("--configurations", default="physical,constant",
                     help="Comma-separated physical,constant,map,bmap")
    run.add_argument("--sweepfinder", default="SweepFinder2")
    run.add_argument("--bcftools", default="bcftools")
    run.add_argument("--grid", type=int, default=1000)
    run.add_argument("--constant-rate", type=float, default=1e-8)
    run.add_argument("--unfolded", action="store_true", help="Require valid INFO/AA")
    run.add_argument("--evo", action="store_true")
    run.add_argument("--seqrepo", type=Path)
    run.add_argument("--model", default="auto",
                     choices=["auto", "evo2_40b", "evo2_7b", "evo2_1b_base"])
    run.add_argument("--allow-rocm", action="store_true")
    run.add_argument("--context", type=int, default=8192)
    run.add_argument("--offsets", default="1024,4096,7168",
                     help="0-based target offsets in oriented context; must be >=1")
    run.add_argument("--evo-stride", type=int, default=1, help="Score every Nth eligible SNP")
    run.add_argument("--window-size", type=int, default=100000,
                     help="Physical bins for exploratory Evo aggregation")
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
    else:
        benchmark.evaluate(args.summary, args.output, args.alpha)


if __name__ == "__main__":
    main()
