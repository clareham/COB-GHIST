#!/usr/bin/env python3
"""Build a placeholder PSMC-fa file from a single haploid consensus FASTA.

CAVEAT (see modules/genetic_history/psmc/main.nf and the README): the
GHIST `*_low_coverage` challenges hand out exactly one haploid consensus
sequence per scenario, with no second haplotype and no heterozygosity
signal at all. Classic PSMC infers coalescence time from heterozygous
sites in a single *diploid* genome, so there is no real het/hom signal to
encode here - every callable window necessarily becomes 'T' (no het).

This script exists so the pipeline is still mechanically runnable end to
end (every 100bp window with only ACGT bases -> 'T', ambiguous/N bases ->
'N'), producing a flat, uninformative PSMC curve. Treat any output from
this path as a structural placeholder, not a validated inference - a real
analysis needs either two haplotypes per individual or actual low-coverage
reads + genotype likelihoods (not provided in this data drop).
"""
import argparse

VALID = set("ACGTacgt")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fasta")
    ap.add_argument("out_psmcfa")
    ap.add_argument("--window", type=int, default=100)
    args = ap.parse_args()

    with open(args.fasta) as fh, open(args.out_psmcfa, "w") as out:
        contig, seq = None, []

        def flush():
            if contig is None:
                return
            out.write(f">{contig}\n")
            s = "".join(seq)
            calls = []
            for i in range(0, len(s), args.window):
                window = s[i:i + args.window]
                calls.append("T" if any(c in VALID for c in window) else "N")
            call_str = "".join(calls)
            for i in range(0, len(call_str), 60):
                out.write(call_str[i:i + 60] + "\n")

        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                flush()
                contig, seq = line[1:].split()[0], []
            else:
                seq.append(line)
        flush()


if __name__ == "__main__":
    main()
