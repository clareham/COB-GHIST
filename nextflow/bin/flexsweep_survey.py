#!/usr/bin/env python3
"""Exploratory FlexSweep 2.0 training and scan, with checked real artifacts.

Source contract: immutable PyPI flexsweep-2.0.tar.gz, SHA256
9b255c48ec7f4d6c6cdadf6db0807e02fa69a8aca29b4c762cc6831e504b8f3b.
https://pypi.org/project/flexsweep/2.0/ (upstream has no Git tags/releases).
Inspected files: flexsweep/{main,simulate_discoal,data,fv,cnn}.py.

Actual CLI names are simulator, fvs-discoal, fvs-vcf, cnn; options use
underscores. The 2.0 cnn Click callback omits its declared center/windows/step
parameters, so use its CNN Python API, not the nonfunctional cnn CLI. The
simulator also passes unsupported compression='gzip' to Polars.write_csv.
Instead, the local simulator driver samples the documented parameter mixture
using Simulator.random_distribution(), writes its own parameter table through
gzip.open(), and calls the public Simulator.simulate_batch() API. Training
uses unmodified CNN.train(): at most 1000 epochs, with upstream validation
accuracy early stopping (patience 5). No framework methods or package files
are modified. The packaged discoal is copied/chmodded locally; no runtime
download or build is permitted.

Input: phased, complete diploid GTs; retain only biallelic SNPs in the BED mask.
REF is assumed ancestral, explicitly configurable only by accepting/rejecting
that assumption. No ancestral inference, imputation, phasing or BGS training.
Contigs require header lengths and are reversibly aliased to chr1, chr2, ...
because upstream infers contigs from numeric filenames. Normalization follows
2.0: neutral simulation bins for training, empirical bins for VCF features.
Mask coverage filters entire predicted windows; boundary/partial windows are
excluded. This is not a statistically calibrated discovery/submission model.
"""

import argparse
import csv
import gzip
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from callable_regions import callable_bases_in_range, load_mask


DEFAULTS = {
    "threads": 4,
    "simulations_per_class": 1000,
    "effective_population_size": 10000,
    "mutation_rate": "5e-9,2e-8",
    "recombination_rate": "1e-10,1e-7,1e-8",
    "selection_coefficients": "0.01,0.05",
    "sweep_time_generations": "0,5000",
    "soft_start_frequency": "0,0.1",
    "incomplete_end_frequency": "0.5,1",
    "locus_length": 1200000,
    "window_size": 100000,
    "statistic_step": 100000,
    "scan_step": 100000,
    "min_accessible_fraction": 0.8,
    "sweep_probability": 0.9,
    "ancestral_mode": "reference",
    "keep_simulations": False,
}

CNN_STATS = (
    "dind,dist_kurtosis,dist_skew,dist_var,h1,h12,h2_h1,haf,hapdaf_o,"
    "hapdaf_s,high_freq,ihs,isafe,k_counts,low_freq,max_fda,nsl,omega_max,"
    "pi,s_ratio,tajima_d,theta_h,theta_w,zns"
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_options(path):
    supplied = json.loads(Path(path).read_text())
    require(isinstance(supplied, dict), "options must be a JSON object")
    require(not (supplied.keys() - DEFAULTS.keys()),
            f"Unknown FlexSweep options: {supplied.keys() - DEFAULTS.keys()}")
    options = DEFAULTS | supplied
    for name in ("threads", "simulations_per_class",
                 "effective_population_size", "locus_length", "window_size",
                 "statistic_step", "scan_step"):
        require(type(options[name]) is int and options[name] > 0,
                f"{name} must be a positive integer")
    require(options["simulations_per_class"] >= 200
            and options["simulations_per_class"] % 20 == 0,
            "simulations_per_class must be >=200 and a multiple of 20 "
            "(upstream 50:50 hard/soft, 10% complete sweep allocation)")
    require(options["locus_length"] % options["statistic_step"] == 0,
            "locus_length must be divisible by statistic_step")
    require(options["window_size"] <= options["locus_length"],
            "window_size must not exceed locus_length")
    require(options["locus_length"] // options["statistic_step"] >= 2,
            "At least two statistic centers are required by the CNN pooling layers")
    require(options["ancestral_mode"] == "reference",
            "Only ancestral_mode=reference is supported: polarize the VCF first "
            "if REF cannot be assumed ancestral")
    require(type(options["keep_simulations"]) is bool, "keep_simulations must be boolean")
    for name in ("min_accessible_fraction", "sweep_probability"):
        require(isinstance(options[name], (int, float)) and 0 < options[name] <= 1,
                f"{name} must be in (0,1]")
    for name, count in (("mutation_rate", 2), ("recombination_rate", 3),
                        ("selection_coefficients", 2), ("sweep_time_generations", 2),
                        ("soft_start_frequency", 2), ("incomplete_end_frequency", 2)):
        values = [float(v) for v in options[name].split(",")]
        require(len(values) == count and all(math.isfinite(v) for v in values),
                f"{name} must contain {count} finite comma-separated numbers")
        require(0 <= values[0] <= values[1], f"Invalid {name} bounds")
        if name in ("mutation_rate", "recombination_rate", "selection_coefficients"):
            require(values[0] > 0, f"{name} must be positive")
        if "frequency" in name:
            require(values[1] <= 1, f"{name} must not exceed 1")
        if count == 3:
            require(values[0] < values[1] and values[2] > 0,
                    "recombination exponential bounds must differ and mean must be positive")
    return options


def coverage(intervals, start, end):
    return callable_bases_in_range(intervals, start, end) / (end - start)


def prepare_vcfs(vcf, mask, directory):
    opener = gzip.open if str(vcf).endswith(".gz") else open
    headers, contigs, handles, counts = [], {}, {}, {}
    samples = None
    last_positions = {}
    try:
        with opener(vcf, "rt") as source:
            for line in source:
                if line.startswith("##contig="):
                    name = re.search(r"(?:<|,)ID=([^,>]+)", line)
                    length = re.search(r"(?:<|,)length=(\d+)", line, re.I)
                    require(name is not None and length is not None,
                            "Every VCF contig header requires ID and length")
                    original = name.group(1)
                    require(original not in contigs, f"Duplicate contig header: {original}")
                    contigs[original] = (f"chr{len(contigs) + 1}", int(length.group(1)))
                elif line.startswith("#"):
                    headers.append(line)
                    if line.startswith("#CHROM"):
                        samples = line.rstrip().split("\t")[9:]
                        require(len(samples) >= 2, "VCF requires at least two diploid samples")
                else:
                    require(samples is not None, "Missing #CHROM VCF header")
                    fields = line.rstrip().split("\t")
                    chrom, position, _, ref, alt = fields[:5]
                    require(chrom in contigs, f"Missing length header for {chrom}")
                    pos = int(position)
                    alias, length = contigs[chrom]
                    require(1 <= pos <= length, f"Variant outside contig length: {chrom}:{pos}")
                    require(pos > last_positions.get(chrom, 0),
                            f"VCF must be sorted with unique positions per contig: {chrom}:{pos}")
                    last_positions[chrom] = pos
                    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                        continue
                    if chrom not in mask or coverage(mask[chrom], pos - 1, pos) == 0:
                        continue
                    require(len(fields) == 9 + len(samples), "Inconsistent VCF sample count")
                    formats = fields[8].split(":")
                    require("GT" in formats, "GT FORMAT is required")
                    gt_index = formats.index("GT")
                    for sample in fields[9:]:
                        values = sample.split(":")
                        require(gt_index < len(values)
                                and re.fullmatch(r"[01]\|[01]", values[gt_index]) is not None,
                                f"Complete phased diploid biallelic GT required: {chrom}:{pos}")
                    if chrom not in handles:
                        target = directory / f"{alias}.vcf"
                        handle = target.open("w")
                        handles[chrom] = handle
                        for header in headers:
                            if header.startswith("#CHROM"):
                                handle.write(f"##contig=<ID={alias},length={length}>\n")
                            handle.write(header)
                        counts[chrom] = 0
                    fields[0] = alias
                    handles[chrom].write("\t".join(fields) + "\n")
                    counts[chrom] += 1
    finally:
        for handle in handles.values():
            handle.close()
    require(bool(counts), "No accessible biallelic SNPs remain")
    for chrom in counts:
        alias, _ = contigs[chrom]
        raw = directory / f"{alias}.vcf"
        compressed = raw.with_suffix(".vcf.gz")
        with compressed.open("wb") as output:
            subprocess.run(["bgzip", "-c", str(raw)], stdout=output, check=True)
        subprocess.run(["tabix", "-p", "vcf", str(compressed)], check=True)
        raw.unlink()
    return samples, contigs, counts


def simulate(options, haplotypes):
    """Own parameter sampling/serialization; leave the upstream simulator unchanged."""
    import numpy as np
    import polars as pl
    from flexsweep import Simulator

    def pair(name):
        return list(map(float, options[name].split(",")))

    mu_min, mu_max = pair("mutation_rate")
    r_min, r_max, r_mean = pair("recombination_rate")
    simulator = Simulator(
        sample_size=haplotypes, demes="demography.yaml", output_folder="simulations",
        mutation_rate={"dist": "uniform", "min": mu_min, "max": mu_max},
        recombination_rate={"dist": "exponential", "min": r_min, "max": r_max, "mean": r_mean},
        locus_length=options["locus_length"], discoal_path=str(Path("discoal").resolve()),
        num_simulations=options["simulations_per_class"], nthreads=options["threads"],
        saf=pair("soft_start_frequency"), eaf=pair("incomplete_end_frequency"),
        s=pair("selection_coefficients"), time=pair("sweep_time_generations"),
    )
    simulator.check_inputs()
    count = options["simulations_per_class"]
    hard = count // 2
    complete = hard // 10
    incomplete = hard - complete
    initial = np.concatenate([
        np.zeros(hard),
        np.random.uniform(*simulator.f_i, complete),
        np.random.uniform(*simulator.f_i, incomplete),
    ])
    final = np.concatenate([
        np.ones(complete), np.random.uniform(*simulator.f_t, incomplete),
        np.ones(complete), np.random.uniform(*simulator.f_t, incomplete),
    ])
    tables = []
    for model in ("neutral", "sweep"):
        mu, recombination = simulator.random_distribution(count)
        selected = model == "sweep"
        tables.append(pl.DataFrame({
            "iter": np.arange(1, count + 1), "mu": mu, "r": recombination,
            "eaf": final if selected else np.zeros(count),
            "saf": initial if selected else np.zeros(count),
            "s": np.random.uniform(*simulator.s, count) if selected else np.zeros(count),
            "t": np.random.uniform(*simulator.time, count) if selected else np.zeros(count),
            "pos": np.full(count, 0.5) if selected else np.zeros(count),
            "model": [model] * count,
        }))
    simulator.parameters = pl.concat(tables, how="vertical")
    with gzip.open("simulations/params.txt.gz", "wb") as handle:
        simulator.parameters.write_csv(handle)
    simulator.simulate_batch()


def child(stage, arguments):
    require(importlib.metadata.version("flexsweep") == "2.0",
            "This adapter requires exactly FlexSweep 2.0")
    import flexsweep as fs

    os.environ["JOBLIB_TEMP_FOLDER"] = str(Path("scratch").resolve())
    if stage == "cli":
        from flexsweep.main import cli
        cli.main(args=arguments, standalone_mode=False)
    elif stage == "simulate":
        simulate(json.loads(Path(arguments[0]).read_text()), int(arguments[1]))
    elif stage in ("train", "predict"):
        options = json.loads(Path(arguments[0]).read_text())
        import numpy as np
        import polars as pl

        features = "simulations/fvs.parquet" if stage == "train" else "empirical/fvs_survey.parquet"
        table = pl.read_parquet(features)
        require(table.height > 0, f"No feature rows in {features}")
        feature_columns = [name for name in table.columns
                           if re.fullmatch(r".+_\d+_\d+", name)]
        centers = range(options["statistic_step"] // 2, options["locus_length"],
                        options["statistic_step"])
        expected_columns = {
            f"{stat}_{options['window_size']}_{center}"
            for stat in CNN_STATS.split(",") for center in centers
        }
        require(set(feature_columns) == expected_columns,
                f"Feature geometry does not match CNN windows/centers in {features}: "
                f"missing={sorted(expected_columns - set(feature_columns))}, "
                f"extra={sorted(set(feature_columns) - expected_columns)}")
        require(np.isfinite(table.select(feature_columns).to_numpy()).all(),
                f"Nonfinite or missing statistic features in {features}")
        if stage == "train":
            classes = set(table["model"].to_list())
            require("neutral" in classes and any(name != "neutral" for name in classes),
                    "Training features must contain both neutral and sweep classes")
            require(table.height == 2 * options["simulations_per_class"],
                    "Upstream feature extraction dropped training simulations; "
                    "inspect logs rather than silently accepting incomplete training")

        runner = fs.CNN(
            train_data="simulations/fvs.parquet" if stage == "train" else None,
            predict_data="empirical/fvs_survey.parquet" if stage == "predict" else None,
            model="training/model.keras" if stage == "predict" else None,
            output_folder="training" if stage == "train" else "prediction",
            center=[0, options["locus_length"]],
            windows=[options["window_size"]],
            step=options["statistic_step"],
        )
        if stage == "train":
            runner.train()
        else:
            prediction = runner.predict()
            require(prediction.height > 0, "CNN produced no predictions")
            prediction.write_csv("predictions.raw.csv")
    else:
        raise ValueError(f"Unknown child stage: {stage}")


def artifact(path):
    require(Path(path).is_file() and Path(path).stat().st_size > 0,
            f"Missing/empty expected FlexSweep artifact: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--options", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    options = read_options(args.options)
    require(args.threads > 0, "--threads must be positive")
    options["threads"] = args.threads
    vcf, mask_path = Path(args.vcf).resolve(), Path(args.mask).resolve()
    helper = str(Path(__file__).resolve())
    root = Path(args.output)
    root.mkdir()
    os.chdir(root)
    for name in ("scratch", "empirical", "training", "prediction"):
        Path(name).mkdir()
    os.environ["TMPDIR"] = str(Path("scratch").resolve())
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    require(importlib.metadata.version("flexsweep") == "2.0", "Requires flexsweep=2.0")
    package = Path(importlib.metadata.distribution("flexsweep").locate_file("flexsweep"))
    binary = package / "data" / "discoal"
    require(binary.is_file(), "FlexSweep package lacks bundled discoal; no runtime build allowed")
    shutil.copyfile(binary, "discoal")
    Path("discoal").chmod(0o755)
    print("WARNING: FlexSweep trains from scratch; CPU training is expensive. "
          "Survey-sized simulations and constant-Ne demography are exploratory; "
          "BGS and fitted demographic histories are not modeled. "
          "Unmodified upstream training allows up to 1000 epochs with early stopping "
          "(validation accuracy, patience 5).", file=sys.stderr)
    mask = load_mask(mask_path)
    samples, contigs, counts = prepare_vcfs(vcf, mask, Path("empirical"))
    Path("options.json").write_text(json.dumps(options, indent=2) + "\n")
    Path("demography.yaml").write_text(
        "time_units: generations\ndemes:\n- name: survey\n  epochs:\n"
        f"  - start_size: {options['effective_population_size']}\n")

    def run(stage, *arguments):
        print("FlexSweep stage:", stage, *arguments, flush=True)
        subprocess.run([sys.executable, helper, "__child", stage, *map(str, arguments)],
                       check=True)

    run("simulate", "options.json", 2 * len(samples))
    artifact("simulations/params.txt.gz")
    for kind in ("neutral", "sweep"):
        paths = list(Path("simulations", kind).glob("*.ms.gz"))
        require(len(paths) == options["simulations_per_class"],
                f"Wrong number of {kind} simulations")
        for path in paths:
            with gzip.open(path, "rt") as handle:
                require(any(line.startswith("segsites:") for line in handle),
                        f"Failed/malformed discoal simulation: {path}")
    geometry = ["--locus_length", options["locus_length"],
                "--windows", options["window_size"], "--step", options["statistic_step"],
                "--stats", CNN_STATS, "--nthreads", args.threads]
    run("cli", "fvs-discoal", "--simulations_path", "simulations", *geometry)
    artifact("simulations/fvs.parquet")
    artifact("simulations/neutral_bins.pickle")
    run("train", "options.json")
    artifact("training/model.keras")
    run("cli", "fvs-vcf", "--vcf_path", "empirical",
        "--suffix", "survey", "--step_vcf", options["scan_step"], *geometry)
    artifact("empirical/fvs_survey.parquet")
    artifact("empirical/empirical_bins_survey.pickle")
    run("predict", "options.json")
    artifact("predictions.raw.csv")
    aliases = {alias: (chrom, length) for chrom, (alias, length) in contigs.items()}
    retained = 0
    with open("predictions.raw.csv") as source, open("predictions.tsv", "w") as output, \
            open("sweeps.bed", "w") as hits:
        reader = csv.DictReader(source)
        columns = ["chrom", "start", "end", "prob_sweep", "prob_neutral",
                   "predicted_model", "accessible_fraction"]
        writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in reader:
            require(row["chr"] in aliases, f"Unexpected predicted contig: {row['chr']}")
            chrom, length = aliases[row["chr"]]
            start, end = int(row["start"]) - 1, int(row["end"])
            if start < 0 or end > length or end - start != options["locus_length"]:
                continue
            fraction = coverage(mask[chrom], start, end)
            if fraction < options["min_accessible_fraction"]:
                continue
            probability, neutral = float(row["prob_sweep"]), float(row["prob_neutral"])
            require(math.isfinite(probability) and 0 <= probability <= 1
                    and math.isfinite(neutral) and 0 <= neutral <= 1
                    and math.isclose(probability + neutral, 1.0, abs_tol=1e-5),
                    "CNN returned invalid probabilities")
            writer.writerow(dict(zip(columns, [chrom, start, end, probability, neutral,
                                               row["predicted_model"], fraction])))
            retained += 1
            if probability >= options["sweep_probability"]:
                hits.write(f"{chrom}\t{start}\t{end}\tFlexSweep\t{round(1000 * probability)}\n")
    require(retained > 0, "No full-length predictions pass accessible-mask coverage; "
            "review locus_length/min_accessible_fraction, not an empty successful scan")
    Path("manifest.json").write_text(json.dumps({
        "flexsweep_version": "2.0", "options": options, "samples": len(samples),
        "release_sdist_sha256": "9b255c48ec7f4d6c6cdadf6db0807e02fa69a8aca29b4c762cc6831e504b8f3b",
        "compatibility_adapters": ["local parameter sampling/gzip serialization with public Simulator API",
                                   "CNN Python API instead of broken Click callback",
                                   "executable local bundled discoal"],
        "training_schedule": "unmodified CNN.train(): maximum 1000 epochs; "
                             "validation accuracy early stopping, patience 5",
        "haplotypes": 2 * len(samples), "retained_predictions": retained,
        "contig_aliases": contigs, "accessible_biallelic_snps": counts,
        "coordinates": "predictions.tsv and sweeps.bed: zero-based half-open",
        "bed_score": "integer round(1000 * prob_sweep); intervals are whole scan windows",
        "training_mix": "50% neutral; 50% sweep, half hard/half soft; 10% complete",
        "normalization": "simulation neutral bins for training; empirical bins for VCF",
        "statistics": CNN_STATS.split(","),
        "limitations": "Exploratory constant-Ne/no-BGS model; REF assumed ancestral; "
                       "masked variants omitted but training simulations are unmasked; "
                       "no deterministic parallel simulation seed; no calibrated FDR",
    }, indent=2) + "\n")
    if not options["keep_simulations"]:
        for kind in ("neutral", "sweep"):
            shutil.rmtree(Path("simulations", kind))
    for path in Path("empirical").glob("*.vcf.gz*"):
        path.unlink()
    shutil.rmtree("scratch")
    Path("discoal").unlink()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "__child":
        child(sys.argv[2], sys.argv[3:])
    else:
        main()
