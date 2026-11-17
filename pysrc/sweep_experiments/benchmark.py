"""Manifest-driven scans and held-out, per-configuration null calibration."""

import csv
import hashlib
import json
import logging
import math
from pathlib import Path
import re
import shutil
import subprocess

from .inputs import (Reference, check_mask, read_frequency, variants, write_bvalues, write_frequency,
                     write_recombination)

LOG = logging.getLogger(__name__)


def file_metadata(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "sha256": digest.hexdigest()}


def run_command(command: list[str], log: Path) -> None:
    LOG.info("Running %s", " ".join(command))
    with log.open("w") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


def scan_peak(path: Path) -> tuple[float, float]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: empty SweepFinder2 scan")
    points = [(float(row["LR"]), float(row["location"])) for row in rows]
    if any(not math.isfinite(v) for point in points for v in point):
        raise ValueError(f"{path}: nonfinite scan")
    return max(points, key=lambda item: item[0])


def evo_peak(path: Path, window_size: int) -> tuple[float, float]:
    # Average strands/offsets per SNP first so missing contexts do not reweight SNPs.
    by_site: dict[int, list[float]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] == "ok":
                by_site.setdefault(int(row["position"]), []).append(abs(float(row["delta_logp"])))
    bins: dict[int, list[float]] = {}
    for position, scores in by_site.items():
        bins.setdefault((position - 1) // window_size, []).append(sum(scores) / len(scores))
    if not bins:
        raise ValueError("No valid Evo scores to aggregate")
    index = max(bins, key=lambda key: sum(bins[key]) / len(bins[key]))
    return sum(bins[index]) / len(bins[index]), index * window_size + (window_size + 1) / 2


def run(args) -> Path:
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    base = manifest_path.parent
    cases = manifest["cases"]
    if not cases:
        raise ValueError("Manifest has no cases")
    ids = [case["id"] for case in cases]
    if len(set(ids)) != len(ids) or any(not re.fullmatch(r"[A-Za-z0-9_-]+", i) for i in ids):
        raise ValueError("Case IDs must be unique, safe directory names")
    evaluation_vcfs = {(base / case["vcf"]).resolve() for case in cases
                       if case.get("split") == "evaluation"}
    calibration_vcfs = {(base / case["vcf"]).resolve() for case in cases
                        if case.get("split") == "calibration"}
    if evaluation_vcfs & calibration_vcfs:
        raise ValueError("Calibration and evaluation must not reuse a VCF")
    if any((base / case["background_vcf"]).resolve() in evaluation_vcfs
           for case in cases if case.get("background_vcf")):
        raise ValueError("An evaluation VCF cannot be an external background VCF")
    configurations = args.configurations.split(",")
    allowed = {"physical", "constant", "map", "bmap"}
    if not configurations or set(configurations) - allowed or len(set(configurations)) != len(configurations):
        raise ValueError(f"Configurations must be unique members of {sorted(allowed)}")
    if args.grid < 2 or args.evo_stride < 1 or args.window_size < 1:
        raise ValueError("Need grid >=2, positive Evo stride and window size")
    if args.evo and args.seqrepo is None:
        raise ValueError("--evo requires a locally populated --seqrepo")
    if args.evo:
        from .evo import windows

        offsets = [int(value) for value in args.offsets.split(",")]
        windows(args.context + 1, 3 * args.context, args.context, offsets)
    executable = shutil.which(args.sweepfinder)
    if executable is None:
        raise FileNotFoundError(f"SweepFinder2 executable not found: {args.sweepfinder}")
    bcftools_version = subprocess.run(
        [args.bcftools, "--version"], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    output = args.output.resolve()
    if output.exists():
        LOG.warning("Reusing output directory %s; generated files will be overwritten", output)
    output.mkdir(parents=True, exist_ok=True)
    provenance = {"manifest": file_metadata(manifest_path), "arguments": {
        key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()
        if key != "func"}, "inputs": {}, "commands": [], "status": "running",
        "sweepfinder_binary": file_metadata(Path(executable)), "bcftools_version": bcftools_version}
    metadata_file = output / "run.json"

    def save() -> None:
        metadata_file.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")

    save()
    predictor = None
    summary = output / "summary.tsv"
    with summary.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[
            "case", "scenario", "split", "configuration", "score", "peak_position",
            "sweep_position", "localization_error_bp"])
        writer.writeheader()
        for case in cases:
            name = case["id"]
            if case.get("split") not in ("calibration", "evaluation", "exploratory"):
                raise ValueError(f"{name}: split must be calibration/evaluation/exploratory")
            scenario = case.get("scenario", "unknown")
            if scenario not in ("neutral", "bgs", "sweep", "bgs_sweep", "unknown"):
                raise ValueError(f"{name}: invalid scenario")
            truth = case.get("sweep_position")
            if scenario in ("neutral", "bgs") and truth is not None:
                raise ValueError(f"{name}: a no-sweep control cannot have sweep_position")
            if truth is not None and (not math.isfinite(float(truth)) or float(truth) < 1):
                raise ValueError(f"{name}: sweep_position must be a positive coordinate")
            directory = output / name
            directory.mkdir(exist_ok=True)

            def source(key: str) -> Path:
                path = (base / case[key]).resolve()
                provenance["inputs"][f"{name}:{key}"] = file_metadata(path)
                return path

            vcf = source("vcf")
            mask = source("mask")
            sites = list(variants(vcf, not args.unfolded, mask, args.bcftools))
            freq = directory / "sites.freq"
            if not sites:
                raise ValueError(f"{name}: no eligible polymorphic SNPs")
            frequency_sites = sites
            if case.get("sf2_frequency"):
                frequency_source = source("sf2_frequency")
                frequency_sites = read_frequency(
                    frequency_source, sites[0].chrom, not args.unfolded)
                check_mask(frequency_sites, mask)
                by_position = {site.position: site for site in frequency_sites}
                if any(site.position not in by_position
                       or (by_position[site.position].x, by_position[site.position].n)
                       != (site.x, site.n) for site in sites):
                    raise ValueError("sf2_frequency must include matching counts for every eligible VCF SNP")
                shutil.copyfile(frequency_source, freq)
            else:
                write_frequency(freq, sites)
            if case.get("background_spectrum") and case.get("background_vcf"):
                raise ValueError("Choose background_spectrum or background_vcf, not both")
            spectrum = directory / "background.spect"
            if case.get("background_spectrum"):
                shutil.copyfile(source("background_spectrum"), spectrum)
            else:
                background_freq = freq
                if case.get("background_vcf"):
                    background_vcf = source("background_vcf")
                    background_mask = source("background_mask")
                    if case["split"] == "evaluation" and background_vcf == vcf:
                        raise ValueError("External background VCF must be independent of evaluation VCF")
                    background_freq = directory / "background.freq"
                    write_frequency(background_freq, list(variants(
                        background_vcf, not args.unfolded, background_mask, args.bcftools)))
                command = [executable, "-f", str(background_freq), str(spectrum)]
                provenance["commands"].append(command)
                save()
                run_command(command, directory / "spectrum.log")
            for configuration in configurations:
                command = [executable, "-l", str(args.grid), str(freq), str(spectrum)]
                if configuration != "physical":
                    recfile = directory / f"{configuration}.recombination.tsv"
                    map_path = source("recombination_map") if configuration in ("map", "bmap") else None
                    write_recombination(recfile, frequency_sites, map_path, args.constant_rate)
                    command[1] = "-lr"
                    command.append(str(recfile))
                if configuration == "bmap":
                    if (not args.unfolded or not case.get("background_spectrum")
                            or not case.get("sf2_frequency")
                            or not any(site.x == site.n for site in frequency_sites)):
                        raise ValueError("B correction requires --unfolded and an externally prepared "
                                         "background_spectrum including fixed differences; the SNP-only "
                                         "converter cannot manufacture outgroup divergence")
                    n1, n2, divergence = (float(case[key]) for key in ("n1", "n2", "divergence_time"))
                    if any(not math.isfinite(value) or value <= 0
                           for value in (n1, n2, divergence)):
                        raise ValueError("B correction requires positive finite n1, n2, divergence_time")
                    if divergence / (2 * n1) < 4:
                        raise ValueError("B correction needs an adequately diverged outgroup; "
                                         "this runner requires T/(2*N1) >= 4")
                    bfile = directory / "bvalues.tsv"
                    write_bvalues(bfile, frequency_sites, source("bmap"))
                    command[1] = "-lrb"
                    command.extend([str(bfile), str(n1), str(n2), str(divergence)])
                scan = directory / f"{configuration}.sf2.out"
                command.append(str(scan))
                provenance["commands"].append(command)
                save()
                run_command(command, directory / f"{configuration}.log")
                score, position = scan_peak(scan)
                row = {"case": name, "scenario": scenario, "split": case["split"],
                       "configuration": configuration, "score": score, "peak_position": position,
                       "sweep_position": truth if truth is not None else "",
                       "localization_error_bp": abs(position - float(truth)) if truth is not None else ""}
                writer.writerow(row)
                handle.flush()
            if args.evo:
                from .evo import EvoPredictor, score_sites

                reference = Reference(args.seqrepo, case["dataset"])
                provenance.setdefault("references", {})[name] = {
                    "assembly": "GRCh38.p14", "accession": reference.accession}
                selected = sites[::args.evo_stride]
                # Fail on reference mismatch before a model download or GPU allocation.
                for site in selected:
                    reference.check(site)
                if predictor is None:
                    predictor = EvoPredictor(args.model, args.allow_rocm, args.context)
                    provenance["evo"] = predictor.metadata()
                    save()
                scores = directory / "evo.tsv"
                score_sites(selected, reference, predictor, scores, args.context, offsets)
                score, position = evo_peak(scores, args.window_size)
                row = {"case": name, "scenario": scenario, "split": case["split"],
                       "configuration": "evo_mean_abs_delta", "score": score,
                       "peak_position": position, "sweep_position": truth if truth is not None else "",
                       "localization_error_bp": abs(position - float(truth)) if truth is not None else ""}
                writer.writerow(row)
                handle.flush()
    provenance["status"] = "complete"
    save()
    return summary


def evaluate(summary: Path, output: Path, alpha: float) -> None:
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    with summary.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = [(row["case"], row["configuration"]) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate case/configuration summary rows")
    if ({row["case"] for row in rows if row["split"] == "calibration"}
            & {row["case"] for row in rows if row["split"] == "evaluation"}):
        raise ValueError("Calibration and evaluation case IDs overlap")
    configurations = sorted({row["configuration"] for row in rows})
    if not configurations:
        raise ValueError("Empty summary")
    report = []
    for configuration in configurations:
        group = [row for row in rows if row["configuration"] == configuration]
        calibration = [row for row in group if row["split"] == "calibration"]
        if not calibration or any(row["scenario"] not in ("neutral", "bgs") for row in calibration):
            raise ValueError(f"{configuration}: calibration must contain only no-sweep controls")
        null = [float(row["score"]) for row in calibration]
        if any(not math.isfinite(value) for value in null):
            raise ValueError("Nonfinite calibration scores")
        if 1 / (len(null) + 1) > alpha:
            raise ValueError(f"{configuration}: need at least {math.ceil(1 / alpha) - 1} "
                             f"null calibration replicates for alpha={alpha}")
        for row in group:
            if row["split"] != "evaluation":
                continue
            if row["scenario"] == "unknown":
                raise ValueError("Unknown truth must use exploratory, not evaluation, split")
            score = float(row["score"])
            if not math.isfinite(score):
                raise ValueError("Nonfinite evaluation score")
            pvalue = (1 + sum(value >= score for value in null)) / (len(null) + 1)
            report.append({**row, "empirical_p": pvalue, "detected": pvalue <= alpha})
    if not report:
        raise ValueError("No evaluation cases")
    output.write_text(json.dumps({
        "alpha": alpha, "unit": "maximum score per case/region",
        "calibration": "pooled no-sweep controls per configuration; >= handles ties conservatively",
        "results": report,
        "rates": [{"configuration": config, "scenario": scenario,
                   "n": len(subset), "detection_rate": sum(r["detected"] for r in subset) / len(subset)}
                  for config in configurations for scenario in ("neutral", "bgs", "sweep", "bgs_sweep")
                  if (subset := [r for r in report if r["configuration"] == config
                                 and r["scenario"] == scenario])]
    }, indent=2, allow_nan=False) + "\n")
