"""GHIST BED3 predictions with separate SF2 calibration and raw gLM evidence."""

from bisect import bisect_left
import csv
from dataclasses import dataclass
import json
import logging
import math
from pathlib import Path

from .calibration import empirical_p
from .inputs import callable_intervals, file_metadata as digest

LOG = logging.getLogger(__name__)
CONFIGURATIONS = {"physical", "constant", "map", "bmap"}
RULES = "https://www.synapse.org/Synapse:syn74803798/wiki/640755"


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_bed(path: Path, chrom: str, intervals: list[tuple[int, int]]) -> None:
    with path.open("w") as handle:
        for start, end in intervals:
            handle.write(f"{chrom.removeprefix('chr')}\t{start}\t{end}\n")


def null_maxima(rows: list[dict[str, str]]) -> dict[str, list[float]]:
    keys = [(row["case"], row["configuration"]) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate case/configuration summary rows")
    result: dict[str, list[float]] = {}
    for row in rows:
        if row["configuration"] not in CONFIGURATIONS or row["split"] != "calibration":
            continue
        if row["scenario"] not in ("neutral", "bgs"):
            raise ValueError("SF2 calibration must contain only neutral/bgs no-sweep controls")
        score = float(row["score"])
        if not math.isfinite(score) or score < 0:
            raise ValueError("SF2 calibration maxima must be finite and nonnegative")
        result.setdefault(row["configuration"], []).append(score)
    return {key: sorted(values) for key, values in result.items()}


def check_independence(source: dict, targets: dict, calibration_cases: set[str]) -> None:
    hashes = []
    for name in calibration_cases:
        hashes.append(source["inputs"][f"{name}:vcf"]["sha256"])
    if len(set(hashes)) != len(hashes):
        raise ValueError("Calibration replicates reuse identical VCF content")
    for name, case in targets["cases"].items():
        if case["split"] != "calibration" and targets["inputs"][f"{name}:vcf"]["sha256"] in hashes:
            raise ValueError(f"{name}: calibration and prediction VCFs must be independent")


@dataclass(frozen=True)
class Point:
    position: float
    clr: float
    start: int
    end: int


def scan_cells(path: Path, start: int, end: int) -> list[Point]:
    rows = read_table(path)
    positions = [float(row["location"]) for row in rows]
    scores = [float(row["LR"]) for row in rows]
    if (not rows or any(not math.isfinite(p) or not start + 1 <= p <= end for p in positions)
            or any(not math.isfinite(s) or s < 0 for s in scores)
            or any(b <= a for a, b in zip(positions, positions[1:]))):
        raise ValueError(f"{path}: expected sorted finite SF2 locations within the scanned span and nonnegative CLR")
    boundaries = [start, *(math.floor((a + b) / 2)
                           for a, b in zip(positions, positions[1:])), end]
    # Sub-base grid spacing can map several locations to one BED base.
    return [Point(p, s, min(math.floor(p) - 1, boundaries[i]),
                  max(math.floor(p), boundaries[i + 1]))
            for i, (p, s) in enumerate(zip(positions, scores))]


def peak_indices(points: list[Point]) -> list[int]:
    result = []
    left = 0
    while left < len(points):
        right = left + 1
        while right < len(points) and points[right].clr == points[left].clr:
            right += 1
        if ((left == 0 or points[left].clr > points[left - 1].clr)
                and (right == len(points) or points[left].clr > points[right].clr)):
            result.append(left)
        left = right
    return result


def clipped_cells(points: list[Point], mask: list[tuple[int, int]]) -> list[tuple[int, int, int]]:
    result = []
    index = 0
    for i, point in enumerate(points):
        while index < len(mask) and mask[index][1] <= point.start:
            index += 1
        j = index
        while j < len(mask) and mask[j][0] < point.end:
            start, end = max(point.start, mask[j][0]), min(point.end, mask[j][1])
            if start < end:
                result.append((start, end, i))
            j += 1
    return result


class GLMEvidence:
    def __init__(self, path: Path | None, chrom: str, bed: Path):
        self.enabled = path is not None
        sites: dict[int, tuple[float, int]] = {}
        self.maxima: dict[int, float] = {}
        if path is not None:
            with path.open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row["chrom"].removeprefix("chr") != chrom.removeprefix("chr"):
                        raise ValueError(f"{path}: gLM chromosome does not match SF2")
                    position = int(row["position"])
                    if position < 1:
                        raise ValueError(f"{path}: invalid gLM position")
                    if row["status"] != "ok":
                        continue
                    value = abs(float(row["delta_logp"]))
                    if not math.isfinite(value):
                        raise ValueError(f"{path}: nonfinite gLM score")
                    total, count = sites.get(position, (0.0, 0))
                    sites[position] = total + value, count + 1
                    self.maxima[position] = max(self.maxima.get(position, 0.0), value)
        self.positions = sorted(sites)
        self.sums = [0.0]
        self.counts = [0]
        for position in self.positions:
            total, count = sites[position]
            self.sums.append(self.sums[-1] + total / count)
            self.counts.append(self.counts[-1] + count)
        write_bed(bed, chrom, [(p - 1, p) for p in self.positions])

    def annotate(self, start: int, end: int) -> dict:
        left, right = bisect_left(self.positions, start + 1), bisect_left(self.positions, end + 1)
        n = right - left
        return {
            "glm_status": "available" if n else "no_scored_sites" if self.enabled else "disabled",
            "glm_snps": n, "glm_rows": self.counts[right] - self.counts[left],
            "glm_mean_abs_delta": (self.sums[right] - self.sums[left]) / n if n else ".",
            "glm_max_abs_delta": max(self.maxima[p] for p in self.positions[left:right]) if n else "."}


def export_regions(run_dir: Path, alpha: float = .05, null_summary: Path | None = None,
                   require_complete: bool = False) -> dict:
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    run_dir = run_dir.resolve()
    provenance = json.loads((run_dir / "run.json").read_text())
    if require_complete and provenance["status"] != "complete":
        raise ValueError("BED re-export requires a complete run")
    if not provenance.get("cases"):
        raise ValueError("Run lacks case coordinate metadata; rerun with the BED-enabled pipeline")
    summary = run_dir / "summary.tsv"
    source_path = null_summary.resolve() if null_summary is not None else summary
    source = provenance
    if source_path != summary:
        source = json.loads((source_path.parent / "run.json").read_text())
        if source["status"] != "complete":
            raise ValueError("External null summary must belong to a completed run")
        for key in ("grid", "unfolded", "constant_rate"):
            if source["arguments"][key] != provenance["arguments"][key]:
                raise ValueError(f"External null run has mismatched {key}")
        if source["sweepfinder_binary"]["sha256"] != provenance["sweepfinder_binary"]["sha256"]:
            raise ValueError("External null run uses a different SweepFinder2 binary")
    null_rows = read_table(source_path)
    null = null_maxima(null_rows)
    calibration_cases = {row["case"] for row in null_rows
                         if row["split"] == "calibration" and row["configuration"] in CONFIGURATIONS}
    check_independence(source, provenance, calibration_cases)
    output = run_dir / "regions"
    output.mkdir(exist_ok=True)
    metadata = {
        "status": "running",
        "format": "GHIST headerless tab-separated BED3; 0-based half-open; no chr prefix",
        "rules": RULES, "alpha": alpha, "null_summary": digest(source_path),
        "run_summary": digest(summary),
        "pvalue": "(1 + count(null_region_maximum >= CLR)) / (1 + n_null)",
        "scope": "per case and SF2 configuration; not corrected across cases or configurations",
        "boundaries": "midpoints between SF2 grid locations, clipped to scanned span and callable mask; not confidence intervals",
        "glm_role": "raw evidence and annotations only; no gLM threshold or combined p-value",
        "cases": {}, "null_counts": {key: len(value) for key, value in null.items()}}
    configurations = provenance["arguments"]["configurations"].split(",")
    metadata["configurations"] = configurations
    evidence_columns = ["chrom", "start0", "end0", "grid_position1", "clr", "empirical_p",
                        "n_null", "status", "significant"]
    region_columns = ["kind", "chrom", "start0", "end0", "peak_position1", "peak_clr",
                      "empirical_p", "n_null", "status", "glm_status", "glm_snps", "glm_rows",
                      "glm_mean_abs_delta", "glm_max_abs_delta"]
    metadata["evidence_columns"] = evidence_columns
    metadata["region_columns"] = region_columns
    metadata_file = output / "metadata.json"
    metadata_file.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    for name, case in provenance["cases"].items():
        directory = output / name
        directory.mkdir(exist_ok=True)
        raw = run_dir / name / "ntv3.tsv" if provenance["arguments"]["ntv3"] else None
        glm = GLMEvidence(raw, case["chrom"], directory / "ntv3.sites.bed")
        mask_path = Path(case["mask"])
        if digest(mask_path)["sha256"] != provenance["inputs"][f"{name}:mask"]["sha256"]:
            raise ValueError(f"{name}: callable mask changed since the scan")
        mask = callable_intervals(mask_path, case["chrom"])
        metadata["cases"][name] = {
            **case, "mask": digest(mask_path), "raw_glm": digest(raw) if raw else None,
            "scans": {}, "prediction_counts": {},
            "files": ["ntv3.sites.bed", *(
                f"{configuration}.{suffix}" for configuration in configurations
                for suffix in ("candidates.bed", "sweeps.bed", "sf2.tsv", "regions.tsv"))]}
        for configuration in configurations:
            path = run_dir / name / f"{configuration}.sf2.out"
            points = scan_cells(path, case["start0"], case["end0"])
            cells = clipped_cells(points, mask)
            maxima = null.get(configuration, [])
            status = ("calibration_only" if case["split"] == "calibration" else
                      "uncalibrated" if not maxima else
                      "insufficient_nulls" if 1 / (len(maxima) + 1) > alpha else "calibrated")
            if status in ("uncalibrated", "insufficient_nulls"):
                LOG.warning("%s/%s: %s; no significant sweep predictions", name, configuration, status)
            probabilities = [empirical_p(point.clr, maxima) if maxima and status != "calibration_only"
                             else None for point in points]
            significant = [status == "calibrated" and p is not None and p <= alpha
                           for p in probabilities]
            evidence = [{"chrom": case["chrom"].removeprefix("chr"), "start0": start, "end0": end,
                         "grid_position1": points[i].position, "clr": points[i].clr,
                         "empirical_p": probabilities[i] if probabilities[i] is not None else ".",
                         "n_null": len(maxima), "status": status, "significant": int(significant[i])}
                        for start, end, i in cells]
            peaks = set(peak_indices(points))
            candidates = [(start, end, i) for start, end, i in cells if i in peaks]
            predictions: list[tuple[int, int, int]] = []
            for start, end, i in cells:
                if not significant[i]:
                    continue
                if predictions and start <= predictions[-1][1]:
                    left, right, peak = predictions[-1]
                    peak = i if points[i].clr > points[peak].clr else peak
                    predictions[-1] = left, max(right, end), peak
                else:
                    predictions.append((start, end, i))
            for kind, intervals in (("candidates", candidates), ("sweeps", predictions)):
                write_bed(directory / f"{configuration}.{kind}.bed", case["chrom"],
                          [(start, end) for start, end, _ in intervals])
            regions = []
            for kind, intervals in (("candidate", candidates), ("prediction", predictions)):
                for start, end, peak in intervals:
                    regions.append({
                        "kind": kind, "chrom": case["chrom"].removeprefix("chr"),
                        "start0": start, "end0": end, "peak_position1": points[peak].position,
                        "peak_clr": points[peak].clr,
                        "empirical_p": probabilities[peak] if probabilities[peak] is not None else ".",
                        "n_null": len(maxima), "status": status, **glm.annotate(start, end)})
            write_table(directory / f"{configuration}.sf2.tsv", evidence_columns, evidence)
            write_table(directory / f"{configuration}.regions.tsv", region_columns, regions)
            metadata["cases"][name]["scans"][configuration] = digest(path)
            metadata["cases"][name]["prediction_counts"][configuration] = len(predictions)
    metadata["status"] = "complete"
    metadata_file.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    result = {"directory": str(output), "metadata": digest(metadata_file), "alpha": alpha}
    if require_complete:
        provenance["region_exports"] = result
        (run_dir / "run.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return result
