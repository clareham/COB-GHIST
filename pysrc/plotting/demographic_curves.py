"""Read and scale demographic curves without assuming a calendar generation length."""

import json
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import median

from .common import number, read_table
from .summary import Estimate


@dataclass
class Curve:
    left: list[float]
    right: list[float]
    size: list[float]
    label: str
    source: Path


def read_psmc(path: Path) -> tuple[list[float], list[float], float | None]:
    """Return RS time/size and TR theta from the same last complete iteration."""
    completed = None
    current: list[tuple[float, float]] | None = None
    theta = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "RD":
                current, theta = [], None
            elif fields[0] == "TR" and current is not None:
                theta = number(fields[1])
                if theta <= 0:
                    raise ValueError(f"{path}: PSMC theta must be positive")
            elif fields[0] == "RS" and current is not None:
                current.append((number(fields[2]), number(fields[3])))
            elif fields[0] == "//" and current:
                completed = current, theta
                current = None
    if completed is None:
        raise ValueError("No complete PSMC iteration containing RS records")
    points, theta = completed
    times, sizes = map(list, zip(*points))
    if any(t < 0 for t in times) or any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("PSMC times must be nonnegative and strictly increasing")
    if any(size <= 0 for size in sizes):
        raise ValueError("PSMC relative sizes must be positive")
    return times, sizes, theta


def mutation_rate(path: Path) -> float | None:
    """Use published assumptions, never species-name or current-config guesses."""
    if "genetic_history" not in path.parts:
        return None
    category_index = path.parts.index("genetic_history")
    challenge = Path(*path.parts[:category_index + 2])
    candidates = sorted((challenge / "msmc2").glob("*.run/assumptions.json"))
    own = path.with_name(path.name + ".assumptions.json")
    if own.is_file():
        candidates.append(own)
    rates = set()
    for candidate in candidates:
        metadata = json.loads(candidate.read_text(encoding="utf-8"))
        if not isinstance(metadata, dict):
            raise ValueError(f"{candidate}: assumptions must be a JSON object")
        options = metadata.get("options", metadata)
        if not isinstance(options, dict):
            raise ValueError(f"{candidate}: options must be a JSON object")
        raw = options.get("mutation_rate")
        if raw is not None:
            rate = number(str(raw))
            if rate <= 0:
                raise ValueError(f"{candidate}: mutation_rate must be positive")
            rates.add(rate)
    if len(rates) > 1:
        raise ValueError(f"{challenge}: conflicting published mutation rates")
    return rates.pop() if rates else None


def scaled_psmc(path: Path, rate: float) -> Curve:
    times, sizes, theta = read_psmc(path)
    if theta is None:
        raise ValueError(f"{path}: final complete iteration has no TR theta")
    metadata_path = path.with_name(path.name + ".assumptions.json")
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    if not isinstance(metadata, dict):
        raise ValueError(f"{metadata_path}: assumptions must be a JSON object")
    # Both pipeline converters use 100-bp bins; newer metadata can record it.
    bin_size = number(str(metadata.get("window_size_bp", 100)))
    if not math.isfinite(rate) or rate <= 0 or bin_size <= 0:
        raise ValueError("Mutation rate and PSMC bin size must be positive")
    n0 = theta / (4 * rate * bin_size)
    left = [2 * n0 * time for time in times]
    return Curve(left, left[1:] + [math.inf], [n0 * size for size in sizes],
                 "PSMC (one diploid sample)", path)


def scaled_msmc(path: Path, rate: float) -> list[Curve]:
    if not math.isfinite(rate) or rate <= 0:
        raise ValueError("Mutation rate must be positive")
    rows = read_table(path)
    left = [number(row["left_time_boundary"]) / rate for row in rows]
    right = [float(row["right_time_boundary"]) / rate for row in rows]
    if any(not (0 <= a < b) for a, b in zip(left, right)):
        raise ValueError(f"{path}: invalid MSMC time intervals")
    if any(a != b for a, b in zip(right, left[1:])):
        raise ValueError(f"{path}: MSMC intervals must be contiguous")
    columns = [key for key in ("lambda", "lambda_00", "lambda_11") if key in rows[0]]
    if not columns:
        raise ValueError("No within-population MSMC2 rate columns")
    curves = []
    for column in columns:
        rates = [number(row[column]) for row in rows]
        if any(value <= 0 for value in rates):
            raise ValueError(f"{path}: nonpositive within-population coalescence rate")
        curves.append(Curve(left, right, [1 / (2 * rate * value) for value in rates],
                            f"MSMC2 {column}", path))
    return curves


def change_estimates(curve: Curve, growth: bool, source: str) -> list[Estimate]:
    """Fit two constant log-size epochs to finite bins; no confidence claim."""
    points = [(a, size) for a, b, size in zip(curve.left, curve.right, curve.size)
              if math.isfinite(b)]
    target = "post_growth_scaling" if growth else "post_decline_fraction"
    method = ("Two-epoch least-squares split of log(Ne) across finite bins, at least "
              "3 bins per epoch; epoch sizes are geometric means. Equal bin weights; "
              "not an event likelihood fit or confidence interval.")
    if len(points) < 6:
        reason = "Insufficient finite bins for a two-epoch curve summary"
        return [Estimate(name, None, units, reason, source, "unavailable")
                for name, units in (("generations", "generations"), (target, "ratio"))]
    logs = [math.log(size) for _, size in points]
    total = sum((x - sum(logs) / len(logs)) ** 2 for x in logs)
    candidates = []
    for index in range(3, len(logs) - 2):
        recent = sum(logs[:index]) / index
        ancient = sum(logs[index:]) / (len(logs) - index)
        error = sum((x - recent) ** 2 for x in logs[:index])
        error += sum((x - ancient) ** 2 for x in logs[index:])
        candidates.append((error, index, recent, ancient))
    error, index, recent, ancient = min(candidates)
    ratio = math.exp(recent - ancient)
    threshold = 1.2 if growth else .8
    direction_matches = ((ratio >= threshold if growth else ratio <= threshold)
                         or math.isclose(ratio, threshold, rel_tol=1e-12))
    if total == 0 or 1 - error / total < 0.25 or not direction_matches:
        reason = "No direction-consistent >=20% size change with >=25% log-size SSE improvement"
        return [Estimate(name, None, units, reason, source, "unavailable")
                for name, units in (("generations", "generations"), (target, "ratio"))]
    return [Estimate("generations", points[index][0], "generations", method, source),
            Estimate(target, ratio, "ratio", method, source)]


def rccr_split(rows: list[dict[str, str]], source: str) -> Estimate:
    """Bracket the first sustained upward half-crossing; do not interpolate steps."""
    values = [number(row["relative_cross_coalescence"]) for row in rows]
    left = [number(row["left_generations"]) for row in rows]
    if (any(value < 0 for value in values) or any(time < 0 for time in left)
            or any(b <= a for a, b in zip(left, left[1:]))):
        raise ValueError("Invalid RCCR values or time order")
    for index in range(1, len(values) - 1):
        if values[index - 1] < 0.5 <= values[index] and values[index + 1] >= 0.5:
            return Estimate("generations_since_split", left[index], "generations",
                            "First upward RCCR=0.5 transition sustained for two bins; "
                            "split-time proxy only, confounded by migration. Bracket is "
                            "adjacent left boundaries, not a confidence interval.",
                            source, interval=(left[index - 1], left[index]))
    return Estimate("generations_since_split", None, "generations",
                    "No sustained upward RCCR=0.5 crossing; split time not identifiable",
                    source, "unavailable")


def recent_size(curve: Curve) -> float:
    return median(curve.size[:3])
