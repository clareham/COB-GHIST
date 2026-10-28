"""Shared result types for challenge overviews and exploratory estimates."""

from dataclasses import dataclass, field

from .common import plt
from matplotlib.figure import Figure


@dataclass
class Estimate:
    target: str
    value: float | int | str | None
    units: str
    method: str
    source: str
    status: str = "exploratory"
    interval: tuple[float, float] | None = None


@dataclass
class ChallengeSummary:
    figure: Figure
    estimates: list[Estimate] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    assignments: list[dict[str, str | int | float]] = field(default_factory=list)


def overview_axes(count: int, title: str):
    figure, axes = plt.subplots(count, 1, figsize=(12, 4 * count),
                               layout="constrained", squeeze=False)
    figure.suptitle(title, fontsize=15)
    return figure, axes[:, 0]
