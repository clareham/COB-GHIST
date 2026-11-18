"""Empirical upper-tail probabilities against independent region maxima."""

from bisect import bisect_left
import math


def empirical_p(score: float, sorted_null: list[float]) -> float:
    if not math.isfinite(score) or not sorted_null or any(
            not math.isfinite(value) for value in (sorted_null[0], sorted_null[-1])):
        raise ValueError("Empirical p-values require a finite score and finite sorted null maxima")
    return (1 + len(sorted_null) - bisect_left(sorted_null, score)) / (len(sorted_null) + 1)
