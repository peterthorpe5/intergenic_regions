"""Sequence-level exact enrichment tests and log-space FDR correction."""

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
from scipy.stats import hypergeom


def enrichment_test(
    *,
    positive_hits: int,
    negative_hits: int,
    positive_total: int,
    negative_total: int,
) -> dict[str, Any]:
    """Test motif presence using one-sided Fisher/hypergeometric enrichment.

    Args:
        positive_hits: Positive sequences containing one or more sites.
        negative_hits: Negative sequences containing one or more sites.
        positive_total: Number of positive sequences.
        negative_total: Number of negative sequences.

    Returns:
        Fractions, fold enrichment, log-space p-value and a Haldane-corrected
        odds ratio with approximate 95% confidence limits. Pseudocounts are
        used only for the effect estimate, never for the significance test.
        Infinite fold enrichment is represented by ``None`` in strict JSON.

    Raises:
        ValueError: Counts or totals are invalid.
    """
    if positive_total < 1 or negative_total < 1:
        raise ValueError("Both sequence totals must be positive")
    if (
        not 0 <= positive_hits <= positive_total
        or not 0 <= negative_hits <= negative_total
    ):
        raise ValueError("Hit counts are outside sequence totals")
    a, b = positive_hits, positive_total - positive_hits
    c, d = negative_hits, negative_total - negative_hits
    log_p = float(
        hypergeom.logsf(
            k=a - 1,
            M=positive_total + negative_total,
            n=a + c,
            N=positive_total,
        )
    )
    if not math.isfinite(log_p):
        # logpmf plus log-sum-exp avoids survival-function underflow.
        from scipy.special import logsumexp

        last = min(positive_total, a + c)
        log_p = float(
            logsumexp(
                hypergeom.logpmf(
                    k=np.arange(a, last + 1),
                    M=positive_total + negative_total,
                    n=a + c,
                    N=positive_total,
                )
            )
        )
    log_p = min(0.0, log_p)
    log_odds = math.log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)))
    se = math.sqrt(sum(1 / (v + 0.5) for v in (a, b, c, d)))
    foreground = a / positive_total
    background = c / negative_total
    return {
        "positive_hits": a,
        "positive_total": positive_total,
        "negative_hits": c,
        "negative_total": negative_total,
        "positive_fraction": foreground,
        "negative_fraction": background,
        "fold_enrichment": foreground / background
        if background
        else None
        if foreground
        else 1.0,
        "odds_ratio_haldane": math.exp(log_odds),
        "odds_ratio_ci_low": math.exp(log_odds - 1.96 * se),
        "odds_ratio_ci_high": math.exp(log_odds + 1.96 * se),
        "p_value": math.exp(log_p),
        "log_p_value": log_p,
        "minus_log10_p": -log_p / math.log(10),
    }


def adjust_fdr(
    *, log_p_values: Sequence[float], family_size: int | None = None
) -> list[dict[str, float]]:
    """Apply Benjamini-Hochberg in log space, including unobserved tests.

    Args:
        log_p_values: Natural logarithms of p-values.
        family_size: Total hypotheses, including unreported p=1 hypotheses.

    Returns:
        Q-values and their logarithms in the original input order.

    Raises:
        ValueError: Logs are invalid or the family is smaller than the list.
    """
    count = len(log_p_values)
    total = count if family_size is None else family_size
    if (
        total < count
        or total < 0
        or any(not math.isfinite(p) or p > 0 for p in log_p_values)
    ):
        raise ValueError("Invalid p-values or multiple-testing family size")
    if not count:
        return []
    order = sorted(range(count), key=lambda i: log_p_values[i])
    adjusted = [0.0] * count
    minimum = 0.0
    for rank_index in range(count - 1, -1, -1):
        original = order[rank_index]
        value = (
            log_p_values[original] + math.log(total) - math.log(rank_index + 1)
        )
        minimum = min(minimum, value)
        adjusted[original] = minimum
    return [
        {
            "q_value": math.exp(value),
            "log_q_value": value,
            "minus_log10_q": -value / math.log(10),
        }
        for value in adjusted
    ]


def kmer_family_size(
    *, lengths: Sequence[int], both_strands: bool = True
) -> int:
    """Count all possible tested k-mers, combining reverse complements.

    Args:
        lengths: Unique k-mer lengths, between 2 and 10 inclusive.
        both_strands: Treat a word and its reverse complement as one test.

    Returns:
        Full hypothesis family, including words absent from both sets.

    Raises:
        ValueError: Lengths are duplicated or outside the supported range.
    """
    if len(set(lengths)) != len(lengths) or any(
        k < 2 or k > 10 for k in lengths
    ):
        raise ValueError("Unique k-mer lengths between 2 and 10 are required")
    return sum(
        (4**k + (4 ** (k // 2) if k % 2 == 0 else 0)) // 2
        if both_strands
        else 4**k
        for k in lengths
    )
