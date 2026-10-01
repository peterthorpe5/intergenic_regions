"""Exact interventional SHAP for linear sequence models in log-odds space."""

import logging
import math
from collections import defaultdict
from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
from scipy.special import expit

LOGGER = logging.getLogger(__name__)
COMPOSITION_FEATURES = ("gc_fraction", "log_length", "ambiguous_fraction")
SHAP_NOTE = (
    "Exact interventional SHAP explains model log-odds, "
    "not enhancer activity. "
    "Correlated k-mers share information; feature contributions are neither "
    "causal effects nor independent biological evidence."
)


def explanation_sample(
    *, labels: NDArray[np.int64], maximum: int, seed: int
) -> NDArray[np.int64]:
    """Select a reproducible, approximately class-proportional SHAP sample.

    Args:
        labels: Binary labels in sequence order.
        maximum: Maximum explained sequences, at least two.
        seed: Sampling seed.

    Returns:
        Sorted indices, retaining both classes when subsampling.

    Raises:
        ValueError: Labels or the limit are invalid.
    """
    if (
        not isinstance(maximum, int)
        or isinstance(maximum, bool)
        or maximum < 2
        or labels.ndim != 1
        or set(labels.tolist()) != {0, 1}
    ):
        raise ValueError("SHAP needs binary labels and at least two samples")
    if len(labels) <= maximum:
        return np.arange(len(labels), dtype=np.int64)
    rng = np.random.default_rng(seed=seed)
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    count = min(
        len(positive),
        max(1, min(maximum - 1, round(maximum * len(positive) / len(labels)))),
    )
    negative_count = min(len(negative), maximum - count)
    count = min(len(positive), maximum - negative_count)
    return np.sort(
        np.concatenate(
            (
                rng.choice(a=positive, size=count, replace=False),
                rng.choice(a=negative, size=negative_count, replace=False),
            )
        )
    ).astype(np.int64)


def linear_shap_values(
    *,
    matrix: Any,
    coefficients: NDArray[np.float64],
    background_mean: NDArray[np.float64],
    intercept: float,
) -> tuple[NDArray[np.float64], float]:
    """Calculate exact linear SHAP values against a training-only background.

    Args:
        matrix: Dense or sparse query features in the model's feature space.
        coefficients: One coefficient per feature.
        background_mean: Mean training features in the same transformed space.
        intercept: Model intercept.

    Returns:
        Contributions and expected log-odds. Their row sums reconstruct logits.

    Raises:
        ValueError: Shapes or numerical values are invalid.
    """
    values = np.asarray(
        matrix.toarray() if hasattr(matrix, "toarray") else matrix,
        dtype=float,
    )
    coefficients = np.asarray(coefficients, dtype=float)
    background_mean = np.asarray(background_mean, dtype=float)
    if (
        values.ndim != 2
        or coefficients.ndim != 1
        or background_mean.shape != coefficients.shape
        or values.shape[1] != len(coefficients)
        or not np.isfinite(values).all()
        or not np.isfinite(coefficients).all()
        or not np.isfinite(background_mean).all()
        or not math.isfinite(intercept)
    ):
        raise ValueError("Invalid linear SHAP feature space or parameters")
    contributions = (values - background_mean) * coefficients
    baseline = float(intercept + background_mean @ coefficients)
    return contributions, baseline


def summarise_shap(
    *,
    contexts: Sequence[dict[str, Any]],
    identifiers: Sequence[str],
    labels: Sequence[int] | None = None,
    max_features: int = 20,
    scope: str = "held_out",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Aggregate fold-specific SHAP and retain exact additive explanations.

    All features contribute to global importance. Only the leading features
    are retained in the long table; the remainder has an explicit additive row.
    Feature blocks bound memory, even for a large sparse word vocabulary.

    Args:
        contexts: Query matrices, coefficients, training background means,
            feature names, intercepts, fold IDs and original row indices.
        identifiers: Complete analysed sequence order.
        labels: Optional binary labels in the same order.
        max_features: Number of leading features shown in per-sequence tables.
        scope: ``held_out`` or ``candidate``; queries remain unvalidated.

    Returns:
        Long explanations, global mean-absolute importance and method summary.

    Raises:
        ValueError: Settings, contexts or sequence indices are invalid.
    """
    if (
        not isinstance(max_features, int)
        or isinstance(max_features, bool)
        or not 1 <= max_features <= 100
        or scope not in {"held_out", "candidate"}
        or len(set(identifiers)) != len(identifiers)
        or (labels is not None and len(labels) != len(identifiers))
    ):
        raise ValueError("Invalid SHAP summary settings or sequence IDs")
    totals: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    seen: set[int] = set()
    contribution_totals: list[NDArray[np.float64]] = []
    for context in contexts:
        names = context["feature_names"]
        matrix = csr_matrix(context["matrix"])
        indices = list(context["indices"])
        if (
            matrix.shape != (len(indices), len(names))
            or len(set(names)) != len(names)
            or any(
                i < 0 or i >= len(identifiers) or i in seen for i in indices
            )
            or len(set(indices)) != len(indices)
        ):
            raise ValueError("Invalid or repeated SHAP context indices")
        seen.update(int(i) for i in indices)
        summed = np.zeros(len(indices), dtype=float)
        for first in range(0, len(names), 256):
            last = min(first + 256, len(names))
            values, _ = linear_shap_values(
                matrix=matrix[:, first:last],
                coefficients=context["coefficients"][first:last],
                background_mean=context["background_mean"][first:last],
                intercept=context["intercept"],
            )
            summed += values.sum(axis=1)
            for offset, name in enumerate(names[first:last]):
                totals[name][0] += float(np.abs(values[:, offset]).sum())
                totals[name][1] += float(values[:, offset].sum())
                totals[name][2] += len(indices)
        contribution_totals.append(summed)
    if not seen:
        raise ValueError("No sequences available for SHAP explanations")
    importance: list[dict[str, Any]] = [
        {
            "feature": name,
            "feature_type": "composition"
            if name in COMPOSITION_FEATURES
            else "kmer",
            "mean_absolute_shap": values[0] / len(seen),
            "mean_signed_shap": values[1] / len(seen),
            "sequences_with_feature_in_model": int(values[2]),
            "explained_sequences": len(seen),
            "scope": scope,
        }
        for name, values in totals.items()
    ]
    importance.sort(
        key=lambda row: (-row["mean_absolute_shap"], row["feature"])
    )
    leading = [row["feature"] for row in importance[:max_features]]
    rows: list[dict[str, Any]] = []
    maximum_error = 0.0
    for context, summed in zip(contexts, contribution_totals, strict=True):
        matrix = csr_matrix(context["matrix"])
        coefficients = context["coefficients"]
        background = context["background_mean"]
        names = {name: i for i, name in enumerate(context["feature_names"])}
        baseline = float(context["intercept"] + background @ coefficients)
        logits = (
            np.asarray(matrix @ coefficients).ravel() + context["intercept"]
        )
        present = [name for name in leading if name in names]
        columns = [names[name] for name in present]
        raw = matrix[:, columns].toarray()
        contributions, _ = linear_shap_values(
            matrix=raw,
            coefficients=coefficients[columns],
            background_mean=background[columns],
            intercept=0.0,
        )
        positions = {name: i for i, name in enumerate(present)}
        for local, original in enumerate(context["indices"]):
            common = {
                "sequence_id": identifiers[int(original)],
                "label": int(labels[int(original)])
                if labels is not None
                else None,
                "fold": context["fold"],
                "base_value": baseline,
                "model_logit": float(logits[local]),
                "signature_score": float(expit(logits[local])),
                "scope": scope,
            }
            subtotal = float(contributions[local].sum())
            for name in leading:
                position = positions.get(name)
                rows.append(
                    {
                        **common,
                        "feature": name,
                        "feature_value": float(raw[local, position])
                        if position is not None
                        else None,
                        "shap_value": float(contributions[local, position])
                        if position is not None
                        else 0.0,
                        "feature_available": position is not None,
                    }
                )
            remainder = float(summed[local] - subtotal)
            rows.append(
                {
                    **common,
                    "feature": "other_features",
                    "feature_value": None,
                    "shap_value": remainder,
                    "feature_available": True,
                }
            )
            maximum_error = max(
                maximum_error,
                abs(float(logits[local]) - baseline - subtotal - remainder),
            )
    rows.sort(key=lambda row: (row["sequence_id"], row["feature"]))
    summary = {
        "status": "completed",
        "method": "exact_interventional_linear_shap",
        "scope": scope,
        "output_scale": "positive_class_log_odds",
        "background": "training-only feature means for each explained model",
        "explained_sequences": len(seen),
        "total_sequences": len(identifiers),
        "subsampled": len(seen) < len(identifiers),
        "features_assessed": len(importance),
        "features_displayed": len(leading),
        "maximum_reconstruction_error": maximum_error,
        "interpretation": SHAP_NOTE,
    }
    LOGGER.info("SHAP explained %d sequences (%s)", len(seen), scope)
    return rows, importance, summary
