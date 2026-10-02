"""Weakly supervised window models with source-level held-out validation."""

import hashlib
import logging
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, hstack
from sklearn.metrics import average_precision_score, roc_auc_score

from intergenic_regions.background import canonical_sequence
from intergenic_regions.explanations import (
    COMPOSITION_FEATURES,
    explanation_sample,
    summarise_shap,
)
from intergenic_regions.learning import (
    composition_matrix,
    fit_classifier,
    fit_vocabulary,
    validation_splits,
    word_matrix,
)
from intergenic_regions.motifs import kmer_counts
from intergenic_regions.statistics import kmer_family_size
from intergenic_regions.windows import RegulatoryWindow

LOGGER = logging.getLogger(__name__)
WINDOW_MODEL_NOTE = (
    "Window labels are inherited from source flanks, not experimentally "
    "labelled enhancer intervals. Scores measure positive-class sequence "
    "resemblance. Validation metrics use one mean score per source at each "
    "scale and do not measure enhancer detection or boundary accuracy. "
    "Peak scores across windows/scales are exploratory, are affected by "
    "search opportunity and have no region-level p-value or FDR."
)


def parent_window_weights(
    *, parent_indices: NDArray[np.int64], parent_labels: NDArray[np.int64]
) -> NDArray[np.float64]:
    """Give each source equal training mass within its supplied class.

    Each class receives half of the total mass. Weight mass equals the number
    of represented parents, making regularisation independent of how many
    overlapping windows a long source contributes.

    Args:
        parent_indices: Window-to-parent indices in the training subset.
        parent_labels: Binary labels indexed by original parent index.

    Returns:
        Positive window weights, summing to the represented parent count.

    Raises:
        ValueError: Shapes, labels or parent indices are invalid.
    """
    if (
        parent_indices.ndim != 1
        or not len(parent_indices)
        or parent_labels.ndim != 1
        or not np.issubdtype(parent_indices.dtype, np.integer)
        or min(parent_indices.tolist()) < 0
        or max(parent_indices.tolist()) >= len(parent_labels)
        or not set(parent_labels.tolist()) <= {0, 1}
    ):
        raise ValueError("Invalid window-to-parent training indices")
    counts = Counter(parent_indices.tolist())
    classes = Counter(int(parent_labels[i]) for i in counts)
    if set(classes) != {0, 1}:
        raise ValueError("Window training parents require both classes")
    return np.asarray(
        [
            len(counts) / (2 * classes[int(parent_labels[i])] * counts[int(i)])
            for i in parent_indices
        ],
        dtype=float,
    )


def fit_window_model(
    *,
    windows: Sequence[RegulatoryWindow],
    parent_groups: Mapping[str, str],
    lengths: Sequence[int] = (4, 5, 6),
    folds: int = 5,
    seed: int = 42,
    max_features: int = 100000,
    regularisation: float = 1.0,
    shap: bool = True,
    shap_max_sequences: int = 1000,
    shap_max_features: int = 20,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
]:
    """Train one scale using grouped sources and parent-balanced windows.

    Splits are made on original parents before expanding to windows. The word
    vocabulary uses training-parent document frequency, composition scaling
    uses training weights and SHAP uses the weighted training background.
    All predictions for a source use a model that excluded its entire group.

    Args:
        windows: Complete windows at exactly one length.
        parent_groups: Group unions from all scales, including shared windows.
        lengths: Short motif-feature lengths, distinct from window lengths.
        folds: Source-level validation folds.
        seed: Deterministic split/classifier seed.
        max_features: Training-only vocabulary limit.
        regularisation: Fixed inverse L2 strength C.
        shap: Produce exact held-out linear SHAP explanations.
        shap_max_sequences: Maximum explained windows, at least two.
        shap_max_features: Leading explanation features, from one to 100.

    Returns:
        Window predictions, source predictions, fold metrics, summary,
        reusable JSON model and optional SHAP tables.

    Raises:
        ValueError: Data, settings or groups cannot support fitting.
    """
    kmer_family_size(lengths=lengths)
    if (
        not windows
        or not lengths
        or len({w.width for w in windows}) != 1
        or len({w.sequence_id for w in windows}) != len(windows)
        or max_features < 1
        or not isinstance(shap_max_features, int)
        or isinstance(shap_max_features, bool)
        or not 1 <= shap_max_features <= 100
    ):
        raise ValueError("Invalid single-scale window model settings")
    width = windows[0].width
    parents = list(dict.fromkeys(w.parent_id for w in windows))
    parent_index = {identifier: i for i, identifier in enumerate(parents)}
    by_parent = {w.parent_id: int(w.label == "positive") for w in windows}
    if any(
        by_parent[w.parent_id] != int(w.label == "positive") for w in windows
    ):
        raise ValueError("Window parent labels are inconsistent")
    parent_labels = np.asarray([by_parent[p] for p in parents], dtype=np.int64)
    if min(Counter(parent_labels.tolist()).get(n, 0) for n in (0, 1)) < 5:
        raise ValueError(
            "Window models require at least five parents per class"
        )
    if any(p not in parent_groups or not parent_groups[p] for p in parents):
        raise ValueError("Every window parent needs an independent group")
    groups = [parent_groups[p] for p in parents]
    splits = validation_splits(
        labels=parent_labels,
        folds=folds,
        seed=seed,
        groups=groups if len(set(groups)) != len(groups) else None,
    )
    indices = np.asarray(
        [parent_index[w.parent_id] for w in windows], dtype=np.int64
    )
    labels = parent_labels[indices]
    explanation_indices = explanation_sample(
        labels=labels, maximum=shap_max_sequences, seed=seed
    )
    counters = [
        kmer_counts(sequence=w.sequence, lengths=lengths) for w in windows
    ]
    parent_words: list[Counter[str]] = [Counter() for _ in parents]
    for counts, index in zip(counters, indices, strict=True):
        parent_words[int(index)].update(counts.keys())
    composition = composition_matrix(sequences=[w.sequence for w in windows])
    scores = np.empty(len(windows), dtype=float)
    baseline_scores = np.empty(len(windows), dtype=float)
    assignments = np.empty(len(windows), dtype=np.int64)
    parent_scores = np.empty(len(parents), dtype=float)
    parent_baseline = np.empty(len(parents), dtype=float)
    parent_folds = np.empty(len(parents), dtype=np.int64)
    metrics: list[dict[str, Any]] = []
    contexts: list[dict[str, Any]] = []

    def fit(
        train: NDArray[np.int64],
    ) -> tuple[Any, Any, Any, dict[str, int], Any, NDArray[np.float64]]:
        parent_train = np.unique(indices[train])
        vocabulary = fit_vocabulary(
            counters=[parent_words[int(i)] for i in parent_train],
            max_features=max_features,
        )
        words = word_matrix(counters=counters, vocabulary=vocabulary)
        weights = parent_window_weights(
            parent_indices=indices[train], parent_labels=parent_labels
        )
        composition_mean = np.average(
            composition[train], axis=0, weights=weights
        )
        constant = np.ptp(composition[train], axis=0) == 0
        composition_mean[constant] = composition[train[0], constant]
        variance = np.average(
            (composition[train] - composition_mean) ** 2,
            axis=0,
            weights=weights,
        )
        composition_scale = np.sqrt(variance)
        composition_scale[constant] = 1.0
        scaled = (composition - composition_mean) / composition_scale
        matrix = hstack([words, csr_matrix(scaled)], format="csr")
        model = fit_classifier(
            matrix=matrix[train],
            labels=labels[train],
            seed=seed,
            regularisation=regularisation,
            sample_weight=weights,
        )
        baseline = fit_classifier(
            matrix=scaled[train],
            labels=labels[train],
            seed=seed,
            regularisation=regularisation,
            sample_weight=weights,
        )
        mean = np.asarray(matrix[train].T @ weights).ravel() / weights.sum()
        return (
            model,
            baseline,
            (composition_mean, composition_scale),
            vocabulary,
            matrix,
            mean,
        )

    for fold, (train_parents, test_parents) in enumerate(splits, start=1):
        train = np.flatnonzero(np.isin(indices, train_parents))
        test = np.flatnonzero(np.isin(indices, test_parents))
        model, baseline, scaling, vocabulary, matrix, mean = fit(train)
        scores[test] = model.predict_proba(X=matrix[test])[:, 1]
        baseline_scores[test] = baseline.predict_proba(
            X=(composition[test] - scaling[0]) / scaling[1]
        )[:, 1]
        assignments[test] = fold
        for parent in test_parents:
            members = np.flatnonzero(indices == parent)
            parent_scores[parent] = scores[members].mean()
            parent_baseline[parent] = baseline_scores[members].mean()
            parent_folds[parent] = fold
        if shap:
            explained = test[np.isin(test, explanation_indices)]
            if len(explained):
                contexts.append(
                    {
                        "matrix": matrix[explained],
                        "coefficients": model.coef_[0].copy(),
                        "intercept": float(model.intercept_[0]),
                        "background_mean": mean,
                        "feature_names": [*vocabulary, *COMPOSITION_FEATURES],
                        "indices": explained,
                        "fold": fold,
                    }
                )
        metrics.append(
            {
                "window_length": width,
                "fold": fold,
                "train_parents": len(train_parents),
                "test_parents": len(test_parents),
                "train_windows": len(train),
                "test_windows": len(test),
                "roc_auc": float(
                    roc_auc_score(
                        parent_labels[test_parents],
                        parent_scores[test_parents],
                    )
                ),
                "average_precision": float(
                    average_precision_score(
                        parent_labels[test_parents],
                        parent_scores[test_parents],
                    )
                ),
                "baseline_roc_auc": float(
                    roc_auc_score(
                        parent_labels[test_parents],
                        parent_baseline[test_parents],
                    )
                ),
                "metric_unit": "one_mean_score_per_parent",
            }
        )
        LOGGER.info(
            "Window scale %d bp: source validation fold %d/%d",
            width,
            fold,
            folds,
        )
    predictions = [
        {
            "sequence_id": w.sequence_id,
            "fold": int(assignments[i]),
            "group": parent_groups[w.parent_id],
            "held_out_signature_score": float(scores[i]),
            "held_out_composition_score": float(baseline_scores[i]),
        }
        for i, w in enumerate(windows)
    ]
    parent_predictions = [
        {
            "parent_id": identifier,
            "window_length": width,
            "label": int(parent_labels[i]),
            "group": groups[i],
            "fold": int(parent_folds[i]),
            "windows": int(np.sum(indices == i)),
            "mean_held_out_signature_score": float(parent_scores[i]),
            "mean_held_out_composition_score": float(parent_baseline[i]),
            "peak_held_out_signature_score": float(scores[indices == i].max()),
        }
        for i, identifier in enumerate(parents)
    ]
    summary: dict[str, Any] = {
        "status": "completed",
        "window_length": width,
        "windows": len(windows),
        "parents": len(parents),
        "positive_parents": int(parent_labels.sum()),
        "negative_parents": int(len(parents) - parent_labels.sum()),
        "independent_groups": len(set(groups)),
        "folds": folds,
        "roc_auc": float(roc_auc_score(parent_labels, parent_scores)),
        "average_precision": float(
            average_precision_score(parent_labels, parent_scores)
        ),
        "baseline_roc_auc": float(
            roc_auc_score(parent_labels, parent_baseline)
        ),
        "baseline_average_precision": float(
            average_precision_score(parent_labels, parent_baseline)
        ),
        "metric_unit": "one_mean_score_per_parent",
        "training_weight_unit": "class_balanced_parents",
        "permutation_p_value": None,
        "interpretation": WINDOW_MODEL_NOTE,
        "shap": {"status": "disabled"},
    }
    explanations: dict[str, Any] = {}
    if shap:
        rows, importance, shap_summary = summarise_shap(
            contexts=contexts,
            identifiers=[w.sequence_id for w in windows],
            labels=labels.tolist(),
            max_features=shap_max_features,
        )
        shap_summary["sampling_unit"] = "windows_with_inherited_parent_labels"
        shap_summary["background"] = (
            "parent-balanced training-window means within each fold"
        )
        summary["shap"] = shap_summary
        explanations.update(
            rows=rows, importance=importance, summary=shap_summary
        )
    model, _, scaling, vocabulary, matrix, mean = fit(
        np.arange(len(windows), dtype=np.int64)
    )
    exported = {
        "schema_version": 1,
        "lengths": list(lengths),
        "window_length": width,
        "vocabulary": vocabulary,
        "coefficients": model.coef_[0].tolist(),
        "intercept": float(model.intercept_[0]),
        "composition_mean": scaling[0].tolist(),
        "composition_scale": scaling[1].tolist(),
        "feature_background_mean": mean.tolist(),
        "training_sequence_sha256": [
            hashlib.sha256(
                canonical_sequence(sequence=w.sequence).encode("ascii")
            ).hexdigest()
            for w in windows
        ],
        "summary": summary,
    }
    return (
        predictions,
        parent_predictions,
        metrics,
        summary,
        exported,
        explanations,
    )
