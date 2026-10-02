"""Interpretable sequence learning with held-out and grouped validation."""

import csv
import hashlib
import logging
import math
from collections import Counter, defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, hstack
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from intergenic_regions.background import (
    canonical_sequence,
    check_sequence_sets,
)
from intergenic_regions.explanations import (
    COMPOSITION_FEATURES,
    explanation_sample,
    summarise_shap,
)
from intergenic_regions.genome import sequence_composition
from intergenic_regions.io import open_text
from intergenic_regions.motifs import kmer_counts
from intergenic_regions.statistics import kmer_family_size

LOGGER = logging.getLogger(__name__)


def read_groups(*, path: Path, identifiers: Sequence[str]) -> list[str]:
    """Read related-sequence groups for leakage-resistant cross-validation.

    Args:
        path: TSV with sequence_id and group columns.
        identifiers: Analysed FASTA IDs in required order.

    Returns:
        A group for every analysed sequence.

    Raises:
        ValueError: IDs duplicate or analysed sequences lack a group.
    """
    groups: dict[str, str] = {}
    with open_text(path=path) as stream:
        reader = csv.DictReader(f=stream, delimiter="\t")
        if not {"sequence_id", "group"} <= set(reader.fieldnames or []):
            raise ValueError("Groups TSV requires sequence_id and group")
        for row in reader:
            identifier, group = row.get("sequence_id"), row.get("group")
            if not identifier or not group or identifier in groups:
                raise ValueError("Malformed or duplicated group assignment")
            groups[identifier] = group
    if set(identifiers) - groups.keys():
        raise ValueError("Every analysed sequence needs an explicit group")
    return [groups[identifier] for identifier in identifiers]


def composition_matrix(*, sequences: Sequence[str]) -> NDArray[np.float64]:
    """Build GC, log-length and ambiguity covariates for the baseline.

    Args:
        sequences: DNA sequences.

    Returns:
        A three-column dense matrix.
    """
    rows: list[list[float]] = []
    for sequence in sequences:
        composition = sequence_composition(sequence=sequence)
        rows.append(
            [
                float(composition["gc_fraction"]),
                math.log1p(len(sequence)),
                float(composition["ambiguous_fraction"]),
            ]
        )
    return np.asarray(rows, dtype=float)


def fit_vocabulary(
    *, counters: Sequence[Counter[str]], max_features: int = 100000
) -> dict[str, int]:
    """Select a vocabulary using training sequences only, without labels.

    Args:
        counters: Training sequence k-mer counts.
        max_features: Maximum number of features.

    Returns:
        Deterministic feature indices, ranked by document frequency.

    Raises:
        ValueError: No valid words exist or the feature limit is invalid.
    """
    if max_features < 1:
        raise ValueError("max_features must be positive")
    frequency: Counter[str] = Counter()
    for counts in counters:
        frequency.update(counts.keys())
    selected = sorted(frequency, key=lambda word: (-frequency[word], word))[
        :max_features
    ]
    if not selected:
        raise ValueError("No valid k-mers for learning")
    return {word: index for index, word in enumerate(selected)}


def word_matrix(
    *, counters: Sequence[Counter[str]], vocabulary: dict[str, int]
) -> csr_matrix:
    """Encode k-mer frequencies as percentages of all valid word windows.

    Args:
        counters: Per-sequence counts, including words outside the vocabulary.
        vocabulary: Training-only vocabulary.

    Returns:
        Sparse features. Novel test words do not alter the denominator.
    """
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for index, counts in enumerate(counters):
        total = sum(counts.values())
        for word, count in counts.items():
            if word in vocabulary and total:
                rows.append(index)
                columns.append(vocabulary[word])
                values.append(100 * count / total)
    return csr_matrix(
        (values, (rows, columns)),
        shape=(len(counters), len(vocabulary)),
        dtype=float,
    )


def validation_splits(
    *,
    labels: NDArray[np.int64],
    folds: int,
    seed: int,
    groups: Sequence[str] | None = None,
) -> list[tuple[NDArray[np.int64], NDArray[np.int64]]]:
    """Create deterministic splits and reject invalid grouped folds.

    Args:
        labels: Binary sequence labels.
        folds: Number of held-out folds.
        seed: Random seed.
        groups: Optional chromosome, family or homology-cluster groups.

    Returns:
        Train/test indices for every fold.

    Raises:
        ValueError: There are insufficient independent observations or a
            fold lacks either class. Grouping is never silently relaxed.
    """
    if (
        folds < 2
        or set(labels.tolist()) != {0, 1}
        or min(np.bincount(labels)) < folds
    ):
        raise ValueError(
            "At least one sequence per class per fold is required"
        )
    if groups is not None and (
        len(groups) != len(labels) or len(set(groups)) < folds
    ):
        raise ValueError("Insufficient or malformed independent groups")
    x = np.zeros((len(labels), 1))
    if groups is None:
        splitter = StratifiedKFold(
            n_splits=folds, shuffle=True, random_state=seed
        )
        splits = list(splitter.split(X=x, y=labels))
    else:
        grouped_splitter = StratifiedGroupKFold(
            n_splits=folds, shuffle=True, random_state=seed
        )
        splits = list(
            grouped_splitter.split(X=x, y=labels, groups=np.asarray(groups))
        )
    for train, test in splits:
        if (
            len(np.unique(labels[train])) != 2
            or len(np.unique(labels[test])) != 2
        ):
            raise ValueError(
                "A validation fold lacks a class; use fewer folds "
                "or revise independent groups"
            )
        if groups is not None and {groups[int(i)] for i in train} & {
            groups[int(i)] for i in test
        }:
            raise ValueError("Validation groups overlap")
    return splits


def fit_classifier(
    *,
    matrix: Any,
    labels: NDArray[np.int64],
    seed: int,
    regularisation: float = 1.0,
    sample_weight: NDArray[np.float64] | None = None,
) -> LogisticRegression:
    """Fit a fixed, regularised, class-balanced logistic sequence model.

    Args:
        matrix: Training feature matrix.
        labels: Training labels.
        seed: Random seed.
        regularisation: Inverse L2 regularisation strength C.
        sample_weight: Optional positive weights including class balancing.
            Supplied weights replace automatic row-level class balancing.

    Returns:
        A converged classifier.

    Raises:
        ValueError: Regularisation is invalid or optimisation fails.
    """
    if regularisation <= 0 or not math.isfinite(regularisation):
        raise ValueError("Regularisation C must be finite and positive")
    if sample_weight is not None and (
        sample_weight.shape != labels.shape
        or not np.isfinite(sample_weight).all()
        or (sample_weight <= 0).any()
    ):
        raise ValueError(
            "Classifier sample weights must be finite and positive"
        )
    model = LogisticRegression(
        C=regularisation,
        class_weight="balanced" if sample_weight is None else None,
        solver="liblinear",
        max_iter=2000,
        random_state=seed,
    )
    model.fit(X=matrix, y=labels, sample_weight=sample_weight)
    if int(np.max(model.n_iter_)) >= model.max_iter:
        raise ValueError("Sequence classifier failed to converge")
    return model


def cross_validate_sequences(
    *,
    counters: Sequence[Counter[str]],
    composition: NDArray[np.float64],
    labels: NDArray[np.int64],
    folds: int,
    seed: int,
    groups: Sequence[str] | None = None,
    max_features: int = 100000,
    regularisation: float = 1.0,
    explanation_contexts: list[dict[str, Any]] | None = None,
    explanation_indices: NDArray[np.int64] | None = None,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.int64],
    list[dict[str, Any]],
    dict[str, list[float]],
]:
    """Generate held-out predictions with all preprocessing inside each fold.

    Args:
        counters: Per-sequence word counts.
        composition: Baseline covariates.
        labels: Binary target labels.
        folds: Validation folds.
        seed: Random seed.
        groups: Optional independent groups.
        max_features: Vocabulary limit fitted separately in each fold.
        regularisation: Fixed model C; no validation-driven tuning.
        explanation_contexts: Optional destination for training-only SHAP
            contexts; omitted for permutation runs.
        explanation_indices: Optional original row indices to explain.

    Returns:
        Sequence-model scores, baseline scores, fold IDs, fold metrics and
        feature coefficients observed across training folds.
    """
    splits = validation_splits(
        labels=labels, folds=folds, seed=seed, groups=groups
    )
    scores = np.empty(len(labels), dtype=float)
    baseline_scores = np.empty(len(labels), dtype=float)
    assignments = np.empty(len(labels), dtype=np.int64)
    metrics: list[dict[str, Any]] = []
    coefficients: dict[str, list[float]] = defaultdict(list)
    for fold, (train, test) in enumerate(splits, start=1):
        vocabulary = fit_vocabulary(
            counters=[counters[int(i)] for i in train],
            max_features=max_features,
        )
        words = word_matrix(counters=counters, vocabulary=vocabulary)
        scaler = StandardScaler().fit(X=composition[train])
        scaled = scaler.transform(X=composition)
        combined = hstack([words, csr_matrix(scaled)], format="csr")
        model = fit_classifier(
            matrix=combined[train],
            labels=labels[train],
            seed=seed,
            regularisation=regularisation,
        )
        baseline = fit_classifier(
            matrix=scaled[train],
            labels=labels[train],
            seed=seed,
            regularisation=regularisation,
        )
        scores[test] = model.predict_proba(X=combined[test])[:, 1]
        baseline_scores[test] = baseline.predict_proba(X=scaled[test])[:, 1]
        assignments[test] = fold
        if explanation_contexts is not None:
            explained = (
                test[np.isin(test, explanation_indices)]
                if explanation_indices is not None
                else test
            )
            if len(explained):
                explanation_contexts.append(
                    {
                        "matrix": combined[explained],
                        "coefficients": model.coef_[0].copy(),
                        "intercept": float(model.intercept_[0]),
                        "background_mean": np.asarray(
                            combined[train].mean(axis=0)
                        ).ravel(),
                        "feature_names": [*vocabulary, *COMPOSITION_FEATURES],
                        "indices": explained,
                        "fold": fold,
                    }
                )
        for word, index in vocabulary.items():
            coefficients[word].append(float(model.coef_[0, index]))
        metrics.append(
            {
                "fold": fold,
                "train_sequences": len(train),
                "test_sequences": len(test),
                "vocabulary_size": len(vocabulary),
                "roc_auc": float(
                    roc_auc_score(y_true=labels[test], y_score=scores[test])
                ),
                "average_precision": float(
                    average_precision_score(
                        y_true=labels[test], y_score=scores[test]
                    )
                ),
                "baseline_roc_auc": float(
                    roc_auc_score(
                        y_true=labels[test], y_score=baseline_scores[test]
                    )
                ),
                "baseline_average_precision": float(
                    average_precision_score(
                        y_true=labels[test], y_score=baseline_scores[test]
                    )
                ),
            }
        )
    return scores, baseline_scores, assignments, metrics, dict(coefficients)


def permute_labels(
    *,
    labels: NDArray[np.int64],
    rng: np.random.Generator,
    groups: Sequence[str] | None = None,
) -> NDArray[np.int64]:
    """Permute labels respecting documented group exchangeability blocks.

    Args:
        labels: Binary target labels.
        rng: Seeded random generator.
        groups: Optional related-sequence groups. Mixed groups permute within
            group; pure groups exchange labels only with equal-sized groups.

    Returns:
        A permutation preserving class totals and pure-group labels.
    """
    if groups is None:
        return rng.permutation(x=labels)
    output = labels.copy()
    members: dict[str, list[int]] = defaultdict(list)
    for index, group in enumerate(groups):
        members[group].append(index)
    pure: dict[int, list[list[int]]] = defaultdict(list)
    for indices in members.values():
        if len(np.unique(labels[indices])) == 1:
            pure[len(indices)].append(indices)
        else:
            output[indices] = rng.permutation(x=labels[indices])
    for block in pure.values():
        shuffled = rng.permutation(
            x=np.asarray([labels[indices[0]] for indices in block])
        )
        for indices, label in zip(block, shuffled, strict=True):
            output[indices] = label
    return output


def fit_sequence_model(
    *,
    positive: dict[str, str],
    negative: dict[str, str],
    lengths: Sequence[int] = (4, 5, 6),
    folds: int = 5,
    seed: int = 42,
    permutations: int = 99,
    groups: Sequence[str] | None = None,
    max_features: int = 100000,
    regularisation: float = 1.0,
    shap: bool = True,
    shap_max_sequences: int = 1000,
    shap_max_features: int = 20,
    explanation_outputs: dict[str, Any] | None = None,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, Any],
]:
    """Learn candidate regulatory signatures and assess predictive evidence.

    Args:
        positive: Foreground sequences.
        negative: Control sequences.
        lengths: Strand-canonical word lengths.
        folds: Held-out validation folds.
        seed: Reproducibility seed.
        permutations: Label-permutation replicates; zero disables the test.
        groups: Optional related-sequence or chromosome groups in FASTA order.
        max_features: Vocabulary limit.
        regularisation: Fixed inverse regularisation strength C.
        shap: Calculate held-out SHAP automatically unless disabled.
        shap_max_sequences: Reproducible explanation sample limit.
        shap_max_features: Leading features in long explanations.
        explanation_outputs: Optional destination for SHAP tables, separate
            from the reusable model and its concise summary.

    Returns:
        Held-out predictions, signature coefficients, fold metrics, summary
        and a reusable, strict-JSON model (no executable pickle data).

    Raises:
        ValueError: Sequence sets, sample counts or settings are invalid.
    """
    check_sequence_sets(positive=positive, negative=negative)
    if min(len(positive), len(negative)) < 5:
        raise ValueError(
            "AI analysis requires at least five independent "
            "sequences per class"
        )
    kmer_family_size(lengths=lengths)
    if not lengths or permutations < 0:
        raise ValueError(
            "AI lengths must be supplied; permutations cannot be negative"
        )
    identifiers = [*positive, *negative]
    sequences = [*positive.values(), *negative.values()]
    labels = np.asarray(
        [1] * len(positive) + [0] * len(negative), dtype=np.int64
    )
    selected_explanations = explanation_sample(
        labels=labels, maximum=shap_max_sequences, seed=seed
    )
    if (
        not isinstance(shap_max_features, int)
        or isinstance(shap_max_features, bool)
        or not 1 <= shap_max_features <= 100
    ):
        raise ValueError("SHAP max_features must be between one and 100")
    contexts: list[dict[str, Any]] = []
    counters = [kmer_counts(sequence=s, lengths=lengths) for s in sequences]
    composition = composition_matrix(sequences=sequences)
    scores, baseline_scores, assignments, metrics, coefficients = (
        cross_validate_sequences(
            counters=counters,
            composition=composition,
            labels=labels,
            folds=folds,
            seed=seed,
            groups=groups,
            max_features=max_features,
            regularisation=regularisation,
            explanation_contexts=contexts if shap else None,
            explanation_indices=selected_explanations,
        )
    )
    average_precision = float(
        average_precision_score(y_true=labels, y_score=scores)
    )
    null_scores: list[float] = []
    rng = np.random.default_rng(seed=seed)
    for replicate in range(permutations):
        permuted = permute_labels(labels=labels, rng=rng, groups=groups)
        if groups is not None and np.array_equal(permuted, labels):
            # An unchanged permutation is valid when exchanges are possible.
            exchangeable = any(
                not np.array_equal(
                    permute_labels(labels=labels, rng=rng, groups=groups),
                    labels,
                )
                for _ in range(20)
            )
            if not exchangeable:
                LOGGER.warning(
                    "No non-trivial grouped label exchanges; "
                    "permutation test unavailable"
                )
                null_scores = []
                break
        try:
            permuted_scores, _, _, _, _ = cross_validate_sequences(
                counters=counters,
                composition=composition,
                labels=permuted,
                folds=folds,
                seed=seed,
                groups=groups,
                max_features=max_features,
                regularisation=regularisation,
            )
        except ValueError as exc:
            raise ValueError(
                "Permutation validation failed; use fewer folds "
                "or disable permutations explicitly"
            ) from exc
        null_scores.append(
            float(
                average_precision_score(
                    y_true=permuted, y_score=permuted_scores
                )
            )
        )
        LOGGER.info("AI permutation %d/%d", replicate + 1, permutations)
    predictions = [
        {
            "sequence_id": identifier,
            "label": int(labels[i]),
            "fold": int(assignments[i]),
            "held_out_signature_score": float(scores[i]),
            "held_out_composition_score": float(baseline_scores[i]),
            "group": groups[i] if groups is not None else identifier,
        }
        for i, identifier in enumerate(identifiers)
    ]
    vocabulary = fit_vocabulary(counters=counters, max_features=max_features)
    words = word_matrix(counters=counters, vocabulary=vocabulary)
    scaler = StandardScaler().fit(X=composition)
    matrix = hstack(
        [words, csr_matrix(scaler.transform(X=composition))], format="csr"
    )
    model = fit_classifier(
        matrix=matrix, labels=labels, seed=seed, regularisation=regularisation
    )
    signatures: list[dict[str, Any]] = []
    for word, index in vocabulary.items():
        observed = coefficients.get(word, [])
        signatures.append(
            {
                "kmer": word,
                "coefficient": float(model.coef_[0, index]),
                "mean_fold_coefficient": float(np.mean(observed))
                if observed
                else 0.0,
                "folds_present": len(observed),
                "positive_fold_fraction": sum(c > 0 for c in observed)
                / len(observed)
                if observed
                else 0.0,
            }
        )
    signatures.sort(key=lambda r: (-r["coefficient"], r["kmer"]))
    summary: dict[str, Any] = {
        "roc_auc": float(roc_auc_score(y_true=labels, y_score=scores)),
        "average_precision": average_precision,
        "baseline_roc_auc": float(
            roc_auc_score(y_true=labels, y_score=baseline_scores)
        ),
        "baseline_average_precision": float(
            average_precision_score(y_true=labels, y_score=baseline_scores)
        ),
        "positive_prevalence": len(positive) / len(labels),
        "permutation_p_value": (
            1 + sum(s >= average_precision for s in null_scores)
        )
        / (1 + len(null_scores))
        if null_scores
        else None,
        "permutation_scores": null_scores,
        "permutations_completed": len(null_scores),
        "grouped_validation": groups is not None,
        "folds": folds,
        "seed": seed,
        "model": (
            "L2 logistic regression on canonical k-mer frequency "
            "plus composition"
        ),
        "interpretation": (
            "Scores measure resemblance to the supplied positive class; "
            "they are not calibrated probabilities of enhancer function."
        ),
        "warnings": [
            "Ungrouped validation does not protect against paralogy or "
            "related loci; supply groups for those data."
        ]
        if groups is None
        else [],
    }
    if shap:
        explanations, importance, shap_summary = summarise_shap(
            contexts=contexts,
            identifiers=identifiers,
            labels=labels.tolist(),
            max_features=shap_max_features,
        )
        summary["shap"] = shap_summary
        if explanation_outputs is not None:
            explanation_outputs.update(
                rows=explanations, importance=importance, summary=shap_summary
            )
    else:
        summary["shap"] = {
            "status": "disabled",
            "reason": "User disabled SHAP",
        }
    exported = {
        "schema_version": 1,
        "lengths": list(lengths),
        "vocabulary": vocabulary,
        "coefficients": model.coef_[0].tolist(),
        "intercept": float(model.intercept_[0]),
        "composition_mean": scaler.mean_.tolist(),
        "composition_scale": scaler.scale_.tolist(),
        "feature_background_mean": np.asarray(matrix.mean(axis=0))
        .ravel()
        .tolist(),
        "training_sequence_sha256": [
            hashlib.sha256(
                canonical_sequence(sequence=s).encode("ascii")
            ).hexdigest()
            for s in sequences
        ],
        "summary": summary,
    }
    return predictions, signatures, metrics, summary, exported


def predict_sequences(
    *, model: dict[str, Any], sequences: dict[str, str]
) -> list[dict[str, Any]]:
    """Score new sequences with a validated, inert JSON model.

    Args:
        model: Exported sequence model.
        sequences: Query FASTA collection.

    Returns:
        Candidate signature scores and a flag for exact training duplicates.

    Raises:
        ValueError: The model schema or numerical parameters are invalid.
    """
    try:
        vocabulary = model["vocabulary"]
        lengths = model["lengths"]
        coefficients = np.asarray(model["coefficients"], dtype=float)
        means = np.asarray(model["composition_mean"], dtype=float)
        scales = np.asarray(model["composition_scale"], dtype=float)
        intercept = float(model["intercept"])
        if (
            model["schema_version"] != 1
            or not isinstance(vocabulary, dict)
            or set(vocabulary.values()) != set(range(len(vocabulary)))
            or any(
                not isinstance(word, str) or not set(word) <= set("ACGT")
                for word in vocabulary
            )
        ):
            raise ValueError("Invalid vocabulary/model schema")
        kmer_family_size(lengths=lengths)
        if not lengths or any(len(word) not in lengths for word in vocabulary):
            raise ValueError("Vocabulary words have unsupported lengths")
        if (
            coefficients.shape != (len(vocabulary) + 3,)
            or means.shape != (3,)
            or scales.shape != (3,)
            or not np.isfinite(coefficients).all()
            or not np.isfinite(means).all()
            or not np.isfinite(scales).all()
            or (scales <= 0).any()
            or not math.isfinite(intercept)
        ):
            raise ValueError(
                "Invalid model coefficients or scaling parameters"
            )
    except (KeyError, TypeError, OverflowError) as exc:
        raise ValueError("Malformed sequence-model JSON") from exc
    width = model.get("window_length")
    if width is not None and (
        not isinstance(width, int)
        or isinstance(width, bool)
        or not 20 <= width <= 100000
        or any(len(sequence) != width for sequence in sequences.values())
    ):
        raise ValueError(
            "Window models require complete sequences at their training length"
        )
    counters = [
        kmer_counts(sequence=s, lengths=lengths) for s in sequences.values()
    ]
    words = word_matrix(counters=counters, vocabulary=vocabulary)
    composition = (
        composition_matrix(sequences=list(sequences.values())) - means
    ) / scales
    logits = (
        np.asarray(words @ coefficients[: len(vocabulary)]).ravel()
        + composition @ coefficients[len(vocabulary) :]
        + intercept
    )
    training = set(model.get("training_sequence_sha256", []))
    return [
        {
            "sequence_id": identifier,
            "candidate_signature_score": float(expit(logits[i])),
            "matches_training_sequence": hashlib.sha256(
                canonical_sequence(sequence=sequence).encode("ascii")
            ).hexdigest()
            in training,
            "interpretation": (
                "positive-class sequence signature; enhancer function "
                "requires evidence"
            ),
        }
        for i, (identifier, sequence) in enumerate(sequences.items())
    ]
