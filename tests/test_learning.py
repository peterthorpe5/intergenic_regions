"""Held-out sequence learning, group integrity, model export and null tests."""

from collections import Counter

import numpy as np
import pytest

from intergenic_regions.io import read_fasta
from intergenic_regions.learning import (
    composition_matrix,
    cross_validate_sequences,
    fit_classifier,
    fit_sequence_model,
    fit_vocabulary,
    permute_labels,
    predict_sequences,
    read_groups,
    validation_splits,
    word_matrix,
)
from intergenic_regions.motifs import kmer_counts


def test_feature_vocabulary_frequency_and_composition():
    counters = [Counter({"AA": 2, "AC": 1}), Counter({"AC": 3, "GG": 1})]
    vocabulary = fit_vocabulary(counters=counters, max_features=2)
    assert vocabulary == {"AC": 0, "AA": 1}
    matrix = word_matrix(counters=counters, vocabulary=vocabulary).toarray()
    assert matrix[0].tolist() == pytest.approx([100 / 3, 200 / 3])
    assert matrix[1].tolist() == [75, 0]
    assert word_matrix(counters=[Counter()], vocabulary=vocabulary).sum() == 0
    composition = composition_matrix(sequences=["ACGN"])
    assert composition[0].tolist() == pytest.approx([2 / 3, np.log(5), 0.25])
    for counters, limit in [([], 2), ([Counter({"A": 1})], 0)]:
        with pytest.raises(ValueError):
            fit_vocabulary(counters=counters, max_features=limit)


def test_group_input_completeness(tmp_path):
    path = tmp_path / "groups.tsv"
    path.write_text("sequence_id\tgroup\na\tfamily1\nb\tfamily1\nc\tfamily2\n")
    assert read_groups(path=path, identifiers=["c", "a"]) == [
        "family2",
        "family1",
    ]
    with pytest.raises(ValueError):
        read_groups(path=path, identifiers=["missing"])
    for text in (
        "bad",
        "sequence_id\tgroup\na\tx\na\ty\n",
        "sequence_id\tgroup\na\t\n",
    ):
        path.write_text(text)
        with pytest.raises(ValueError):
            read_groups(path=path, identifiers=["a"])


def test_grouped_validation_no_family_crosses_folds():
    labels = np.asarray([1, 0] * 10, dtype=np.int64)
    groups = [f"family{i // 2}" for i in range(20)]
    splits = validation_splits(labels=labels, folds=5, seed=42, groups=groups)
    assert len(splits) == 5
    assert sorted(
        np.concatenate([test for _, test in splits]).tolist()
    ) == list(range(20))
    for train, test in splits:
        assert not {groups[i] for i in train} & {groups[i] for i in test}
    for kwargs in (
        {"folds": 1},
        {"folds": 11},
        {"groups": ["single"] * 20},
        {"groups": ["bad"]},
    ):
        options = dict(labels=labels, folds=5, seed=42)
        options.update(kwargs)
        with pytest.raises(ValueError):
            validation_splits(**options)
    with pytest.raises(ValueError, match="lacks a class"):
        validation_splits(
            labels=np.asarray([1] * 5 + [0] * 5),
            folds=2,
            seed=42,
            groups=["p"] * 5 + ["n"] * 5,
        )


def test_permutations_preserve_exchangeability_blocks():
    labels = np.asarray([1, 1, 0, 0, 1, 0], dtype=np.int64)
    groups = ["a", "a", "b", "b", "c", "c"]
    permuted = permute_labels(
        labels=labels, groups=groups, rng=np.random.default_rng(seed=1)
    )
    assert permuted.sum() == labels.sum()
    assert permuted[0] == permuted[1] and permuted[2] == permuted[3]
    assert set(
        permute_labels(labels=labels, rng=np.random.default_rng(seed=3))
    ) == {0, 1}


def test_classifier_errors():
    with pytest.raises(ValueError):
        fit_classifier(
            matrix=np.ones((4, 2)),
            labels=np.asarray([0, 1, 0, 1]),
            seed=1,
            regularisation=0,
        )


def test_training_only_vocabulary_and_fold_scaling(dataset, monkeypatch):
    positive = read_fasta(path=dataset["positive"])
    negative = read_fasta(path=dataset["negative"])
    sequences = [*positive.values(), *negative.values()]
    counters = [kmer_counts(sequence=s, lengths=[4]) for s in sequences]
    observed_training_sizes = []
    original = fit_vocabulary

    def monitored(*, counters, max_features):
        observed_training_sizes.append(len(counters))
        return original(counters=counters, max_features=max_features)

    monkeypatch.setattr(
        "intergenic_regions.learning.fit_vocabulary", monitored
    )
    scores, baseline, folds, metrics, coefficients = cross_validate_sequences(
        counters=counters,
        composition=composition_matrix(sequences=sequences),
        labels=np.asarray([1] * 15 + [0] * 15),
        folds=3,
        seed=42,
        max_features=100,
    )
    assert observed_training_sizes == [20, 20, 20]
    assert len(scores) == len(baseline) == len(folds) == 30
    assert len(metrics) == 3 and coefficients
    assert np.isfinite(scores).all()


@pytest.fixture
def fitted(dataset):
    positive = read_fasta(path=dataset["positive"])
    negative = read_fasta(path=dataset["negative"])
    return fit_sequence_model(
        positive=positive,
        negative=negative,
        lengths=[4, 5, 6],
        folds=3,
        permutations=3,
        seed=42,
        max_features=1000,
    )


def test_planted_signature_and_reusable_model(dataset, fitted):
    predictions, signatures, metrics, summary, model = fitted
    assert len(predictions) == 30 and len(metrics) == 3
    assert summary["roc_auc"] > 0.85
    assert summary["average_precision"] > summary["baseline_average_precision"]
    assert summary["permutation_p_value"] >= 1 / 4
    assert summary["permutations_completed"] == 3
    assert any(
        r["kmer"] == "CACGTG" and r["coefficient"] > 0 for r in signatures
    )
    positive = read_fasta(path=dataset["positive"])
    results = predict_sequences(
        model=model, sequences={"new": "ACGT" * 40, "training": positive["p0"]}
    )
    assert 0 <= results[0]["candidate_signature_score"] <= 1
    assert not results[0]["matches_training_sequence"]
    assert results[1]["matches_training_sequence"]
    assert all("enhancer" in r["interpretation"] for r in results)


def test_grouped_learning_and_disabled_permutations(dataset):
    positive = read_fasta(path=dataset["positive"])
    negative = read_fasta(path=dataset["negative"])
    groups = [f"family{i}" for i in range(15)] * 2
    _, _, _, summary, _ = fit_sequence_model(
        positive=positive,
        negative=negative,
        groups=groups,
        folds=3,
        permutations=0,
        lengths=[4],
        max_features=100,
    )
    assert summary["grouped_validation"]
    assert summary["permutation_p_value"] is None
    assert summary["warnings"] == []


@pytest.mark.parametrize(
    "changes",
    [
        {"schema_version": 2},
        {"coefficients": [1]},
        {"composition_scale": [0, 1, 1]},
        {"composition_mean": [float("nan"), 0, 0]},
        {"intercept": float("inf")},
        {"lengths": []},
        {"vocabulary": {"bad": 0}},
        {"vocabulary": {"AA": 1}},
    ],
)
def test_model_validation(fitted, changes):
    model = dict(fitted[-1])
    model.update(changes)
    with pytest.raises(ValueError):
        predict_sequences(model=model, sequences={"x": "AAAA"})
    with pytest.raises(ValueError):
        predict_sequences(model={}, sequences={"x": "AAAA"})


def test_small_or_invalid_learning_data(dataset):
    positive = read_fasta(path=dataset["positive"])
    negative = read_fasta(path=dataset["negative"])
    with pytest.raises(ValueError, match="five independent"):
        fit_sequence_model(
            positive=dict(list(positive.items())[:4]), negative=negative
        )
    with pytest.raises(ValueError):
        fit_sequence_model(
            positive=positive, negative=negative, lengths=[], permutations=0
        )
    with pytest.raises(ValueError):
        fit_sequence_model(
            positive=positive, negative=negative, permutations=-1
        )
