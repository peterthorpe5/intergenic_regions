"""Exact SHAP, held-out backgrounds and official plotting regressions."""

from collections import defaultdict

import numpy as np
import pytest
import shap
from scipy.sparse import csr_matrix
from scipy.special import expit

from intergenic_regions.explanations import (
    explanation_sample,
    linear_shap_values,
    summarise_shap,
)
from intergenic_regions.io import read_fasta
from intergenic_regions.learning import fit_sequence_model
from intergenic_regions.shap_reporting import plot_shap, shap_explanation
from intergenic_regions.workflows import learning_outputs


@pytest.mark.parametrize("positive_count", [1, 2, 10, 19])
def test_reproducible_stratified_subsampling(positive_count):
    labels = np.asarray([1] * positive_count + [0] * (20 - positive_count))
    sample = explanation_sample(labels=labels, maximum=7, seed=3)
    assert len(sample) == 7
    assert set(labels[sample]) == {0, 1}
    assert sample.tolist() == sorted(set(sample.tolist()))
    np.testing.assert_array_equal(
        sample, explanation_sample(labels=labels, maximum=7, seed=3)
    )
    assert len(explanation_sample(labels=labels, maximum=20, seed=3)) == 20


@pytest.mark.parametrize(
    "labels,maximum",
    [
        ([1, 1], 10),
        ([0, 1], 1),
        ([0, 1], True),
        ([[0, 1]], 5),
    ],
)
def test_sample_invalid(labels, maximum):
    with pytest.raises(ValueError):
        explanation_sample(labels=np.asarray(labels), maximum=maximum, seed=1)


def test_exact_values_agree_with_official_linear_explainer():
    background = np.asarray([[1, 2, 3], [4, 2, 1], [0, 4, 7]], dtype=float)
    query = np.asarray([[2, 4, 1], [1, 0, 3]], dtype=float)
    coefficients = np.asarray([0.3, -0.7, 0.1])
    intercept = -0.4
    values, baseline = linear_shap_values(
        matrix=csr_matrix(query),
        coefficients=coefficients,
        background_mean=background.mean(axis=0),
        intercept=intercept,
    )
    official = shap.LinearExplainer(
        model=(coefficients, intercept),
        masker=shap.maskers.Independent(
            data=background, max_samples=len(background)
        ),
    )(query)
    np.testing.assert_allclose(values, official.values, atol=1e-12)
    np.testing.assert_allclose(baseline, official.base_values, atol=1e-12)
    np.testing.assert_allclose(
        values.sum(axis=1) + baseline,
        query @ coefficients + intercept,
    )
    dense, _ = linear_shap_values(
        matrix=query,
        coefficients=coefficients,
        background_mean=background.mean(axis=0),
        intercept=intercept,
    )
    np.testing.assert_array_equal(dense, values)


@pytest.mark.parametrize(
    "changes",
    [
        {"matrix": [1, 2]},
        {"coefficients": np.ones((1, 2))},
        {"background_mean": np.zeros(3)},
        {"matrix": [[np.nan, 1]]},
        {"coefficients": [np.inf, 1]},
        {"background_mean": [0, np.nan]},
        {"intercept": np.inf},
    ],
)
def test_linear_validation(changes):
    options = dict(
        matrix=[[1, 2]],
        coefficients=np.ones(2),
        background_mean=np.zeros(2),
        intercept=0,
    )
    with pytest.raises(ValueError):
        linear_shap_values(**{**options, **changes})


@pytest.fixture
def contexts():
    return [
        dict(
            matrix=csr_matrix([[1, 2], [3, 0]]),
            coefficients=np.asarray([1.0, -0.5]),
            background_mean=np.asarray([0.5, 1]),
            intercept=0.2,
            feature_names=["AAAA", "gc_fraction"],
            indices=[0, 1],
            fold=1,
        ),
        dict(
            matrix=csr_matrix([[2, 4]]),
            coefficients=np.asarray([-0.2, 0.5]),
            background_mean=np.asarray([1.0, 2]),
            intercept=-0.3,
            feature_names=["CCCC", "gc_fraction"],
            indices=[2],
            fold=2,
        ),
    ]


def test_fold_specific_features_and_additive_remainder(contexts):
    rows, importance, summary = summarise_shap(
        contexts=contexts,
        identifiers=["a", "b", "c"],
        labels=[1, 0, 1],
        max_features=1,
    )
    assert summary["features_assessed"] == 3
    assert summary["scope"] == "held_out"
    assert summary["output_scale"] == "positive_class_log_odds"
    assert importance[0]["feature"] == "AAAA"
    assert importance[0]["mean_absolute_shap"] == pytest.approx(1.0)
    unavailable = next(
        r for r in rows if r["sequence_id"] == "c" and r["feature"] == "AAAA"
    )
    assert unavailable["feature_value"] is None
    assert unavailable["shap_value"] == 0
    assert unavailable["feature_available"] is False
    explained = shap_explanation(rows=rows)
    np.testing.assert_allclose(
        explained.values.sum(axis=1) + explained.base_values, [0.2, 3.2, 1.3]
    )
    assert all("enhancer" in summary["interpretation"] for _ in [0])


@pytest.mark.parametrize(
    "changes",
    [
        {"max_features": 0},
        {"max_features": True},
        {"scope": "validated"},
        {"identifiers": ["a", "a", "c"]},
        {"labels": [1]},
        {"contexts": []},
    ],
)
def test_summary_invalid(contexts, changes):
    options = dict(contexts=contexts, identifiers=["a", "b", "c"])
    with pytest.raises(ValueError):
        summarise_shap(**{**options, **changes})


def test_context_shape_and_repeated_indices(contexts):
    for changes in (
        {"indices": [0, 0]},
        {"indices": [1, 3]},
        {"feature_names": ["A", "A"]},
        {"indices": [0]},
    ):
        invalid = [{**contexts[0], **changes}, contexts[1]]
        with pytest.raises(ValueError):
            summarise_shap(contexts=invalid, identifiers=["a", "b", "c"])
    with pytest.raises(ValueError):
        summarise_shap(
            contexts=[contexts[0], contexts[0]], identifiers=["a", "b", "c"]
        )


def test_explanation_rejects_bad_rows(contexts):
    rows, _, _ = summarise_shap(contexts=contexts, identifiers=["a", "b", "c"])
    for invalid in (
        [],
        [*rows, rows[0]],
        rows[1:],
        [{**r, "shap_value": np.nan} for r in rows],
        [{**r, "model_logit": 100} for r in rows],
    ):
        with pytest.raises(ValueError):
            shap_explanation(rows=invalid)
    with pytest.raises(ValueError):
        plot_shap(rows=rows, importance=[], directory=None)


def test_official_plots_are_png_and_pdf(contexts, tmp_path):
    rows, importance, _ = summarise_shap(
        contexts=contexts,
        identifiers=["a", "b", "c"],
        max_features=3,
    )
    paths = plot_shap(rows=rows, importance=importance, directory=tmp_path)
    assert {p.stem for p in paths} == {
        "shap_global_bar",
        "shap_beeswarm",
        "shap_heatmap",
        "shap_waterfall_positive",
        "shap_waterfall_negative",
        "shap_dependence",
    }
    for path in paths:
        assert path.read_bytes().startswith(b"\x89PNG")
        assert path.with_suffix(".pdf").read_bytes().startswith(b"%PDF")


def test_held_out_explanations_reconstruct_predictions(dataset, monkeypatch):
    positive = read_fasta(path=dataset["positive"])
    negative = read_fasta(path=dataset["negative"])
    outputs = {}
    from intergenic_regions import learning

    original = learning.summarise_shap
    observed = []

    def monitored(**kwargs):
        observed.extend(kwargs["contexts"])
        return original(**kwargs)

    monkeypatch.setattr(learning, "summarise_shap", monitored)
    predictions, _, _, summary, model = fit_sequence_model(
        positive=positive,
        negative=negative,
        lengths=[4],
        folds=3,
        permutations=0,
        shap_max_sequences=12,
        shap_max_features=3,
        explanation_outputs=outputs,
    )
    totals = defaultdict(float)
    first = {}
    for row in outputs["rows"]:
        totals[row["sequence_id"]] += row["shap_value"]
        first[row["sequence_id"]] = row
    assert len(totals) == 12
    assert summary["shap"]["subsampled"] is True
    for prediction in predictions:
        identifier = prediction["sequence_id"]
        if identifier in totals:
            assert expit(
                totals[identifier] + first[identifier]["base_value"]
            ) == pytest.approx(prediction["held_out_signature_score"])
            assert first[identifier]["fold"] == prediction["fold"]
    # Training-scaled composition columns have zero means within each fold.
    for context in observed:
        np.testing.assert_allclose(
            context["background_mean"][-3:], 0, atol=1e-12
        )
    assert len(model["feature_background_mean"]) == len(model["coefficients"])
    _, _, _, disabled, _ = fit_sequence_model(
        positive=positive,
        negative=negative,
        lengths=[4],
        folds=3,
        permutations=0,
        shap=False,
    )
    assert disabled["shap"]["status"] == "disabled"
    for changes in (
        {"shap_max_features": True},
        {"shap_max_features": 0},
        {"shap_max_sequences": 1},
    ):
        with pytest.raises(ValueError):
            fit_sequence_model(
                positive=positive, negative=negative, permutations=0, **changes
            )


def test_missing_shap_plots_keep_model_and_numerical_values(
    dataset,
    tmp_path,
    monkeypatch,
):
    def unavailable(**kwargs):
        raise ImportError("No shap")

    monkeypatch.setattr(
        "intergenic_regions.shap_reporting.plot_shap", unavailable
    )
    summary, _ = learning_outputs(
        directory=tmp_path,
        positive=read_fasta(path=dataset["positive"]),
        negative=read_fasta(path=dataset["negative"]),
        settings={"lengths": [6], "folds": 3, "permutations": 0},
    )
    assert summary["shap"]["plot_status"] == "unavailable"
    assert (tmp_path / "model.json").is_file()
    assert (tmp_path / "shap_values.tsv").is_file()
