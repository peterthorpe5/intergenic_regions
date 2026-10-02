"""Parent-balanced, held-out localisation and integrated regional outputs."""

import csv
import json
import warnings
from collections import defaultdict
from dataclasses import replace

import numpy as np
import pytest

from intergenic_regions.cli import build_parser, main, region_options
from intergenic_regions.io import read_fasta
from intergenic_regions.learning import fit_classifier, predict_sequences
from intergenic_regions.regional_analysis import regional_outputs
from intergenic_regions.window_learning import (
    fit_window_model,
    parent_window_weights,
)
from intergenic_regions.window_reporting import plot_regional_results
from intergenic_regions.windows import assign_window_groups, make_windows
from intergenic_regions.workflows import (
    enrichment_outputs,
    enrichment_workflow,
    pipeline_workflow,
)


def read_tsv(path):
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_constant_window_covariates_have_exact_safe_scaling(windows):
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        *_, model, _ = fit_window_model(
            windows=windows,
            parent_groups=assign_window_groups(windows=windows),
            lengths=[6],
            folds=3,
            shap=False,
        )
    assert model["composition_mean"][1] == pytest.approx(np.log(61))
    assert model["composition_scale"][1:] == [1.0, 1.0]


def test_local_reference_overlap_and_single_midpoint_plot(dataset, tmp_path):
    from intergenic_regions.references import import_reference

    bed = tmp_path / "reference.bed"
    bed.write_text("p0\t40\t70\n")
    reference = tmp_path / "reference"
    import_reference(
        output=reference,
        source="custom",
        local_bed=bed,
        organism="Test_species",
        assembly="test_build",
    )
    output = tmp_path / "reference_pipeline"
    pipeline_workflow(
        genome_path=dataset["genome"],
        annotation_path=dataset["annotation"],
        positive_genes=dataset["positive_genes"],
        negative_genes=dataset["negative_genes"],
        output=output,
        use_ai=False,
        regional_settings={
            "lengths": [120],
            "step": 120,
            "max_locus_plots": 1,
        },
        evidence_settings={
            "references": [reference],
            "organism": "Test_species",
            "assembly": "test_build",
        },
    )
    rows = read_tsv(output / "regions/window_scores.tsv")
    assert (
        next(r for r in rows if r["parent_id"] == "p0|upstream")[
            "reference_support_count"
        ]
        == "1"
    )
    assert (output / "regions/figures/region_locus_01.png").is_file()


@pytest.fixture
def windows(dataset):
    items, _ = make_windows(
        positive=read_fasta(path=dataset["positive"]),
        negative=read_fasta(path=dataset["negative"]),
        lengths=[60],
        step=30,
    )
    return items


def test_parent_weights_equalise_long_flanks_and_classes():
    indices = np.asarray([0, 0, 1, 2, 2, 2], dtype=np.int64)
    labels = np.asarray([1, 1, 0], dtype=np.int64)
    weights = parent_window_weights(
        parent_indices=indices, parent_labels=labels
    )
    assert weights.sum() == pytest.approx(3)
    assert [weights[indices == i].sum() for i in range(3)] == pytest.approx(
        [0.75, 0.75, 1.5]
    )
    assert (
        weights[labels[indices] == 1].sum()
        == weights[labels[indices] == 0].sum()
    )
    model = fit_classifier(
        matrix=np.asarray([[0], [0], [1], [2], [2], [2]]),
        labels=labels[indices],
        sample_weight=weights,
        seed=42,
    )
    assert model.class_weight is None
    assert np.isfinite(model.coef_).all()
    for weights in (
        np.asarray([0.0] * 6),
        np.ones(5),
        np.asarray([float("nan")] * 6),
    ):
        with pytest.raises(ValueError, match="weights"):
            fit_classifier(
                matrix=np.ones((6, 2)),
                labels=labels[indices],
                sample_weight=weights,
                seed=42,
            )


@pytest.mark.parametrize(
    "indices,labels",
    [
        ([], [1, 0]),
        ([[0]], [1, 0]),
        ([0.5], [1, 0]),
        ([-1], [1, 0]),
        ([2], [1, 0]),
        ([0], [[1, 0]]),
        ([0, 1], [2, 0]),
        ([0], [1, 0]),
    ],
)
def test_invalid_parent_weights(indices, labels):
    with pytest.raises(ValueError):
        parent_window_weights(
            parent_indices=np.asarray(indices),
            parent_labels=np.asarray(labels),
        )


def test_window_model_holds_out_sources_and_families_and_vocabulary(
    windows, monkeypatch
):
    parents = {w.parent_id: w.parent_id[1:] for w in windows}
    groups = assign_window_groups(windows=windows, parent_groups=parents)
    training_units = []
    from intergenic_regions import window_learning

    original = window_learning.fit_vocabulary

    def monitored(*, counters, max_features):
        training_units.append(len(counters))
        return original(counters=counters, max_features=max_features)

    monkeypatch.setattr(window_learning, "fit_vocabulary", monitored)
    predictions, parent_rows, metrics, summary, model, explanations = (
        fit_window_model(
            windows=windows,
            parent_groups=groups,
            lengths=[6],
            folds=3,
            shap_max_sequences=20,
            shap_max_features=5,
        )
    )
    assert training_units == [20, 20, 20, 30]
    assert len(predictions) == 90 and len(parent_rows) == 30
    assert len(metrics) == 3
    folds_by_group = defaultdict(set)
    for prediction in predictions:
        folds_by_group[prediction["group"]].add(prediction["fold"])
    assert all(len(folds) == 1 for folds in folds_by_group.values())
    assert summary["metric_unit"] == "one_mean_score_per_parent"
    assert summary["training_weight_unit"] == "class_balanced_parents"
    assert summary["permutation_p_value"] is None
    assert "not experimentally" in summary["interpretation"]
    by_parent = defaultdict(list)
    for item, prediction in zip(windows, predictions, strict=True):
        by_parent[item.parent_id].append(
            prediction["held_out_signature_score"]
        )
    for record in parent_rows:
        assert record["mean_held_out_signature_score"] == pytest.approx(
            np.mean(by_parent[record["parent_id"]])
        )
    assert model["window_length"] == 60
    scores = predict_sequences(model=model, sequences={"new": "ACGT" * 15})
    assert 0 <= scores[0]["candidate_signature_score"] <= 1
    for width in (True, "60", 19, 100001):
        with pytest.raises(ValueError, match="training length"):
            predict_sequences(
                model={**model, "window_length": width},
                sequences={"new": "ACGT" * 15},
            )
    with pytest.raises(ValueError, match="training length"):
        predict_sequences(model=model, sequences={"new": "ACGT" * 30})
    assert summary["shap"]["background"].startswith("parent-balanced")
    assert summary["shap"]["maximum_reconstruction_error"] < 1e-10
    reconstruction = defaultdict(lambda: {"sum": 0, "base": 0, "score": 0})
    for record in explanations["rows"]:
        local = reconstruction[record["sequence_id"]]
        local["sum"] += record["shap_value"]
        local.update(
            base=record["base_value"], score=record["signature_score"]
        )
    by_id = {r["sequence_id"]: r for r in predictions}
    for identifier, record in reconstruction.items():
        assert record["sum"] + record["base"] == pytest.approx(
            np.log(record["score"] / (1 - record["score"]))
        )
        assert record["score"] == pytest.approx(
            by_id[identifier]["held_out_signature_score"]
        )


def test_window_model_without_shap_and_shared_windows_block_validation(
    windows,
):
    *_, summary, model, explanations = fit_window_model(
        windows=windows,
        parent_groups=assign_window_groups(windows=windows),
        lengths=[6],
        folds=3,
        shap=False,
    )
    assert summary["shap"]["status"] == "disabled" and explanations == {}
    assert model["summary"] is summary
    shared = [replace(w, sequence="A" * 60) for w in windows]
    with pytest.raises(ValueError, match="independent groups"):
        fit_window_model(
            windows=shared,
            parent_groups=assign_window_groups(windows=shared),
            lengths=[6],
            folds=3,
        )


@pytest.mark.parametrize(
    "changes",
    [
        {"windows": []},
        {"lengths": []},
        {"max_features": 0},
        {"shap_max_features": 0},
        {"shap_max_features": True},
        {"parent_groups": {}},
        {"folds": 1},
        {"shap_max_sequences": 1},
    ],
)
def test_bad_window_learning_inputs(windows, changes):
    options = dict(
        windows=windows,
        parent_groups=assign_window_groups(windows=windows),
        lengths=[6],
        folds=3,
    )
    options.update(changes)
    with pytest.raises(ValueError):
        fit_window_model(**options)


def test_window_labels_duplicate_rows_lengths_and_small_sample_validation(
    windows,
):
    groups = assign_window_groups(windows=windows)
    for items in (
        [*windows, windows[0]],
        [
            replace(windows[0], sequence="A" * 20, sequence_end=20),
            *windows[1:],
        ],
        [replace(windows[0], label="negative"), *windows[1:]],
        windows[:3],
    ):
        with pytest.raises(ValueError):
            fit_window_model(
                windows=items, parent_groups=groups, lengths=[6], folds=3
            )


def test_cli_defaults_regions_and_named_configuration():
    options = build_parser().parse_args(
        [
            "enrich",
            "--positive-fasta",
            "p.fa",
            "--negative-fasta",
            "n.fa",
            "--output-dir",
            "out",
        ]
    )
    assert options.no_regions is False
    assert region_options(args=options)["lengths"] == [
        100,
        200,
        400,
        800,
        1600,
    ]
    assert region_options(args=options)["step"] == 50


def test_integrated_coordinates_evidence_shap_exports_and_portable_graphs(
    dataset, tmp_path
):
    output = tmp_path / "pipeline"
    summary = pipeline_workflow(
        genome_path=dataset["genome"],
        annotation_path=dataset["annotation"],
        positive_genes=dataset["positive_genes"],
        negative_genes=dataset["negative_genes"],
        output=output,
        motif_path=dataset["motifs"],
        learning_settings={
            "lengths": [6],
            "folds": 3,
            "permutations": 0,
            "shap_max_sequences": 10,
            "shap_max_features": 5,
        },
        regional_settings={
            "lengths": [60, 120, 200],
            "step": 30,
            "score_threshold": 0.0,
            "max_locus_plots": 1,
        },
        evidence_settings={
            "accessibility_bed": dataset["peaks"],
            "evidence_tsv": dataset["evidence"],
        },
    )
    regional = summary["multiscale"]
    assert (
        regional["status"] == "completed" and regional["total_windows"] == 120
    )
    assert [s["status"] for s in regional["scales"]] == [
        "completed",
        "completed",
        "not_estimable",
    ]
    rows = read_tsv(output / "regions/window_scores.tsv")
    p0 = [
        r
        for r in rows
        if r["parent_id"] == "p0|upstream" and r["window_length"] == "60"
    ]
    assert [int(r["accessibility_overlap_bp"]) for r in p0] == [20, 30, 10]
    assert all(
        r["gene_id"] == "p0"
        and "reporter_assay" in r["linked_gene_evidence_json"]
        for r in p0
    )
    assert all(float(r["distance_to_gene_start"]) < 0 for r in rows)
    assert all(
        0 <= int(r["genomic_start"]) < int(r["genomic_end"]) <= 120
        for r in rows
    )
    candidates = read_tsv(output / "regions/candidate_regions.tsv")
    assert len(candidates) == 30
    assert all(
        r["enhancer_status"] == "unvalidated_candidate"
        and r["region_q_value"] == ""
        for r in candidates
    )
    assert "functional boundaries" in candidates[0]["uncertainty"]
    assert (
        len(read_fasta(path=output / "regions/candidate_regions.fasta")) == 30
    )
    assert (
        len(
            (output / "regions/candidate_regions.bed").read_text().splitlines()
        )
        == 30
    )
    assert (output / "regions/scale_60/shap_values.tsv").is_file()
    assert (
        len(list((output / "regions/scale_60/figures").glob("shap_*.png")))
        == 6
    )
    report = (output / "report.html").read_text()
    assert "Multi-scale regulatory region report" in report
    assert "region_score_distance" in report
    assert "region_motif_distance" in report
    assert "data:application/pdf;base64," in report
    detail = (output / "regions/report.html").read_text()
    assert "No region-level" in detail and "60 bp SHAP report" in detail
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["settings"]["regional_search"]["lengths"] == [60, 120, 200]


@pytest.fixture
def region_input(dataset, tmp_path):
    positive = read_fasta(path=dataset["positive"])
    negative = read_fasta(path=dataset["negative"])
    native = tmp_path / "motifs"
    enrichment_outputs(
        directory=native,
        positive=positive,
        negative=negative,
        motif_path=dataset["motifs"],
    )
    return dict(
        directory=tmp_path / "regions",
        positive=positive,
        negative=negative,
        motif_directory=native,
        lengths=[60],
        step=60,
        max_locus_plots=0,
        learning_settings={"lengths": [6], "folds": 3, "shap": False},
    )


def test_regional_disabled_unavailable_and_small_source_statuses(
    region_input, monkeypatch
):
    assert (
        regional_outputs(**region_input, enabled=False)["status"] == "disabled"
    )
    assert (
        regional_outputs(**region_input, use_ai=False)["scales"][0]["status"]
        == "disabled"
    )

    def missing(**kwargs):
        raise ImportError("missing sklearn")

    monkeypatch.setattr(
        "intergenic_regions.window_learning.fit_window_model", missing
    )
    assert (
        regional_outputs(**region_input)["scales"][0]["status"]
        == "unavailable"
    )
    assert all(
        r["held_out_signature_score"] == ""
        for r in read_tsv(region_input["directory"] / "window_scores.tsv")
    )


def test_regional_too_short_groups_and_invalid_inputs(region_input, tmp_path):
    options = {**region_input, "lengths": [200]}
    assert regional_outputs(**options)["status"] == "no_complete_windows"
    options = {
        **region_input,
        "positive": dict(list(region_input["positive"].items())[:4]),
    }
    assert (
        regional_outputs(**options)["scales"][0]["status"] == "not_estimable"
    )
    for change in (
        {"score_threshold": -1},
        {"motif_q_value": 0},
        {"max_motifs": 0},
        {"max_locus_plots": 21},
        {"position_bin_width": 0},
        {"learning_settings": {"folds": 1}},
    ):
        with pytest.raises(ValueError):
            regional_outputs(**{**region_input, **change})
    groups = tmp_path / "groups.tsv"
    groups.write_text(
        "sequence_id\tgroup\n"
        + "".join(
            f"{p}\t{p}\n"
            for p in [*region_input["positive"], *region_input["negative"]]
        )
    )
    assert (
        regional_outputs(**region_input, groups_path=groups)["scales"][0][
            "status"
        ]
        == "completed"
    )
    for text in (
        "bad",
        "sequence_id\tgroup\np0\tx\np0\ty\n",
        "sequence_id\tgroup\np0\t\n",
        "sequence_id\tgroup\np0\tx\n",
    ):
        groups.write_text(text)
        with pytest.raises(ValueError):
            regional_outputs(**region_input, groups_path=groups)


def test_numeric_window_shap_retained_when_official_plots_unavailable(
    region_input, monkeypatch
):
    def unavailable(**kwargs):
        raise ImportError("missing shap")

    monkeypatch.setattr(
        "intergenic_regions.shap_reporting.plot_shap", unavailable
    )
    result = regional_outputs(
        **{
            **region_input,
            "learning_settings": {"lengths": [6], "folds": 3, "shap": True},
        }
    )
    assert result["scales"][0]["shap"]["plot_status"] == "unavailable"
    assert (region_input["directory"] / "scale_60/shap_values.tsv").is_file()


def test_enrich_no_regions_no_ml_and_transactional_limits(dataset, tmp_path):
    kwargs = dict(
        positive_path=dataset["positive"],
        negative_path=dataset["negative"],
        output=tmp_path / "analysis",
        learning_settings={"lengths": [6], "folds": 3, "permutations": 0},
    )
    summary = enrichment_workflow(**kwargs, use_regions=False, use_ai=False)
    assert summary["multiscale"]["status"] == "disabled"
    assert not (kwargs["output"] / "regions/window_scores.tsv").exists()
    output = tmp_path / "bad"
    result = main(
        argv=[
            "enrich",
            "--positive-fasta",
            str(dataset["positive"]),
            "--negative-fasta",
            str(dataset["negative"]),
            "--region-max-windows",
            "1",
            "--no-ml",
            "--output-dir",
            str(output),
        ]
    )
    assert result == 2 and not output.exists()


def test_plot_empty_and_invalid_limit(tmp_path):
    paths = plot_regional_results(
        rows=[],
        audit=[],
        candidates=[],
        profiles=[],
        scales=[],
        directory=tmp_path,
        max_loci=0,
    )
    assert len(paths) == 2
    assert all(p.is_file() and p.with_suffix(".pdf").is_file() for p in paths)
    with pytest.raises(ValueError):
        plot_regional_results(
            rows=[],
            audit=[],
            candidates=[],
            profiles=[],
            scales=[],
            directory=tmp_path,
            max_loci=21,
        )
