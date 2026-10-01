"""Automatic ML, explicit uncertainty and transparent evidence rankings."""

import json

import pytest

from intergenic_regions.cli import build_parser
from intergenic_regions.prioritisation import prioritise_candidates
from intergenic_regions.reporting import write_report
from intergenic_regions.workflows import (
    automatic_learning,
    enrichment_workflow,
    pipeline_workflow,
)


def test_default_cli_enables_learning():
    parser = build_parser()
    options = [
        "enrich",
        "--positive-fasta",
        "p.fa",
        "--negative-fasta",
        "n.fa",
        "--output-dir",
        "out",
    ]
    assert parser.parse_args(args=options).ai is True
    assert parser.parse_args(args=[*options, "--no-ml"]).ai is False
    with pytest.raises(SystemExit):
        parser.parse_args(args=[*options, "--ai", "--no-ml"])


def test_automatic_learning_completed(dataset, tmp_path):
    result = enrichment_workflow(
        positive_path=dataset["positive"],
        negative_path=dataset["negative"],
        motif_path=dataset["motifs"],
        output=tmp_path / "analysis",
        learning_settings={"folds": 3, "permutations": 0, "lengths": [6]},
    )
    assert result["ai"]["status"] == "completed"
    assert result["ai"]["roc_auc"] > 0.8
    report = (tmp_path / "analysis" / "report.html").read_text()
    assert "Held-out model validation" in report
    assert "data:image/png;base64," in report
    assert "candidate_priorities.tsv" in report
    assert "unvalidated_candidate" in report
    assert (tmp_path / "analysis" / "ai" / "model.json").is_file()


def test_pipeline_runs_learning_without_opt_in(dataset, tmp_path):
    result = pipeline_workflow(
        genome_path=dataset["genome"],
        annotation_path=dataset["annotation"],
        positive_genes=dataset["positive_genes"],
        negative_genes=dataset["negative_genes"],
        output=tmp_path / "pipeline",
        motif_path=dataset["motifs"],
        learning_settings={"folds": 3, "permutations": 0, "lengths": [6]},
    )
    assert result["ai"]["status"] == "completed"
    assert result["ai"]["grouped_validation"] is True
    assert (
        "unvalidated_candidate"
        in (tmp_path / "pipeline" / "candidate_priorities.tsv").read_text()
    )


def test_small_dataset_preserves_motifs_and_reports_unavailable_model(
    tmp_path,
):
    positive = tmp_path / "p.fa"
    negative = tmp_path / "n.fa"
    positive.write_text(">p\nACGTCACGTG\n")
    negative.write_text(">n\nGCTTTAGTGC\n")
    result = enrichment_workflow(
        positive_path=positive,
        negative_path=negative,
        output=tmp_path / "small",
        settings={"lengths": [4]},
    )
    assert result["ai"]["status"] == "not_estimable"
    assert "at least five" in result["ai"]["reason"]
    assert (tmp_path / "small" / "motif_enrichment.tsv").is_file()
    assert (
        "unvalidated_candidate"
        in (tmp_path / "small" / "candidate_priorities.tsv").read_text()
    )
    assert not (tmp_path / "small" / "ai" / "model.json").exists()


def test_learning_disabled_missing_dependency_and_bad_settings(
    tmp_path, monkeypatch
):
    settings = dict(
        directory=tmp_path / "disabled",
        positive={"p": "ACGT"},
        negative={"n": "CCTA"},
    )
    result, predictions = automatic_learning(**settings, enabled=False)
    assert result["status"] == "disabled" and predictions == []

    def missing(**kwargs):
        raise ImportError("No sklearn")

    monkeypatch.setattr(
        "intergenic_regions.workflows.learning_outputs", missing
    )
    result, predictions = automatic_learning(**settings)
    assert result["status"] == "unavailable" and predictions == []
    with pytest.raises(ValueError, match="two validation folds"):
        automatic_learning(**settings, settings={"folds": 1})

    def invalid(**kwargs):
        raise ValueError("Malformed group TSV")

    monkeypatch.setattr(
        "intergenic_regions.workflows.learning_outputs", invalid
    )
    with pytest.raises(ValueError, match="Malformed group"):
        automatic_learning(**settings)


def test_prioritisation_separates_sequence_scores_from_evidence():
    rows = [
        {"sequence_id": "sequence", "held_out_signature_score": 0.99},
        {
            "sequence_id": "open",
            "held_out_signature_score": 0.3,
            "accessibility_overlap_bp": 10,
        },
        {
            "sequence_id": "both",
            "held_out_signature_score": 0.2,
            "accessibility_overlap_bp": 10,
            "reference_support_count": 1,
        },
        {"sequence_id": "annotation", "enhancer_annotation_overlap_bp": 3},
        {
            "sequence_id": "linked",
            "linked_gene_evidence_json": json.dumps(
                [{"value": "positive"}, {"value": "unknown"}]
            ),
        },
        {
            "sequence_id": "contradicted",
            "linked_gene_evidence_json": json.dumps(
                [{"value": "positive"}, {"value": "negative"}]
            ),
        },
    ]
    result = prioritise_candidates(rows=rows)
    assert [r["sequence_id"] for r in result] == [
        "both",
        "open",
        "annotation",
        "linked",
        "sequence",
        "contradicted",
    ]
    assert all(r["enhancer_status"] == "unvalidated_candidate" for r in result)
    assert result[3]["uninterpreted_gene_evidence_count"] == 1
    assert result[-1]["priority_tier"] == 0
    assert "contradictory" in result[-1]["uncertainty"]
    assert "priority_rank" not in rows[0]
    assert prioritise_candidates(rows=[]) == []


@pytest.mark.parametrize(
    "rows",
    [
        [{}],
        [{"sequence_id": ""}],
        [{"sequence_id": "x"}, {"sequence_id": "x"}],
        [{"sequence_id": "x", "held_out_signature_score": float("nan")}],
        [{"sequence_id": "x", "held_out_signature_score": 1.1}],
        [{"sequence_id": "x", "held_out_signature_score": "bad"}],
        [{"sequence_id": "x", "linked_gene_evidence_json": "bad"}],
        [{"sequence_id": "x", "linked_gene_evidence_json": "{}"}],
        [{"sequence_id": "x", "linked_gene_evidence_json": "[{}]"}],
    ],
)
def test_prioritisation_rejects_malformed_records(rows):
    with pytest.raises(ValueError):
        prioritise_candidates(rows=rows)


def test_dashboard_safe_links_previews_and_interactivity(tmp_path):
    path = tmp_path / "report.html"
    write_report(
        path=path,
        title="<script>unsafe</script>",
        summary={"ai": {"status": "not_estimable", "reason": "too small"}},
        tables={"Candidates": [{"gene_id": f"g{i}"} for i in range(501)]},
        links={"Full table": "candidates.tsv"},
    )
    report = path.read_text()
    assert "&lt;script&gt;unsafe&lt;/script&gt;" in report
    assert "500 preview rows of 501" in report
    assert "g500</td>" not in report
    assert "data-sort" in report and "type='search'" in report
    assert "ML: too small" in report
    assert "fetch(" not in report
    for href in ("javascript:alert(1)", "//external.example/path"):
        with pytest.raises(ValueError, match="relative paths"):
            write_report(
                path=path, title="bad", summary={}, links={"bad": href}
            )
