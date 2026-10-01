"""Complete workflows, CLI failures and portable scientific reports."""

import json
import subprocess
import sys

import matplotlib.pyplot as plt
import pytest

from intergenic_regions.cli import (
    build_parser,
    dispatch,
    evidence_options,
    extraction_options,
    learning_options,
    main,
    motif_options,
)
from intergenic_regions.io import read_fasta
from intergenic_regions.models import Region
from intergenic_regions.reporting import (
    plot_composition,
    plot_logos,
    plot_motif_results,
    save_figure,
    write_report,
)
from intergenic_regions.workflows import (
    attach_evidence,
    enrichment_workflow,
    evidence_input_paths,
    extract_workflow,
    learning_workflow,
    pipeline_workflow,
    read_region_table,
    region_rows,
    serialise_settings,
    write_regions,
)


def test_region_files_roundtrip_and_coordinates(dataset, tmp_path):
    output = tmp_path / "extract"
    summary = extract_workflow(
        genome_path=dataset["genome"],
        annotation_path=dataset["annotation"],
        output=output,
        identifiers_path=dataset["positive_genes"],
        extraction_settings={"length": 100},
    )
    assert summary["retained_regions"] == 15
    assert len(read_fasta(path=output / "regions.fasta")) == 15
    first = read_region_table(path=output / "regions.tsv")[0]
    assert (first.start, first.end) == (20, 120)
    assert (output / "regions.bed").read_text().splitlines()[0].split("\t")[
        1:3
    ] == ["20", "120"]
    assert (output / "regions.gff3").read_text().splitlines()[1].split("\t")[
        3:5
    ] == ["21", "120"]
    assert "not_supplied" in (output / "evidence.tsv").read_text()
    assert json.loads((output / "manifest.json").read_text())["inputs"][0][
        "sha256"
    ]
    assert "Intergenic extraction" in (output / "report.html").read_text()


def test_empty_extraction_audit_and_region_helpers(tmp_path):
    region = Region(
        gene_id="g",
        contig="x",
        start=0,
        end=0,
        strand="+",
        direction="upstream",
        sequence="",
        status="no_intergenic_space",
        stop_reason="contig_boundary",
        available_length=0,
    )
    assert region_rows(regions=[region])[0]["length"] == 0
    directory = tmp_path / "empty"
    assert (
        write_regions(directory=directory, regions=[region])[
            "retained_regions"
        ]
        == 0
    )
    assert (directory / "regions.bed").read_text() == ""
    assert (
        read_region_table(path=directory / "regions.tsv")[0].status
        == "no_intergenic_space"
    )
    assert serialise_settings(
        settings={"path": directory, "paths": [directory], "enabled": True}
    )["path"] == str(directory.resolve())
    assert evidence_input_paths(
        settings={
            "references": [directory],
            "enhancer_bed": tmp_path / "e.bed",
        }
    ) == [tmp_path / "e.bed", directory / "reference.json"]


@pytest.mark.parametrize(
    "text",
    [
        "bad",
        "gene_id\tcontig\tstart\tend\tstrand\tdirection\tstatus\ng\tx\t-1\t2\t+\tupstream\tretained\n",
        "gene_id\tcontig\tstart\tend\tstrand\tdirection\tstatus\ng\tx\t0\t0\t+\tupstream\tretained\n",
        "gene_id\tcontig\tstart\tend\tstrand\tdirection\tstatus\ng\tx\t0\t2\t+\tupstream\tretained\ng\tx\t0\t2\t+\tupstream\tretained\n",
    ],
)
def test_region_table_validation(tmp_path, text):
    path = tmp_path / "bad.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_region_table(path=path)


def test_evidence_mismatch_and_required_reference_metadata(tmp_path):
    region = Region(
        gene_id="g",
        contig="x",
        start=0,
        end=5,
        strand="+",
        direction="upstream",
        sequence="AAAAA",
        status="retained",
        stop_reason="contig_boundary",
        available_length=5,
    )
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("wrong\t0\t5\n")
    with pytest.raises(ValueError, match="contig"):
        attach_evidence(
            directory=tmp_path, regions=[region], accessibility_bed=peaks
        )
    with pytest.raises(ValueError, match="organism"):
        attach_evidence(
            directory=tmp_path, regions=[region], references=[tmp_path]
        )


def test_native_workflow_reports_and_graphics(dataset, tmp_path):
    output = tmp_path / "motifs"
    summary = enrichment_workflow(
        positive_path=dataset["positive"],
        negative_path=dataset["negative"],
        output=output,
        motif_path=dataset["motifs"],
        settings={"lengths": [6]},
    )
    assert summary["significant_q_0_05"] > 0
    assert "Gbox" in (output / "motif_enrichment.tsv").read_text()
    assert "data:image/png;base64," in (output / "report.html").read_text()
    for name in (
        "background_composition",
        "motif_prevalence",
        "motif_enrichment",
        "motif_logos",
    ):
        assert (
            (output / "figures" / f"{name}.png")
            .read_bytes()
            .startswith(b"\x89PNG")
        )
        assert (
            (output / "figures" / f"{name}.pdf")
            .read_bytes()
            .startswith(b"%PDF")
        )
    assert not plt.get_fignums()


def test_learning_workflow_candidates_and_reports(dataset, tmp_path):
    output = tmp_path / "ai"
    summary = learning_workflow(
        positive_path=dataset["positive"],
        negative_path=dataset["negative"],
        output=output,
        settings={"folds": 3, "lengths": [4, 6], "permutations": 2},
        candidates_path=dataset["positive"],
    )
    assert summary["roc_auc"] > 0.8
    assert (output / "model.json").is_file()
    assert (
        "matches_training_sequence"
        in (output / "candidate_scores.tsv").read_text()
    )
    assert (output / "figures" / "ai_permutation.pdf").is_file()


def test_integrated_pipeline_optional_evidence_and_ai(dataset, tmp_path):
    output = tmp_path / "integrated"
    summary = pipeline_workflow(
        genome_path=dataset["genome"],
        annotation_path=dataset["annotation"],
        positive_genes=dataset["positive_genes"],
        negative_genes=dataset["negative_genes"],
        output=output,
        motif_path=dataset["motifs"],
        use_ai=True,
        learning_settings={"folds": 3, "lengths": [4, 6], "permutations": 0},
        matching_settings={"max_gc_difference": 0.25},
        evidence_settings={
            "accessibility_bed": dataset["peaks"],
            "enhancer_bed": dataset["peaks"],
            "evidence_tsv": dataset["evidence"],
        },
    )
    assert summary["positive_regions"] == summary["negative_regions"] == 15
    assert summary["ai"]["grouped_validation"]
    assert summary["optional_evidence_supplied"]
    assert "reporter_assay" in (output / "candidate_summary.tsv").read_text()
    assert (output / "background_matching.tsv").read_text().count("\n") == 16
    assert (output / "motifs" / "report.html").is_file()
    assert (output / "ai" / "report.html").is_file()
    with pytest.raises(ValueError, match="overlap"):
        pipeline_workflow(
            genome_path=dataset["genome"],
            annotation_path=dataset["annotation"],
            positive_genes=dataset["positive_genes"],
            negative_genes=dataset["positive_genes"],
            output=tmp_path / "bad",
        )
    with pytest.raises(ValueError):
        pipeline_workflow(
            genome_path=dataset["genome"],
            annotation_path=dataset["annotation"],
            positive_genes=dataset["positive_genes"],
            negative_genes=dataset["negative_genes"],
            output=tmp_path / "bad",
            group_by="bad",
        )


def test_portable_report_escaping_and_empty_plots(tmp_path):
    path = tmp_path / "report.html"
    write_report(
        path=path,
        title="<script>alert(1)</script>",
        summary={"x": "<&>"},
        tables={
            "Empty": [],
            "Values": [{"x": None, "y": 0.01, "z": "<script>"}],
        },
        notes=["<note>"],
    )
    report = path.read_text()
    assert "<script>" not in report and "&lt;script&gt;" in report
    assert "No records" in report
    assert plot_motif_results(rows=[], directory=tmp_path) == []
    assert plot_logos(rows=[], motifs=[], directory=tmp_path) is None
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])
    assert save_figure(
        figure=figure, directory=tmp_path, name="unit"
    ).is_file()
    assert plot_composition(
        positive={"p": "ACGT"}, negative={"n": "GGGT"}, directory=tmp_path
    ).is_file()


def test_cli_parser_settings_and_validation(tmp_path):
    parser = build_parser()
    args = parser.parse_args(
        args=[
            "pipeline",
            "--genome",
            "g.fa",
            "--annotation",
            "a.gff",
            "--positive-genes",
            "p.txt",
            "--negative-genes",
            "n.txt",
            "--output-dir",
            "out",
            "--full-gap",
            "--ai",
            "--forward-only",
            "--ai-kmer-lengths",
            "4",
            "6",
        ]
    )
    assert extraction_options(args=args)["length"] is None
    assert evidence_options(args=args)["references"] == []
    assert motif_options(args=args)["both_strands"] is False
    assert learning_options(args=args)["lengths"] == [4, 6]
    with pytest.raises(SystemExit):
        parser.parse_args(args=["extract", "--output-dir", "out"])
    with pytest.raises(ValueError):
        dispatch(
            args=type("Args", (), {"command": "bad", "output_dir": tmp_path})()
        )


def test_cli_commands_and_summaries(dataset, tmp_path, capsys):
    assert main(argv=["references"]) == 0
    assert "screen-human-enhancers" in capsys.readouterr().out
    reference = tmp_path / "ref"
    assert (
        main(
            argv=[
                "import-reference",
                "--source",
                "custom",
                "--bed",
                str(dataset["peaks"]),
                "--organism",
                "test_species",
                "--assembly",
                "test_build",
                "--output-dir",
                str(reference),
            ]
        )
        == 0
    )
    extract = tmp_path / "cli-extract"
    assert (
        main(
            argv=[
                "extract",
                "--genome",
                str(dataset["genome"]),
                "--annotation",
                str(dataset["annotation"]),
                "--genes",
                str(dataset["positive_genes"]),
                "--reference-dir",
                str(reference),
                "--organism",
                "test_species",
                "--assembly",
                "test_build",
                "--output-dir",
                str(extract),
            ]
        )
        == 0
    )
    assert (extract / "reference_overlaps.tsv").read_text().count("\n") == 16
    assert (
        main(
            argv=[
                "annotate",
                "--regions-tsv",
                str(extract / "regions.tsv"),
                "--evidence-tsv",
                str(dataset["evidence"]),
                "--output-dir",
                str(tmp_path / "annotate"),
            ]
        )
        == 0
    )
    assert (
        main(
            argv=[
                "homer",
                "--positive-fasta",
                str(dataset["positive"]),
                "--negative-fasta",
                str(dataset["negative"]),
                "--output-dir",
                str(tmp_path / "homer"),
                "--dry-run",
            ]
        )
        == 0
    )
    ai = tmp_path / "cli-ai"
    assert (
        main(
            argv=[
                "ai",
                "--positive-fasta",
                str(dataset["positive"]),
                "--negative-fasta",
                str(dataset["negative"]),
                "--folds",
                "3",
                "--permutations",
                "0",
                "--ai-kmer-lengths",
                "4",
                "--output-dir",
                str(ai),
            ]
        )
        == 0
    )
    assert (
        main(
            argv=[
                "predict",
                "--model",
                str(ai / "model.json"),
                "--fasta",
                str(dataset["positive"]),
                "--output-dir",
                str(tmp_path / "predict"),
            ]
        )
        == 0
    )
    assert (
        main(
            argv=[
                "enrich",
                "--positive-fasta",
                str(dataset["positive"]),
                "--negative-fasta",
                str(dataset["negative"]),
                "--motifs",
                str(dataset["motifs"]),
                "--output-dir",
                str(tmp_path / "cli-enrich"),
            ]
        )
        == 0
    )


@pytest.mark.integration
def test_subprocess_end_to_end_pipeline_and_failed_input(dataset, tmp_path):
    output = tmp_path / "subprocess-output"
    command = [
        sys.executable,
        "-m",
        "intergenic_regions",
        "pipeline",
        "--genome",
        str(dataset["genome"]),
        "--annotation",
        str(dataset["annotation"]),
        "--positive-genes",
        str(dataset["positive_genes"]),
        "--negative-genes",
        str(dataset["negative_genes"]),
        "--motifs",
        str(dataset["motifs"]),
        "--kmer-lengths",
        "6",
        "--ai",
        "--folds",
        "3",
        "--permutations",
        "0",
        "--ai-kmer-lengths",
        "4",
        "6",
        "--output-dir",
        str(output),
    ]
    completed = subprocess.run(
        args=command, capture_output=True, text=True, check=False, timeout=60
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["positive_regions"] == 15
    failed = subprocess.run(
        args=command, capture_output=True, text=True, check=False, timeout=60
    )
    assert failed.returncode == 2
    assert "Output already exists" in failed.stderr
    failed = subprocess.run(
        args=[
            sys.executable,
            "-m",
            "intergenic_regions",
            "extract",
            "--genome",
            "missing",
            "--annotation",
            "missing",
            "--output-dir",
            str(tmp_path / "missing"),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=20,
    )
    assert failed.returncode == 2
    assert not (tmp_path / "missing").exists()


def test_cli_error_interrupt_dependency_and_logging(tmp_path, monkeypatch):
    def missing_dependency(*, args):
        raise ImportError("missing analysis")

    monkeypatch.setattr("intergenic_regions.cli.dispatch", missing_dependency)
    assert main(argv=["references"]) == 2

    def interrupted(*, args):
        raise KeyboardInterrupt

    monkeypatch.setattr("intergenic_regions.cli.dispatch", interrupted)
    assert main(argv=["references"]) == 130
    assert (
        main(
            argv=[
                "references",
                "--log-file",
                str(tmp_path / "missing" / "log"),
            ]
        )
        == 2
    )
    monkeypatch.setattr(
        "intergenic_regions.cli.dispatch", lambda **kwargs: {"ok": True}
    )
    assert (
        main(
            argv=[
                "references",
                "--verbose",
                "--log-file",
                str(tmp_path / "log"),
            ]
        )
        == 0
    )
