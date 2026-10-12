import argparse
import csv
import json
import runpy
from pathlib import Path

import pytest

from intergenic_regions.cli import build_parser, main
from intergenic_regions.rna_cli import (
    add_rna_start_options,
    dispatch_rna,
    rna_start_options,
)
from intergenic_regions.rna_inputs import read_rna_samples
from intergenic_regions.rna_workflows import (
    rna_de_workflow,
    rna_input_paths,
    rna_mapping_workflow,
    rna_pipeline_workflow,
    rna_start_outputs,
    rna_starts_workflow,
)


def read_rows(*, path):
    with path.open() as stream:
        return list(csv.DictReader(f=stream, delimiter="\t"))


def pipeline_inputs(*, directory, output):
    return dict(
        samples_path=directory / "samples.tsv",
        genome_path=directory / "genome.fasta",
        annotation_path=directory / "genes.gff3",
        output=output,
        condition="stimulated",
        reference="reference",
        positive_genes=directory / "positive_genes.txt",
        negative_genes=directory / "negative_genes.txt",
        max_loci=0,
        pipeline_settings={
            "motif_path": directory / "motifs.tsv",
            "use_ai": False,
            "use_regions": False,
            "extraction_settings": {"length": None},
        },
    )


def test_independent_calls_recover_both_strands_planted_truth(
    rna_dataset, tmp_path
):
    samples = read_rna_samples(path=rna_dataset / "samples.tsv")
    identifiers = (rna_dataset / "positive_genes.txt").read_text().splitlines()
    data = rna_start_outputs(
        directory=tmp_path / "starts",
        samples=samples,
        genome_path=rna_dataset / "genome.fasta",
        annotation_path=rna_dataset / "genes.gff3",
        identifiers=identifiers,
        condition="stimulated",
        reference="reference",
    )
    truth = {
        (r["gene_id"], r["sample_id"]): int(r["planted_start"])
        for r in read_rows(path=rna_dataset / "planted_starts.tsv")
    }
    assert all(
        r["candidate_start"] == truth[(r["gene_id"], r["sample_id"])]
        for r in data["calls"]
    )
    assert len(data["calls"]) == 72
    assert data["summary"]["supported_gene_conditions"] == 24
    assert data["summary"]["edge_p_values"] is None
    assert all(r["promoter_usage_p_value"] is None for r in data["shifts"])
    assert rna_input_paths(samples=samples)


@pytest.mark.parametrize(
    "settings",
    [
        {"max_genes": 1},
        {"max_profile_rows": 1},
        {"max_profile_rows": 0},
        {"profile_bin_width": 0},
        {"condition": "absent"},
        {"reference": "reference"},
        {"identifiers": ["absent"]},
        {"identifiers": []},
        {"identifiers": ["positive_001", "positive_001"]},
        {"samples": []},
        {"window": 1},
    ],
)
def test_start_caps_and_inputs_fail_before_publication(
    rna_dataset, tmp_path, settings
):
    samples = read_rna_samples(path=rna_dataset / "samples.tsv")
    with pytest.raises(ValueError):
        rna_start_outputs(
            **{
                "directory": tmp_path / "invalid",
                "samples": samples,
                "genome_path": rna_dataset / "genome.fasta",
                "annotation_path": rna_dataset / "genes.gff3",
                **settings,
            }
        )


def test_curated_trimmed_lists_are_primary_without_de(rna_dataset, tmp_path):
    positive, negative = tmp_path / "p.txt", tmp_path / "n.txt"
    positive.write_text("positive_001\npositive_002\npositive_003\n")
    negative.write_text("negative_001\nnegative_002\nnegative_003\n")
    output = tmp_path / "curated"
    args = pipeline_inputs(directory=rna_dataset, output=output)
    result = rna_pipeline_workflow(
        **{
            **args,
            "positive_genes": positive,
            "negative_genes": negative,
            "reference": None,
        }
    )
    assert result["rna"]["gene_sets"] == "user_supplied"
    assert result["rna"]["de"]["status"] == "not_supplied"
    assert result["rna"]["searched_genes"] == 6
    assert (output / "expression/positive_genes.txt").read_text() == (
        positive.read_text()
    )
    assert (output / "expression/negative_genes.txt").read_text() == (
        negative.read_text()
    )
    assert not (output / "expression/de_results.tsv").exists()
    assert (output / "rna/report.html").is_file()
    report = (output / "report.html").read_text()
    assert "How to interpret this figure" in report
    assert "RNA evidence report" in report
    assert "unvalidated" in report.lower()
    assert result["rna"]["applied_rna_anchors"] == 6
    priorities = read_rows(path=output / "candidate_priorities.tsv")
    assert len(priorities) == 6
    assert all(
        r["enhancer_status"] != "validated_enhancer" for r in priorities
    )


def test_imported_de_context_cannot_override_curated_labels(
    rna_dataset, tmp_path
):
    imported = tmp_path / "opposite.tsv"
    positive = (rna_dataset / "positive_genes.txt").read_text().splitlines()
    negative = (rna_dataset / "negative_genes.txt").read_text().splitlines()
    imported.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\n"
        + "".join(f"{g}\t100\t0\t0.9\n" for g in positive)
        + "".join(f"{g}\t100\t4\t0.0001\n" for g in negative)
    )
    output = tmp_path / "de_context"
    result = rna_pipeline_workflow(
        **pipeline_inputs(directory=rna_dataset, output=output),
        de_path=imported,
    )
    assert result["rna"]["de"]["status"] == "imported"
    assert result["rna"]["gene_sets"] == "user_supplied"
    rows = read_rows(path=output / "candidate_priorities.tsv")
    assert {r["gene_id"] for r in rows if r["label"] == "positive"} == set(
        positive
    )
    assert all(r["de_padj"] == "0.9" for r in rows if r["label"] == "positive")
    assert {
        r["reason"]
        for r in read_rows(path=output / "expression/gene_set_audit.tsv")
    } == {"user_supplied_gene_set"}


def test_de_sets_are_generated_only_when_lists_are_absent(
    rna_dataset, tmp_path
):
    output = tmp_path / "de_selection"
    args = pipeline_inputs(directory=rna_dataset, output=output)
    result = rna_pipeline_workflow(
        **{
            **args,
            "positive_genes": None,
            "negative_genes": None,
            "counts_path": rna_dataset / "gene_counts.tsv",
        }
    )
    assert result["rna"]["de"]["status"] == "completed"
    assert result["rna"]["gene_sets"] == "DE_threshold_screening_sets"
    assert result["rna"]["positive_gene_set"] == 12
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["settings"]["command"] == "rna-pipeline"
    assert all(".tmp" not in r["path"] for r in manifest["inputs"])


@pytest.mark.parametrize(
    "changes",
    [
        {"negative_genes": None},
        {"condition": "absent"},
        {"reference": "stimulated"},
        {"reference": None, "run_de": True},
        {"positive_genes": None, "negative_genes": None, "reference": None},
    ],
)
def test_pipeline_invalid_inputs_leave_no_partial_results(
    rna_dataset, tmp_path, changes
):
    output = tmp_path / "invalid"
    args = pipeline_inputs(directory=rna_dataset, output=output)
    with pytest.raises(ValueError):
        rna_pipeline_workflow(**{**args, **changes})
    assert not output.exists()


def test_standalone_starts_de_and_mapping_workflows(rna_dataset, tmp_path):
    starts = rna_starts_workflow(
        samples_path=rna_dataset / "samples.tsv",
        genome_path=rna_dataset / "genome.fasta",
        annotation_path=rna_dataset / "genes.gff3",
        identifiers_path=rna_dataset / "positive_genes.txt",
        output=tmp_path / "starts",
        max_loci=1,
        settings={"condition": "stimulated", "reference": "reference"},
    )
    assert starts["searched_genes"] == 12
    assert (tmp_path / "starts/figures/rna_locus_01.pdf").is_file()
    de = rna_de_workflow(
        samples_path=rna_dataset / "samples.tsv",
        counts_path=rna_dataset / "gene_counts.tsv",
        output=tmp_path / "de",
        condition="stimulated",
        reference="reference",
    )
    assert de["positive_genes"] == 12
    assert (tmp_path / "de/figures/rna_de_volcano.pdf").is_file()
    plan = rna_mapping_workflow(
        samples_path=rna_dataset / "raw_samples.tsv",
        genome_path=rna_dataset / "genome.fasta",
        output=tmp_path / "mapping",
        dry_run=True,
    )
    assert plan["status"] == "dry_run"
    assert plan["working_directory"] == str((tmp_path / "mapping").resolve())
    assert not (tmp_path / "mapping/star_index").exists()
    with pytest.raises(FileExistsError):
        rna_mapping_workflow(
            samples_path=rna_dataset / "raw_samples.tsv",
            genome_path=rna_dataset / "genome.fasta",
            output=tmp_path / "mapping",
            dry_run=True,
        )


def test_rna_cli_public_helpers_and_all_commands(rna_dataset, tmp_path):
    parser = argparse.ArgumentParser()
    add_rna_start_options(parser=parser)
    settings = rna_start_options(args=parser.parse_args([]))
    assert settings["exclude_duplicates"] is True
    parser = build_parser()
    args = parser.parse_args(
        [
            "rna-pipeline",
            "--samples",
            str(rna_dataset / "samples.tsv"),
            "--genome",
            str(rna_dataset / "genome.fasta"),
            "--annotation",
            str(rna_dataset / "genes.gff3"),
            "--condition",
            "stimulated",
            "--run-de",
            "--output-dir",
            str(tmp_path / "unused"),
        ]
    )
    assert args.ai is True and not args.no_shap and not args.no_regions
    assert args.run_de is True
    with pytest.raises(ValueError):
        dispatch_rna(args=argparse.Namespace(command="rna-unknown"))
    common = [
        "--samples",
        str(rna_dataset / "samples.tsv"),
        "--genome",
        str(rna_dataset / "genome.fasta"),
        "--annotation",
        str(rna_dataset / "genes.gff3"),
    ]
    assert (
        main(
            argv=[
                "rna-starts",
                *common,
                "--genes",
                str(rna_dataset / "positive_genes.txt"),
                "--rna-max-locus-plots",
                "0",
                "--output-dir",
                str(tmp_path / "cli_starts"),
            ]
        )
        == 0
    )
    assert (
        main(
            argv=[
                "rna-map",
                "--samples",
                str(rna_dataset / "raw_samples.tsv"),
                "--genome",
                str(rna_dataset / "genome.fasta"),
                "--dry-run",
                "--output-dir",
                str(tmp_path / "cli_mapping"),
            ]
        )
        == 0
    )
    assert (
        main(
            argv=[
                "rna-de",
                "--samples",
                str(rna_dataset / "samples.tsv"),
                "--counts-tsv",
                str(rna_dataset / "gene_counts.tsv"),
                "--condition",
                "stimulated",
                "--reference-condition",
                "reference",
                "--output-dir",
                str(tmp_path / "cli_de"),
            ]
        )
        == 0
    )
    assert (
        main(
            argv=[
                "rna-pipeline",
                *common,
                "--positive-genes",
                str(rna_dataset / "positive_genes.txt"),
                "--negative-genes",
                str(rna_dataset / "negative_genes.txt"),
                "--condition",
                "stimulated",
                "--motifs",
                str(rna_dataset / "motifs.tsv"),
                "--no-ml",
                "--no-regions",
                "--rna-max-locus-plots",
                "0",
                "--output-dir",
                str(tmp_path / "cli_pipeline"),
            ]
        )
        == 0
    )


def test_rna_demo_is_reproducible_bounded_and_explicitly_synthetic(tmp_path):
    module = runpy.run_path(
        path_name=str(
            Path(__file__).resolve().parents[1] / "examples/make_rna_demo.py"
        )
    )
    result = module["make_rna_demo"](directory=tmp_path / "a", count=2)
    module["make_rna_demo"](directory=tmp_path / "b", count=2)
    assert result["synthetic"] is True
    assert (tmp_path / "a/genome.fasta").read_bytes() == (
        tmp_path / "b/genome.fasta"
    ).read_bytes()
    with pytest.raises(ValueError):
        module["make_rna_demo"](directory=tmp_path / "bad", count=1)
    assert (
        module["main"](
            argv=["--output-dir", str(tmp_path / "cli"), "--count", "2"]
        )
        == 0
    )
    assert module["main"](argv=["--output-dir", str(tmp_path / "cli")]) == 2
