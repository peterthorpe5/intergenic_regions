import gzip
import json
import shutil
import subprocess
from dataclasses import replace
from types import SimpleNamespace

import pytest

from intergenic_regions.annotation import read_annotation
from intergenic_regions.rna_inputs import (
    RnaSample,
    read_gene_exons,
    read_rna_samples,
)
from intergenic_regions.rna_tools import (
    count_rna_genes,
    inspect_fastq,
    map_rna_reads,
    plan_rna_mapping,
)
from intergenic_regions.rna_workflows import rna_pipeline_workflow


def test_fastq_plain_gzip_and_bounded_check(tmp_path):
    text = "@r\nACGTRN\n+\nIIIIII\n"
    path = tmp_path / "reads.fastq"
    path.write_text(text + "broken remainder\n")
    assert (
        inspect_fastq(path=path, max_records=1)["maximum_sampled_read_length"]
        == 6
    )
    with pytest.raises(ValueError):
        inspect_fastq(path=path)
    gz = tmp_path / "reads.fastq.gz"
    with gzip.open(filename=gz, mode="wt") as stream:
        stream.write(text)
    assert inspect_fastq(path=gz)["sampled_records"] == 1


@pytest.mark.parametrize(
    "text",
    [
        "",
        "r\nACG\n+\nIII\n",
        "@r\nACG\n-\nIII\n",
        "@r\nACG\n+\nII\n",
        "@r\nACG\n+\n",
        "@r\nXYZ\n+\nIII\n",
        "@r\nA\n+\n \n",
    ],
)
def test_invalid_fastq(tmp_path, text):
    path = tmp_path / "reads.fq"
    path.write_text(text)
    with pytest.raises(ValueError):
        inspect_fastq(path=path)
    with pytest.raises(ValueError):
        inspect_fastq(path=path, max_records=0)


def test_mapping_dry_run_safe_paths_and_bounded_star_index(
    rna_dataset, tmp_path
):
    samples = read_rna_samples(path=rna_dataset / "raw_samples.tsv")[:1]
    suspicious = tmp_path / "reads with spaces;$(touch unsafe).fastq.gz"
    with gzip.open(filename=suspicious, mode="wb") as stream:
        stream.write(samples[0].read1.read_bytes())
    samples = [replace(samples[0], read1=suspicious)]
    plan = plan_rna_mapping(
        samples=samples,
        genome_path=rna_dataset / "genome.fasta",
        directory=tmp_path / "plan",
        threads=2,
    )
    assert plan["read_links"][0]["source"] == str(suspicious.resolve())
    argv = plan["commands"][1]["argv"]
    assert "gzip" in argv and "-cd" in argv
    assert all("$(touch" not in arg for arg in argv)
    assert argv[argv.index("--readFilesIn") + 1].startswith("reads/")
    index_argv = plan["commands"][0]["argv"]
    assert index_argv[index_argv.index("--sjdbOverhang") + 1] == "0"
    assert int(index_argv[index_argv.index("--genomeSAindexNbases") + 1]) < 14
    mapped, summary = map_rna_reads(
        samples=samples,
        genome_path=rna_dataset / "genome.fasta",
        directory=tmp_path / "dry",
        dry_run=True,
    )
    assert mapped == samples and summary["status"] == "dry_run"
    assert not (tmp_path / "unsafe").exists()
    assert not (tmp_path / "dry/reads").exists()


@pytest.mark.parametrize(
    "settings",
    [
        {"samples": []},
        {"threads": 0},
        {"sort_ram": 1},
        {"index_ram": 1},
        {"read_length": 59},
        {"read_length": True},
    ],
)
def test_mapping_limits(rna_dataset, tmp_path, settings):
    samples = read_rna_samples(path=rna_dataset / "raw_samples.tsv")[:1]
    with pytest.raises(ValueError):
        plan_rna_mapping(
            **{
                "samples": samples,
                "genome_path": rna_dataset / "genome.fasta",
                "directory": tmp_path,
                **settings,
            }
        )


def test_existing_star_index_identity_and_paired_compression(
    rna_dataset, tmp_path
):
    samples = read_rna_samples(path=rna_dataset / "raw_samples.tsv")[:1]
    index = tmp_path / "index"
    index.mkdir()
    with pytest.raises(ValueError):
        plan_rna_mapping(
            samples=samples,
            genome_path=rna_dataset / "genome.fasta",
            directory=tmp_path,
            star_index=index,
        )
    for name in ("Genome", "SA", "SAindex", "genomeParameters.txt"):
        (index / name).touch()
    plan = plan_rna_mapping(
        samples=samples,
        genome_path=rna_dataset / "genome.fasta",
        directory=tmp_path,
        star_index=index,
    )
    assert len(plan["commands"]) == 2
    assert (
        plan["reference_verification"] == "supplied_index_identity_unverified"
    )
    (index / "intergenic_reference.json").write_text(
        '{"genome_sha256":"wrong"}'
    )
    with pytest.raises(ValueError, match="fingerprint"):
        plan_rna_mapping(
            samples=samples,
            genome_path=rna_dataset / "genome.fasta",
            directory=tmp_path,
            star_index=index,
        )
    paired = replace(samples[0], layout="paired", read2=tmp_path / "r2.gz")
    with gzip.open(filename=paired.read2, mode="wt") as stream:
        stream.write("@r\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="compression"):
        plan_rna_mapping(
            samples=[paired],
            genome_path=rna_dataset / "genome.fasta",
            directory=tmp_path,
        )


def test_tool_failure_timeout_and_missing_output_are_explicit(
    rna_dataset, tmp_path, monkeypatch
):
    import intergenic_regions.rna_tools as module

    samples = read_rna_samples(path=rna_dataset / "raw_samples.tsv")[:1]
    args = dict(samples=samples, genome_path=rna_dataset / "genome.fasta")
    monkeypatch.setattr(module.shutil, "which", lambda exe: None)
    with pytest.raises(FileNotFoundError):
        map_rna_reads(**args, directory=tmp_path / "absent")
    monkeypatch.setattr(module.shutil, "which", lambda exe: "/fake/tool")
    monkeypatch.setattr(
        module.subprocess, "run", lambda **kw: SimpleNamespace(returncode=1)
    )
    with pytest.raises(RuntimeError, match="failed"):
        map_rna_reads(**args, directory=tmp_path / "failed")

    def expire(**kwargs):
        raise subprocess.TimeoutExpired(cmd=kwargs["args"], timeout=1)

    monkeypatch.setattr(module.subprocess, "run", expire)
    with pytest.raises(TimeoutError):
        map_rna_reads(**args, directory=tmp_path / "timeout")
    with pytest.raises(ValueError, match="timeout"):
        map_rna_reads(**args, directory=tmp_path / "invalid", timeout=0)
    monkeypatch.setattr(
        module.subprocess, "run", lambda **kw: SimpleNamespace(returncode=0)
    )
    with pytest.raises(RuntimeError, match="complete indexed BAM"):
        map_rna_reads(**args, directory=tmp_path / "incomplete")


def test_counting_input_validation(rna_dataset, tmp_path):
    genes = read_annotation(path=rna_dataset / "genes.gff3")
    samples = [
        RnaSample(
            sample_id="a",
            condition="c",
            replicate_id="1",
            strandedness="forward",
            layout="single",
            bam=tmp_path / "a.bam",
        )
    ]
    exons = read_gene_exons(path=rna_dataset / "genes.gff3", genes=genes)
    args = dict(
        samples=samples, genes=genes, exons=exons, directory=tmp_path / "c"
    )
    for changes in (
        {"samples": []},
        {"exons": {}},
        {"threads": 0},
        {"exons": {"absent": [(0, 10)]}},
        {"exons": {genes[0].gene_id: [(0, 999999)]}},
        {},
    ):
        with pytest.raises(ValueError):
            count_rna_genes(**{**args, **changes})


@pytest.mark.integration
@pytest.mark.skipif(
    not all(
        shutil.which(name) for name in ("STAR", "samtools", "featureCounts")
    ),
    reason="Optional external RNA tools are not installed",
)
def test_actual_raw_reads_mapping_counting_de_with_curated_labels(
    rna_dataset, tmp_path
):
    output = tmp_path / "real_tools"
    summary = rna_pipeline_workflow(
        samples_path=rna_dataset / "raw_samples.tsv",
        genome_path=rna_dataset / "genome.fasta",
        annotation_path=rna_dataset / "genes.gff3",
        positive_genes=rna_dataset / "positive_genes.txt",
        negative_genes=rna_dataset / "negative_genes.txt",
        condition="stimulated",
        reference="reference",
        output=output,
        map_reads=True,
        run_de=True,
        max_loci=0,
        mapping_settings={
            "threads": 2,
            "sort_ram": 100000000,
            "index_ram": 200000000,
            "timeout": 120,
        },
        counting_settings={"threads": 2},
        pipeline_settings={
            "motif_path": rna_dataset / "motifs.tsv",
            "use_ai": False,
            "use_regions": False,
        },
    )
    assert summary["rna"]["mapping"]["status"] == "completed"
    assert summary["rna"]["de"]["status"] == "completed"
    assert summary["rna"]["gene_sets"] == "user_supplied"
    assert summary["rna"]["applied_rna_anchors"] == 24
    assert not (output / "mapping/reads").exists()
    mapped = read_rna_samples(path=output / "mapping/mapped_samples.tsv")
    assert all(sample.bam.is_file() for sample in mapped)
    plan = json.loads(
        (output / "expression/counting/counting_plan.json").read_text()
    )
    assert plan["tool_versions"]["featureCounts"] != "version_unavailable"
    counts = json.loads((output / "rna/summary.json").read_text())
    assert counts["coverage_divisors"]["reference_1"] != 1
