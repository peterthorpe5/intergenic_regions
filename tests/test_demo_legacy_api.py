"""Public API, reproducible demo and preserved legacy boundary regressions."""

import runpy
from pathlib import Path

import pytest

from intergenic_regions import (
    GeneIndex,
    Genome,
    extract_regions,
    read_annotation,
    reverse_complement,
)
from intergenic_regions.io import read_fasta, write_fasta


def test_variable_length_synthetic_demo_keeps_both_strand_sequences(tmp_path):
    module = runpy.run_path(path_name=str(ROOT / "examples/make_demo.py"))
    path = tmp_path / "multiscale"
    result = module["make_demo"](
        output=path,
        sequences_per_class=6,
        flank_length=1000,
        module_lengths=[120, 240, 480],
        seed=31,
    )
    assert result["synthetic_module_lengths"] == [120, 240, 480]
    expected = {
        **read_fasta(path=path / "positive.fasta"),
        **read_fasta(path=path / "negative.fasta"),
    }
    genes = read_annotation(path=path / "genes.gff3")
    with Genome(path=path / "genome.fasta") as genome:
        regions = extract_regions(
            index=GeneIndex(genes=genes, lengths=genome.lengths),
            genome=genome,
            identifiers=list(expected),
            length=None,
        )
    assert {r.gene_id: r.sequence for r in regions} == expected
    assert all(len(r.sequence) == 1000 for r in regions)
    assert (path / "synthetic_modules.tsv").is_file()
    assert "GCACTG" in expected["negative_000"]
    for settings in (
        {"flank_length": 119},
        {"module_lengths": [201]},
        {"module_lengths": [True]},
    ):
        with pytest.raises(ValueError):
            module["make_demo"](output=tmp_path / "invalid", **settings)
    assert (
        module["main"](
            argv=[
                "--output-dir",
                str(tmp_path / "cli_multiscale"),
                "--flank-length",
                "500",
                "--module-lengths",
                "120",
                "240",
            ]
        )
        == 0
    )


ROOT = Path(__file__).resolve().parents[1]


def test_original_fixture_boundaries_are_strictly_intergenic(tmp_path):
    genes = read_annotation(
        path=ROOT / "tests/inputs/tests_gene_indexing_simplified.txt",
        annotation_format="tsv",
    )
    records = read_fasta(path=ROOT / "tests/inputs/genome_simplified.fasta")
    # Legacy fixture ends before one annotated gene. Preserve its bases and
    # pad only the missing tail; production validation must remain strict.
    with pytest.raises(ValueError, match="outside genome"):
        GeneIndex(genes=genes, lengths={k: len(v) for k, v in records.items()})
    padded = tmp_path / "complete_fixture.fasta"
    write_fasta(
        path=padded,
        records={
            contig: sequence
            + "N"
            * max(
                0,
                max((g.end for g in genes if g.contig == contig), default=0)
                - len(sequence),
            )
            for contig, sequence in records.items()
        },
    )
    with Genome(path=padded) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        regions = extract_regions(index=index, genome=genome, length=5)
        by_id = {r.gene_id: r for r in regions}
        minus = by_id["GPLIN_000000300"]
        # The next gene begins at one-based 55: that base cannot be returned.
        assert (minus.start, minus.end) == (50, 54)
        assert minus.sequence == reverse_complement(
            sequence=genome.fetch(contig=minus.contig, start=50, end=54)
        )
        assert by_id["GPLIN_000000100"].sequence == "ATGT"
        for region in regions:
            for gene in genes:
                if gene.contig == region.contig:
                    assert not (
                        region.start < gene.end and region.end > gene.start
                    )


def test_demo_generator_matches_extracted_sequences(tmp_path):
    module = runpy.run_path(path_name=str(ROOT / "examples/make_demo.py"))
    path = tmp_path / "demo"
    result = module["make_demo"](output=path, sequences_per_class=6, seed=13)
    assert result["synthetic"] is True
    expected = {
        **read_fasta(path=path / "positive.fasta"),
        **read_fasta(path=path / "negative.fasta"),
    }
    genes = read_annotation(path=path / "genes.gff3")
    with Genome(path=path / "genome.fasta") as genome:
        regions = extract_regions(
            index=GeneIndex(genes=genes, lengths=genome.lengths),
            genome=genome,
            identifiers=list(expected),
            length=1000,
        )
    assert {r.gene_id: r.sequence for r in regions} == expected
    other = tmp_path / "other"
    module["make_demo"](output=other, sequences_per_class=6, seed=13)
    assert (other / "genome.fasta").read_bytes() == (
        path / "genome.fasta"
    ).read_bytes()
    with pytest.raises(ValueError):
        module["make_demo"](output=tmp_path / "bad", sequences_per_class=4)
    with pytest.raises(FileExistsError):
        module["make_demo"](output=path)
    assert module["main"](argv=["--output-dir", str(tmp_path / "cli")]) == 0
    assert module["main"](argv=["--output-dir", str(path)]) == 2
    assert (
        module["main"](
            argv=[
                "--output-dir",
                str(tmp_path / "bad"),
                "--sequences-per-class",
                "1",
            ]
        )
        == 2
    )
