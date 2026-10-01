"""Exact strand/boundary tests and an independent per-base interval oracle."""

import gzip

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from intergenic_regions.extraction import (
    GeneIndex,
    extract_regions,
    merge_intervals,
)
from intergenic_regions.genome import (
    Genome,
    normalise_genome,
    reverse_complement,
    sequence_composition,
)
from intergenic_regions.io import write_fasta
from intergenic_regions.models import Gene


def test_reverse_complement_and_composition():
    assert (
        reverse_complement(sequence="AaCGTRYSWKMBDHVN") == "NBDHVKMWSRYACGtT"
    )
    assert reverse_complement(sequence="") == ""
    with pytest.raises(ValueError):
        reverse_complement(sequence="AZ")
    assert sequence_composition(sequence="GcAN") == {
        "length": 4,
        "gc_fraction": 2 / 3,
        "ambiguous_fraction": 0.25,
    }
    assert sequence_composition(sequence="")["ambiguous_fraction"] == 0
    assert sequence_composition(sequence="NN")["gc_fraction"] == 0


def test_normalise_irregular_fasta_preserves_case_and_bases(tmp_path):
    source = tmp_path / "irregular.fa"
    source.write_text(">x description\n" + "Ac" * 41 + "\n\nGT\n>y\nTT\nT\n")
    destination = tmp_path / "normalised.fa"
    normalise_genome(source=source, destination=destination)
    with Genome(path=source) as genome:
        assert genome.lengths == {"x": 84, "y": 3}
        assert genome.fetch(contig="x", start=0, end=84) == "Ac" * 41 + "GT"
    assert not (source.parent / "irregular.fa.fai").exists()
    with pytest.raises(ValueError, match="must differ"):
        normalise_genome(source=source, destination=source)


@pytest.mark.parametrize(
    "text",
    [
        "",
        "AAA",
        ">\nAAA",
        ">x\nAAA\n>x\nTTT",
        ">x\nAZZ",
        ">x\n>y\nTTT",
        ">x\nAAA\n>y\n",
    ],
)
def test_normalise_rejects_bad_genomes(tmp_path, text):
    source = tmp_path / "bad.fa"
    source.write_text(text)
    with pytest.raises(ValueError):
        normalise_genome(source=source, destination=tmp_path / "out.fa")


def test_genome_plain_gzip_bounds_and_cleanup(tmp_path):
    path = tmp_path / "genome.fa"
    write_fasta(path=path, records={"x": "ACgtRYN" * 20})
    accessor = Genome(path=path)
    with pytest.raises(RuntimeError):
        _ = accessor.lengths
    with pytest.raises(RuntimeError):
        accessor.fetch(contig="x", start=0, end=1)
    with accessor as genome:
        assert genome.lengths == {"x": 140}
        assert genome.fetch(contig="x", start=1, end=8) == "CgtRYNA"
        for contig, start, end in [
            ("missing", 0, 1),
            ("x", -1, 2),
            ("x", 2, 1),
            ("x", 0, 141),
        ]:
            with pytest.raises(ValueError):
                genome.fetch(contig=contig, start=start, end=end)
        assert genome.fetch(contig="x", start=3, end=3) == ""
    assert not (tmp_path / "genome.fa.fai").exists()
    accessor.close()
    compressed = tmp_path / "genome.fa.gz"
    with gzip.open(filename=compressed, mode="wb") as stream:
        stream.write(path.read_bytes())
    with Genome(path=compressed) as genome:
        assert genome.fetch(contig="x", start=0, end=4) == "ACgt"


@pytest.mark.parametrize(
    "text",
    [
        ">x\nAXX\n",
        ">x\nAAA\n>x\nTTT\n",
        "",
        ">x\n>y\nAAA\n",
    ],
)
def test_invalid_genomes(tmp_path, text):
    path = tmp_path / "bad.fa"
    path.write_text(text)
    with pytest.raises((ValueError, RuntimeError)):
        with Genome(path=path) as genome:
            genome.fetch(contig="x", start=0, end=3)


def test_merge_intervals():
    assert merge_intervals(intervals=[(8, 10), (1, 5), (4, 8), (2, 3)]) == [
        (1, 10)
    ]
    assert merge_intervals(intervals=[]) == []
    with pytest.raises(ValueError):
        merge_intervals(intervals=[(3, 3)])


def test_unindexed_target_cannot_be_queried():
    gene = Gene(gene_id="g", contig="x", start=5, end=10, strand="+")
    index = GeneIndex(genes=[], lengths={"x": 20})
    with pytest.raises(ValueError, match="not in this index"):
        index.free_flank(gene=gene, direction="upstream")


@given(
    spans=st.lists(
        st.tuples(
            st.integers(min_value=0, max_value=80),
            st.integers(min_value=1, max_value=20),
        ),
        min_size=1,
        max_size=20,
    ),
    strand=st.sampled_from(["+", "-"]),
    direction=st.sampled_from(["upstream", "downstream"]),
)
@settings(max_examples=200, deadline=None)
def test_interval_query_matches_independent_base_oracle(
    spans, strand, direction
):
    genes = [
        Gene(
            gene_id=f"g{i}",
            contig="x",
            start=start,
            end=start + width,
            strand=strand if i == 0 else ".",
        )
        for i, (start, width) in enumerate(spans)
    ]
    index = GeneIndex(genes=list(reversed(genes)), lengths={"x": 100})
    target = genes[0]
    occupied = {base for gene in genes for base in range(gene.start, gene.end)}
    left = (strand == "+") == (direction == "upstream")
    anchor = target.start if left else target.end
    if left:
        boundary = anchor
        while boundary > 0 and boundary - 1 not in occupied:
            boundary -= 1
        expected = (boundary, anchor)
    else:
        boundary = anchor
        while boundary < 100 and boundary not in occupied:
            boundary += 1
        expected = (anchor, boundary)
    assert index.free_flank(gene=target, direction=direction)[:2] == expected


def test_index_validation_and_unknown_strand():
    gene = Gene(gene_id="g", contig="x", start=0, end=5, strand=".")
    with pytest.raises(ValueError):
        GeneIndex(genes=[gene, gene], lengths={"x": 10})
    with pytest.raises(ValueError):
        GeneIndex(genes=[gene], lengths={"x": 3})
    with pytest.raises(ValueError):
        GeneIndex(genes=[gene], lengths={"y": 10})
    index = GeneIndex(genes=[gene], lengths={"x": 10})
    with pytest.raises(ValueError):
        index.free_flank(gene=gene, direction="upstream")
    with pytest.raises(ValueError):
        index.free_flank(gene=gene, direction="bad")


def test_exact_flanks_offsets_strands_and_exclusions(tmp_path):
    sequence = "ACGT" * 30
    path = tmp_path / "genome.fa"
    write_fasta(path=path, records={"x": sequence})
    genes = [
        Gene(gene_id="plus", contig="x", start=10, end=20, strand="+"),
        Gene(gene_id="minus", contig="x", start=30, end=40, strand="-"),
        Gene(gene_id="nested", contig="x", start=12, end=15, strand="+"),
        Gene(gene_id="unknown", contig="x", start=50, end=55, strand="."),
    ]
    with Genome(path=path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        regions = extract_regions(
            index=index,
            genome=genome,
            identifiers=["plus", "minus"],
            direction="both",
            length=7,
        )
        assert [(r.start, r.end) for r in regions] == [
            (3, 10),
            (20, 27),
            (40, 47),
            (23, 30),
        ]
        assert regions[2].sequence == reverse_complement(
            sequence=sequence[40:47]
        )
        assert regions[3].sequence == reverse_complement(
            sequence=sequence[23:30]
        )
        full = extract_regions(
            index=index, genome=genome, identifiers=["minus"], length=None
        )[0]
        assert (full.start, full.end) == (40, 50)
        offset = extract_regions(
            index=index,
            genome=genome,
            identifiers=["minus"],
            length=3,
            offset=4,
        )[0]
        assert (offset.start, offset.end) == (44, 47)
        offset = extract_regions(
            index=index, genome=genome, identifiers=["minus"], offset=11
        )[0]
        assert offset.status == "offset_exceeds_gap"
        assert (
            extract_regions(
                index=index, genome=genome, identifiers=["nested"]
            )[0].status
            == "no_intergenic_space"
        )
        assert (
            extract_regions(
                index=index, genome=genome, identifiers=["unknown"]
            )[0].status
            == "unknown_strand"
        )
        assert (
            extract_regions(
                index=index, genome=genome, identifiers=["plus"], min_length=11
            )[0].status
            == "below_minimum_length"
        )
        assert len(extract_regions(index=index, genome=genome)) == 4
        for kwargs in (
            {"direction": "bad"},
            {"length": 0},
            {"offset": -1},
            {"min_length": 0},
            {"max_ambiguous_fraction": 2},
            {"identifiers": ["missing"]},
            {"identifiers": ["plus", "plus"]},
        ):
            with pytest.raises(ValueError):
                extract_regions(index=index, genome=genome, **kwargs)


def test_soft_masking_and_terminal_bases(tmp_path):
    path = tmp_path / "genome.fa"
    write_fasta(path=path, records={"x": "aaNCAAAGTn"})
    gene = Gene(gene_id="g", contig="x", start=4, end=7, strand="+")
    with Genome(path=path) as genome:
        index = GeneIndex(genes=[gene], lengths=genome.lengths)
        regions = extract_regions(
            index=index, genome=genome, length=None, direction="both"
        )
        assert [r.sequence for r in regions] == ["AANC", "GTN"]
        assert (
            extract_regions(index=index, genome=genome, mask_lowercase=True)[
                0
            ].sequence
            == "NNNC"
        )
        assert (
            extract_regions(
                index=index, genome=genome, max_ambiguous_fraction=0.2
            )[0].status
            == "excess_ambiguity"
        )
