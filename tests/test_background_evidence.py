"""Shared-promoter controls, composition matching and optional evidence."""

from dataclasses import replace

import pytest

from intergenic_regions.background import (
    canonical_sequence,
    check_sequence_sets,
    collapse_overlaps,
    match_background,
    sanitise_regions,
)
from intergenic_regions.evidence import (
    PeakIndex,
    annotate_evidence,
    read_bed,
    read_functional_evidence,
)
from intergenic_regions.models import Region


def test_invalid_peak_query():
    index = PeakIndex(intervals={})
    for start, end in ((-1, 5), (5, 4)):
        with pytest.raises(ValueError, match="query interval"):
            index.overlap_bases(contig="x", start=start, end=end)


@pytest.fixture
def regions():
    return [
        Region(
            gene_id=f"g{i}",
            contig="x",
            start=start,
            end=start + len(sequence),
            strand="+",
            direction="upstream",
            sequence=sequence,
            status="retained",
            stop_reason="gene_boundary",
            available_length=len(sequence),
        )
        for i, (start, sequence) in enumerate(
            [
                (0, "AAAACCC"),
                (4, "AACCCCC"),
                (20, "ACGTGAA"),
                (40, "GGGGTAA"),
                (60, "CGTAAAA"),
            ]
        )
    ]


def test_sequence_sets():
    assert canonical_sequence(sequence="TTT") == "AAA"
    check_sequence_sets(positive={"p": "ACG"}, negative={"n": "AAA"})
    for positive, negative in [
        ({}, {"n": "AAA"}),
        ({"x": "ACG"}, {"x": "AAA"}),
        ({"p": "ACG"}, {"n": "CGT"}),
        ({"p1": "ACG", "p2": "CGT"}, {"n": "AAA"}),
    ]:
        with pytest.raises(ValueError):
            check_sequence_sets(positive=positive, negative=negative)


def test_overlap_collapse_and_cross_set_leakage(regions):
    collapsed, audit = collapse_overlaps(
        regions=regions
        + [
            replace(
                regions[0], gene_id="excluded", status="below_minimum_length"
            )
        ],
        label="positive",
    )
    assert len(collapsed) == 4
    assert audit[0]["reason"] == "overlapping_region"
    positive, negative, audit = sanitise_regions(
        positive=regions[:1],
        negative=regions[1:]
        + [replace(regions[0], gene_id="dup", start=100, end=107)],
    )
    assert len(positive) == 1 and len(negative) == 3
    assert {row["reason"] for row in audit} == {
        "duplicate_sequence",
        "overlaps_positive_region",
    }
    with pytest.raises(ValueError):
        sanitise_regions(positive=regions[:1], negative=regions[:1])
    positive, negative, audit = sanitise_regions(
        positive=[
            regions[0],
            replace(regions[0], gene_id="identical", start=200, end=207),
        ],
        negative=regions[2:3],
    )
    assert len(positive) == len(negative) == 1
    assert audit[0]["reason"] == "duplicate_sequence"


def test_background_matching_is_deterministic(regions):
    selected, pairs = match_background(
        positive=regions[:1], negative=regions[1:], max_gc_difference=0.5
    )
    assert selected[0].sequence_id == regions[2].sequence_id
    assert pairs[0]["length_ratio"] == 1
    assert (
        match_background(
            positive=regions[:1],
            negative=list(reversed(regions[1:])),
            max_gc_difference=0.5,
        )[1]
        == pairs
    )
    for kwargs in (
        {"ratio": 0},
        {"max_gc_difference": -1},
        {"max_length_ratio": 0.5},
        {"ratio": 10},
    ):
        with pytest.raises(ValueError):
            match_background(
                positive=regions[:1], negative=regions[2:], **kwargs
            )
    with pytest.raises(ValueError):
        match_background(positive=[], negative=regions)
    with pytest.raises(ValueError):
        match_background(
            positive=regions[:1], negative=regions[1:2], max_gc_difference=0
        )


def test_peak_union_and_overlap_boundaries(tmp_path):
    path = tmp_path / "peaks.bed"
    path.write_text(
        "track name=x\nbrowser position=x\n#comment\n"
        "x\t2\t5\tpeak1\nx\t4\t8\tpeak2\nx\t10\t12\tpeak3\n"
    )
    index = read_bed(path=path)
    assert index.overlap_bases(contig="x", start=0, end=12) == 8
    assert index.overlap_bases(contig="x", start=8, end=10) == 0
    assert index.overlap_bases(contig="unknown", start=0, end=10) == 0
    assert index.overlap_bases(contig="x", start=5, end=11) == 4


@pytest.mark.parametrize(
    "text", ["bad", "x\ta\tb", "x\t-1\t3", "x\t2\t2", "\t0\t5"]
)
def test_peak_failures(tmp_path, text):
    path = tmp_path / "bad.bed"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_bed(path=path)


def test_optional_functional_and_interval_evidence(tmp_path, regions):
    path = tmp_path / "evidence.tsv"
    path.write_text(
        "gene_id\tevidence_type\tvalue\tsource\ng0\treporter\tpositive\texperiment\n"
    )
    functional = read_functional_evidence(path=path)
    index = PeakIndex(intervals={"x": [(2, 5), (3, 6)]})
    rows = annotate_evidence(
        regions=regions,
        accessibility=index,
        enhancers=index,
        functional=functional,
    )
    assert rows[0]["accessibility_overlap_bp"] == 4
    assert rows[0]["linked_gene_evidence_count"] == 1
    assert rows[0]["enhancer_annotation_overlap_fraction"] == 4 / 7
    assert rows[2]["evidence_status"] == "no_overlap_in_supplied_data"
    absent = annotate_evidence(
        regions=[regions[0], replace(regions[1], status="excluded")]
    )
    assert absent[0]["evidence_status"] == "not_supplied"
    assert absent[0]["accessibility_overlap_bp"] is None
    assert len(absent) == 1
    for text in (
        "bad",
        "gene_id\tevidence_type\tvalue\tsource\ng\t\tx\tx\n",
        "gene_id\tevidence_type\tvalue\tsource\ng\tx\tx\tx\textra\n",
    ):
        path.write_text(text)
        with pytest.raises(ValueError):
            read_functional_evidence(path=path)
