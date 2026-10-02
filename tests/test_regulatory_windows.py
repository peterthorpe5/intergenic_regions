"""Variable scales, direction-aware coordinates and explicit uncertainty."""

import random
from dataclasses import asdict, replace

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from intergenic_regions.extraction import GeneIndex, extract_regions
from intergenic_regions.genome import Genome, reverse_complement
from intergenic_regions.io import write_fasta
from intergenic_regions.models import Gene, Region
from intergenic_regions.regional_analysis import (
    window_motif_counts,
    window_position_profiles,
)
from intergenic_regions.windows import (
    RegulatoryWindow,
    assign_window_groups,
    make_windows,
    merge_window_candidates,
)


def window(**changes):
    values = dict(
        parent_id="p",
        label="positive",
        sequence_start=0,
        sequence_end=20,
        sequence="ACGT" * 5,
    )
    values.update(changes)
    return RegulatoryWindow(**values)


def row(item, score=0.9):
    record = asdict(item)
    record.pop("sequence")
    record.update(
        sequence_id=item.sequence_id,
        window_length=item.width,
        held_out_signature_score=score,
        motif_sites_per_kb=10.0,
    )
    return record


def test_multiscale_windows_complete_end_aligned_and_audited():
    windows, audit = make_windows(
        positive={"p": "ACGT" * 25 + "AAA"},
        negative={"n": "GGTA" * 12},
        lengths=[20, 50, 100, 200],
        step=20,
    )
    assert len({w.sequence_id for w in windows}) == len(windows)
    assert all(w.width == len(w.sequence) for w in windows)
    assert {
        w.sequence_start
        for w in windows
        if w.parent_id == "p" and w.width == 100
    } == {0, 3}
    assert max(w.sequence_end for w in windows if w.parent_id == "p") == 103
    assert max(w.sequence_end for w in windows if w.parent_id == "n") == 48
    assert all(
        w.genomic_start is None and w.distance_to_gene_start is None
        for w in windows
    )
    assert (
        next(
            r
            for r in audit
            if r["parent_id"] == "n" and r["window_length"] == 50
        )["status"]
        == "parent_too_short"
    )
    assert sum(r["windows"] for r in audit) == len(windows)


@given(
    length=st.integers(20, 300),
    width=st.integers(20, 150),
    step=st.integers(1, 20),
)
@settings(max_examples=35)
def test_sliding_windows_cover_every_base_without_truncation(
    length, width, step
):
    windows, _ = make_windows(
        positive={"p": "A" * length},
        negative={"n": "C" * length},
        lengths=[width],
        step=step,
    )
    selected = [w for w in windows if w.parent_id == "p"]
    if width > length:
        assert selected == []
    else:
        covered = {
            i
            for w in selected
            for i in range(w.sequence_start, w.sequence_end)
        }
        assert covered == set(range(length))
        assert all(
            w.sequence_end - w.sequence_start == width for w in selected
        )


def test_original_direction_aware_gene_boundaries_and_reverse_complement(
    tmp_path,
):
    sequence = "".join(random.Random(8).choices("ACGT", k=380))
    path = tmp_path / "genome.fa"
    write_fasta(path=path, records={"chr": sequence})
    genes = [
        Gene(gene_id="left", contig="chr", start=0, end=17, strand="-"),
        Gene(gene_id="p", contig="chr", start=120, end=150, strand="+"),
        Gene(gene_id="middle", contig="chr", start=150, end=180, strand="+"),
        Gene(gene_id="n", contig="chr", start=180, end=210, strand="-"),
        Gene(gene_id="right", contig="chr", start=330, end=380, strand="+"),
    ]
    with Genome(path=path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        regions = extract_regions(
            index=index,
            genome=genome,
            identifiers=["p", "n"],
            length=None,
            offset=7,
        )
    assert [(r.start, r.end) for r in regions] == [(17, 113), (217, 330)]
    windows, audit = make_windows(
        positive={regions[0].sequence_id: regions[0].sequence},
        negative={regions[1].sequence_id: regions[1].sequence},
        regions=regions,
        genes=genes,
        lengths=[20, 50, 100, 200],
        step=13,
    )
    for item in windows:
        source = regions[0] if item.label == "positive" else regions[1]
        assert (
            source.start <= item.genomic_start < item.genomic_end <= source.end
        )
        genomic = sequence[item.genomic_start : item.genomic_end]
        expected = (
            genomic
            if item.strand == "+"
            else reverse_complement(sequence=genomic)
        )
        assert item.sequence == expected
        assert not any(
            item.genomic_start < g.end and g.start < item.genomic_end
            for g in genes
        )
        assert item.distance_to_gene_start < 0
    assert all(r["windows"] == 0 for r in audit if r["window_length"] == 200)
    negative_first = next(
        w for w in windows if w.label == "negative" and w.width == 20
    )
    assert (negative_first.genomic_start, negative_first.genomic_end) == (
        310,
        330,
    )
    assert negative_first.distance_to_gene_start == pytest.approx(-110.5)


@pytest.mark.parametrize(
    "changes",
    [
        {"lengths": []},
        {"lengths": [19]},
        {"lengths": [100001]},
        {"lengths": [20, 20]},
        {"lengths": [True]},
        {"lengths": [20.0]},
        {"step": 0},
        {"step": 21},
        {"step": True},
        {"max_windows": 0},
        {"max_windows": False},
        {"max_windows": 1},
    ],
)
def test_window_generation_rejects_bad_settings_and_over_allocation(changes):
    options = dict(
        positive={"p": "A" * 100},
        negative={"n": "C" * 100},
        lengths=[20],
        step=10,
    )
    options.update(changes)
    with pytest.raises(ValueError):
        make_windows(**options)


@pytest.mark.parametrize(
    "changes",
    [
        {"parent_id": ""},
        {"parent_id": "has space"},
        {"label": "unknown"},
        {"strand": "?"},
        {"sequence_start": -1},
        {"sequence_end": 0},
        {"sequence": "A"},
        {"genomic_start": 0},
        {"contig": "chr"},
        {
            "genomic_start": -1,
            "genomic_end": 19,
            "contig": "chr",
            "gene_id": "g",
            "strand": "+",
            "direction": "upstream",
        },
        {
            "genomic_start": 0,
            "genomic_end": 21,
            "contig": "chr",
            "gene_id": "g",
            "strand": "+",
            "direction": "upstream",
        },
        {"genomic_start": 0, "genomic_end": 20, "contig": "chr"},
        {"distance_to_gene_start": 1},
    ],
)
def test_window_model_rejects_inconsistent_metadata(changes):
    with pytest.raises(ValueError):
        window(**changes)


def test_window_coordinates_and_parent_metadata_validation():
    source = Region(
        gene_id="p",
        contig="chr",
        start=100,
        end=120,
        strand="+",
        direction="upstream",
        sequence="A" * 20,
        status="retained",
        stop_reason="gene",
        available_length=20,
    )
    kwargs = dict(
        positive={source.sequence_id: source.sequence},
        negative={"n|upstream": "C" * 20},
        lengths=[20],
        step=20,
    )
    for sources in ([source], [source, source], []):
        with pytest.raises(ValueError, match="flank metadata"):
            make_windows(**kwargs, regions=sources)
    other = replace(source, gene_id="n", sequence="C" * 20)
    for bad in (
        replace(source, status="excluded"),
        replace(source, strand="."),
        replace(source, sequence="G" * 20),
        replace(source, end=121),
    ):
        with pytest.raises(ValueError, match="metadata does not match"):
            make_windows(**kwargs, regions=[bad, other])
    with pytest.raises(ValueError, match="metadata disagree"):
        make_windows(
            **kwargs,
            regions=[source, other],
            genes=[
                Gene(
                    gene_id="p",
                    contig="different",
                    start=120,
                    end=140,
                    strand="+",
                )
            ],
        )
    genomic = window(
        contig="chr",
        gene_id="p",
        genomic_start=100,
        genomic_end=120,
        strand="+",
        direction="upstream",
    )
    with pytest.raises(ValueError, match="distance"):
        replace(genomic, distance_to_gene_start=float("inf"))


def test_groups_union_reverse_complement_windows_and_user_families():
    sequence = "AAAACGCGTTTTAGACGACT"
    items = [
        window(parent_id="a", sequence=sequence),
        window(
            parent_id="b",
            label="negative",
            sequence=reverse_complement(sequence=sequence),
        ),
        window(parent_id="c", sequence="C" * 20),
    ]
    assert assign_window_groups(windows=items) == {
        "a": "a",
        "b": "a",
        "c": "c",
    }
    assert set(
        assign_window_groups(
            windows=items,
            parent_groups={"a": "family", "b": "separate", "c": "family"},
        ).values()
    ) == {"a"}
    assert assign_window_groups(windows=[]) == {}
    for groups in ({}, {"a": "family", "b": "", "c": "other"}):
        with pytest.raises(ValueError, match="non-empty group"):
            assign_window_groups(windows=items, parent_groups=groups)


def test_window_unions_variable_length_provenance_and_no_claimed_fdr():
    items = [
        row(window(sequence_end=100, sequence="A" * 100)),
        row(
            window(sequence_start=80, sequence_end=180, sequence="A" * 100),
            0.8,
        ),
        row(window(sequence_start=180, sequence_end=200), 0.9),
        row(window(parent_id="n", label="negative"), 0.3),
        row(window(parent_id="unscored"), None),
    ]
    candidates = merge_window_candidates(rows=items)
    assert [(r["sequence_start"], r["sequence_end"]) for r in candidates] == [
        (0, 180),
        (180, 200),
    ]
    assert candidates[0]["supporting_windows"] == 2
    assert candidates[0]["length"] == 180
    assert candidates[0]["region_q_value"] is None
    assert "not_functional_boundary" in candidates[0]["boundary_status"]
    assert candidates[0]["genomic_start"] is None
    assert "boundary_status" not in items[0]
    genomic_rows = [
        row(
            window(
                contig="chr",
                gene_id="g",
                strand="-",
                direction="upstream",
                genomic_start=200,
                genomic_end=220,
            )
        ),
        row(
            window(
                contig="chr",
                gene_id="g",
                strand="-",
                direction="upstream",
                sequence_start=10,
                sequence_end=30,
                genomic_start=190,
                genomic_end=210,
            )
        ),
    ]
    assert (
        merge_window_candidates(rows=genomic_rows)[0]["genomic_start"] == 190
    )


@pytest.mark.parametrize("threshold", [-1, 2, float("nan")])
def test_bad_candidate_threshold(threshold):
    with pytest.raises(ValueError):
        merge_window_candidates(rows=[], score_threshold=threshold)


@pytest.mark.parametrize(
    "changes",
    [
        {"held_out_signature_score": "bad"},
        {"held_out_signature_score": float("nan")},
        {"held_out_signature_score": 1.2},
        {"parent_id": ""},
        {"sequence_start": -1},
        {"sequence_end": 0},
    ],
)
def test_bad_scored_window(changes):
    record = row(window())
    record.update(changes)
    with pytest.raises(ValueError):
        merge_window_candidates(rows=[record])


def test_motif_density_physical_sites_containment_and_diversity():
    items = [window(), window(sequence_start=10, sequence_end=30)]
    sites = [
        {
            "sequence_id": "p",
            "motif_id": "A",
            "start": 0,
            "end": 6,
            "strand": "+",
        },
        {
            "sequence_id": "p",
            "motif_id": "A",
            "start": 0,
            "end": 6,
            "strand": "-",
        },
        {"sequence_id": "p", "motif_id": "B", "start": 14, "end": 20},
        {"sequence_id": "p", "motif_id": "A", "start": 18, "end": 24},
        {"sequence_id": "p", "motif_id": "ignored", "start": 0, "end": 6},
    ]
    result = window_motif_counts(
        windows=items, sites=sites, motif_ids=["A", "B"]
    )
    assert result[items[0].sequence_id] == {
        "motif_sites": 2,
        "motif_diversity": 2,
        "motif_sites_per_kb": 100.0,
    }
    assert result[items[1].sequence_id]["motif_sites"] == 2
    assert (
        window_motif_counts(windows=items, sites=[], motif_ids=[])[
            items[0].sequence_id
        ]["motif_sites"]
        == 0
    )
    with pytest.raises(ValueError):
        window_motif_counts(windows=[], sites=[], motif_ids=["A", "A"])
    with pytest.raises(ValueError):
        window_motif_counts(
            windows=[],
            sites=[
                {"sequence_id": "p", "motif_id": "A", "start": -1, "end": 2}
            ],
            motif_ids=["A"],
        )


def test_position_profiles_average_sources_not_correlated_windows():
    items = [
        row(window(), 0.1),
        row(window(sequence_start=10, sequence_end=30), 0.3),
        row(window(parent_id="other"), 0.8),
    ]
    profile = window_position_profiles(rows=items)[0]
    assert profile["parents"] == profile["scored_parents"] == 2
    assert profile["mean_held_out_signature_score"] == pytest.approx(0.5)
    assert profile["coordinate_system"] == "oriented_sequence_offset"
    for record in items:
        record["held_out_signature_score"] = None
        record["distance_to_gene_start"] = -10.5
    profile = window_position_profiles(rows=items)[0]
    assert profile["mean_held_out_signature_score"] is None
    assert (profile["bin_start"], profile["bin_end"]) == (-100, 0)
    assert profile["coordinate_system"] == "annotated_gene_5prime_base"
    assert window_position_profiles(rows=items, limit=5, bin_width=5) == []
    items[0]["distance_to_gene_start"] = None
    with pytest.raises(ValueError, match="mix coordinate"):
        window_position_profiles(rows=items)


@pytest.mark.parametrize(
    "changes",
    [
        {"bin_width": 0},
        {"bin_width": True},
        {"limit": 10},
        {"limit": True},
        {"bin_width": 1, "limit": 1000},
    ],
)
def test_invalid_profile_bins(changes):
    with pytest.raises(ValueError):
        window_position_profiles(rows=[], **changes)
