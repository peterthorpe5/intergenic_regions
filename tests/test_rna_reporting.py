from dataclasses import replace

import pytest

from intergenic_regions.models import Region
from intergenic_regions.rna_reporting import (
    plot_rna_results,
    rna_motif_profiles,
)


def profile_inputs():
    regions = [
        Region(
            gene_id=identifier,
            contig="chr1",
            start=left,
            end=left + 10,
            strand=strand,
            direction="upstream",
            sequence="ACGTACGTAC",
            status="retained",
            stop_reason="gene_boundary",
            available_length=10,
        )
        for identifier, left, strand in [("p", 100, "+"), ("n", 200, "-")]
    ]
    return dict(
        regions=regions,
        anchors=[
            {
                "sequence_id": r.sequence_id,
                "applied_anchor": anchor,
                "anchor_source": source,
            }
            for r, anchor, source in zip(
                regions,
                [110, 199],
                ["rna_coverage_candidate", "annotation_fallback"],
                strict=True,
            )
        ],
        sites=[
            {
                "sequence_id": r.sequence_id,
                "motif_id": "m",
                "start": 0,
                "end": 3,
                "strand": strand,
            }
            for r in regions
            for strand in ("+", "-")
        ],
        motifs=[
            {
                "motif_id": "m",
                "consensus": "AAA",
                "q_value": 0.01,
                "positive_fraction": 1,
                "negative_fraction": 0,
            }
        ],
        labels={
            regions[0].sequence_id: "positive",
            regions[1].sequence_id: "negative",
        },
        bin_width=5,
        upstream=10,
    )


def test_motif_distance_actual_sites_denominators_and_anchor_sources():
    result = rna_motif_profiles(**profile_inputs())
    assert len(result) == 4
    assert sum(r["site_count"] for r in result) == 2
    assert all(r["eligible_motif_centres"] == 4 for r in result)
    assert all(r["site_count"] == 1 for r in result if r["bin_start"] == -10)
    assert {r["anchor_source"] for r in result} == {
        "rna_coverage_candidate",
        "annotation_fallback",
    }
    assert all(r["position_test_p_value"] is None for r in result)
    args = profile_inputs()
    args["motifs"][0]["q_value"] = 0.8
    assert rna_motif_profiles(**args) == []


@pytest.mark.parametrize(
    "settings",
    [
        {"bin_width": 0},
        {"upstream": 1},
        {"max_motifs": 0},
        {"q_value": 0},
        {"upstream": 10000, "bin_width": 1},
    ],
)
def test_motif_profile_bounds(settings):
    with pytest.raises(ValueError):
        rna_motif_profiles(**{**profile_inputs(), **settings})


def test_motif_profiles_reject_conflicting_sites_and_missing_metadata():
    args = profile_inputs()
    for changes in (
        {"anchors": []},
        {"anchors": args["anchors"] * 2},
        {"regions": args["regions"] * 2},
        {"labels": {}},
        {"sites": [{**args["sites"][0], "end": 11}]},
        {"motifs": [{**args["motifs"][0], "consensus": ""}]},
    ):
        with pytest.raises(ValueError):
            rna_motif_profiles(**{**args, **changes})
    excluded = replace(
        args["regions"][0], status="unusable_rna_anchor", sequence=""
    )
    result = rna_motif_profiles(**{**args, "regions": [excluded]})
    assert result == []


def test_plot_empty_inputs_and_limits(tmp_path):
    args = dict(
        profiles=[],
        starts=[],
        calls=[],
        shifts=[],
        de_rows=[],
        motif_profiles=[],
        labels={},
        directory=tmp_path,
    )
    assert plot_rna_results(**args) == []
    for changes in (
        {"max_loci": -1},
        {"max_loci": True},
        {"size_factors": {"a": 0}},
    ):
        with pytest.raises(ValueError):
            plot_rna_results(**{**args, **changes})
