import random

import pytest

from intergenic_regions.extraction import GeneIndex, extract_regions
from intergenic_regions.genome import Genome, reverse_complement
from intergenic_regions.io import write_fasta
from intergenic_regions.models import Gene
from intergenic_regions.rna_starts import (
    compare_rna_starts,
    consensus_rna_starts,
    detect_coverage_start,
    extract_rna_upstream,
    rna_search_interval,
)


def gene(*, identifier, start, end, strand="+"):
    return Gene(
        gene_id=identifier, contig="chr1", start=start, end=end, strand=strand
    )


def edge(*, position=100, replicate="1", condition="treated", strand="+"):
    return {
        "gene_id": "g",
        "contig": "chr1",
        "strand": strand,
        "annotation_anchor": 200,
        "condition": condition,
        "replicate_id": replicate,
        "strandedness": "forward",
        "status": "candidate",
        "candidate_start": position,
        "neighbourhood_start": position - 20,
        "neighbourhood_end": position + 21,
        "fold_rise": 21.0,
    }


@pytest.mark.parametrize("strand", ["+", "-"])
def test_directional_coverage_edge(strand):
    oriented = [0] * 80 + [10] * 90
    result = detect_coverage_start(
        depth=oriented if strand == "+" else oriented[::-1],
        genomic_start=1000,
        strand=strand,
    )
    assert result["candidate_start"] == (1080 if strand == "+" else 1089)
    assert result["status"] == "candidate"
    assert result["fold_rise"] == 21
    assert result["neighbourhood_end"] - result["neighbourhood_start"] == 41
    assert "p_value" not in result


@pytest.mark.parametrize(
    "depth,status",
    [
        ([], "insufficient_search_space"),
        ([0] * 90, "not_detected"),
        ([10] * 90, "not_detected"),
        ([0] * 40 + [10] * 5 + [0] * 80, "not_detected"),
    ],
)
def test_depth_is_not_start_significance(depth, status):
    result = detect_coverage_start(depth=depth, genomic_start=0, strand="+")
    assert result["status"] == status
    assert result["candidate_start"] is None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"strand": "."},
        {"genomic_start": -1},
        {"window": 1},
        {"window": True},
        {"window": 1001},
        {"min_depth": 0},
        {"min_fold_rise": 1},
        {"min_covered_fraction": 0},
        {"pseudocount": 0},
        {"depth": [float("nan")]},
        {"depth": [-1]},
    ],
)
def test_invalid_edge_settings(kwargs):
    with pytest.raises(ValueError):
        detect_coverage_start(
            **{"depth": [0] * 80, "genomic_start": 0, "strand": "+", **kwargs}
        )


def test_multiple_edges_are_explicit():
    result = detect_coverage_start(
        depth=[0] * 80 + [10] * 80 + [0] * 80 + [20] * 80,
        genomic_start=0,
        strand="+",
    )
    assert result["candidate_start"] == 240
    assert result["alternative_edges"] == 1


def test_replicate_consensus_does_not_chain_or_miss_dense_clusters():
    calls = [
        edge(position=p, replicate=str(i))
        for i, p in enumerate([100, 130, 145, 155])
    ]
    result = consensus_rna_starts(
        calls=calls, agreement_tolerance=40, min_support_fraction=0.75
    )[0]
    assert result["supporting_replicates"] == 3
    assert result["candidate_start"] == 145
    assert result["replicate_start_max"] - result["replicate_start_min"] == 25
    assert result["status"] == "replicate_supported_candidate"
    chained = consensus_rna_starts(
        calls=[
            edge(position=p, replicate=str(i))
            for i, p in enumerate([100, 130, 160])
        ],
        min_support_fraction=1,
        agreement_tolerance=40,
    )[0]
    assert chained["status"] == "discordant_replicates"
    assert chained["supporting_replicates"] == 2


def test_singletons_missing_edges_and_unstranded_uncertainty():
    single = consensus_rna_starts(calls=[edge()])[0]
    assert single["status"] == "provisional_single_replicate"
    missing = {
        **edge(replicate="2"),
        "status": "not_detected",
        "candidate_start": None,
    }
    result = consensus_rna_starts(calls=[edge(), missing])[0]
    assert result["status"] == "insufficient_replicate_support"
    assert consensus_rna_starts(calls=[missing])[0]["status"] == "not_detected"
    ambiguous = consensus_rna_starts(
        calls=[{**edge(), "strandedness": "unstranded"}, edge(replicate="2")]
    )[0]
    assert ambiguous["strand_evidence"] == "unstranded_ambiguous"
    assert ambiguous["tss_status"] == "unvalidated_coverage_edge"


@pytest.mark.parametrize(
    "changes",
    [
        {"candidate_start": None},
        {"candidate_start": True},
        {"candidate_start": -1},
        {"candidate_start": 1.5},
        {"neighbourhood_start": 101},
        {"neighbourhood_end": 100},
        {"fold_rise": float("nan")},
        {"strand": "."},
    ],
)
def test_consensus_invalid_candidate_records(changes):
    with pytest.raises(ValueError):
        consensus_rna_starts(calls=[{**edge(), **changes}])


def test_consensus_invalid_observations_and_limits():
    for calls in (
        [edge(), edge()],
        [edge(), {**edge(replicate="2"), "contig": "chr2"}],
    ):
        with pytest.raises(ValueError):
            consensus_rna_starts(calls=calls)
    for settings in (
        {"min_replicates": 1},
        {"agreement_tolerance": -1},
        {"min_support_fraction": 0.49},
    ):
        with pytest.raises(ValueError):
            consensus_rna_starts(calls=[], **settings)
    assert consensus_rna_starts(calls=[]) == []


def test_start_shifts_are_oriented_and_descriptive():
    calls = [
        edge(condition=condition, position=pos, strand="-", replicate=rep)
        for condition, pos in [("ref", 110), ("treated", 100)]
        for rep in ("1", "2")
    ]
    starts = consensus_rna_starts(calls=calls)
    result = compare_rna_starts(
        starts=starts, condition="treated", reference="ref"
    )[0]
    assert result["oriented_shift_bp"] == 10
    assert result["neighbourhoods_overlap"] is True
    assert result["promoter_usage_p_value"] is None
    starts[0]["status"] = "not_detected"
    assert (
        compare_rna_starts(
            starts=starts, condition="treated", reference="ref"
        )[0]["oriented_shift_bp"]
        is None
    )
    for kwargs in (
        {"condition": "treated", "reference": "treated"},
        {"condition": "absent", "reference": "ref"},
    ):
        with pytest.raises(ValueError):
            compare_rna_starts(starts=starts, **kwargs)
    with pytest.raises(ValueError):
        compare_rna_starts(
            starts=[*starts, starts[0]], condition="treated", reference="ref"
        )


@pytest.fixture
def safe_genome(*, tmp_path):
    genes = [
        gene(identifier="left", start=100, end=160, strand="-"),
        gene(identifier="p", start=400, end=500),
        gene(identifier="middle", start=600, end=660),
        gene(identifier="n", start=800, end=900, strand="-"),
        gene(identifier="right", start=1050, end=1100),
    ]
    rng = random.Random(91)
    path = tmp_path / "genome.fa"
    write_fasta(
        path=path, records={"chr1": "".join(rng.choices("ACGT", k=1300))}
    )
    return path, genes


def test_search_stops_at_opposite_strand_neighbours_and_internal_genes(
    safe_genome,
):
    path, genes = safe_genome
    with Genome(path=path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        assert rna_search_interval(index=index, gene=index.genes["p"]) == (
            160,
            500,
            "gene_span_fallback",
        )
        assert rna_search_interval(
            index=index, gene=index.genes["n"], exons=[(800, 820), (870, 900)]
        ) == (870, 1050, "first_exonic_union")
        with pytest.raises(ValueError):
            rna_search_interval(
                index=index, gene=index.genes["p"], exons=[(390, 410)]
            )
        with pytest.raises(ValueError):
            rna_search_interval(index=index, gene=index.genes["p"], upstream=0)
        inner = [
            *genes,
            gene(identifier="nestedp", start=430, end=440),
            gene(identifier="nestedn", start=860, end=875),
        ]
        index = GeneIndex(genes=inner, lengths=genome.lengths)
        assert (
            rna_search_interval(index=index, gene=index.genes["p"])[1] == 430
        )
        assert (
            rna_search_interval(index=index, gene=index.genes["n"])[0] == 875
        )
        overlapping = GeneIndex(
            genes=[*genes, gene(identifier="over", start=390, end=420)],
            lengths=genome.lengths,
        )
        first, last, source = rna_search_interval(
            index=overlapping, gene=overlapping.genes["p"]
        )
        assert first == last and source == "overlapping_gene"


def anchors(
    *, positive=350, negative=950, status="replicate_supported_candidate"
):
    return [
        {
            "gene_id": identifier,
            "condition": "treated",
            "contig": "chr1",
            "strand": strand,
            "candidate_start": position,
            "status": status,
            "supporting_replicates": 2,
            "strand_evidence": "strand_specific",
            "tss_status": "unvalidated_coverage_edge",
        }
        for identifier, strand, position in [
            ("p", "+", positive),
            ("n", "-", negative),
        ]
    ]


def test_rna_anchor_excludes_transcribed_base_and_reverse_complements(
    safe_genome,
):
    path, genes = safe_genome
    with Genome(path=path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        regions, audit = extract_rna_upstream(
            index=index,
            genome=genome,
            starts=anchors(),
            condition="treated",
            identifiers=["p", "n"],
            length=None,
        )
        assert [(r.start, r.end) for r in regions] == [(160, 350), (951, 1050)]
        assert regions[1].sequence == reverse_complement(
            sequence=genome.fetch(contig="chr1", start=951, end=1050)
        )
        assert all(
            r["anchor_source"] == "rna_coverage_candidate" for r in audit
        )
        for region in regions:
            assert not any(
                region.start < g.end and region.end > g.start for g in genes
            )
        shifted, _ = extract_rna_upstream(
            index=index,
            genome=genome,
            starts=anchors(),
            condition="treated",
            identifiers=["p", "n"],
            length=20,
            offset=5,
        )
        assert [(r.start, r.end) for r in shifted] == [(325, 345), (956, 976)]
        assert index.genes == {g.gene_id: g for g in genes}


@pytest.mark.parametrize(
    "positive,negative,reason",
    [
        (410, 880, "rna_start_inside_annotated_target"),
        (150, 1060, "rna_start_outside_contiguous_safe_gap"),
    ],
)
def test_rna_start_cannot_authorise_genic_extraction(
    safe_genome, positive, negative, reason
):
    path, genes = safe_genome
    with Genome(path=path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        args = dict(
            index=index,
            genome=genome,
            starts=anchors(positive=positive, negative=negative),
            condition="treated",
            identifiers=["p", "n"],
            length=None,
        )
        regions, audit = extract_rna_upstream(**args)
        assert regions == extract_regions(
            index=index, genome=genome, identifiers=["p", "n"], length=None
        )
        assert {r["anchor_reason"] for r in audit} == {reason}
        excluded, audit = extract_rna_upstream(**args, fallback="exclude")
        assert all(r.status == "unusable_rna_anchor" for r in excluded)
        assert all(r["applied_anchor"] is None for r in audit)


def test_rna_anchor_provisional_unknown_and_input_validation(safe_genome):
    path, genes = safe_genome
    with Genome(path=path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        args = dict(
            index=index,
            genome=genome,
            starts=anchors(status="provisional_single_replicate"),
            condition="treated",
            identifiers=["p", "n"],
        )
        _, audit = extract_rna_upstream(**args)
        assert all(r["anchor_source"] == "annotation_fallback" for r in audit)
        _, audit = extract_rna_upstream(**args, allow_provisional=True)
        assert all(
            r["anchor_source"] == "rna_coverage_candidate" for r in audit
        )
        for updates in (
            {"fallback": "guess"},
            {"condition": "absent"},
            {"starts": [*anchors(), anchors()[0]]},
            {"starts": anchors(positive=-1)},
            {"starts": [{**anchors()[0], "strand": "-"}]},
        ):
            with pytest.raises(ValueError):
                extract_rna_upstream(**{**args, **updates})
