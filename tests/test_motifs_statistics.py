"""Exact statistical references, finite-state PWM oracle and motif parsers."""

import itertools
import math

import numpy as np
import pytest
from scipy.stats import false_discovery_control, fisher_exact

from intergenic_regions.motifs import (
    Motif,
    analyse_motifs,
    kmer_counts,
    pattern_motif,
    pooled_background,
    prepare_pwm,
    read_motifs,
    scan_pattern,
    scan_pwm,
)
from intergenic_regions.statistics import (
    adjust_fdr,
    enrichment_test,
    kmer_family_size,
)


@pytest.mark.parametrize(
    "a,c,np_,nn",
    [(5, 1, 10, 10), (0, 0, 5, 7), (10, 0, 10, 10), (3, 6, 20, 40)],
)
def test_enrichment_matches_fisher(a, c, np_, nn):
    result = enrichment_test(
        positive_hits=a, negative_hits=c, positive_total=np_, negative_total=nn
    )
    expected = fisher_exact(
        table=[[a, np_ - a], [c, nn - c]], alternative="greater"
    ).pvalue
    assert result["p_value"] == pytest.approx(expected)
    assert (
        result["odds_ratio_ci_low"]
        < result["odds_ratio_haldane"]
        < result["odds_ratio_ci_high"]
    )
    assert (
        result["fold_enrichment"] is None
        if a and not c
        else result["fold_enrichment"] is not None
    )


def test_logspace_fdr_and_extreme_p_values():
    p_values = np.asarray([0.01, 0.04, 0.001, 0.9])
    actual = adjust_fdr(log_p_values=np.log(p_values).tolist())
    assert [r["q_value"] for r in actual] == pytest.approx(
        false_discovery_control(ps=p_values)
    )
    extended = false_discovery_control(
        ps=np.concatenate([p_values, np.ones(6)])
    )[:4]
    assert [
        r["q_value"]
        for r in adjust_fdr(
            log_p_values=np.log(p_values).tolist(), family_size=10
        )
    ] == pytest.approx(extended)
    extreme = enrichment_test(
        positive_hits=10000,
        negative_hits=0,
        positive_total=10000,
        negative_total=10000,
    )
    assert math.isfinite(extreme["log_p_value"])
    assert extreme["minus_log10_p"] > 1000
    assert (
        adjust_fdr(log_p_values=[-10000], family_size=100)[0]["log_q_value"]
        < -9000
    )
    assert adjust_fdr(log_p_values=[]) == []


@pytest.mark.parametrize(
    "kwargs",
    [
        {"positive_hits": -1},
        {"negative_hits": 11},
        {"positive_total": 0},
        {"negative_total": 0},
    ],
)
def test_statistical_failures(kwargs):
    options = dict(
        positive_hits=1, negative_hits=1, positive_total=10, negative_total=10
    )
    options.update(kwargs)
    with pytest.raises(ValueError):
        enrichment_test(**options)


@pytest.mark.parametrize(
    "logs,family",
    [
        ([0.1], 1),
        ([float("nan")], 1),
        ([float("-inf")], 1),
        ([-1, -2], 1),
        ([], -1),
    ],
)
def test_fdr_failures(logs, family):
    with pytest.raises(ValueError):
        adjust_fdr(log_p_values=logs, family_size=family)


def test_kmer_family_and_counts():
    assert kmer_family_size(lengths=[2]) == 10
    assert kmer_family_size(lengths=[3]) == 32
    assert kmer_family_size(lengths=[2, 3], both_strands=False) == 80
    for lengths in ([1], [11], [4, 4]):
        with pytest.raises(ValueError):
            kmer_family_size(lengths=lengths)
    assert kmer_counts(sequence="AACGNTT", lengths=[2]) == {
        "AA": 2,
        "AC": 1,
        "CG": 1,
    }
    assert kmer_counts(sequence="ACN", lengths=[2], both_strands=False) == {
        "AC": 1
    }


def test_pattern_matching_overlaps_reverse_and_ambiguity():
    assert (
        len(scan_pattern(sequence="AAAA", pattern="AAA", both_strands=False))
        == 2
    )
    hits = scan_pattern(sequence="ACGNNNCGT", pattern="ACG")
    assert {(r["start"], r["site_strand"]) for r in hits} == {
        (0, "+"),
        (6, "-"),
    }
    assert scan_pattern(sequence="NNN", pattern="NNN") == []
    assert len(scan_pattern(sequence="CACGTG", pattern="CACGTG")) == 1
    with pytest.raises(ValueError):
        pattern_motif(motif_id="bad", pattern="AZ")


@pytest.mark.parametrize(
    "matrix",
    [
        (),
        ((1, 0),),
        ((1, -1, 1, 0),),
        ((0.5, 0.5, 0.5, 0.5),),
        ((float("nan"), 0, 0, 1),),
    ],
)
def test_motif_validation(matrix):
    with pytest.raises(ValueError):
        Motif(motif_id="m", name="m", matrix=matrix)
    with pytest.raises(ValueError):
        Motif(motif_id="bad id", name="m", matrix=((1, 0, 0, 0),))


def test_pwm_distribution_matches_enumerated_word_oracle():
    motif = Motif(
        motif_id="m",
        name="m",
        matrix=(
            (0.7, 0.1, 0.1, 0.1),
            (0.1, 0.4, 0.4, 0.1),
            (0.1, 0.1, 0.1, 0.7),
        ),
    )
    background = np.asarray([0.3, 0.2, 0.2, 0.3])
    scanner = prepare_pwm(
        motif=motif, background=background, site_p_value=0.15
    )
    scores = []
    probabilities = []
    for indices in itertools.product(range(4), repeat=3):
        scores.append(
            sum(
                int(scanner.weights[i, base]) for i, base in enumerate(indices)
            )
        )
        probabilities.append(math.prod(background[base] for base in indices))
    for score in set(scores):
        expected = sum(
            p for s, p in zip(scores, probabilities, strict=True) if s >= score
        )
        assert scanner.tail[score - scanner.minimum] == pytest.approx(expected)
    for word in ("AAT", "ACT", "AGT", "TTT", "NNN"):
        hits = scan_pwm(sequence=word, scanner=scanner, both_strands=False)
        if "N" in word:
            assert not hits
        else:
            score = sum(
                int(scanner.weights[i, "ACGT".index(base)])
                for i, base in enumerate(word)
            )
            assert bool(hits) == (score >= scanner.threshold)
    assert scan_pwm(sequence="A", scanner=scanner) == []
    strict = prepare_pwm(
        motif=motif, background=background, site_p_value=1e-10
    )
    assert not scan_pwm(sequence="ACTACT", scanner=strict)


def test_pwm_reverse_palindrome_and_settings():
    motif = pattern_motif(motif_id="pal", pattern="ACGT")
    background = pooled_background(sequences=["AAAACCCC", "TGGGNN"])
    assert np.allclose(background, background[::-1])
    scanner = prepare_pwm(
        motif=motif, background=background, site_p_value=0.02
    )
    assert scanner.palindrome
    assert len(scan_pwm(sequence="ACGT", scanner=scanner)) == 1
    non_palindrome = prepare_pwm(
        motif=pattern_motif(motif_id="x", pattern="AAC"),
        background=np.full(4, 0.25),
        site_p_value=0.02,
    )
    assert (
        scan_pwm(sequence="GTT", scanner=non_palindrome)[0]["site_strand"]
        == "-"
    )
    for kwargs in (
        {"site_p_value": 0},
        {"resolution": 0},
        {"pseudocount": 0},
        {"background": np.asarray([0.3, 0.2, 0.3, 0.2])},
    ):
        options = dict(motif=motif, background=background)
        options.update(kwargs)
        with pytest.raises(ValueError):
            prepare_pwm(**options)


def test_all_motif_formats(tmp_path):
    meme = tmp_path / "motifs.meme"
    meme.write_text(
        "MEME version 4\nALPHABET= ACGT\nMOTIF x X\n"
        "letter-probability matrix: alength= 4 w= 2\n1 0 0 0\n0 1 0 0\n"
    )
    assert read_motifs(path=meme)[0].matrix == ((1, 0, 0, 0), (0, 1, 0, 0))
    jaspar = tmp_path / "motifs.jaspar"
    jaspar.write_text(">x X\nA [ 10 0 ]\nC [ 0 10 ]\nG [ 0 0 ]\nT [ 0 0 ]\n")
    assert (
        read_motifs(path=jaspar)[0].matrix == read_motifs(path=meme)[0].matrix
    )
    iupac = tmp_path / "motifs.tsv"
    iupac.write_text("motif_id\tpattern\nx\tACY\n")
    assert read_motifs(path=iupac)[0].pattern == "ACY"


@pytest.mark.parametrize(
    "text,format_name",
    [
        ("", "auto"),
        ("bad", "auto"),
        ("x\n", "jaspar"),
        (">x\nA [ 1 ]", "jaspar"),
        (">x\nA [ 1 ]\nC [ 0 1 ]\nG [ 0 ]\nT [ 0 ]", "jaspar"),
        (">x\nA [ -1 ]\nC [ 0 ]\nG [ 0 ]\nT [ 0 ]", "jaspar"),
        ("MEME version 4\nALPHABET= ACGU", "meme"),
        ("MEME version 4\nMOTIF x\nbad", "meme"),
        (
            "MEME version 4\nMOTIF x\nletter-probability matrix: alength=4",
            "meme",
        ),
        (
            "MEME version 4\nMOTIF x\nletter-probability matrix: w=2\n1 0 0 0",
            "meme",
        ),
        ("motif_id\tpattern\nx\tAC\nx\tGT\n", "iupac"),
        ("bad", "unknown"),
        ("bad", "meme"),
    ],
)
def test_motif_parser_failures(tmp_path, text, format_name):
    path = tmp_path / "bad"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_motifs(path=path, motif_format=format_name)


def test_planted_motif_and_full_testing_family(dataset):
    from intergenic_regions.io import read_fasta

    positive = read_fasta(path=dataset["positive"])
    negative = read_fasta(path=dataset["negative"])
    motif = pattern_motif(motif_id="Gbox", pattern="CACGTG")
    rows, sites, summary = analyse_motifs(
        positive=positive, negative=negative, motifs=[motif], kmer_lengths=[6]
    )
    row = next(r for r in rows if r["motif_id"] == "Gbox")
    assert row["positive_hits"] == 15
    assert row["q_value"] < 0.05
    assert summary["testing_family_size"] == kmer_family_size(lengths=[6]) + 1
    assert any(site["motif_id"] == "Gbox" for site in sites)
    rows, _, _ = analyse_motifs(
        positive=positive,
        negative=negative,
        motifs=[Motif(motif_id="pwm", name="PWM", matrix=motif.matrix)],
        site_p_value=0.001,
    )
    assert rows[0]["positive_hits"] == 15
    with pytest.raises(ValueError):
        analyse_motifs(positive=positive, negative=negative)
    with pytest.raises(ValueError):
        analyse_motifs(
            positive=positive, negative=negative, motifs=[motif, motif]
        )
