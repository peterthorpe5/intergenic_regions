"""Genome-wide substitution screens, coordinates and opportunity profiles."""

import csv
import json

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from intergenic_regions.cli import build_parser, dispatch, scan_options
from intergenic_regions.genome import reverse_complement
from intergenic_regions.io import write_fasta, write_tsv
from intergenic_regions.motifs import IUPAC
from intergenic_regions.scan_reporting import plot_scan, scan_report
from intergenic_regions.scanning import (
    ScanMotif,
    consensus_sites,
    genome_scan_outputs,
    read_scan_motifs,
)
from intergenic_regions.workflows import (
    genome_scan_workflow,
    pipeline_workflow,
    scan_targets,
)


def records(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


@given(
    sequence=st.text(alphabet="ACGTN", min_size=0, max_size=45),
    pattern=st.text(alphabet="ACGTRYN", min_size=1, max_size=8),
    requested=st.integers(min_value=0, max_value=7),
    both=st.booleans(),
)
@settings(max_examples=100, deadline=None)
def test_vectorised_sites_agree_with_brute_force(
    sequence,
    pattern,
    requested,
    both,
):
    limit = min(requested, len(pattern) - 1)
    expected = []
    reverse = reverse_complement(sequence=pattern)
    for start in range(len(sequence) - len(pattern) + 1):
        word = sequence[start : start + len(pattern)]
        if "N" in word:
            continue
        forward = sum(
            base not in IUPAC[letter]
            for base, letter in zip(word, pattern, strict=True)
        )
        backward = sum(
            base not in IUPAC[letter]
            for base, letter in zip(word, reverse, strict=True)
        )
        orientations = [
            strand
            for distance, strand in ((forward, "+"), (backward, "-"))
            if distance <= limit and (strand == "+" or both)
        ]
        if orientations:
            expected.append(
                dict(
                    start=start,
                    end=start + len(pattern),
                    site_strand="."
                    if len(orientations) == 2
                    else orientations[0],
                    mismatches=min(forward, backward) if both else forward,
                )
            )
    assert (
        consensus_sites(
            sequence=sequence,
            pattern=pattern,
            max_mismatches=limit,
            both_strands=both,
        )
        == expected
    )


def test_iupac_ambiguity_overlaps_and_strand_merging():
    assert [
        r["start"]
        for r in consensus_sites(
            sequence="aaaa", pattern="AA", both_strands=False
        )
    ] == [0, 1, 2]
    assert consensus_sites(sequence="ANG", pattern="ANG") == []
    assert consensus_sites(sequence="AGG", pattern="ARG")[0]["mismatches"] == 0
    assert (
        consensus_sites(sequence="CACGTG", pattern="CACGTG")[0]["site_strand"]
        == "."
    )
    assert (
        consensus_sites(sequence="CGT", pattern="ACG")[0]["site_strand"] == "-"
    )
    assert consensus_sites(sequence="AC", pattern="ACG") == []


@pytest.mark.parametrize(
    "options",
    [
        {"max_mismatches": -1},
        {"max_mismatches": 3},
        {"max_mismatches": True},
        {"pattern": "X"},
        {"sequence": "ACX"},
    ],
)
def test_invalid_consensus_screen(options):
    with pytest.raises(ValueError):
        consensus_sites(**{**dict(sequence="ACG", pattern="ACG"), **options})


def test_scan_target_validation():
    for options in (
        dict(motif_id="", pattern="ACG"),
        dict(motif_id="x", pattern="acg"),
        dict(motif_id="x", pattern="ACG", source_q_value=np.nan),
        dict(motif_id="x", pattern="ACG", source_q_value=1.1),
    ):
        with pytest.raises(ValueError):
            ScanMotif(**options)


def test_enrichment_selection_and_explicit_pwm_consensus(tmp_path):
    path = tmp_path / "enrichment.tsv"
    rows = [
        dict(
            motif_id="b",
            consensus="ACG",
            kind="kmer",
            q_value=0.02,
            positive_fraction=0.8,
            negative_fraction=0.1,
        ),
        dict(
            motif_id="a",
            consensus="CGT",
            kind="pwm",
            q_value=0.01,
            positive_fraction=0.9,
            negative_fraction=0.2,
        ),
        dict(
            motif_id="depleted",
            consensus="TTT",
            kind="kmer",
            q_value=0.001,
            positive_fraction=0.1,
            negative_fraction=0.8,
        ),
        dict(
            motif_id="weak",
            consensus="AAA",
            kind="kmer",
            q_value=0.1,
            positive_fraction=0.5,
            negative_fraction=0.1,
        ),
    ]
    write_tsv(path=path, rows=rows, fields=tuple(rows[0]))
    targets = read_scan_motifs(enrichment_path=path, max_motifs=1)
    assert [m.motif_id for m in targets] == ["a"]
    assert targets[0].source_kind == "pwm_consensus"
    assert targets[0].source_q_value == 0.01
    assert [
        m.motif_id
        for m in read_scan_motifs(enrichment_path=path, motif_ids=["b"])
    ] == ["b"]
    assert read_scan_motifs(enrichment_path=path, q_threshold=0) == []
    for options in (
        {"motif_ids": ["absent"]},
        {"motif_ids": ["a", "a"]},
        {"max_motifs": True},
        {"q_threshold": np.nan},
        {"q_threshold": -1},
    ):
        with pytest.raises(ValueError):
            read_scan_motifs(enrichment_path=path, **options)
    targets, options = scan_targets(
        enrichment_path=path, settings={"max_motifs": 1, "max_mismatches": 1}
    )
    assert len(targets) == 1 and options == {"max_mismatches": 1}
    with pytest.raises(ValueError):
        read_scan_motifs()
    with pytest.raises(ValueError):
        read_scan_motifs(enrichment_path=path, motif_path=path)
    path.write_text("wrong\tcolumns\n")
    with pytest.raises(ValueError, match="columns"):
        read_scan_motifs(enrichment_path=path)
    for change in (
        {"q_value": "bad"},
        {"negative_fraction": np.inf},
        {"positive_fraction": None},
    ):
        write_tsv(
            path=path, rows=[{**rows[0], **change}], fields=tuple(rows[0])
        )
        with pytest.raises(ValueError):
            read_scan_motifs(enrichment_path=path)
    write_tsv(path=path, rows=[rows[0], rows[0]], fields=tuple(rows[0]))
    with pytest.raises(ValueError, match="Duplicate"):
        read_scan_motifs(enrichment_path=path)
    jaspar = tmp_path / "motif.jaspar"
    jaspar.write_text(">J demo\nA [9 1]\nC [1 9]\nG [0 0]\nT [0 0]\n")
    target = read_scan_motifs(motif_path=jaspar)[0]
    assert target.pattern == "AC" and target.source_kind == "pwm_consensus"


@pytest.fixture
def scan_data(tmp_path):
    bases = list("T" * 120)
    for start, word in (
        (5, "ACG"),
        (22, "ACG"),
        (32, "CGT"),
        (44, "ACG"),
        (61, "ATG"),
        (90, "ACG"),
        (100, "acg"),
        (110, "ANG"),
    ):
        bases[start : start + len(word)] = word
    genome = tmp_path / "genome.fa"
    annotation = tmp_path / "genes.gff3"
    write_fasta(
        path=genome, records={"chr": "".join(bases), "no_genes": "ACG"}
    )
    annotation.write_text(
        "##gff-version 3\n"
        "chr\ttest\tgene\t71\t80\t.\t-\t.\tID=n\n"
        "chr\ttest\tgene\t21\t26\t.\t+\t.\tID=p\n"
        "chr\ttest\tgene\t41\t50\t.\t.\t.\tID=unknown\n"
    )
    positive = tmp_path / "positive.txt"
    negative = tmp_path / "negative.txt"
    positive.write_text("p\n")
    negative.write_text("n\n")
    return dict(
        genome_path=genome,
        annotation_path=annotation,
        positive_genes=positive,
        negative_genes=negative,
    )


def test_chunk_boundaries_match_full_scan(scan_data, tmp_path):
    from intergenic_regions.genome import Genome

    targets = [
        ScanMotif(motif_id="x", pattern="ACG"),
        ScanMotif(motif_id="y", pattern="TACGT"),
    ]
    for chunk in (1, 7, 31, 500):
        directory = tmp_path / f"chunk{chunk}"
        result = genome_scan_outputs(
            directory=directory,
            motifs=targets,
            max_mismatches=1,
            chunk_size=chunk,
            genome_path=scan_data["genome_path"],
        )
        actual = records(directory / "genome_motif_sites.tsv")
        expected = []
        with Genome(path=scan_data["genome_path"]) as genome:
            for contig, length in genome.lengths.items():
                for motif in targets:
                    for row in consensus_sites(
                        sequence=genome.fetch(
                            contig=contig, start=0, end=length
                        ),
                        pattern=motif.pattern,
                        max_mismatches=1,
                    ):
                        expected.append(
                            (
                                contig,
                                motif.motif_id,
                                row["start"],
                                row["site_strand"],
                                row["mismatches"],
                            )
                        )
        observed = [
            (
                r["contig"],
                r["motif_id"],
                int(r["start"]),
                r["site_strand"],
                int(r["mismatches"]),
            )
            for r in actual
        ]
        assert sorted(observed) == sorted(expected)
        assert result["total_sites"] == len(expected)
        assert all(r["nearest_gene_id"] == "" for r in actual)
        assert "annotation_not_supplied" in {r["context"] for r in actual}


def test_gene_orientation_intergenic_mask_and_density(scan_data, tmp_path):
    output = tmp_path / "scan"
    result = genome_scan_outputs(
        directory=output,
        motifs=[ScanMotif(motif_id="x", pattern="ACG")],
        max_mismatches=1,
        intergenic_only=True,
        mask_lowercase=True,
        chunk_size=9,
        upstream=50,
        downstream=50,
        bin_width=10,
        **scan_data,
    )
    sites = records(output / "genome_motif_sites.tsv")
    by_start = {int(r["start"]): r for r in sites if r["contig"] == "chr"}
    assert 22 not in by_start and 44 not in by_start and 100 not in by_start
    assert 110 not in by_start
    assert by_start[5]["distance_to_gene_start"] == "-14.0"
    assert by_start[32]["site_strand"] == "-"
    assert by_start[32]["distance_to_gene_start"] == "13.0"
    assert by_start[61]["distance_to_gene_start"] == "17.0"
    assert by_start[61]["mismatches"] == "1"
    assert by_start[90]["distance_to_gene_start"] == "-12.0"
    assert by_start[90]["gene_start_proxy"] == "79"
    assert by_start[90]["gene_cohort"] == "negative"
    assert result["genic_overlap_sites"] == 0
    assert all(
        r["enhancer_status"] == "unvalidated_sequence_match" for r in sites
    )
    profiles = records(output / "distance_profiles.tsv")
    assert sum(int(r["site_count"]) for r in profiles) == len(by_start)
    for row in profiles:
        count, denominator = (
            int(row["site_count"]),
            int(row["eligible_windows"]),
        )
        assert count <= denominator
        if denominator:
            assert float(row["hits_per_million_windows"]) == pytest.approx(
                1e6 * count / denominator
            )
    assert "TSS proxy" in result["distance_anchor"]
    assert len(list((output / "figures").glob("*.png"))) == 5


def test_nearest_ties_are_audited_and_excluded(tmp_path):
    genome = tmp_path / "genome.fa"
    bases = list("G" * 100)
    bases[49:51] = "AT"
    write_fasta(path=genome, records={"chr": "".join(bases)})
    annotation = tmp_path / "genes.gff3"
    annotation.write_text(
        "chr\tt\tgene\t21\t30\t.\t+\t.\tID=a\n"
        "chr\tt\tgene\t71\t80\t.\t-\t.\tID=b\n"
    )
    result = genome_scan_outputs(
        directory=tmp_path / "tie",
        genome_path=genome,
        annotation_path=annotation,
        motifs=[ScanMotif(motif_id="x", pattern="AT")],
        upstream=100,
        downstream=100,
        bin_width=25,
    )
    assert result["total_sites"] == result["tied_nearest_gene_sites"] == 1
    assert (
        records(tmp_path / "tie/genome_motif_sites.tsv")[0][
            "nearest_gene_tied"
        ]
        == "True"
    )
    assert (
        sum(
            int(r["site_count"])
            for r in records(tmp_path / "tie/distance_profiles.tsv")
        )
        == 0
    )
    annotation.write_text(
        "chr\tt\tgene\t21\t30\t.\t+\t.\tID=a\n"
        "chr\tt\tgene\t21\t40\t.\t+\t.\tID=b\n"
    )
    duplicate = genome_scan_outputs(
        directory=tmp_path / "duplicate",
        genome_path=genome,
        annotation_path=annotation,
        motifs=[ScanMotif(motif_id="x", pattern="AT")],
    )
    assert duplicate["tied_nearest_gene_sites"] == 1
    assert all(
        int(r["eligible_windows"]) == 0
        for r in records(tmp_path / "duplicate/distance_profiles.tsv")
    )


@pytest.mark.parametrize(
    "changes",
    [
        {"chunk_size": 0},
        {"max_hits": True},
        {"max_mismatches": 3},
        {"upstream": -1},
        {"upstream": 0, "downstream": 0},
        {"bin_width": 0},
        {"upstream": 5000, "bin_width": 1},
        {"max_mismatches": True},
        {"annotation_path": None, "intergenic_only": True},
        {"motifs": [ScanMotif(motif_id="x", pattern="ACG")] * 2},
    ],
)
def test_defensive_scan_settings(scan_data, tmp_path, changes):
    options = {
        **scan_data,
        "directory": tmp_path / "bad",
        "motifs": [ScanMotif(motif_id="x", pattern="ACG")],
    }
    with pytest.raises(ValueError):
        genome_scan_outputs(**{**options, **changes})


def test_bad_cohorts_and_atomic_hit_cap(scan_data, tmp_path):
    scan_data["positive_genes"].write_text("absent\n")
    with pytest.raises(ValueError, match="unknown"):
        genome_scan_outputs(
            directory=tmp_path / "unknown",
            **scan_data,
            motifs=[ScanMotif(motif_id="x", pattern="ACG")],
        )
    scan_data["positive_genes"].write_text("n\n")
    with pytest.raises(ValueError, match="overlap"):
        genome_scan_outputs(
            directory=tmp_path / "overlap", **scan_data, motifs=[]
        )
    motifs = tmp_path / "motifs.tsv"
    motifs.write_text("motif_id\tpattern\nx\tACG\n")
    output = tmp_path / "capped"
    with pytest.raises(ValueError, match="max_hits"):
        genome_scan_workflow(
            genome_path=scan_data["genome_path"],
            motif_path=motifs,
            output=output,
            settings={"max_hits": 1, "chunk_size": 7},
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".capped.*"))


def test_no_selected_motifs_and_empty_plot_paths(scan_data, tmp_path):
    summary = genome_scan_outputs(
        directory=tmp_path / "empty",
        genome_path=scan_data["genome_path"],
        motifs=[],
    )
    assert summary["status"] == "no_selected_motifs"
    assert summary["total_sites"] == 0
    assert records(tmp_path / "empty/genome_motif_sites.tsv") == []
    assert (
        plot_scan(profiles=[], burden=[], genes=[], directory=tmp_path) == []
    )
    scan_report(
        directory=tmp_path,
        summary=summary,
        sites=[],
        profiles=[],
        burden=[],
        genes=[],
    )
    assert "No records" in (tmp_path / "report.html").read_text()


def test_cli_scan_and_pipeline_integration(dataset, tmp_path):
    parser = build_parser()
    args = parser.parse_args(
        [
            "scan-genome",
            "--genome",
            str(dataset["genome"]),
            "--annotation",
            str(dataset["annotation"]),
            "--motifs",
            str(dataset["motifs"]),
            "--scan-max-motifs",
            "1",
            "--max-mismatches",
            "1",
            "--scan-intergenic-only",
            "--output-dir",
            str(tmp_path / "cli"),
        ]
    )
    assert scan_options(args=args)["max_mismatches"] == 1
    summary = dispatch(args=args)
    assert summary["total_sites"] > 0
    assert summary["genic_overlap_sites"] == 0
    assert json.loads((tmp_path / "cli/manifest.json").read_text())["inputs"]
    summary = pipeline_workflow(
        genome_path=dataset["genome"],
        annotation_path=dataset["annotation"],
        positive_genes=dataset["positive_genes"],
        negative_genes=dataset["negative_genes"],
        motif_path=dataset["motifs"],
        output=tmp_path / "pipeline",
        use_ai=False,
        scan_genome=True,
        scan_settings={
            "max_motifs": 1,
            "max_mismatches": 1,
            "intergenic_only": True,
        },
    )
    assert summary["scan"]["selected_motifs"] == 1
    assert summary["scan"]["total_sites"] > 0
    report = (tmp_path / "pipeline/report.html").read_text()
    assert "Whole-genome scan report" in report
    assert "genome_distance_density" in report
    assert "data:application/pdf;base64," in report
    assert "Download PNG" in report
