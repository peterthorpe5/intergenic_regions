"""Generate explicitly synthetic RNA-start and curated-gene fixtures."""

import argparse
import logging
import random
from pathlib import Path
from typing import Any

from intergenic_regions.genome import reverse_complement
from intergenic_regions.io import (
    output_bundle,
    write_fasta,
    write_json,
    write_tsv,
)

LOGGER = logging.getLogger(__name__)


def make_rna_demo(
    *,
    directory: Path,
    count: int = 12,
    replicates: int = 3,
    flank_length: int = 1200,
    seed: int = 17,
    write_reads: bool = False,
) -> dict[str, Any]:
    """Create strand-balanced synthetic coverage, counts and curated gene sets.

    Args:
        directory: New destination; existing paths are refused.
        count: Genes per curated class, including both transcript strands.
        replicates: Independent simulated observations per condition.
        flank_length: Safe upstream gap before the annotated gene boundary.
        seed: Reproducible sequence generator seed.
        write_reads: Also export small splice-aware FASTQ mapping fixtures.

    Returns:
        Fixture metadata, with planted starts and motif/control definitions.

    Raises:
        ValueError: Counts, flank length or seed are invalid.
    """
    if (
        isinstance(count, bool)
        or not 2 <= count <= 100
        or isinstance(replicates, bool)
        or not 1 <= replicates <= 20
        or isinstance(flank_length, bool)
        or not 600 <= flank_length <= 10000
        or isinstance(seed, bool)
        or not isinstance(seed, int)
        or seed < 0
    ):
        raise ValueError("Invalid synthetic RNA demo settings")
    rng = random.Random(seed)
    genome: dict[str, str] = {}
    annotation: list[str] = ["##gff-version 3\n"]
    targets: list[dict[str, Any]] = []
    genes: dict[str, str] = {}
    truth: list[dict[str, Any]] = []
    positive: list[str] = []
    negative: list[str] = []
    size = flank_length + 780
    for label, identifiers in (("positive", positive), ("negative", negative)):
        for number in range(count):
            gene_id = f"{label}_{number + 1:03d}"
            contig = f"locus_{gene_id}"
            strand = "+" if number % 2 == 0 else "-"
            start = flank_length + 80 if strand == "+" else 200
            end = start + 500
            blocker_start, blocker_end = (
                (0, 80)
                if strand == "+"
                else (end + flank_length, end + flank_length + 80)
            )
            symbols = list(rng.choices(population="ACGT", k=size))
            if strand == "+":
                symbols[start + 100 : start + 102] = "GT"
                symbols[end - 102 : end - 100] = "AG"
            else:
                symbols[start + 100 : start + 102] = "CT"
                symbols[end - 102 : end - 100] = "AC"
            motif = "CACGTG" if label == "positive" else "GCACTG"
            for distance in (190, 220, 250, 280, 310, 340):
                if strand == "+":
                    first = start - distance
                    symbols[first : first + 6] = motif
                else:
                    first = end + distance - 6
                    symbols[first : first + 6] = reverse_complement(
                        sequence=motif
                    )
            genome[contig] = "".join(symbols)
            exons = [(start, start + 100), (end - 100, end)]
            for exon_number, (left, right) in enumerate(exons, start=1):
                annotation.append(
                    f"{contig}\tdemo\texon\t{left + 1}\t{right}\t.\t"
                    f"{strand}\t.\t"
                    f"ID={gene_id}.exon{exon_number};Parent={gene_id}.t1\n"
                )
            annotation += [
                f"{contig}\tdemo\tmRNA\t{start + 1}\t{end}\t.\t{strand}\t.\t"
                f"ID={gene_id}.t1;Parent={gene_id}\n",
                f"{contig}\tdemo\tgene\t{start + 1}\t{end}\t.\t{strand}\t.\t"
                f"ID={gene_id}\n",
                f"{contig}\tdemo\tgene\t{blocker_start + 1}\t{blocker_end}"
                "\t.\t"
                f"+\t.\tID=blocker_{gene_id}\n",
                f"{contig}\tdemo\texon\t{blocker_start + 1}\t{blocker_end}"
                "\t.\t"
                f"+\t.\tParent=blocker_{gene_id}\n",
            ]
            identifiers.append(gene_id)
            genes[gene_id] = label
            genes[f"blocker_{gene_id}"] = "background"
            targets.append(
                {
                    "gene_id": gene_id,
                    "contig": contig,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "label": label,
                }
            )
    with output_bundle(path=directory) as stage:
        write_fasta(path=stage / "genome.fasta", records=genome)
        (stage / "genes.gff3").write_text(
            data="".join(annotation), encoding="utf-8"
        )
        (stage / "positive_genes.txt").write_text(
            data="".join(f"{gene}\n" for gene in positive), encoding="utf-8"
        )
        (stage / "negative_genes.txt").write_text(
            data="".join(f"{gene}\n" for gene in negative), encoding="utf-8"
        )
        (stage / "motifs.tsv").write_text(
            data="motif_id\tpattern\nsynthetic_G_box\tCACGTG\n",
            encoding="utf-8",
        )
        manifests: list[dict[str, Any]] = []
        raw_manifests: list[dict[str, Any]] = []
        count_rows: dict[str, dict[str, Any]] = {
            gene: {"gene_id": gene} for gene in genes
        }
        for condition in ("reference", "stimulated"):
            for replicate in range(replicates):
                sample_id = f"{condition}_{replicate + 1}"
                jitter = (replicate % 3 - 1) * 2
                tracks: dict[str, list[tuple[str, int, int, float]]] = {
                    "+": [],
                    "-": [],
                }
                reads: list[str] = []
                for target in targets:
                    extension = (
                        150 if condition == "stimulated" else 110
                    ) + jitter
                    gene_id, start, end = (
                        target["gene_id"],
                        target["start"],
                        target["end"],
                    )
                    strand, contig = target["strand"], target["contig"]
                    anchor = (
                        start - extension
                        if strand == "+"
                        else end - 1 + extension
                    )
                    depth = (
                        (45 if condition == "stimulated" else 15)
                        if target["label"] == "positive"
                        else 20
                    )
                    depth += replicate
                    spans = (
                        [(anchor, start + 100), (end - 100, end)]
                        if strand == "+"
                        else [(start, start + 100), (end - 100, anchor + 1)]
                    )
                    tracks[strand] += [
                        (contig, left, right, float(depth))
                        for left, right in spans
                    ]
                    truth.append(
                        {
                            "gene_id": gene_id,
                            "sample_id": sample_id,
                            "condition": condition,
                            "strand": strand,
                            "planted_start": anchor,
                            "status": "synthetic_fixture_not_experimental_TSS",
                        }
                    )
                    if write_reads:
                        pieces = [
                            genome[contig][left:right] for left, right in spans
                        ]
                        transcript = (
                            "".join(pieces)
                            if strand == "+"
                            else reverse_complement(sequence="".join(pieces))
                        )
                        stride = (
                            2
                            if condition == "stimulated"
                            and target["label"] == "positive"
                            else 3
                        )
                        keep_probability = rng.uniform(a=0.6, b=0.95)
                        for offset in range(0, len(transcript) - 59, stride):
                            if rng.random() > keep_probability:
                                continue
                            sequence = transcript[offset : offset + 60]
                            reads.append(
                                f"@{sample_id}_{gene_id}_{offset}\n"
                                f"{sequence}\n+\n{'I' * 60}\n"
                            )
                for strand, name in (("+", "plus"), ("-", "minus")):
                    path = stage / f"{sample_id}.{name}.bedgraph"
                    path.write_text(
                        data="".join(
                            f"{contig}\t{left}\t{right}\t{value}\n"
                            for contig, left, right, value in sorted(
                                tracks[strand]
                            )
                        ),
                        encoding="utf-8",
                    )
                manifests.append(
                    {
                        "sample_id": sample_id,
                        "condition": condition,
                        "replicate_id": f"bio{replicate + 1}",
                        "strandedness": "forward",
                        "layout": "coverage",
                        "coverage_plus": f"{sample_id}.plus.bedgraph",
                        "coverage_minus": f"{sample_id}.minus.bedgraph",
                        "size_factor": 1,
                    }
                )
                if write_reads:
                    (stage / f"{sample_id}.fastq").write_text(
                        data="".join(reads), encoding="utf-8"
                    )
                    raw_manifests.append(
                        {
                            "sample_id": sample_id,
                            "condition": condition,
                            "replicate_id": f"bio{replicate + 1}",
                            "strandedness": "forward",
                            "layout": "single",
                            "read1": f"{sample_id}.fastq",
                        }
                    )
                for gene, label in genes.items():
                    baseline = 100 + sum(ord(char) for char in gene) % 80
                    dispersion = rng.gammavariate(alpha=60, beta=1 / 60)
                    scale = (
                        3.5
                        if label == "positive" and condition == "stimulated"
                        else 1
                    )
                    count_rows[gene][sample_id] = max(
                        1, int(baseline * dispersion * scale)
                    )
        write_tsv(
            path=stage / "samples.tsv",
            rows=manifests,
            fields=tuple(manifests[0]),
        )
        if raw_manifests:
            write_tsv(
                path=stage / "raw_samples.tsv",
                rows=raw_manifests,
                fields=tuple(raw_manifests[0]),
            )
        write_tsv(
            path=stage / "gene_counts.tsv",
            rows=list(count_rows.values()),
            fields=tuple(next(iter(count_rows.values()))),
        )
        write_tsv(
            path=stage / "planted_starts.tsv",
            rows=truth,
            fields=tuple(truth[0]),
        )
        summary = {
            "synthetic": True,
            "validation_scope": "software_fixture_not_biological_performance",
            "genes_per_curated_class": count,
            "biological_replicates_per_condition": replicates,
            "flank_length": flank_length,
            "seed": seed,
            "raw_reads_exported": write_reads,
            "positive_motif": "CACGTG",
            "control_pattern": "GCACTG",
            "coverage_tracks": (
                "idealised_piecewise_depth_not_recomputed_from_exported_FASTQ"
            ),
            "RNA_extension_reference": 110,
            "RNA_extension_stimulated": 150,
        }
        write_json(path=stage / "fixture.json", data=summary)
    LOGGER.info("Created explicit synthetic RNA fixture: %s", directory)
    return summary


def main(*, argv: list[str] | None = None) -> int:
    """Generate a synthetic fixture from named command-line options.

    Args:
        argv: Optional named-option arguments for tests.

    Returns:
        Zero on completion; two for invalid inputs or an existing destination.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--flank-length", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--write-reads", action="store_true")
    args = parser.parse_args(args=argv)
    logging.basicConfig(level=logging.INFO)
    try:
        make_rna_demo(
            directory=args.output_dir,
            count=args.count,
            replicates=args.replicates,
            flank_length=args.flank_length,
            seed=args.seed,
            write_reads=args.write_reads,
        )
    except (ValueError, OSError) as exc:
        LOGGER.error("%s", exc)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
