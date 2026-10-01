"""Generate a reproducible synthetic genome and regulatory analysis inputs."""

import argparse
import logging
import random
from pathlib import Path
from typing import Any
from urllib.parse import quote

from intergenic_regions.genome import reverse_complement
from intergenic_regions.io import (
    output_bundle,
    write_fasta,
    write_json,
    write_tsv,
)

LOGGER = logging.getLogger(__name__)


def make_demo(
    *, output: Path, sequences_per_class: int = 20, seed: int = 17
) -> dict[str, Any]:
    """Write a planted dataset with both strands and explicit gene blockers.

    Args:
        output: New destination directory.
        sequences_per_class: Positive and negative observation count.
        seed: Deterministic random seed.

    Returns:
        Synthetic dataset metadata.

    Raises:
        ValueError: Fewer than five sequences per class are requested.
        FileExistsError: The destination already exists.
    """
    if sequences_per_class < 5:
        raise ValueError("Demo needs at least five sequences per class")
    rng = random.Random(seed)
    genome: dict[str, str] = {}
    records: dict[str, dict[str, str]] = {"positive": {}, "negative": {}}
    annotations: list[str] = []
    groups: list[dict[str, str]] = []
    peaks: list[str] = []
    assays: list[dict[str, str]] = []
    for label in records:
        for index in range(sequences_per_class):
            identifier = f"{label}_{index:03d}"
            contig = f"contig_{identifier}"
            promoter = "".join(rng.choices("ACGT", k=240))
            if label == "positive":
                promoter = promoter[:80] + "CACGTGCACGTG" + promoter[92:]
            strand = "+" if index % 2 == 0 else "-"
            start, end = (280, 310) if strand == "+" else (40, 70)
            genome[contig] = (
                "A" * 40 + promoter + "G" * 30 + "T" * 20
                if strand == "+"
                else "A" * 40
                + "G" * 30
                + reverse_complement(sequence=promoter)
                + "T" * 20
            )
            for gene_id, left, right, gene_strand in (
                (identifier, start, end, strand),
                (f"{identifier}_left_blocker", 0, 40, "-"),
                (f"{identifier}_right_blocker", 310, 330, "+"),
            ):
                annotations.append(
                    f"{contig}\tdemo\tgene\t{left + 1}\t{right}\t.\t"
                    f"{gene_strand}\t.\tID={quote(gene_id, safe='._-')}\n"
                )
            records[label][identifier] = promoter
            groups.append(
                {
                    "sequence_id": f"{identifier}|upstream",
                    "group": contig,
                }
            )
            if label == "positive" and index < sequences_per_class // 2:
                left = 60 if strand == "+" else 120
                peaks.append(f"{contig}\t{left}\t{left + 60}\tdemo_peak\n")
            if label == "positive" and index == 0:
                assays.append(
                    {
                        "gene_id": identifier,
                        "evidence_type": "synthetic_assay",
                        "value": "positive",
                        "source": "planted demonstration",
                    }
                )
    summary = {
        "synthetic": True,
        "seed": seed,
        "sequences_per_class": sequences_per_class,
        "planted_pattern": "CACGTG",
        "promoter_length": 240,
        "interpretation": (
            "Synthetic software demonstration; no biological validation"
        ),
    }
    with output_bundle(path=output) as stage:
        write_fasta(path=stage / "genome.fasta", records=genome)
        (stage / "genes.gff3").write_text(
            "##gff-version 3\n" + "".join(reversed(annotations)),
            encoding="utf-8",
        )
        for label, collection in records.items():
            (stage / f"{label}_genes.txt").write_text(
                "\n".join(collection) + "\n",
                encoding="utf-8",
            )
            write_fasta(path=stage / f"{label}.fasta", records=collection)
        write_tsv(
            path=stage / "groups.tsv",
            rows=groups,
            fields=("sequence_id", "group"),
        )
        write_tsv(
            path=stage / "motifs.tsv",
            rows=[
                {
                    "motif_id": "Gbox",
                    "pattern": "CACGTG",
                    "name": "G box",
                }
            ],
            fields=("motif_id", "pattern", "name"),
        )
        write_tsv(
            path=stage / "functional.tsv",
            rows=assays,
            fields=("gene_id", "evidence_type", "value", "source"),
        )
        (stage / "accessibility.bed").write_text(
            "".join(peaks),
            encoding="utf-8",
        )
        write_json(path=stage / "dataset.json", data=summary)
    LOGGER.info("Synthetic inputs written to %s", output)
    return summary


def main(*, argv: list[str] | None = None) -> int:
    """Run the named-option demo generator.

    Args:
        argv: Optional command arguments.

    Returns:
        Zero on success or two on invalid input/output.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sequences-per-class", type=int, default=20)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args(args=argv)
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s"
    )
    try:
        make_demo(
            output=args.output_dir,
            sequences_per_class=args.sequences_per_class,
            seed=args.seed,
        )
    except (ValueError, OSError) as exc:
        LOGGER.error("%s", exc)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
