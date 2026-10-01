"""Small deterministic datasets for biological and complete-workflow tests."""

import random
from pathlib import Path

import pytest

from intergenic_regions.io import write_fasta


@pytest.fixture
def dataset(*, tmp_path: Path) -> dict[str, Path]:
    """Create a planted, composition-matched regulatory sequence dataset."""
    rng = random.Random(17)
    sequences: dict[str, str] = {}
    annotation: list[str] = []
    positive: dict[str, str] = {}
    negative: dict[str, str] = {}
    for label, records in (("p", positive), ("n", negative)):
        for index in range(15):
            identifier = f"{label}{index}"
            promoter = "".join(rng.choices("ACGT", k=120))
            if label == "p":
                promoter = promoter[:45] + "CACGTGCACGTG" + promoter[57:]
            sequences[identifier] = promoter + "A" * 20 + "T" * 10
            records[identifier] = promoter
            annotation.append(
                f"{identifier}\ttest\tgene\t121\t140\t.\t+\t.\tID={identifier}\n"
            )
    paths = {
        key: tmp_path / filename
        for key, filename in (
            ("genome", "genome.fasta"),
            ("annotation", "genes.gff3"),
            ("positive", "positive.fasta"),
            ("negative", "negative.fasta"),
            ("positive_genes", "positive.txt"),
            ("negative_genes", "negative.txt"),
            ("motifs", "motifs.tsv"),
            ("peaks", "peaks.bed"),
            ("evidence", "functional.tsv"),
        )
    }
    write_fasta(path=paths["genome"], records=sequences)
    write_fasta(path=paths["positive"], records=positive)
    write_fasta(path=paths["negative"], records=negative)
    paths["annotation"].write_text(
        "##gff-version 3\n" + "".join(reversed(annotation))
    )
    paths["positive_genes"].write_text("\n".join(positive) + "\n")
    paths["negative_genes"].write_text("\n".join(negative) + "\n")
    paths["motifs"].write_text(
        "motif_id\tpattern\tname\nGbox\tCACGTG\tG box\nAbsent\tAAAAAAA\tabsent\n"
    )
    paths["peaks"].write_text("p0\t40\t70\tpeak\n")
    paths["evidence"].write_text(
        "gene_id\tevidence_type\tvalue\tsource\np0\treporter_assay\tpositive\tuser experiment\n"
    )
    return paths
