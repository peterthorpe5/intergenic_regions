"""Logarithmic flank queries over the union of all annotated gene spans."""

import logging
from bisect import bisect_left
from collections import defaultdict
from collections.abc import Sequence

from intergenic_regions.genome import (
    Genome,
    reverse_complement,
    sequence_composition,
)
from intergenic_regions.models import Gene, Region

LOGGER = logging.getLogger(__name__)


def merge_intervals(
    *, intervals: Sequence[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Merge overlapping and touching intervals, independent of input order.

    Args:
        intervals: Non-empty, zero-based half-open spans.

    Returns:
        Sorted disjoint spans.

    Raises:
        ValueError: An input span is invalid.
    """
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError("Invalid interval")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


class GeneIndex:
    """Sorted merged genic spans; both strands always block extraction."""

    def __init__(
        self, *, genes: Sequence[Gene], lengths: dict[str, int]
    ) -> None:
        """Validate genes against genome lengths and build interval indices.

        Args:
            genes: All genes, including genes outside a target list.
            lengths: Genome contig lengths.

        Raises:
            ValueError: IDs repeat or a gene falls outside the genome.
        """
        self.genes: dict[str, Gene] = {}
        self.lengths = dict(lengths)
        by_contig: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for gene in genes:
            if gene.gene_id in self.genes:
                raise ValueError(f"Duplicate gene ID: {gene.gene_id}")
            if gene.contig not in lengths or gene.end > lengths[gene.contig]:
                raise ValueError(
                    f"Gene outside genome: {gene.gene_id} ({gene.contig})"
                )
            self.genes[gene.gene_id] = gene
            by_contig[gene.contig].append((gene.start, gene.end))
        self.spans = {
            contig: merge_intervals(intervals=spans)
            for contig, spans in by_contig.items()
        }
        self.starts = {
            contig: [s for s, _ in spans]
            for contig, spans in self.spans.items()
        }

    def free_flank(
        self, *, gene: Gene, direction: str
    ) -> tuple[int, int, str]:
        """Return the contiguous intergenic interval touching a gene boundary.

        Args:
            gene: An indexed target gene.
            direction: ``upstream`` or ``downstream``.

        Returns:
            Start, end and limiting boundary type. An occupied boundary
            produces a zero-length interval; blockers are never jumped.

        Raises:
            ValueError: Direction or target strand is invalid.
        """
        if self.genes.get(gene.gene_id) != gene:
            raise ValueError("Target gene is not in this index")
        if direction not in {"upstream", "downstream"} or gene.strand == ".":
            raise ValueError(
                "A known strand and valid flank direction are required"
            )
        left = (gene.strand == "+") == (direction == "upstream")
        anchor = gene.start if left else gene.end
        spans = self.spans[gene.contig]
        position = bisect_left(a=self.starts[gene.contig], x=anchor)
        previous_end = spans[position - 1][1] if position else 0
        if previous_end > anchor:
            return anchor, anchor, "overlapping_gene"
        if left:
            reason = "gene_boundary" if position else "contig_boundary"
            return previous_end, anchor, reason
        next_start = (
            spans[position][0]
            if position < len(spans)
            else self.lengths[gene.contig]
        )
        reason = (
            "gene_boundary" if position < len(spans) else "contig_boundary"
        )
        return anchor, next_start, reason


def extract_regions(
    *,
    index: GeneIndex,
    genome: Genome,
    identifiers: Sequence[str] | None = None,
    direction: str = "upstream",
    length: int | None = 1000,
    offset: int = 0,
    min_length: int = 1,
    max_ambiguous_fraction: float = 1.0,
    mask_lowercase: bool = False,
) -> list[Region]:
    """Extract selected flanks without including any annotated genic bases.

    Args:
        index: Index of the full annotation.
        genome: Open indexed genome.
        identifiers: Exact gene IDs, or all genes when omitted.
        direction: ``upstream``, ``downstream`` or ``both``.
        length: Maximum bases to return; ``None`` returns the full gap.
        offset: Bases to skip nearest the gene, within the same free flank.
        min_length: Minimum retained length, inclusive.
        max_ambiguous_fraction: Maximum fraction of non-ACGT bases.
        mask_lowercase: Replace soft-masked bases with N.

    Returns:
        Regions including audited exclusions, in requested gene order.

    Raises:
        ValueError: Parameters or target identifiers are invalid.
    """
    if direction not in {"upstream", "downstream", "both"}:
        raise ValueError("Invalid extraction direction")
    if (length is not None and length < 1) or offset < 0 or min_length < 1:
        raise ValueError(
            "Length/minimum must be positive; offset non-negative"
        )
    if not 0 <= max_ambiguous_fraction <= 1:
        raise ValueError("Ambiguous fraction must lie between zero and one")
    targets = (
        list(identifiers) if identifiers is not None else list(index.genes)
    )
    if len(set(targets)) != len(targets):
        raise ValueError("Target identifiers contain duplicates")
    missing = set(targets) - index.genes.keys()
    if missing:
        raise ValueError(
            f"Unknown gene identifiers: {' '.join(sorted(missing)[:20])}"
        )
    directions = (
        ("upstream", "downstream") if direction == "both" else (direction,)
    )
    regions: list[Region] = []
    for identifier in targets:
        gene = index.genes[identifier]
        for side in directions:
            if gene.strand == ".":
                regions.append(
                    Region(
                        gene_id=identifier,
                        contig=gene.contig,
                        start=gene.start,
                        end=gene.start,
                        strand=".",
                        direction=side,
                        sequence="",
                        status="unknown_strand",
                        stop_reason="unknown_strand",
                        available_length=0,
                    )
                )
                continue
            start, end, reason = index.free_flank(gene=gene, direction=side)
            available = end - start
            skip = min(offset, available)
            requested = available if length is None else length
            keep = min(requested, available - skip)
            left = (gene.strand == "+") == (side == "upstream")
            if left:
                end -= skip
                start = end - keep
            else:
                start += skip
                end = start + keep
            sequence = genome.fetch(contig=gene.contig, start=start, end=end)
            if mask_lowercase:
                sequence = "".join("N" if c.islower() else c for c in sequence)
            if gene.strand == "-":
                sequence = reverse_complement(sequence=sequence)
            status = "retained"
            if not keep:
                status = (
                    "no_intergenic_space"
                    if not available
                    else "offset_exceeds_gap"
                )
            elif keep < min_length:
                status = "below_minimum_length"
            elif (
                sequence_composition(sequence=sequence)["ambiguous_fraction"]
                > max_ambiguous_fraction
            ):
                status = "excess_ambiguity"
            regions.append(
                Region(
                    gene_id=identifier,
                    contig=gene.contig,
                    start=start,
                    end=end,
                    strand=gene.strand,
                    direction=side,
                    sequence=sequence.upper() if status == "retained" else "",
                    status=status,
                    stop_reason=reason,
                    available_length=available,
                )
            )
    LOGGER.info(
        "Retained %d of %d requested flanks",
        sum(r.status == "retained" for r in regions),
        len(regions),
    )
    return regions
