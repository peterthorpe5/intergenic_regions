"""Replicate-aware, descriptive coverage edges with safe intergenic anchors."""

import logging
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import replace
from statistics import median
from typing import Any

from intergenic_regions.extraction import GeneIndex, extract_regions
from intergenic_regions.genome import (
    Genome,
    reverse_complement,
    sequence_composition,
)
from intergenic_regions.models import Gene, Region

LOGGER = logging.getLogger(__name__)


def rna_search_interval(
    *,
    index: GeneIndex,
    gene: Gene,
    exons: Sequence[tuple[int, int]] = (),
    upstream: int = 2000,
    downstream: int = 300,
) -> tuple[int, int, str]:
    """Find a safe 5-prime search interval including the first annotated exon.

    Args:
        index: Full annotation boundary index.
        gene: Exact indexed gene.
        exons: Merged exonic spans for this gene.
        upstream: Maximum intergenic bases before the annotated boundary.
        downstream: Maximum bases into the 5-prime exon or gene-span fallback.

    Returns:
        Genomic start/end and exon-source status. An overlapping 5-prime
        boundary or another gene inside the search stops the interval.

    Raises:
        ValueError: Limits, exon spans or target strand are invalid.
    """
    if (
        isinstance(upstream, bool)
        or not 1 <= upstream <= 100000
        or isinstance(downstream, bool)
        or not 1 <= downstream <= 10000
    ):
        raise ValueError("RNA upstream/downstream limits are invalid")
    free_start, free_end, reason = index.free_flank(
        gene=gene, direction="upstream"
    )
    if reason == "overlapping_gene":
        return gene.start, gene.start, "overlapping_gene"
    if any(
        left < gene.start or right > gene.end or right <= left
        for left, right in exons
    ):
        raise ValueError("RNA exon is outside its target gene")
    source = "first_exonic_union" if exons else "gene_span_fallback"
    if gene.strand == "+":
        first_exon = min(exons) if exons else (gene.start, gene.end)
        start = max(free_start, gene.start - upstream)
        end = min(gene.start + downstream, first_exon[1])
        positions = index._gene_start_values[gene.contig]
        first = bisect_left(a=positions, x=gene.start)
        last = bisect_left(a=positions, x=end)
        for other_start, other_end, identifier in index._gene_starts[
            gene.contig
        ][first:last]:
            if identifier != gene.gene_id and other_end > gene.start:
                end = min(end, max(gene.start, other_start))
    else:
        first_exon = (
            max(exons, key=lambda span: span[1])
            if exons
            else (gene.start, gene.end)
        )
        end = min(free_end, gene.end + upstream)
        start = max(gene.end - downstream, first_exon[0])
        positions = index._gene_end_values[gene.contig]
        first = bisect_right(a=positions, x=start)
        last = bisect_right(a=positions, x=gene.end)
        for other_end, other_start, identifier in index._gene_ends[
            gene.contig
        ][first:last]:
            if identifier != gene.gene_id and other_start < gene.end:
                start = max(start, min(gene.end, other_end))
    return start, max(start, end), source


def detect_coverage_start(
    *,
    depth: Sequence[float],
    genomic_start: int,
    strand: str,
    window: int = 20,
    min_depth: float = 2.0,
    min_fold_rise: float = 4.0,
    min_covered_fraction: float = 0.6,
    pseudocount: float = 0.5,
) -> dict[str, Any]:
    """Locate a sustained coverage rise in transcriptional orientation.

    Args:
        depth: Finite, non-negative genomic-order aligned-base depths.
        genomic_start: Zero-based first genomic base of depth.
        strand: Transcript strand.
        window: Full upstream block and two sustained downstream blocks.
        min_depth: Minimum raw mean depth in both downstream blocks.
        min_fold_rise: Minimum pseudocount-adjusted downstream/upstream ratio.
        min_covered_fraction: Covered-base fraction in both downstream blocks.
        pseudocount: Positive depth pseudocount stabilising zero background.

    Returns:
        Best descriptive edge, a resolution neighbourhood, alternative-edge
        count, effect estimates and explicit status. No significance p-value
        or calibrated TSS probability is assigned to correlated coverage bases.

    Raises:
        ValueError: Depth, coordinates or settings are invalid.
    """
    if (
        strand not in {"+", "-"}
        or genomic_start < 0
        or isinstance(window, bool)
        or not 2 <= window <= 1000
        or not math.isfinite(min_depth)
        or min_depth <= 0
        or not math.isfinite(min_fold_rise)
        or min_fold_rise <= 1
        or not math.isfinite(min_covered_fraction)
        or not 0 < min_covered_fraction <= 1
        or not math.isfinite(pseudocount)
        or pseudocount <= 0
        or any(not math.isfinite(value) or value < 0 for value in depth)
    ):
        raise ValueError("Invalid RNA coverage-start input or settings")
    empty = {
        "status": "not_detected",
        "candidate_start": None,
        "neighbourhood_start": None,
        "neighbourhood_end": None,
        "upstream_mean": None,
        "downstream_mean": None,
        "fold_rise": None,
        "covered_fraction": None,
        "alternative_edges": 0,
        "resolution_bp": window,
    }
    if len(depth) < 3 * window:
        return {**empty, "status": "insufficient_search_space"}
    oriented = list(depth) if strand == "+" else list(reversed(depth))
    sums = [0.0]
    covered = [0]
    for value in oriented:
        sums.append(sums[-1] + value)
        covered.append(covered[-1] + int(value > 0))
    candidates: list[dict[str, Any]] = []
    for edge in range(window, len(depth) - 2 * window + 1):
        before = (sums[edge] - sums[edge - window]) / window
        first = (sums[edge + window] - sums[edge]) / window
        second = (sums[edge + 2 * window] - sums[edge + window]) / window
        after = min(first, second)
        fraction = (
            min(
                covered[edge + window] - covered[edge],
                covered[edge + 2 * window] - covered[edge + window],
            )
            / window
        )
        ratio = (after + pseudocount) / (before + pseudocount)
        if (
            after >= min_depth
            and ratio >= min_fold_rise
            and fraction >= min_covered_fraction
        ):
            candidates.append(
                {
                    "edge": edge,
                    "upstream_mean": before,
                    "downstream_mean": after,
                    "fold_rise": ratio,
                    "covered_fraction": fraction,
                    "edge_score": math.log2(ratio) * fraction,
                }
            )
    if not candidates:
        return empty
    best = max(candidates, key=lambda row: (row["edge_score"], -row["edge"]))
    clusters: list[list[int]] = []
    for row in candidates:
        if not clusters or row["edge"] - clusters[-1][-1] > window:
            clusters.append([])
        clusters[-1].append(row["edge"])
    position = (
        genomic_start + best["edge"]
        if strand == "+"
        else (genomic_start + len(depth) - 1 - best["edge"])
    )
    return {
        **empty,
        **{
            name: best[name]
            for name in (
                "upstream_mean",
                "downstream_mean",
                "fold_rise",
                "covered_fraction",
            )
        },
        "status": "candidate",
        "candidate_start": position,
        "neighbourhood_start": max(genomic_start, position - window),
        "neighbourhood_end": min(
            genomic_start + len(depth), position + window + 1
        ),
        "alternative_edges": len(clusters) - 1,
    }


def consensus_rna_starts(
    *,
    calls: Sequence[Mapping[str, Any]],
    min_replicates: int = 2,
    agreement_tolerance: int = 40,
    min_support_fraction: float = 0.5,
) -> list[dict[str, Any]]:
    """Summarise independently detected edges within each gene/condition.

    Args:
        calls: Per-biological-observation rows, including non-detections.
        min_replicates: Required independent observations for supported status.
        agreement_tolerance: Maximum coordinate span of an agreeing cluster.
        min_support_fraction: Required fraction of all condition observations.

    Returns:
        Consensus positions with replicate counts/ranges and descriptive
        uncertainty intervals. These intervals are not statistical confidence
        intervals. Pooling never promotes singleton/discordant calls.

    Raises:
        ValueError: Settings or biological observation rows are invalid.
    """
    if (
        isinstance(min_replicates, bool)
        or min_replicates < 2
        or isinstance(agreement_tolerance, bool)
        or agreement_tolerance < 0
        or not math.isfinite(min_support_fraction)
        or not 0.5 <= min_support_fraction <= 1
    ):
        raise ValueError("Invalid RNA replicate agreement settings")
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    seen: set[tuple[str, str, str]] = set()
    for row in calls:
        key = (row["gene_id"], row["condition"], row["replicate_id"])
        if key in seen or row["strand"] not in {"+", "-", "."}:
            raise ValueError(
                "RNA calls repeat a biological observation or strand"
            )
        seen.add(key)
        if row.get("candidate_start") is not None and (
            isinstance(row["candidate_start"], bool)
            or not isinstance(row["candidate_start"], int)
            or row["candidate_start"] < 0
        ):
            raise ValueError(
                "RNA candidate coordinates must be non-negative integers"
            )
        if row["status"] == "candidate":
            position = row.get("candidate_start")
            left, right = (
                row.get(name)
                for name in ("neighbourhood_start", "neighbourhood_end")
            )
            rise = row.get("fold_rise")
            if (
                position is None
                or left is None
                or right is None
                or row["strand"] == "."
                or any(
                    isinstance(value, bool) or not isinstance(value, int)
                    for value in (left, right)
                )
                or not 0 <= left <= position < right
                or not isinstance(rise, (float, int))
                or isinstance(rise, bool)
                or not math.isfinite(rise)
                or rise <= 1
            ):
                raise ValueError("RNA candidate edge metadata is invalid")
        grouped[(row["gene_id"], row["condition"])].append(row)
    summaries: list[dict[str, Any]] = []
    for (gene_id, condition), rows in sorted(grouped.items()):
        if (
            len(
                {
                    (r["contig"], r["strand"], r["annotation_anchor"])
                    for r in rows
                }
            )
            != 1
        ):
            raise ValueError(
                "Conflicting RNA gene metadata between replicates"
            )
        detected = sorted(
            [r for r in rows if r["status"] == "candidate"],
            key=lambda r: r["candidate_start"],
        )
        clusters: list[list[Mapping[str, Any]]] = []
        left = 0
        for right, row in enumerate(detected):
            while (
                row["candidate_start"] - detected[left]["candidate_start"]
                > agreement_tolerance
            ):
                left += 1
            clusters.append(detected[left : right + 1])
        cluster = (
            max(
                clusters,
                key=lambda group: (
                    len(group),
                    sum(r["fold_rise"] for r in group),
                    -abs(
                        float(median([r["candidate_start"] for r in group]))
                        - rows[0]["annotation_anchor"]
                    ),
                ),
            )
            if clusters
            else []
        )
        required = max(
            min_replicates, math.ceil(len(rows) * min_support_fraction)
        )
        if len(cluster) >= required:
            status = "replicate_supported_candidate"
        elif len(rows) == 1 and cluster:
            status = "provisional_single_replicate"
        elif len(detected) >= required and len(cluster) < required:
            status = "discordant_replicates"
        elif cluster:
            status = "insufficient_replicate_support"
        else:
            status = "not_detected"
        position = (
            int(median([r["candidate_start"] for r in cluster]))
            if cluster
            else None
        )
        stranded = all(r["strandedness"] != "unstranded" for r in rows)
        summaries.append(
            {
                "gene_id": gene_id,
                "condition": condition,
                "contig": rows[0]["contig"],
                "strand": rows[0]["strand"],
                "annotation_anchor": rows[0]["annotation_anchor"],
                "candidate_start": position,
                "status": status,
                "biological_observations": len(rows),
                "detected_replicates": len(detected),
                "supporting_replicates": len(cluster),
                "required_replicates": required,
                "neighbourhood_start": min(
                    r["neighbourhood_start"] for r in cluster
                )
                if cluster
                else None,
                "neighbourhood_end": max(
                    r["neighbourhood_end"] for r in cluster
                )
                if cluster
                else None,
                "replicate_start_min": min(
                    r["candidate_start"] for r in cluster
                )
                if cluster
                else None,
                "replicate_start_max": max(
                    r["candidate_start"] for r in cluster
                )
                if cluster
                else None,
                "mean_fold_rise": sum(r["fold_rise"] for r in cluster)
                / len(cluster)
                if cluster
                else None,
                "strand_evidence": "strand_specific"
                if stranded
                else "unstranded_ambiguous",
                "tss_status": "unvalidated_coverage_edge",
                "interval_status": (
                    "resolution_and_replicate_range_not_confidence_interval"
                ),
            }
        )
    return summaries


def compare_rna_starts(
    *, starts: Sequence[Mapping[str, Any]], condition: str, reference: str
) -> list[dict[str, Any]]:
    """Compare supported edges without inventing promoter-usage tests.

    Args:
        starts: Consensus starts for multiple conditions.
        condition: Numerator/target condition.
        reference: Reference condition.

    Returns:
        Signed shift in transcriptional orientation and interval-overlap audit.
        Positive shifts are downstream. No differential promoter-usage p-value
        is inferred from differences in detection coordinates.

    Raises:
        ValueError: Conditions are equal, absent, or starts repeat.
    """
    available = {row["condition"] for row in starts}
    if condition == reference or not {condition, reference} <= available:
        raise ValueError("RNA start comparison needs two existing conditions")
    keyed: dict[tuple[str, str], Mapping[str, Any]] = {}
    for row in starts:
        key = (row["gene_id"], row["condition"])
        if key in keyed:
            raise ValueError("RNA consensus rows repeat a gene/condition")
        keyed[key] = row
    output: list[dict[str, Any]] = []
    for gene_id in sorted({row["gene_id"] for row in starts}):
        target, baseline = (
            keyed.get((gene_id, condition)),
            keyed.get((gene_id, reference)),
        )
        estimable = bool(
            target
            and baseline
            and all(
                r["status"] == "replicate_supported_candidate"
                for r in (target, baseline)
            )
        )
        strand = (target or baseline or {})["strand"]
        output.append(
            {
                "gene_id": gene_id,
                "condition": condition,
                "reference": reference,
                "condition_start": target.get("candidate_start")
                if target
                else None,
                "reference_start": baseline.get("candidate_start")
                if baseline
                else None,
                "oriented_shift_bp": (
                    target["candidate_start"] - baseline["candidate_start"]
                )
                * (1 if strand == "+" else -1)
                if estimable and target is not None and baseline is not None
                else None,
                "neighbourhoods_overlap": max(
                    target["neighbourhood_start"],
                    baseline["neighbourhood_start"],
                )
                < min(
                    target["neighbourhood_end"], baseline["neighbourhood_end"]
                )
                if estimable and target is not None and baseline is not None
                else None,
                "status": "descriptive_candidate_shift"
                if estimable
                else "not_estimable",
                "promoter_usage_p_value": None,
            }
        )
    return output


def extract_rna_upstream(
    *,
    index: GeneIndex,
    genome: Genome,
    starts: Sequence[Mapping[str, Any]],
    condition: str,
    identifiers: Sequence[str],
    length: int | None = 1000,
    offset: int = 0,
    min_length: int = 1,
    max_ambiguous_fraction: float = 1.0,
    mask_lowercase: bool = False,
    fallback: str = "annotation",
    allow_provisional: bool = False,
) -> tuple[list[Region], list[dict[str, Any]]]:
    """Refine promoters using supported edges entirely within the safe gap.

    Args:
        index: Original full annotation; it is never rewritten.
        genome: Open genome.
        starts: Consensus edge records.
        condition: Condition used identically for positive and negative genes.
        identifiers: Exact target IDs.
        length: Maximum flank length, or full gap when None.
        offset: Bases skipped nearest the applied anchor.
        min_length: Minimum returned length.
        max_ambiguous_fraction: Sequence ambiguity filter.
        mask_lowercase: Mask lowercase bases before extraction.
        fallback: Audited annotation fallback or exclusion for unusable starts.
        allow_provisional: Permit a singleton edge; it stays unvalidated.

    Returns:
        Original-compatible regions and an anchor audit. Starts inside any
        annotated gene, including the target, never authorise genic extraction.

    Raises:
        ValueError: Inputs, conditions or extraction settings are invalid.
    """
    if fallback not in {"annotation", "exclude"}:
        raise ValueError("RNA anchor fallback must be annotation or exclude")
    keyed: dict[str, Mapping[str, Any]] = {}
    for start_row in starts:
        if start_row["condition"] == condition:
            if start_row["gene_id"] in keyed:
                raise ValueError(
                    "RNA starts repeat a gene within the anchor condition"
                )
            keyed[start_row["gene_id"]] = start_row
    if not keyed:
        raise ValueError("No RNA starts for the selected anchor condition")
    original = extract_regions(
        index=index,
        genome=genome,
        identifiers=identifiers,
        direction="upstream",
        length=length,
        offset=offset,
        min_length=min_length,
        max_ambiguous_fraction=max_ambiguous_fraction,
        mask_lowercase=mask_lowercase,
    )
    output: list[Region] = []
    audit: list[dict[str, Any]] = []
    for base in original:
        gene = index.genes[base.gene_id]
        annotation_anchor = gene.start if gene.strand != "-" else gene.end - 1
        row = keyed.get(gene.gene_id)
        candidate = row.get("candidate_start") if row else None
        eligible = bool(
            row
            and (
                row["status"] == "replicate_supported_candidate"
                or (
                    allow_provisional
                    and row["status"] == "provisional_single_replicate"
                )
            )
        )
        reason = "unsupported_or_missing_rna_start"
        if row and (
            row["contig"] != gene.contig or row["strand"] != gene.strand
        ):
            raise ValueError(
                "RNA anchor metadata differs from target annotation"
            )
        if candidate is not None and (
            isinstance(candidate, bool)
            or not isinstance(candidate, int)
            or not 0 <= candidate < genome.lengths[gene.contig]
        ):
            raise ValueError("RNA anchor falls outside the genome")
        free_start = free_end = gene.start
        if gene.strand != ".":
            free_start, free_end, _ = index.free_flank(
                gene=gene, direction="upstream"
            )
        if eligible and candidate is not None and gene.strand != ".":
            boundary = candidate if gene.strand == "+" else candidate + 1
            if (
                gene.start <= candidate < gene.end
                and candidate != annotation_anchor
            ):
                reason = "rna_start_inside_annotated_target"
                eligible = False
            elif not free_start <= boundary <= free_end:
                reason = "rna_start_outside_contiguous_safe_gap"
                eligible = False
            else:
                reason = "rna_supported_safe_anchor"
        else:
            eligible = False
        region = base
        anchor_source = "annotation_fallback"
        applied: int | None = annotation_anchor
        if eligible:
            assert candidate is not None
            start, end = free_start, free_end
            if gene.strand == "+":
                end = candidate
            else:
                start = candidate + 1
            available = end - start
            skip = min(offset, available)
            keep = min(
                available - skip, available if length is None else length
            )
            if gene.strand == "+":
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
            region = replace(
                base,
                start=start,
                end=end,
                available_length=available,
                sequence=sequence.upper() if status == "retained" else "",
                status=status,
            )
            anchor_source, applied = "rna_coverage_candidate", candidate
        elif fallback == "exclude":
            region = replace(base, sequence="", status="unusable_rna_anchor")
            anchor_source, applied = "excluded", None
        output.append(region)
        audit.append(
            {
                "sequence_id": region.sequence_id,
                "gene_id": gene.gene_id,
                "condition": condition,
                "contig": gene.contig,
                "strand": gene.strand,
                "annotation_anchor": annotation_anchor,
                "rna_candidate_start": candidate,
                "applied_anchor": applied,
                "anchor_source": anchor_source,
                "anchor_reason": reason,
                "rna_status": row.get("status") if row else "missing",
                "supporting_replicates": row.get("supporting_replicates")
                if row
                else 0,
                "strand_evidence": row.get("strand_evidence")
                if row
                else "missing",
                "neighbourhood_start": row.get("neighbourhood_start")
                if row
                else None,
                "neighbourhood_end": row.get("neighbourhood_end")
                if row
                else None,
                "extraction_status": region.status,
                "tss_status": "unvalidated_coverage_edge",
            }
        )
    LOGGER.info(
        "Applied %d safe RNA anchors to %d genes",
        sum(r["anchor_source"] == "rna_coverage_candidate" for r in audit),
        len(audit),
    )
    return output, audit
