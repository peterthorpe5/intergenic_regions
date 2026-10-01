"""Prevent shared-region leakage and optionally match control composition."""

from bisect import bisect_left
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import replace
from typing import Any

from intergenic_regions.extraction import merge_intervals
from intergenic_regions.genome import reverse_complement, sequence_composition
from intergenic_regions.models import Region


def canonical_sequence(*, sequence: str) -> str:
    """Return a strand-independent sequence representation.

    Args:
        sequence: DNA sequence.

    Returns:
        Lexicographically smaller strand, in upper case.
    """
    upper = sequence.upper()
    return min(upper, reverse_complement(sequence=upper))


def check_sequence_sets(
    *, positive: dict[str, str], negative: dict[str, str]
) -> None:
    """Require non-empty, disjoint collections without duplicated sequences.

    Args:
        positive: Positive FASTA collection.
        negative: Negative FASTA collection.

    Raises:
        ValueError: IDs or exact/reverse-complement sequences repeat.
    """
    if not positive or not negative:
        raise ValueError(
            "Both positive and negative sets must contain sequences"
        )
    if positive.keys() & negative.keys():
        raise ValueError("Positive and negative sequence identifiers overlap")
    seen: set[str] = set()
    for sequence in [*positive.values(), *negative.values()]:
        key = canonical_sequence(sequence=sequence)
        if key in seen:
            raise ValueError(
                "Duplicate/shared sequence (including reverse complement); "
                "deduplicate before analysis"
            )
        seen.add(key)


def collapse_overlaps(
    *, regions: Sequence[Region], label: str
) -> tuple[list[Region], list[dict[str, Any]]]:
    """Keep one longest representative per connected component of overlaps.

    Args:
        regions: Extracted regions; exclusions are ignored.
        label: Positive or negative group for audit records.

    Returns:
        Retained representatives and exclusion audit rows.
    """
    by_contig: dict[str, list[Region]] = defaultdict(list)
    for region in regions:
        if region.status == "retained":
            by_contig[region.contig].append(region)
    retained: list[Region] = []
    audit: list[dict[str, Any]] = []
    for items in by_contig.values():
        components: list[list[Region]] = []
        maximum = -1
        for region in sorted(
            items, key=lambda r: (r.start, r.end, r.sequence_id)
        ):
            if not components or region.start >= maximum:
                components.append([region])
                maximum = region.end
            else:
                components[-1].append(region)
                maximum = max(maximum, region.end)
        for component in components:
            representative = min(
                component, key=lambda r: (-len(r.sequence), r.sequence_id)
            )
            retained.append(representative)
            for region in component:
                if region != representative:
                    audit.append(
                        {
                            "sequence_id": region.sequence_id,
                            "label": label,
                            "reason": "overlapping_region",
                            "representative": representative.sequence_id,
                        }
                    )
    return retained, audit


def sanitise_regions(
    *, positive: Sequence[Region], negative: Sequence[Region]
) -> tuple[list[Region], list[Region], list[dict[str, Any]]]:
    """Remove shared genomic intervals and duplicated DNA before analysis.

    Args:
        positive: Foreground extracted regions.
        negative: Candidate controls.

    Returns:
        Independent foreground/control representatives and an exclusion audit.
        Foreground regions take precedence over overlapping controls.

    Raises:
        ValueError: No foreground or controls survive filtering.
    """
    positives, audit = collapse_overlaps(regions=positive, label="positive")
    negatives, negative_audit = collapse_overlaps(
        regions=negative, label="negative"
    )
    audit.extend(negative_audit)
    occupied: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for region in positives:
        occupied[region.contig].append((region.start, region.end))
    occupied = {
        c: merge_intervals(intervals=spans) for c, spans in occupied.items()
    }
    starts = {c: [s for s, _ in spans] for c, spans in occupied.items()}
    seen: dict[str, str] = {}
    groups: list[list[Region]] = [[], []]
    for group_index, items in enumerate((positives, negatives)):
        for region in sorted(items, key=lambda r: r.sequence_id):
            reason, representative = "", ""
            if group_index and region.contig in occupied:
                spans = occupied[region.contig]
                position = (
                    bisect_left(a=starts[region.contig], x=region.end) - 1
                )
                if position >= 0 and spans[position][1] > region.start:
                    reason = "overlaps_positive_region"
            key = canonical_sequence(sequence=region.sequence)
            if not reason and key in seen:
                reason, representative = "duplicate_sequence", seen[key]
            if reason:
                audit.append(
                    {
                        "sequence_id": region.sequence_id,
                        "label": "negative" if group_index else "positive",
                        "reason": reason,
                        "representative": representative,
                    }
                )
            else:
                seen[key] = region.sequence_id
                groups[group_index].append(region)
    if not groups[0] or not groups[1]:
        raise ValueError(
            "No positive or negative regions after leakage filtering"
        )
    return groups[0], groups[1], audit


def match_background(
    *,
    positive: Sequence[Region],
    negative: Sequence[Region],
    ratio: int = 1,
    max_gc_difference: float = 0.1,
    max_length_ratio: float = 1.5,
) -> tuple[list[Region], list[dict[str, Any]]]:
    """Select controls deterministically by length and GC, without replacement.

    Args:
        positive: Foreground representatives.
        negative: Independent candidate controls.
        ratio: Number of controls per foreground region.
        max_gc_difference: Largest permitted absolute GC difference.
        max_length_ratio: Largest permitted longer/shorter length ratio.

    Returns:
        Selected controls and per-pair matching diagnostics.

    Raises:
        ValueError: Parameters are invalid or complete matching is impossible.
    """
    if ratio < 1 or not 0 <= max_gc_difference <= 1 or max_length_ratio < 1:
        raise ValueError("Invalid background matching parameters")
    if not positive or len(negative) < len(positive) * ratio:
        raise ValueError(
            "Insufficient controls for complete background matching"
        )
    candidates = {r.sequence_id: r for r in negative}
    composition = {
        key: sequence_composition(sequence=r.sequence)
        for key, r in candidates.items()
    }
    selected: list[Region] = []
    audit: list[dict[str, Any]] = []
    # The most constrained foreground is allocated first.
    eligible: dict[str, list[tuple[float, str, float, float]]] = {}
    for region in positive:
        gc = float(
            sequence_composition(sequence=region.sequence)["gc_fraction"]
        )
        eligible[region.sequence_id] = []
        for key, candidate in candidates.items():
            gc_difference = abs(gc - float(composition[key]["gc_fraction"]))
            length_ratio = max(
                len(region.sequence), len(candidate.sequence)
            ) / min(len(region.sequence), len(candidate.sequence))
            if (
                gc_difference <= max_gc_difference
                and length_ratio <= max_length_ratio
            ):
                eligible[region.sequence_id].append(
                    (
                        gc_difference + length_ratio - 1,
                        key,
                        gc_difference,
                        length_ratio,
                    )
                )
    for region in sorted(
        positive, key=lambda r: (len(eligible[r.sequence_id]), r.sequence_id)
    ):
        options = sorted(
            o for o in eligible[region.sequence_id] if o[1] in candidates
        )
        if len(options) < ratio:
            raise ValueError(
                f"Cannot match controls for {region.sequence_id}; "
                "supply more controls or relax tolerances"
            )
        for _, key, gc_difference, length_ratio in options[:ratio]:
            selected.append(replace(candidates.pop(key)))
            audit.append(
                {
                    "positive_id": region.sequence_id,
                    "negative_id": key,
                    "gc_difference": gc_difference,
                    "length_ratio": length_ratio,
                }
            )
    return selected, audit
