"""Variable-length regulatory search windows within existing safe flanks."""

import hashlib
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any
from urllib.parse import quote

from intergenic_regions.background import (
    canonical_sequence,
    check_sequence_sets,
)
from intergenic_regions.models import Gene, Region

DEFAULT_WINDOW_LENGTHS = (100, 200, 400, 800, 1600)


@dataclass(frozen=True, slots=True, kw_only=True)
class RegulatoryWindow:
    """A search interval, without any assertion of enhancer function.

    Attributes:
        parent_id: Source FASTA/flank identifier.
        label: Supplied parent class, ``positive`` or ``negative``.
        sequence_start: Inclusive offset in transcription-oriented sequence.
        sequence_end: Exclusive offset in transcription-oriented sequence.
        sequence: Already strand-oriented window sequence.
        gene_id: Optional target gene identifier.
        contig: Optional genomic contig.
        genomic_start: Optional zero-based inclusive genomic coordinate.
        genomic_end: Optional zero-based exclusive genomic coordinate.
        strand: Target strand; ``.`` for coordinate-free FASTA.
        direction: Source flank direction, if known.
        distance_to_gene_start: Signed midpoint distance from the annotated
            gene's 5-prime base; upstream is negative. This is a TSS proxy.
    """

    parent_id: str
    label: str
    sequence_start: int
    sequence_end: int
    sequence: str
    gene_id: str | None = None
    contig: str | None = None
    genomic_start: int | None = None
    genomic_end: int | None = None
    strand: str = "."
    direction: str | None = None
    distance_to_gene_start: float | None = None

    def __post_init__(self) -> None:
        """Reject inconsistent offsets, coordinates and labels."""
        if (
            not self.parent_id
            or any(c.isspace() for c in self.parent_id)
            or self.label not in {"positive", "negative"}
            or self.strand not in {"+", "-", "."}
            or self.sequence_start < 0
            or self.sequence_end <= self.sequence_start
            or len(self.sequence) != self.sequence_end - self.sequence_start
        ):
            raise ValueError("Invalid regulatory window")
        genomic = self.genomic_start is not None
        if genomic != (self.genomic_end is not None) or genomic != (
            self.contig is not None
        ):
            raise ValueError("Window genomic coordinates must be complete")
        if genomic and (
            self.genomic_start < 0  # type: ignore[operator]
            or self.genomic_end - self.genomic_start  # type: ignore[operator]
            != len(self.sequence)
            or not self.gene_id
            or self.strand == "."
            or self.direction not in {"upstream", "downstream"}
        ):
            raise ValueError("Window genomic span or target is inconsistent")
        if self.distance_to_gene_start is not None and (
            not genomic or not math.isfinite(self.distance_to_gene_start)
        ):
            raise ValueError("Window gene-start distance requires coordinates")

    @property
    def sequence_id(self) -> str:
        """Return a deterministic identifier distinct from the source flank."""
        return (
            f"{quote(self.parent_id, safe='._-')}|window:"
            f"{self.sequence_start}-{self.sequence_end}"
        )

    @property
    def width(self) -> int:
        """Return the complete search-window length in bases."""
        return self.sequence_end - self.sequence_start


def make_windows(
    *,
    positive: dict[str, str],
    negative: dict[str, str],
    lengths: Sequence[int] = DEFAULT_WINDOW_LENGTHS,
    step: int = 50,
    max_windows: int = 50000,
    regions: Sequence[Region] | None = None,
    genes: Sequence[Gene] | None = None,
) -> tuple[list[RegulatoryWindow], list[dict[str, Any]]]:
    """Tile each source at several lengths, never padding or crossing it.

    A final end-aligned window covers the flank's last bases even when the
    regular step does not reach them. Too-short flanks are audited per scale;
    truncated windows are never misrepresented as full-length observations.

    Args:
        positive: Foreground parent sequences.
        negative: Control parent sequences.
        lengths: Unique search scales, from 20 to 100000 bp.
        step: Offset step in bases, no larger than the shortest window.
        max_windows: Hard limit, checked before allocating any windows.
        regions: Optional matching retained genomic flanks.
        genes: Optional full gene records for annotated-start distances.

    Returns:
        Complete windows and one availability audit row per parent and scale.

    Raises:
        ValueError: Inputs, limits or coordinate metadata are inconsistent.
    """
    check_sequence_sets(positive=positive, negative=negative)
    if (
        not lengths
        or len(set(lengths)) != len(lengths)
        or any(
            not isinstance(n, int)
            or isinstance(n, bool)
            or not 20 <= n <= 100000
            for n in lengths
        )
        or not isinstance(step, int)
        or isinstance(step, bool)
        or not 1 <= step <= min(lengths)
        or not isinstance(max_windows, int)
        or isinstance(max_windows, bool)
        or max_windows < 1
    ):
        raise ValueError(
            "Invalid window lengths, step or maximum window count"
        )
    sources = {**positive, **negative}
    metadata = {r.sequence_id: r for r in regions or ()}
    gene_by_id = {g.gene_id: g for g in genes or ()}
    if regions is not None and (
        len(metadata) != len(regions) or set(metadata) != set(sources)
    ):
        raise ValueError("Every window parent needs unique flank metadata")
    for identifier, source_region in metadata.items():
        if (
            source_region.status != "retained"
            or source_region.strand not in {"+", "-"}
            or source_region.sequence != sources[identifier]
            or source_region.end - source_region.start
            != len(source_region.sequence)
        ):
            raise ValueError(
                "Window parent metadata does not match its sequence"
            )
    audit: list[dict[str, Any]] = []
    total = 0
    for identifier, sequence in sources.items():
        for width in sorted(lengths):
            remainder = len(sequence) - width
            count = (
                remainder // step + 1 + int(remainder % step != 0)
                if remainder >= 0
                else 0
            )
            total += count
            audit.append(
                {
                    "parent_id": identifier,
                    "label": "positive"
                    if identifier in positive
                    else "negative",
                    "parent_length": len(sequence),
                    "window_length": width,
                    "windows": count,
                    "status": "available" if count else "parent_too_short",
                }
            )
    if total > max_windows:
        raise ValueError(
            f"Multi-scale search needs {total} windows, exceeding "
            f"max_windows={max_windows}; increase the step, restrict inputs "
            "or explicitly raise the limit"
        )
    windows: list[RegulatoryWindow] = []
    for identifier, sequence in sources.items():
        region = metadata.get(identifier)
        gene = gene_by_id.get(region.gene_id) if region else None
        if (
            gene is not None
            and region is not None
            and (gene.contig != region.contig or gene.strand != region.strand)
        ):
            raise ValueError("Window gene and flank metadata disagree")
        for width in sorted(lengths):
            last = len(sequence) - width
            if last < 0:
                continue
            starts = list(range(0, last + 1, step))
            if starts[-1] != last:
                starts.append(last)
            for first in starts:
                end = first + width
                left, right = (
                    (region.end - end, region.end - first)
                    if region is not None and region.strand == "-"
                    else (region.start + first, region.start + end)
                    if region is not None
                    else (None, None)
                )
                distance = None
                if gene is not None:
                    assert left is not None and right is not None
                    midpoint = (left + right - 1) / 2
                    anchor = gene.start if gene.strand == "+" else gene.end - 1
                    distance = (midpoint - anchor) * (
                        1 if gene.strand == "+" else -1
                    )
                windows.append(
                    RegulatoryWindow(
                        parent_id=identifier,
                        label="positive"
                        if identifier in positive
                        else "negative",
                        sequence_start=first,
                        sequence_end=end,
                        sequence=sequence[first:end],
                        gene_id=region.gene_id if region else None,
                        contig=region.contig if region else None,
                        genomic_start=left,
                        genomic_end=right,
                        strand=region.strand if region else ".",
                        direction=region.direction if region else None,
                        distance_to_gene_start=distance,
                    )
                )
    return windows, audit


def assign_window_groups(
    *,
    windows: Sequence[RegulatoryWindow],
    parent_groups: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Group source flanks and exact shared windows before validation.

    Group unions use sequence equality, including reverse complements, and
    optional user family/contig assignments. Labels never affect grouping.
    Partial homology still requires suitable user-supplied groups.

    Args:
        windows: All scales, before filtering scales for model training.
        parent_groups: Optional complete parent-to-independent-group mapping.

    Returns:
        Deterministic group assignments for every represented parent.

    Raises:
        ValueError: A supplied parent assignment is missing or empty.
    """
    parents = sorted({w.parent_id for w in windows})
    links = {identifier: identifier for identifier in parents}

    def root(identifier: str) -> str:
        while links[identifier] != identifier:
            links[identifier] = links[links[identifier]]
            identifier = links[identifier]
        return identifier

    def union(first: str, second: str) -> None:
        a, b = sorted((root(first), root(second)))
        links[b] = a

    assigned: dict[str, str] = {}
    if parent_groups is not None:
        for identifier in parents:
            group = parent_groups.get(identifier)
            if not isinstance(group, str) or not group.strip():
                raise ValueError("Every window parent needs a non-empty group")
            if group in assigned:
                union(identifier, assigned[group])
            assigned[group] = identifier
    sequences: dict[bytes, str] = {}
    for window in windows:
        digest = hashlib.sha256(
            canonical_sequence(sequence=window.sequence).encode("ascii")
        ).digest()
        if digest in sequences:
            union(window.parent_id, sequences[digest])
        sequences[digest] = window.parent_id
    return {identifier: root(identifier) for identifier in parents}


def merge_window_candidates(
    *, rows: Sequence[Mapping[str, Any]], score_threshold: float = 0.75
) -> list[dict[str, Any]]:
    """Join overlapping supported windows into variable-length hypotheses.

    The union describes search support, not a measured enhancer boundary.
    Different-width models are not assumed to have comparable calibration.
    The reported maximum is an exploratory ranking statistic, without a
    region-level p-value or FDR claim.

    Args:
        rows: Window records with held_out_signature_score and parent offsets.
        score_threshold: Fixed display/triage threshold, never tuned on labels.

    Returns:
        Unranked interval unions with peak-window and scale provenance.

    Raises:
        ValueError: Scores, spans or the threshold are malformed.
    """
    if not math.isfinite(score_threshold) or not 0 <= score_threshold <= 1:
        raise ValueError("Window score threshold must lie in [0, 1]")
    selected: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for source in rows:
        row = dict(source)
        score = row.get("held_out_signature_score")
        if score is None:
            continue
        if (
            not isinstance(score, (int, float))
            or not math.isfinite(score)
            or not 0 <= score <= 1
            or not row.get("parent_id")
            or row["sequence_start"] < 0
            or row["sequence_end"] <= row["sequence_start"]
        ):
            raise ValueError("Invalid scored regulatory window")
        if score >= score_threshold:
            selected[row["parent_id"]].append(row)
    result: list[dict[str, Any]] = []
    for parent, items in sorted(selected.items()):
        components: list[list[dict[str, Any]]] = []
        right = -1
        for row in sorted(
            items, key=lambda r: (r["sequence_start"], r["sequence_end"])
        ):
            if not components or row["sequence_start"] >= right:
                components.append([row])
                right = row["sequence_end"]
            else:
                components[-1].append(row)
                right = max(right, row["sequence_end"])
        for members in components:
            peak = min(
                members,
                key=lambda r: (
                    -r["held_out_signature_score"],
                    r["sequence_id"],
                ),
            )
            first = min(r["sequence_start"] for r in members)
            end = max(r["sequence_end"] for r in members)
            genomic = peak.get("genomic_start") is not None
            result.append(
                {
                    "sequence_id": (
                        f"{quote(parent, safe='._-')}|candidate:{first}-{end}"
                    ),
                    "parent_id": parent,
                    "gene_id": peak.get("gene_id"),
                    "label": peak["label"],
                    "contig": peak.get("contig"),
                    "strand": peak.get("strand", "."),
                    "direction": peak.get("direction"),
                    "sequence_start": first,
                    "sequence_end": end,
                    "length": end - first,
                    "genomic_start": min(r["genomic_start"] for r in members)
                    if genomic
                    else None,
                    "genomic_end": max(r["genomic_end"] for r in members)
                    if genomic
                    else None,
                    "held_out_signature_score": peak[
                        "held_out_signature_score"
                    ],
                    "peak_window_id": peak["sequence_id"],
                    "peak_window_length": peak["window_length"],
                    "supporting_windows": len(members),
                    "supporting_lengths": ",".join(
                        str(n)
                        for n in sorted({r["window_length"] for r in members})
                    ),
                    "boundary_status": (
                        "union_of_search_windows_not_functional_boundary"
                    ),
                    "region_q_value": None,
                }
            )
    return result
