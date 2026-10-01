"""Optional genomic and functional evidence, independent of model fitting."""

import csv
from bisect import bisect_left
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from intergenic_regions.extraction import merge_intervals
from intergenic_regions.io import open_text
from intergenic_regions.models import Region


class PeakIndex:
    """Index merged accessibility or enhancer BED annotations."""

    def __init__(self, *, intervals: dict[str, list[tuple[int, int]]]) -> None:
        """Merge peaks so overlapping annotations cannot inflate coverage.

        Args:
            intervals: Half-open BED intervals grouped by contig.
        """
        self.intervals = {
            c: merge_intervals(intervals=spans)
            for c, spans in intervals.items()
        }
        self.ends = {
            c: [e for _, e in spans] for c, spans in self.intervals.items()
        }

    def overlap_bases(self, *, contig: str, start: int, end: int) -> int:
        """Measure union overlap in base pairs.

        Args:
            contig: Record identifier.
            start: Zero-based start.
            end: Exclusive end.

        Returns:
            Number of overlapping bases, without double counting.

        Raises:
            ValueError: Coordinates are negative or reversed.
        """
        if start < 0 or end < start:
            raise ValueError("Invalid evidence query interval")
        spans = self.intervals.get(contig, [])
        position = bisect_left(a=self.ends.get(contig, []), x=start + 1)
        total = 0
        for left, right in spans[position:]:
            if left >= end:
                break
            total += max(0, min(end, right) - max(start, left))
        return total


def read_bed(*, path: Path) -> PeakIndex:
    """Read BED3 or wider peak files, ignoring track/browser/comment lines.

    Args:
        path: BED with zero-based, half-open coordinates.

    Returns:
        A merged genomic peak index.

    Raises:
        ValueError: A record has malformed columns or coordinates.
    """
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(path=path) as stream:
        for number, raw in enumerate(stream, start=1):
            if not raw.strip() or raw.startswith(("#", "track ", "browser ")):
                continue
            fields = raw.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{number}: BED requires at least three TSV columns"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{number}: invalid BED coordinates"
                ) from exc
            if start < 0 or end <= start or not fields[0]:
                raise ValueError(f"{path}:{number}: invalid BED interval")
            intervals[fields[0]].append((start, end))
    return PeakIndex(intervals=dict(intervals))


def read_functional_evidence(*, path: Path) -> dict[str, list[dict[str, str]]]:
    """Read optional user-supplied functional results keyed by gene ID.

    Args:
        path: TSV with gene_id, evidence_type, value and source columns.

    Returns:
        Evidence records grouped by exact gene ID.

    Raises:
        ValueError: Required columns or record fields are missing.
    """
    records: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open_text(path=path) as stream:
        reader = csv.DictReader(f=stream, delimiter="\t")
        required = {"gene_id", "evidence_type", "value", "source"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(
                "Evidence TSV requires gene_id/evidence_type/value/source"
            )
        for row in reader:
            if any(not row.get(key) for key in required) or None in row:
                raise ValueError(
                    "Evidence table contains malformed or empty fields"
                )
            records[row["gene_id"]].append(
                {key: row[key] for key in sorted(required)}
            )
    return dict(records)


def annotate_evidence(
    *,
    regions: Sequence[Region],
    accessibility: PeakIndex | None = None,
    enhancers: PeakIndex | None = None,
    functional: dict[str, list[dict[str, str]]] | None = None,
) -> list[dict[str, Any]]:
    """Attach evidence without asserting enhancer function or changing labels.

    Args:
        regions: Retained genomic regions.
        accessibility: Optional accessibility BED index.
        enhancers: Optional user-supplied enhancer annotation BED index.
        functional: Optional linked gene-level results, kept distinct from
            direct interval evidence.

    Returns:
        Per-region overlap coverage and supplied functional evidence.
    """
    import json

    records: list[dict[str, Any]] = []
    for region in regions:
        if region.status != "retained":
            continue
        accessibility_bp = (
            accessibility.overlap_bases(
                contig=region.contig, start=region.start, end=region.end
            )
            if accessibility
            else None
        )
        enhancer_bp = (
            enhancers.overlap_bases(
                contig=region.contig, start=region.start, end=region.end
            )
            if enhancers
            else None
        )
        linked = (functional or {}).get(region.gene_id, [])
        records.append(
            {
                "sequence_id": region.sequence_id,
                "gene_id": region.gene_id,
                "contig": region.contig,
                "start": region.start,
                "end": region.end,
                "accessibility_overlap_bp": accessibility_bp,
                "accessibility_overlap_fraction": accessibility_bp
                / (region.end - region.start)
                if accessibility_bp is not None
                else None,
                "enhancer_annotation_overlap_bp": enhancer_bp,
                "enhancer_annotation_overlap_fraction": enhancer_bp
                / (region.end - region.start)
                if enhancer_bp is not None
                else None,
                "linked_gene_evidence_count": len(linked),
                "linked_gene_evidence_json": json.dumps(
                    obj=linked, sort_keys=True
                ),
                "evidence_status": "supplied_evidence"
                if linked or accessibility_bp or enhancer_bp
                else "no_overlap_in_supplied_data"
                if accessibility or enhancers or functional is not None
                else "not_supplied",
            }
        )
    return records
