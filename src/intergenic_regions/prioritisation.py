"""Transparent candidate ranking; sequence scores never validate enhancers."""

import json
import math
from collections.abc import Mapping, Sequence
from typing import Any


def prioritise_candidates(
    *, rows: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    """Rank by supplied evidence tier, then held-out sequence signature.

    Accessibility and regulatory annotations provide contextual support.
    Positive and negative gene-level assays are counted separately; neither
    establishes the function of the extracted interval. Unknown assay values
    remain uninterpreted. Missing data are never converted to negative assays.

    Args:
        rows: Sequence records with optional evidence and signature scores.

    Returns:
        Copied records with deterministic ranks, evidence tiers and explicit
        uncertainty. The ordering is a triage rule, not a calibrated model.

    Raises:
        ValueError: Identifiers repeat or scores/evidence are malformed.
    """
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for original in rows:
        row = dict(original)
        identifier = row.get("sequence_id")
        if (
            not isinstance(identifier, str)
            or not identifier
            or identifier in seen
        ):
            raise ValueError("Candidate identifiers must be unique tokens")
        seen.add(identifier)
        score = row.get("held_out_signature_score")
        if score is not None and (
            not isinstance(score, (int, float))
            or not math.isfinite(score)
            or not 0 <= score <= 1
        ):
            raise ValueError("Candidate signature scores must lie in [0, 1]")
        accessible = (row.get("accessibility_overlap_bp") or 0) > 0
        annotated = (row.get("enhancer_annotation_overlap_bp") or 0) > 0 or (
            row.get("reference_support_count") or 0
        ) > 0
        try:
            linked = json.loads(row.get("linked_gene_evidence_json", "[]"))
            if not isinstance(linked, list):
                raise ValueError("Linked evidence must be a JSON list")
            values = [
                str(record["value"]).strip().lower() for record in linked
            ]
        except (TypeError, KeyError, json.JSONDecodeError) as exc:
            raise ValueError("Malformed linked gene evidence") from exc
        supports = sum(
            v in {"positive", "supported", "1", "true"} for v in values
        )
        contradicts = sum(
            v in {"negative", "not_supported", "0", "false"} for v in values
        )
        tier = (
            3
            if accessible and annotated
            else 2
            if accessible or annotated
            else 1
            if supports and not contradicts
            else 0
        )
        description = (
            "accessibility_and_annotation"
            if tier == 3
            else "accessibility_support"
            if accessible
            else "annotation_overlap"
            if annotated
            else "linked_gene_support"
            if tier == 1
            else "sequence_only"
        )
        row.update(
            {
                "priority_tier": tier,
                "support_category": description,
                "supporting_gene_evidence_count": supports,
                "contradicting_gene_evidence_count": contradicts,
                "uninterpreted_gene_evidence_count": len(values)
                - supports
                - contradicts,
                "enhancer_status": "unvalidated_candidate",
                "uncertainty": (
                    "Gene-level evidence is contradictory; "
                    "inspect assay context. "
                    if contradicts
                    else ""
                )
                + (
                    "Supplied context supports prioritisation; "
                    "interval-specific "
                    "function and target-gene linkage remain unestablished."
                    if tier
                    else "Sequence-only hypothesis; regulatory activity and "
                    "enhancer function are unestablished."
                ),
            }
        )
        output.append(row)
    output.sort(
        key=lambda r: (
            -r["priority_tier"],
            r["contradicting_gene_evidence_count"] > 0,
            -(r.get("held_out_signature_score") or 0),
            r["sequence_id"],
        )
    )
    for rank, row in enumerate(output, start=1):
        row["priority_rank"] = rank
    return output
