"""Multi-scale motif localisation, weak-label learning and region reporting."""

import csv
import logging
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from intergenic_regions.genome import sequence_composition
from intergenic_regions.io import write_fasta, write_json, write_tsv
from intergenic_regions.models import Gene, Region
from intergenic_regions.prioritisation import prioritise_candidates
from intergenic_regions.reporting import write_report
from intergenic_regions.windows import (
    DEFAULT_WINDOW_LENGTHS,
    RegulatoryWindow,
    assign_window_groups,
    make_windows,
    merge_window_candidates,
)

LOGGER = logging.getLogger(__name__)
REGION_NOTES = [
    "Search lengths are configurable analysis scales, not universal enhancer "
    "lengths. Too-short flanks are audited; windows are never padded or "
    "extended across the original gene-boundary exclusions.",
    "Labels are inherited from source flanks. Per-scale models use "
    "parent-balanced training weights and hold out complete parent/group "
    "units, including exact shared windows across scales. Partial homology "
    "still needs appropriate family/contig groups.",
    "Validation measures source-class resemblance using one mean prediction "
    "per parent at each scale, not enhancer detection. Short word features "
    "do not explicitly model transcription-factor spacing or cooperativity.",
    "Candidate intervals join overlapping above-threshold windows. Their "
    "boundaries, peak scores and support across scales are exploratory; "
    "longer gaps provide more search opportunities. Scores from different "
    "scales are not assumed to have comparable calibration. No region-level "
    "p-value, FDR or validated enhancer length is assigned.",
    "Motif density uses positive-enriched motifs selected on the complete "
    "dataset. It is descriptive localisation, separate from the training-only "
    "word vocabulary. Related motifs and overlapping sites are not "
    "independent evidence. Motif q-values remain source-level tests.",
    "Accessibility and other evidence are re-queried at window/candidate "
    "coordinates, remain optional and do not establish enhancer function or "
    "target-gene linkage. Gene-start distances use annotation 5-prime bases, "
    "which are TSS proxies; coordinate-free FASTA uses sequence offsets.",
    "Strictly intergenic flanks cover only the extracted part of the genome; "
    "intragenic, more distal or other-side regulatory elements can be missed. "
    "Increase --length or use --full-gap to explore larger safe flanks.",
]


def window_motif_counts(
    *,
    windows: Sequence[RegulatoryWindow],
    sites: Sequence[Mapping[str, Any]],
    motif_ids: Sequence[str],
) -> dict[str, dict[str, Any]]:
    """Count fully contained physical motif sites in each search window.

    Args:
        windows: Strand-oriented search intervals.
        sites: Motif sites in source-sequence offsets.
        motif_ids: Selected positive-enriched motifs for descriptive plots.

    Returns:
        Site count, motif diversity and sites per kb of full window span.

    Raises:
        ValueError: Site coordinates or selected motif IDs are malformed.
    """
    if len(set(motif_ids)) != len(motif_ids):
        raise ValueError("Window motif identifiers must be unique")
    selected = set(motif_ids)
    physical: dict[str, set[tuple[int, int, str]]] = defaultdict(set)
    for row in sites:
        if row["motif_id"] not in selected:
            continue
        first, end = int(row["start"]), int(row["end"])
        if first < 0 or end <= first:
            raise ValueError("Invalid source motif site coordinates")
        physical[row["sequence_id"]].add((first, end, row["motif_id"]))
    ordered = {parent: sorted(items) for parent, items in physical.items()}
    offsets = {
        parent: [s[0] for s in items] for parent, items in ordered.items()
    }
    result: dict[str, dict[str, Any]] = {}
    for window in windows:
        candidates = ordered.get(window.parent_id, [])
        starts = offsets.get(window.parent_id, [])
        first = bisect_left(starts, window.sequence_start)
        last = bisect_right(starts, window.sequence_end - 1)
        contained = [
            s for s in candidates[first:last] if s[1] <= window.sequence_end
        ]
        result[window.sequence_id] = {
            "motif_sites": len(contained),
            "motif_diversity": len({s[2] for s in contained}),
            "motif_sites_per_kb": 1000 * len(contained) / window.width,
        }
    return result


def window_position_profiles(
    *,
    rows: Sequence[Mapping[str, Any]],
    bin_width: int = 100,
    limit: int = 5000,
) -> list[dict[str, Any]]:
    """Summarise position by averaging within parent before between parents.

    Args:
        rows: Window records with density, optional score and signed distance.
        bin_width: Position bin size in bp.
        limit: Absolute plot/profile distance limit; full windows are retained.

    Returns:
        Per-class, per-scale bins with parent coverage and mean measurements.
        Bins use annotated-start distances when supplied, otherwise oriented
        sequence-midpoint offsets. Unestimable scores remain None.

    Raises:
        ValueError: Bin settings or position coordinate systems are invalid.
    """
    if (
        not isinstance(bin_width, int)
        or isinstance(bin_width, bool)
        or not isinstance(limit, int)
        or isinstance(limit, bool)
        or bin_width < 1
        or limit < bin_width
        or 2 * math.ceil(limit / bin_width) > 500
    ):
        raise ValueError("Position profiles require 1 to 500 finite bins")
    has_distance = [r.get("distance_to_gene_start") is not None for r in rows]
    if any(has_distance) and not all(has_distance):
        raise ValueError(
            "Window position profiles cannot mix coordinate systems"
        )
    coordinate = (
        "annotated_gene_5prime_base"
        if any(has_distance)
        else "oriented_sequence_offset"
    )
    grouped: dict[tuple[str, int, int], dict[str, list[Mapping[str, Any]]]] = (
        defaultdict(lambda: defaultdict(list))
    )
    for row in rows:
        position = row.get("distance_to_gene_start")
        if position is None:
            position = (row["sequence_start"] + row["sequence_end"] - 1) / 2
        if not -limit <= position < limit:
            continue
        first = math.floor(position / bin_width) * bin_width
        grouped[(row["label"], row["window_length"], first)][
            row["parent_id"]
        ].append(row)
    result: list[dict[str, Any]] = []
    for (label, width, first), sources in sorted(grouped.items()):
        source_scores: list[float] = []
        source_density: list[float] = []
        for items in sources.values():
            scores = [
                r["held_out_signature_score"]
                for r in items
                if r.get("held_out_signature_score") is not None
            ]
            if scores:
                source_scores.append(sum(scores) / len(scores))
            source_density.append(
                sum(r["motif_sites_per_kb"] for r in items) / len(items)
            )
        result.append(
            {
                "label": label,
                "window_length": width,
                "bin_start": max(first, -limit),
                "bin_end": min(first + bin_width, limit),
                "coordinate_system": coordinate,
                "parents": len(sources),
                "scored_parents": len(source_scores),
                "mean_held_out_signature_score": sum(source_scores)
                / len(source_scores)
                if source_scores
                else None,
                "mean_motif_sites_per_kb": sum(source_density)
                / len(source_density),
            }
        )
    return result


@dataclass(frozen=True, slots=True, kw_only=True)
class _NamedRegion(Region):
    identifier: str

    @property
    def sequence_id(self) -> str:
        return self.identifier


def _coordinate_regions(rows: Sequence[Mapping[str, Any]]) -> list[Region]:
    return [
        _NamedRegion(
            identifier=row["sequence_id"],
            gene_id=row["gene_id"],
            contig=row["contig"],
            start=row["genomic_start"],
            end=row["genomic_end"],
            strand=row["strand"],
            direction=row["direction"],
            sequence="",
            status="retained",
            stop_reason="within_source_intergenic_flank",
            available_length=row["genomic_end"] - row["genomic_start"],
        )
        for row in rows
        if row.get("genomic_start") is not None
    ]


def _attach_context(
    directory: Path, rows: list[dict[str, Any]], settings: dict[str, Any]
) -> list[dict[str, Any]]:
    from intergenic_regions.workflows import attach_evidence

    regions = _coordinate_regions(rows)
    if not regions:
        return rows
    evidence, references = attach_evidence(
        directory=directory, regions=regions, **settings
    )
    by_id = {r["sequence_id"]: r for r in evidence}
    names: dict[str, set[str]] = defaultdict(set)
    for reference in references:
        if reference["overlap_bp"] > 0:
            names[reference["sequence_id"]].add(reference["reference_name"])
    for row in rows:
        context = by_id.get(row["sequence_id"], {})
        row.update(context)
        row["reference_support_count"] = len(names[row["sequence_id"]])
    return rows


def _selected_sites(
    motif_directory: Path,
    parents: dict[str, str],
    maximum: int,
    q_value: float,
    both_strands: bool,
) -> tuple[list[str], list[dict[str, Any]]]:
    from intergenic_regions.motifs import scan_pattern

    with (motif_directory / "motif_enrichment.tsv").open(
        encoding="utf-8"
    ) as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    selected = [
        r
        for r in rows
        if float(r["q_value"]) <= q_value
        and float(r["positive_fraction"]) > float(r["negative_fraction"])
    ][:maximum]
    identifiers = [r["motif_id"] for r in selected]
    with (motif_directory / "motif_sites.tsv").open(
        encoding="utf-8"
    ) as stream:
        sites = [
            r
            for r in csv.DictReader(stream, delimiter="\t")
            if r["motif_id"] in identifiers
        ]
    for motif in selected:
        if motif["kind"] == "kmer":
            for identifier, sequence in parents.items():
                sites.extend(
                    {
                        "sequence_id": identifier,
                        "motif_id": motif["motif_id"],
                        **hit,
                    }
                    for hit in scan_pattern(
                        sequence=sequence,
                        pattern=motif["consensus"],
                        both_strands=both_strands,
                    )
                )
    return identifiers, sites


def _read_parent_groups(
    path: Path, identifiers: Sequence[str]
) -> dict[str, str]:
    groups: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"sequence_id", "group"} <= set(reader.fieldnames or ()):
            raise ValueError("Groups TSV requires sequence_id and group")
        for row in reader:
            identifier, group = row.get("sequence_id"), row.get("group")
            if not identifier or not group or identifier in groups:
                raise ValueError("Malformed or duplicated group assignment")
            groups[identifier] = group
    if set(identifiers) - groups.keys():
        raise ValueError("Every analysed sequence needs an explicit group")
    return groups


def regional_outputs(
    *,
    directory: Path,
    positive: dict[str, str],
    negative: dict[str, str],
    motif_directory: Path,
    enabled: bool = True,
    use_ai: bool = True,
    lengths: Sequence[int] = DEFAULT_WINDOW_LENGTHS,
    step: int = 50,
    max_windows: int = 50000,
    score_threshold: float = 0.75,
    max_motifs: int = 20,
    motif_q_value: float = 0.05,
    both_strands: bool = True,
    position_bin_width: int = 100,
    position_limit: int = 5000,
    max_locus_plots: int = 6,
    learning_settings: dict[str, Any] | None = None,
    groups_path: Path | None = None,
    parent_groups: Mapping[str, str] | None = None,
    regions: Sequence[Region] | None = None,
    genes: Sequence[Gene] | None = None,
    evidence_settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write multi-scale search, separate models and region hypotheses.

    Args:
        directory: Regional result directory within an atomic output bundle.
        positive: Source foreground sequences.
        negative: Source control sequences.
        motif_directory: Native enrichment directory with known-motif sites.
        enabled: Explicit False disables this layer.
        use_ai: Follow the parent workflow's automatic ML setting.
        lengths: Regulatory search-window scales, not motif feature lengths.
        step: Sliding offset step.
        max_windows: Hard allocation limit across all scales.
        score_threshold: Fixed exploratory candidate display threshold.
        max_motifs: Maximum enriched motifs for descriptive density.
        motif_q_value: Source-level enrichment cutoff for density motifs.
        both_strands: Match native enrichment's strand convention.
        position_bin_width: Position profile bin width.
        position_limit: Absolute profile distance limit, only for summaries.
        max_locus_plots: Maximum individual locus figures, zero to 20.
        learning_settings: Existing word-model settings; permutations remain in
            the whole-flank model and are not repeated for exploratory scales.
        groups_path: Optional parent-level group TSV, taking precedence.
        parent_groups: Optional automatic parent-to-gene/contig group mapping.
        regions: Original retained genomic flanks for exact coordinate mapping.
        genes: Full annotation for signed gene-start distances.
        evidence_settings: Optional assay/reference inputs, re-queried locally.

    Returns:
        Availability, scale-specific model status and explicit uncertainty.

    Raises:
        ValueError: Settings, groups, evidence or input metadata are malformed.
    """
    directory.mkdir(parents=True, exist_ok=True)
    if not enabled:
        result: dict[str, Any] = {
            "status": "disabled",
            "reason": "Explicitly disabled by user",
        }
        write_json(path=directory / "summary.json", data=result)
        return result
    if (
        not math.isfinite(score_threshold)
        or not 0 <= score_threshold <= 1
        or not math.isfinite(motif_q_value)
        or not 0 < motif_q_value <= 1
        or not isinstance(max_motifs, int)
        or isinstance(max_motifs, bool)
        or max_motifs < 1
        or not isinstance(max_locus_plots, int)
        or isinstance(max_locus_plots, bool)
        or not 0 <= max_locus_plots <= 20
    ):
        raise ValueError("Invalid regulatory search threshold or plot limits")
    window_position_profiles(
        rows=[], bin_width=position_bin_width, limit=position_limit
    )
    windows, audit = make_windows(
        positive=positive,
        negative=negative,
        lengths=lengths,
        step=step,
        max_windows=max_windows,
        regions=regions,
        genes=genes,
    )
    if groups_path is not None:
        identifiers = [*positive, *negative]
        parent_groups = _read_parent_groups(groups_path, identifiers)
    groups = assign_window_groups(windows=windows, parent_groups=parent_groups)
    parents = {**positive, **negative}
    motif_ids, sites = _selected_sites(
        motif_directory, parents, max_motifs, motif_q_value, both_strands
    )
    densities = window_motif_counts(
        windows=windows, sites=sites, motif_ids=motif_ids
    )
    rows: list[dict[str, Any]] = []
    for window in windows:
        row = asdict(window)
        row.pop("sequence")
        row.update(
            sequence_id=window.sequence_id,
            window_length=window.width,
            held_out_signature_score=None,
            held_out_composition_score=None,
            fold=None,
            group=groups[window.parent_id],
            **sequence_composition(sequence=window.sequence),
            **densities[window.sequence_id],
        )
        rows.append(row)
    by_id = {r["sequence_id"]: r for r in rows}
    scale_summaries: list[dict[str, Any]] = []
    parent_predictions: list[dict[str, Any]] = []
    fold_metrics: list[dict[str, Any]] = []
    settings = dict(learning_settings or {})
    settings.pop("permutations", None)
    if (
        not isinstance(settings.get("folds", 5), int)
        or isinstance(settings.get("folds", 5), bool)
        or settings.get("folds", 5) < 2
    ):
        raise ValueError("At least two window validation folds are required")
    for width in sorted(lengths):
        subset = [w for w in windows if w.width == width]
        scale_directory = directory / f"scale_{width}"
        scale_directory.mkdir(parents=True, exist_ok=True)
        summary: dict[str, Any] = {
            "window_length": width,
            "windows": len(subset),
            "parents": len({w.parent_id for w in subset}),
            "status": "disabled" if not use_ai else "not_estimable",
            "reason": "ML explicitly disabled"
            if not use_ai
            else "No complete windows at this scale",
        }
        if use_ai and subset:
            try:
                from intergenic_regions.window_learning import fit_window_model

                (
                    predictions,
                    parent_rows,
                    metrics,
                    summary,
                    model,
                    explanations,
                ) = fit_window_model(
                    windows=subset, parent_groups=groups, **settings
                )
                for prediction in predictions:
                    by_id[prediction["sequence_id"]].update(prediction)
                parent_predictions.extend(parent_rows)
                fold_metrics.extend(metrics)
                write_json(path=scale_directory / "model.json", data=model)
                if explanations:
                    from intergenic_regions.shap_reporting import plot_shap

                    for name, records in (
                        ("shap_values", explanations["rows"]),
                        ("shap_importance", explanations["importance"]),
                    ):
                        write_tsv(
                            path=scale_directory / f"{name}.tsv",
                            rows=records,
                            fields=tuple(records[0]),
                        )
                    try:
                        images = plot_shap(
                            rows=explanations["rows"],
                            importance=explanations["importance"],
                            directory=scale_directory / "figures",
                        )
                        summary["shap"]["plot_status"] = "completed"
                    except ImportError:
                        images = []
                        summary["shap"]["plot_status"] = "unavailable"
                    write_report(
                        path=scale_directory / "report.html",
                        title=f"{width} bp search-window model",
                        summary={"ai": summary},
                        images=images,
                        tables={
                            "Held-out window SHAP importance": explanations[
                                "importance"
                            ]
                        },
                        notes=REGION_NOTES,
                        links={
                            "Window SHAP explanations TSV": "shap_values.tsv",
                            "Model JSON": "model.json",
                        },
                    )
            except ImportError:
                summary["status"] = "unavailable"
                summary["reason"] = (
                    "Install intergenic-regions[analysis] for window ML"
                )
            except ValueError as exc:
                allowed = (
                    "Window models require",
                    "At least one sequence per class per fold",
                    "Insufficient or malformed independent groups",
                    "A validation fold lacks a class",
                    "No valid k-mers for learning",
                )
                if not str(exc).startswith(allowed):
                    raise
                summary["status"] = "not_estimable"
                summary["reason"] = str(exc)
        LOGGER.info(
            "Regulatory window scale %d bp: %s", width, summary["status"]
        )
        write_json(path=scale_directory / "summary.json", data=summary)
        scale_summaries.append(summary)
    rows = _attach_context(directory, rows, evidence_settings or {})
    candidates = merge_window_candidates(
        rows=rows, score_threshold=score_threshold
    )
    candidate_directory = directory / "candidate_evidence"
    candidate_directory.mkdir(parents=True, exist_ok=True)
    candidates = _attach_context(
        candidate_directory, candidates, evidence_settings or {}
    )
    candidates = prioritise_candidates(rows=candidates)
    for candidate in candidates:
        candidate["uncertainty"] += (
            " Search-window union; functional boundaries and enhancer "
            "length are unestablished. Peak scores have no region-level FDR."
        )
    profiles = window_position_profiles(
        rows=rows, bin_width=position_bin_width, limit=position_limit
    )
    tables = {
        "window_availability": audit,
        "window_scores": rows,
        "parent_predictions": parent_predictions,
        "fold_metrics": fold_metrics,
        "candidate_regions": candidates,
        "position_profiles": profiles,
    }
    for name, records in tables.items():
        write_tsv(
            path=directory / f"{name}.tsv",
            rows=records,
            fields=tuple(records[0]) if records else ("sequence_id",),
        )
    write_fasta(
        path=directory / "candidate_regions.fasta",
        records={
            r["sequence_id"]: parents[r["parent_id"]][
                r["sequence_start"] : r["sequence_end"]
            ]
            for r in candidates
        },
    )
    for name, records in (
        ("windows", rows),
        ("candidate_regions", candidates),
    ):
        with (directory / f"{name}.bed").open(
            mode="w", encoding="utf-8"
        ) as stream:
            for row in records:
                if row["genomic_start"] is not None:
                    score = round(
                        1000 * (row.get("held_out_signature_score") or 0)
                    )
                    stream.write(
                        f"{row['contig']}\t{row['genomic_start']}\t{row['genomic_end']}\t{row['sequence_id']}\t{score}\t{row['strand']}\n"
                    )
    result = {
        "status": "completed" if windows else "no_complete_windows",
        "window_lengths": sorted(lengths),
        "step": step,
        "total_windows": len(windows),
        "parents_with_windows": len(groups),
        "independent_groups": len(set(groups.values())),
        "scales": scale_summaries,
        "candidate_regions": len(candidates),
        "score_threshold": score_threshold,
        "selected_motifs": motif_ids,
        "motif_selection_q_value": motif_q_value,
        "position_bin_width": position_bin_width,
        "position_limit": position_limit,
        "profiled_windows": sum(
            1
            for r in rows
            if -position_limit
            <= (
                r["distance_to_gene_start"]
                if r["distance_to_gene_start"] is not None
                else (r["sequence_start"] + r["sequence_end"] - 1) / 2
            )
            < position_limit
        ),
        "region_q_values": "not_assigned_exploratory_search",
        "interpretation": REGION_NOTES,
    }
    from intergenic_regions.window_reporting import plot_regional_results

    images = plot_regional_results(
        rows=rows,
        audit=audit,
        candidates=candidates,
        profiles=profiles,
        scales=scale_summaries,
        directory=directory / "figures",
        max_loci=max_locus_plots,
    )
    write_json(path=directory / "summary.json", data=result)
    write_report(
        path=directory / "report.html",
        title="Multi-scale regulatory region search",
        summary={"multiscale": result},
        images=images,
        notes=REGION_NOTES,
        tables={
            "Candidate window unions": candidates,
            "Scale availability": audit,
            "Window scores": rows,
        },
        links={
            "All window scores TSV": "window_scores.tsv",
            "Candidate regions TSV": "candidate_regions.tsv",
            "Candidate sequences FASTA": "candidate_regions.fasta",
            "Genomic candidate BED": "candidate_regions.bed",
            "Source-level validation TSV": "parent_predictions.tsv",
            "Position profiles TSV": "position_profiles.tsv",
            **{
                f"{s['window_length']} bp SHAP report": (
                    f"scale_{s['window_length']}/report.html"
                )
                for s in scale_summaries
                if (
                    directory / f"scale_{s['window_length']}" / "report.html"
                ).is_file()
            },
        },
    )
    return result
