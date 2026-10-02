"""Complete extraction, enrichment, AI and evidence-reporting workflows."""

import csv
import json
import logging
from collections import Counter
from collections.abc import Sequence
from dataclasses import asdict
from pathlib import Path
from typing import Any
from urllib.parse import quote

from intergenic_regions.annotation import read_annotation
from intergenic_regions.background import match_background, sanitise_regions
from intergenic_regions.evidence import (
    annotate_evidence,
    read_bed,
    read_functional_evidence,
)
from intergenic_regions.extraction import GeneIndex, extract_regions
from intergenic_regions.genome import Genome, sequence_composition
from intergenic_regions.io import (
    output_bundle,
    provenance,
    read_fasta,
    read_identifiers,
    write_fasta,
    write_json,
    write_tsv,
)
from intergenic_regions.models import Region
from intergenic_regions.prioritisation import prioritise_candidates
from intergenic_regions.references import query_references
from intergenic_regions.reporting import write_report

LOGGER = logging.getLogger(__name__)

REGION_FIELDS = (
    "sequence_id",
    "gene_id",
    "contig",
    "start",
    "end",
    "strand",
    "direction",
    "status",
    "stop_reason",
    "available_length",
    "length",
    "gc_fraction",
    "ambiguous_fraction",
)
AUDIT_FIELDS = ("sequence_id", "label", "reason", "representative")
EVIDENCE_FIELDS = (
    "sequence_id",
    "gene_id",
    "contig",
    "start",
    "end",
    "accessibility_overlap_bp",
    "accessibility_overlap_fraction",
    "enhancer_annotation_overlap_bp",
    "enhancer_annotation_overlap_fraction",
    "linked_gene_evidence_count",
    "linked_gene_evidence_json",
    "evidence_status",
)
REFERENCE_FIELDS = (
    "sequence_id",
    "gene_id",
    "reference_source",
    "reference_name",
    "reference_assembly",
    "reference_sha256",
    "evidence_type",
    "overlap_bp",
    "overlap_fraction",
)
MOTIF_FIELDS = (
    "motif_id",
    "name",
    "kind",
    "consensus",
    "positive_hits",
    "positive_total",
    "negative_hits",
    "negative_total",
    "positive_fraction",
    "negative_fraction",
    "fold_enrichment",
    "odds_ratio_haldane",
    "odds_ratio_ci_low",
    "odds_ratio_ci_high",
    "p_value",
    "log_p_value",
    "minus_log10_p",
    "q_value",
    "log_q_value",
    "minus_log10_q",
)
SITE_FIELDS = (
    "motif_id",
    "sequence_id",
    "label",
    "start",
    "end",
    "site_strand",
    "binned_score",
    "site_p_value",
)


def region_rows(*, regions: Sequence[Region]) -> list[dict[str, Any]]:
    """Convert regions to audited coordinate/composition records.

    Args:
        regions: Retained regions and exclusions.

    Returns:
        Records matching REGION_FIELDS; sequence bases are in FASTA only.
    """
    rows: list[dict[str, Any]] = []
    for region in regions:
        row = asdict(obj=region)
        row.pop("sequence")
        row["sequence_id"] = region.sequence_id
        row.update(sequence_composition(sequence=region.sequence))
        rows.append(row)
    return rows


def write_regions(
    *, directory: Path, regions: Sequence[Region]
) -> dict[str, Any]:
    """Write FASTA, BED, GFF3 and audit records with consistent coordinates.

    Args:
        directory: Existing destination directory.
        regions: Retained regions and exclusions.

    Returns:
        Counts by extraction status.
    """
    directory.mkdir(parents=True, exist_ok=True)
    retained = [r for r in regions if r.status == "retained"]
    write_fasta(
        path=directory / "regions.fasta",
        records={r.sequence_id: r.sequence for r in retained},
    )
    write_tsv(
        path=directory / "regions.tsv",
        rows=region_rows(regions=regions),
        fields=REGION_FIELDS,
    )
    with (directory / "regions.bed").open(mode="w", encoding="utf-8") as bed:
        with (directory / "regions.gff3").open(
            mode="w", encoding="utf-8"
        ) as gff:
            gff.write("##gff-version 3\n")
            for region in retained:
                bed.write(
                    f"{region.contig}\t{region.start}\t{region.end}\t"
                    f"{region.sequence_id}\t0\t{region.strand}\n"
                )
                attributes = (
                    f"ID={quote(region.sequence_id, safe='._-')};"
                    f"gene_id={quote(region.gene_id, safe='._-')};"
                    f"direction={region.direction}"
                )
                gff.write(
                    f"{region.contig}\tintergenic-regions\tintergenic_region\t"
                    f"{region.start + 1}\t{region.end}\t.\t"
                    f"{region.strand}\t.\t{attributes}\n"
                )
    summary = {
        "requested_regions": len(regions),
        "retained_regions": len(retained),
        "status_counts": dict(Counter(r.status for r in regions)),
        "coordinates": (
            "BED/TSV zero-based half-open; GFF3 one-based inclusive"
        ),
    }
    write_json(path=directory / "summary.json", data=summary)
    return summary


def read_region_table(*, path: Path) -> list[Region]:
    """Read a generated audit TSV for subsequent evidence annotation.

    Args:
        path: regions.tsv from an extraction workflow.

    Returns:
        Regions without sequence bases, sufficient for interval queries.

    Raises:
        ValueError: Required fields or intervals are malformed.
    """
    regions: list[Region] = []
    with path.open(mode="r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(f=stream, delimiter="\t")
        required = {
            "gene_id",
            "contig",
            "start",
            "end",
            "strand",
            "direction",
            "status",
        }
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Region TSV lacks required extraction columns")
        seen: set[str] = set()
        for row in reader:
            start, end = int(row["start"]), int(row["end"])
            if (
                start < 0
                or end < start
                or row["strand"] not in {"+", "-", "."}
                or row["direction"] not in {"upstream", "downstream"}
                or (row["status"] == "retained" and end == start)
            ):
                raise ValueError("Invalid interval or strand in region TSV")
            region = Region(
                gene_id=row["gene_id"],
                contig=row["contig"],
                start=start,
                end=end,
                strand=row["strand"],
                direction=row["direction"],
                sequence="",
                status=row["status"],
                stop_reason=row.get("stop_reason", ""),
                available_length=int(row.get("available_length", end - start)),
            )
            if region.sequence_id in seen:
                raise ValueError("Duplicate region in TSV")
            seen.add(region.sequence_id)
            regions.append(region)
    return regions


def attach_evidence(
    *,
    directory: Path,
    regions: Sequence[Region],
    accessibility_bed: Path | None = None,
    enhancer_bed: Path | None = None,
    evidence_tsv: Path | None = None,
    references: Sequence[Path] = (),
    organism: str | None = None,
    assembly: str | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Write optional assay, functional and database evidence tables.

    Args:
        directory: Output directory.
        regions: Extracted regions.
        accessibility_bed: Optional accessibility peaks.
        enhancer_bed: Optional enhancer annotations.
        evidence_tsv: Optional linked gene-level experimental results.
        references: Imported regulatory-reference directories.
        organism: Required query organism when references are used.
        assembly: Required query assembly when references are used.

    Returns:
        Assay/functional evidence and reference overlap rows.

    Raises:
        ValueError: Contig names or reference organism/build do not match.
    """
    accessibility = (
        read_bed(path=accessibility_bed) if accessibility_bed else None
    )
    enhancers = read_bed(path=enhancer_bed) if enhancer_bed else None
    contigs = {r.contig for r in regions if r.status == "retained"}
    for index in (accessibility, enhancers):
        if (
            index is not None
            and index.intervals
            and contigs
            and not contigs & index.intervals.keys()
        ):
            raise ValueError(
                "Evidence and regions have no shared contig identifiers"
            )
    functional = (
        read_functional_evidence(path=evidence_tsv) if evidence_tsv else None
    )
    evidence = annotate_evidence(
        regions=regions,
        accessibility=accessibility,
        enhancers=enhancers,
        functional=functional,
    )
    reference_rows: list[dict[str, Any]] = []
    if references:
        if not organism or not assembly:
            raise ValueError(
                "Reference queries require --organism and --assembly"
            )
        reference_rows = query_references(
            regions=regions,
            reference_directories=references,
            organism=organism,
            assembly=assembly,
        )
    write_tsv(
        path=directory / "evidence.tsv", rows=evidence, fields=EVIDENCE_FIELDS
    )
    write_tsv(
        path=directory / "reference_overlaps.tsv",
        rows=reference_rows,
        fields=REFERENCE_FIELDS,
    )
    return evidence, reference_rows


def extract_workflow(
    *,
    genome_path: Path,
    annotation_path: Path,
    output: Path,
    identifiers_path: Path | None = None,
    annotation_format: str = "auto",
    extraction_settings: dict[str, Any] | None = None,
    evidence_settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run complete indexed extraction into an atomic output bundle.

    Args:
        genome_path: Genome FASTA.
        annotation_path: Full gene annotation.
        output: New result directory.
        identifiers_path: Optional exact gene list.
        annotation_format: Annotation format selection.
        extraction_settings: Keyword arguments for extract_regions.
        evidence_settings: Optional keyword arguments for attach_evidence.

    Returns:
        Extraction summary.
    """
    settings = extraction_settings or {}
    genes = read_annotation(
        path=annotation_path, annotation_format=annotation_format
    )
    identifiers = (
        read_identifiers(path=identifiers_path) if identifiers_path else None
    )
    with Genome(path=genome_path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        regions = extract_regions(
            index=index, genome=genome, identifiers=identifiers, **settings
        )
    with output_bundle(path=output) as stage:
        summary = write_regions(directory=stage, regions=regions)
        evidence, references = attach_evidence(
            directory=stage, regions=regions, **(evidence_settings or {})
        )
        inputs = [genome_path, annotation_path] + (
            [identifiers_path] if identifiers_path else []
        )
        write_json(
            path=stage / "manifest.json",
            data=provenance(
                inputs=inputs
                + evidence_input_paths(settings=evidence_settings or {}),
                settings={
                    "annotation_format": annotation_format,
                    **settings,
                    **serialise_settings(settings=evidence_settings or {}),
                },
            ),
        )
        write_report(
            path=stage / "report.html",
            title="Intergenic extraction",
            summary=summary,
            tables={
                "Region audit": region_rows(regions=regions),
                "Optional evidence": evidence,
                "Reference overlaps": references,
            },
            notes=[
                "All blockers come from the full annotation. Sequences stop "
                "before any annotated gene, regardless of strand. Missing "
                "annotations cannot be detected from the genome sequence "
                "alone."
            ],
        )
    return summary


def serialise_settings(*, settings: dict[str, Any]) -> dict[str, Any]:
    """Convert workflow paths and sequences of paths to JSON settings.

    Args:
        settings: Workflow keyword arguments.

    Returns:
        JSON-compatible settings preserving explicit paths.
    """
    return {
        key: str(value.resolve())
        if isinstance(value, Path)
        else [
            str(item.resolve()) if isinstance(item, Path) else item
            for item in value
        ]
        if isinstance(value, (tuple, list))
        else value
        for key, value in settings.items()
    }


def evidence_input_paths(*, settings: dict[str, Any]) -> list[Path]:
    """Collect evidence files for reproducibility hashing.

    Args:
        settings: Evidence workflow arguments.

    Returns:
        Assay/functional paths and regulatory bundle metadata files.
    """
    paths = [
        settings[key]
        for key in ("accessibility_bed", "enhancer_bed", "evidence_tsv")
        if settings.get(key)
    ]
    paths += [
        directory / "reference.json"
        for directory in settings.get("references", [])
    ]
    return paths


def enrichment_outputs(
    *,
    directory: Path,
    positive: dict[str, str],
    negative: dict[str, str],
    motif_path: Path | None = None,
    motif_format: str = "auto",
    lengths: Sequence[int] = (),
    both_strands: bool = True,
    site_p_value: float = 1e-4,
) -> dict[str, Any]:
    """Write native enrichment tables, motif sites, graphics and HTML.

    Args:
        directory: Destination directory.
        positive: Foreground FASTA records.
        negative: Negative FASTA records.
        motif_path: Optional known-motif file.
        motif_format: Known-motif format.
        lengths: Exact k-mer lengths; defaults to six when no motifs supplied.
        both_strands: Scan/combine both DNA strands.
        site_p_value: PWM per-site threshold.

    Returns:
        Native analysis summary.
    """
    from intergenic_regions.motifs import analyse_motifs, read_motifs
    from intergenic_regions.reporting import (
        plot_composition,
        plot_logos,
        plot_motif_results,
    )

    directory.mkdir(parents=True, exist_ok=True)
    motifs = (
        read_motifs(path=motif_path, motif_format=motif_format)
        if motif_path
        else []
    )
    selected_lengths = lengths or (() if motifs else (6,))
    rows, sites, summary = analyse_motifs(
        positive=positive,
        negative=negative,
        motifs=motifs,
        kmer_lengths=selected_lengths,
        both_strands=both_strands,
        site_p_value=site_p_value,
    )
    write_tsv(
        path=directory / "motif_enrichment.tsv", rows=rows, fields=MOTIF_FIELDS
    )
    write_tsv(
        path=directory / "motif_sites.tsv", rows=sites, fields=SITE_FIELDS
    )
    summary["kmer_lengths"] = list(selected_lengths)
    write_json(path=directory / "summary.json", data=summary)
    images = [
        plot_composition(
            positive=positive,
            negative=negative,
            directory=directory / "figures",
        )
    ]
    images.extend(
        plot_motif_results(rows=rows, directory=directory / "figures")
    )
    logo = plot_logos(
        rows=rows, motifs=motifs, directory=directory / "figures"
    )
    if logo:
        images.append(logo)
    write_report(
        path=directory / "report.html",
        title="Motif enrichment",
        summary=summary,
        tables={"Motif enrichment": rows},
        images=images,
        notes=[
            "The unit of testing is sequence presence, not the number of "
            "overlapping sites. BH correction includes all possible tested "
            "words, including unobserved k-mers.",
            "PWM site p-values use a discretised zero-order background; "
            "enrichment q-values apply to sequence-level tests. Exact "
            "IUPAC patterns are not assigned site p-values.",
            "Check background composition and biological independence. "
            "Motif enrichment supports a regulatory hypothesis; it does "
            "not establish enhancer function.",
        ],
    )
    return summary


def enrichment_workflow(
    *,
    positive_path: Path,
    negative_path: Path,
    output: Path,
    motif_path: Path | None = None,
    settings: dict[str, Any] | None = None,
    mask_lowercase: bool = False,
    use_ai: bool = True,
    learning_settings: dict[str, Any] | None = None,
    groups_path: Path | None = None,
    use_regions: bool = True,
    regional_settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run native enrichment on supplied sequence sets.

    Args:
        positive_path: Foreground FASTA.
        negative_path: Negative FASTA.
        output: New result directory.
        motif_path: Optional known motifs.
        settings: Keyword arguments for enrichment_outputs.
        mask_lowercase: Mask lowercase DNA before scanning.
        use_ai: Run ML by default, independently of motif significance.
        learning_settings: Fixed model and validation settings.
        groups_path: Optional homology/chromosome group assignments.
        use_regions: Run multi-scale regulatory search by default.
        regional_settings: Candidate window and position-profile settings.

    Returns:
        Analysis summary.
    """
    positive = read_fasta(path=positive_path, mask_lowercase=mask_lowercase)
    negative = read_fasta(path=negative_path, mask_lowercase=mask_lowercase)
    with output_bundle(path=output) as stage:
        summary = enrichment_outputs(
            directory=stage,
            positive=positive,
            negative=negative,
            motif_path=motif_path,
            **(settings or {}),
        )
        ai_summary, predictions = automatic_learning(
            directory=stage / "ai",
            positive=positive,
            negative=negative,
            enabled=use_ai,
            settings=learning_settings,
            groups_path=groups_path,
        )
        prediction_by_id = {row["sequence_id"]: row for row in predictions}
        candidate_rows = prioritise_candidates(
            rows=[
                {
                    "sequence_id": identifier,
                    "label": label,
                    **sequence_composition(sequence=sequence),
                    "held_out_signature_score": prediction_by_id.get(
                        identifier, {}
                    ).get("held_out_signature_score"),
                }
                for label, records in (
                    ("positive", positive),
                    ("negative", negative),
                )
                for identifier, sequence in records.items()
            ]
        )
        write_tsv(
            path=stage / "candidate_priorities.tsv",
            rows=candidate_rows,
            fields=tuple(candidate_rows[0])
            if candidate_rows
            else ("sequence_id", "priority_rank", "enhancer_status"),
        )
        summary["ai"] = ai_summary
        from intergenic_regions.regional_analysis import regional_outputs

        summary["multiscale"] = regional_outputs(
            directory=stage / "regions",
            positive=positive,
            negative=negative,
            motif_directory=stage,
            enabled=use_regions,
            use_ai=use_ai,
            learning_settings=learning_settings,
            groups_path=groups_path,
            both_strands=(settings or {}).get("both_strands", True),
            **(regional_settings or {}),
        )
        write_json(path=stage / "summary.json", data=summary)
        analysis_report(
            directory=stage,
            summary=summary,
            candidates=candidate_rows,
            motif_directory=stage,
        )
        inputs = (
            [positive_path, negative_path]
            + ([motif_path] if motif_path else [])
            + ([groups_path] if groups_path else [])
        )
        write_json(
            path=stage / "manifest.json",
            data=provenance(
                inputs=inputs,
                settings={
                    "mask_lowercase": mask_lowercase,
                    "use_ai": use_ai,
                    "learning": learning_settings or {},
                    "use_regions": use_regions,
                    "regional_search": regional_settings or {},
                    **(settings or {}),
                },
            ),
        )
    return summary


def automatic_learning(
    *,
    directory: Path,
    positive: dict[str, str],
    negative: dict[str, str],
    enabled: bool = True,
    settings: dict[str, Any] | None = None,
    groups_path: Path | None = None,
    groups: Sequence[str] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Run ML and report unsupported data without inventing scores.

    Args:
        directory: Model result directory.
        positive: Foreground sequences.
        negative: Negative sequences.
        enabled: False only when ML is explicitly disabled.
        settings: Model settings passed unchanged to learning_outputs.
        groups_path: Optional supplied independent-group TSV.
        groups: Optional automatic gene or contig grouping.

    Returns:
        Explicit model status and held-out predictions when estimable.

    Raises:
        ValueError: User settings or supplied grouping files are invalid.
    """
    if settings is not None and settings.get("folds", 5) < 2:
        raise ValueError("At least two validation folds are required")
    if not enabled:
        result = {
            "status": "disabled",
            "reason": "Explicitly disabled by user",
        }
    else:
        try:
            summary, predictions = learning_outputs(
                directory=directory,
                positive=positive,
                negative=negative,
                settings=settings,
                groups_path=groups_path,
                groups=groups,
            )
            result = {"status": "completed", **summary}
            write_json(path=directory / "summary.json", data=result)
            return result, predictions
        except ImportError:
            result = {
                "status": "unavailable",
                "reason": "Install intergenic-regions[analysis] to run ML",
            }
        except ValueError as exc:
            allowed = (
                "AI analysis requires at least five",
                "At least one sequence per class per fold",
                "Insufficient or malformed independent groups",
                "A validation fold lacks a class",
                "No valid k-mers for learning",
                "Permutation validation failed",
            )
            if not str(exc).startswith(allowed):
                raise
            result = {"status": "not_estimable", "reason": str(exc)}
    LOGGER.warning("ML status %s: %s", result["status"], result["reason"])
    directory.mkdir(parents=True, exist_ok=True)
    write_json(path=directory / "summary.json", data=result)
    write_report(
        path=directory / "report.html",
        title="Sequence model status",
        summary=result,
        notes=[result["reason"]],
    )
    return result, []


def analysis_report(
    *,
    directory: Path,
    summary: dict[str, Any],
    candidates: Sequence[dict[str, Any]],
    motif_directory: Path,
    extra_tables: dict[str, Sequence[dict[str, Any]]] | None = None,
) -> None:
    """Write a single offline dashboard combining motifs, ML and evidence.

    Args:
        directory: Root output directory.
        summary: Integrated analysis metrics and model status.
        candidates: Transparent ranked candidate records.
        motif_directory: Native motif results directory.
        extra_tables: Additional extraction/evidence audit tables.
    """
    with (motif_directory / "motif_enrichment.tsv").open(
        mode="r", encoding="utf-8", newline=""
    ) as stream:
        motifs = list(csv.DictReader(f=stream, delimiter="\t"))
    metrics: list[dict[str, Any]] = []
    fold_path = directory / "ai" / "fold_metrics.tsv"
    if fold_path.is_file():
        with fold_path.open(mode="r", encoding="utf-8", newline="") as stream:
            metrics = list(csv.DictReader(f=stream, delimiter="\t"))
    images = sorted((motif_directory / "figures").glob("*.png"))
    images += sorted((directory / "ai" / "figures").glob("*.png"))
    images += sorted((directory / "genome_scan" / "figures").glob("*.png"))
    images += sorted((directory / "regions" / "figures").glob("*.png"))
    shap_rows: list[dict[str, Any]] = []
    shap_path = directory / "ai" / "shap_importance.tsv"
    if shap_path.is_file():
        with shap_path.open(mode="r", encoding="utf-8", newline="") as stream:
            shap_rows = list(csv.DictReader(f=stream, delimiter="\t"))
    links = {
        "Motif enrichment TSV": str(
            (motif_directory / "motif_enrichment.tsv").relative_to(directory)
        ),
        "Ranked candidates TSV": "candidate_priorities.tsv",
        "Model report": "ai/report.html",
        "Run manifest": "manifest.json",
    }
    if shap_rows:
        links["SHAP feature importance TSV"] = "ai/shap_importance.tsv"
        links["SHAP sequence explanations TSV"] = "ai/shap_values.tsv"
    if (directory / "genome_scan" / "report.html").is_file():
        links["Whole-genome scan report"] = "genome_scan/report.html"
        links["Genome motif sites TSV"] = "genome_scan/genome_motif_sites.tsv"
        links["Gene-start distance profiles TSV"] = (
            "genome_scan/distance_profiles.tsv"
        )
    region_tables: dict[str, list[dict[str, Any]]] = {}
    if (directory / "regions" / "report.html").is_file():
        links["Multi-scale regulatory region report"] = "regions/report.html"
        links["Multi-scale window scores TSV"] = "regions/window_scores.tsv"
        links["Candidate window unions TSV"] = "regions/candidate_regions.tsv"
        links["Candidate window unions BED"] = "regions/candidate_regions.bed"
        for name, filename in (
            ("Multi-scale candidate window unions", "candidate_regions.tsv"),
            ("Multi-scale availability audit", "window_availability.tsv"),
        ):
            with (directory / "regions" / filename).open(
                encoding="utf-8"
            ) as stream:
                region_tables[name] = list(
                    csv.DictReader(stream, delimiter="\t")
                )
    write_report(
        path=directory / "report.html",
        title="Intergenic regulatory sequence analysis",
        summary=summary,
        tables={
            "Ranked candidate regions": [
                {
                    key: row.get(key)
                    for key in (
                        "priority_rank",
                        "sequence_id",
                        "label",
                        "held_out_signature_score",
                        "support_category",
                        "accessibility_overlap_fraction",
                        "reference_support_count",
                        "enhancer_status",
                    )
                }
                for row in candidates
            ],
            "Motif enrichment": [
                {
                    key: row.get(key)
                    for key in (
                        "motif_id",
                        "consensus",
                        "q_value",
                        "p_value",
                        "positive_fraction",
                        "negative_fraction",
                        "fold_enrichment",
                        "kind",
                    )
                }
                for row in motifs
            ],
            "Held-out model validation": metrics,
            "Held-out SHAP feature importance": shap_rows,
            **region_tables,
            **(extra_tables or {}),
        },
        images=images,
        links=links,
        notes=[
            "Motifs and ML use separate evidence: enrichment tests sequence "
            "presence; held-out ML measures positive-class resemblance.",
            "Priorities sort by supplied contextual evidence, then held-out "
            "signature score. This ordering is a transparent triage rule, "
            "not a calibrated probability of enhancer activity.",
            "Sequence-only candidates remain unvalidated. Accessibility and "
            "annotation overlaps support context, while gene-level assays "
            "do not establish the function of this extracted interval.",
            "FDR includes the full tested motif family. Inspect GC/length "
            "balance, independent groups and the composition-only baseline.",
            "SHAP explains held-out model log-odds using each fold's training "
            "background. Correlated words can share or redistribute "
            "attribution; these explanations are not causal evidence.",
            "Multi-scale windows are searched automatically within safe "
            "flanks. Their labels are inherited from parents. Window unions "
            "are exploratory regulatory hypotheses, with unestablished "
            "enhancer boundaries and no region-level FDR. Inspect "
            "availability "
            "and source-level validation at every scale.",
        ],
    )


def learning_outputs(
    *,
    directory: Path,
    positive: dict[str, str],
    negative: dict[str, str],
    settings: dict[str, Any] | None = None,
    groups_path: Path | None = None,
    groups: Sequence[str] | None = None,
    candidates_path: Path | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Write held-out AI analysis, coefficients, reusable model and graphics.

    Args:
        directory: Destination directory.
        positive: Foreground FASTA records.
        negative: Negative FASTA records.
        settings: Keyword arguments for fit_sequence_model.
        groups_path: Optional related-sequence group TSV.
        groups: Optional automatically derived groups.
        candidates_path: Optional new unlabelled candidate FASTA.

    Returns:
        Summary and held-out sequence prediction rows.
    """
    from intergenic_regions.learning import (
        fit_sequence_model,
        predict_sequences,
        read_groups,
    )
    from intergenic_regions.reporting import plot_learning

    directory.mkdir(parents=True, exist_ok=True)
    selected_groups = (
        read_groups(path=groups_path, identifiers=[*positive, *negative])
        if groups_path
        else groups
    )
    explanations: dict[str, Any] = {}
    predictions, signatures, metrics, summary, model = fit_sequence_model(
        positive=positive,
        negative=negative,
        groups=selected_groups,
        explanation_outputs=explanations,
        **(settings or {}),
    )
    write_tsv(
        path=directory / "held_out_predictions.tsv",
        rows=predictions,
        fields=tuple(predictions[0]),
    )
    write_tsv(
        path=directory / "sequence_signatures.tsv",
        rows=signatures,
        fields=tuple(signatures[0]),
    )
    write_tsv(
        path=directory / "fold_metrics.tsv",
        rows=metrics,
        fields=tuple(metrics[0]),
    )
    write_json(path=directory / "model.json", data=model)
    candidates: list[dict[str, Any]] = []
    if candidates_path:
        candidates = predict_sequences(
            model=model, sequences=read_fasta(path=candidates_path)
        )
        write_tsv(
            path=directory / "candidate_scores.tsv",
            rows=candidates,
            fields=tuple(candidates[0]),
        )
    images = plot_learning(
        predictions=predictions,
        signatures=signatures,
        summary=summary,
        directory=directory / "figures",
    )
    importance = explanations.get("importance", [])
    if explanations:
        from intergenic_regions.shap_reporting import plot_shap

        write_tsv(
            path=directory / "shap_values.tsv",
            rows=explanations["rows"],
            fields=tuple(explanations["rows"][0]),
        )
        write_tsv(
            path=directory / "shap_importance.tsv",
            rows=importance,
            fields=tuple(importance[0]),
        )
        try:
            images.extend(
                plot_shap(
                    rows=explanations["rows"],
                    importance=importance,
                    directory=directory / "figures",
                )
            )
            summary["shap"]["plot_status"] = "completed"
        except ImportError:
            summary["shap"]["plot_status"] = "unavailable"
            summary["shap"]["plot_reason"] = (
                "Install intergenic-regions[analysis] for official SHAP plots"
            )
            LOGGER.warning(summary["shap"]["plot_reason"])
    write_json(path=directory / "summary.json", data=summary)
    write_report(
        path=directory / "report.html",
        title="AI sequence signatures",
        summary=summary,
        tables={
            "Held-out predictions": predictions,
            "Exploratory sequence signatures": signatures,
            "New candidate scores": candidates,
            "Held-out SHAP feature importance": importance,
        },
        images=images,
        links={
            "Held-out predictions TSV": "held_out_predictions.tsv",
            "Reusable model JSON": "model.json",
            **(
                {
                    "SHAP feature importance TSV": "shap_importance.tsv",
                    "SHAP sequence explanations TSV": "shap_values.tsv",
                }
                if explanations
                else {}
            ),
        },
        notes=[
            "Vocabulary selection, composition scaling and classifier "
            "fitting occur separately inside each validation fold. "
            "Reported validation scores never come from the model "
            "fitted to the full training set.",
            "The GC/length/ambiguity baseline helps identify "
            "composition-driven separation. Model coefficients are "
            "exploratory signatures, with fold stability rather than "
            "per-feature significance p-values.",
            "Sequence scores measure resemblance to the supplied positive "
            "class. Chromatin overlap and user-supplied functional results "
            "are separate evidence, and do not change training labels.",
            "Interventional SHAP contributions are additive positive-class "
            "log-odds, not probability changes. Each held-out fold uses its "
            "own training-only background. The displayed 'Other features' "
            "term preserves all omitted contributions; the global bar uses "
            "every model feature. Feature correlations affect attribution "
            "and no explanation establishes enhancer function.",
        ],
    )
    return summary, predictions


def learning_workflow(
    *,
    positive_path: Path,
    negative_path: Path,
    output: Path,
    settings: dict[str, Any] | None = None,
    groups_path: Path | None = None,
    candidates_path: Path | None = None,
) -> dict[str, Any]:
    """Run a complete standalone AI sequence workflow.

    Args:
        positive_path: Positive FASTA.
        negative_path: Negative FASTA.
        output: New result directory.
        settings: Sequence-model settings.
        groups_path: Optional related-sequence groups.
        candidates_path: Optional new candidates.

    Returns:
        Held-out predictive summary.
    """
    with output_bundle(path=output) as stage:
        summary, _ = learning_outputs(
            directory=stage,
            positive=read_fasta(path=positive_path),
            negative=read_fasta(path=negative_path),
            settings=settings,
            groups_path=groups_path,
            candidates_path=candidates_path,
        )
        inputs = [positive_path, negative_path] + [
            path for path in (groups_path, candidates_path) if path is not None
        ]
        write_json(
            path=stage / "manifest.json",
            data=provenance(inputs=inputs, settings=settings or {}),
        )
    return summary


def pipeline_workflow(
    *,
    genome_path: Path,
    annotation_path: Path,
    positive_genes: Path,
    negative_genes: Path,
    output: Path,
    annotation_format: str = "auto",
    extraction_settings: dict[str, Any] | None = None,
    motif_path: Path | None = None,
    enrichment_settings: dict[str, Any] | None = None,
    use_ai: bool = True,
    learning_settings: dict[str, Any] | None = None,
    groups_path: Path | None = None,
    group_by: str = "gene",
    matching_settings: dict[str, Any] | None = None,
    evidence_settings: dict[str, Any] | None = None,
    scan_genome: bool = False,
    scan_settings: dict[str, Any] | None = None,
    use_regions: bool = True,
    regional_settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run gene-list extraction, leakage filtering, enrichment and optional AI.

    Args:
        genome_path: Genome FASTA.
        annotation_path: Full gene annotation.
        positive_genes: Exact foreground gene IDs.
        negative_genes: Exact negative/candidate-control IDs.
        output: New result directory.
        annotation_format: Annotation format.
        extraction_settings: Flank settings.
        motif_path: Optional known-motif database.
        enrichment_settings: Native motif settings.
        use_ai: Enable supervised signature learning.
        learning_settings: Fixed AI settings.
        groups_path: Optional explicit independent-sequence groups.
        group_by: Default AI groups: ``gene`` or ``contig``.
        matching_settings: Optional GC/length control-matching settings.
        evidence_settings: Optional experimental/reference evidence.
        scan_genome: Screen positively enriched consensuses across the genome.
        scan_settings: Screen selection, substitution and profile limits.
        use_regions: Run multi-scale candidate search by default.
        regional_settings: Candidate window and position-profile settings.

    Returns:
        Integrated workflow summary.

    Raises:
        ValueError: Lists overlap or no independent data remain.
    """
    if group_by not in {"gene", "contig"}:
        raise ValueError("AI group_by must be gene or contig")
    positives = read_identifiers(path=positive_genes)
    negatives = read_identifiers(path=negative_genes)
    if set(positives) & set(negatives):
        raise ValueError("Positive and negative gene lists overlap")
    genes = read_annotation(
        path=annotation_path, annotation_format=annotation_format
    )
    with Genome(path=genome_path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        positive_regions = extract_regions(
            index=index,
            genome=genome,
            identifiers=positives,
            **(extraction_settings or {}),
        )
        negative_regions = extract_regions(
            index=index,
            genome=genome,
            identifiers=negatives,
            **(extraction_settings or {}),
        )
    positive, negative, audit = sanitise_regions(
        positive=positive_regions, negative=negative_regions
    )
    matching: list[dict[str, Any]] = []
    if matching_settings is not None:
        negative, matching = match_background(
            positive=positive, negative=negative, **matching_settings
        )
    positive_records = {r.sequence_id: r.sequence for r in positive}
    negative_records = {r.sequence_id: r.sequence for r in negative}
    with output_bundle(path=output) as stage:
        write_regions(
            directory=stage / "positive_extraction", regions=positive_regions
        )
        write_regions(
            directory=stage / "negative_extraction", regions=negative_regions
        )
        write_fasta(path=stage / "positive.fasta", records=positive_records)
        write_fasta(path=stage / "negative.fasta", records=negative_records)
        write_tsv(
            path=stage / "background_audit.tsv",
            rows=audit,
            fields=AUDIT_FIELDS,
        )
        write_tsv(
            path=stage / "background_matching.tsv",
            rows=matching,
            fields=(
                "positive_id",
                "negative_id",
                "gc_difference",
                "length_ratio",
            ),
        )
        retained = [*positive, *negative]
        evidence, references = attach_evidence(
            directory=stage, regions=retained, **(evidence_settings or {})
        )
        native = enrichment_outputs(
            directory=stage / "motifs",
            positive=positive_records,
            negative=negative_records,
            motif_path=motif_path,
            **(enrichment_settings or {}),
        )
        groups = [
            r.gene_id if group_by == "gene" else r.contig for r in retained
        ]
        ai_summary, predictions = automatic_learning(
            directory=stage / "ai",
            positive=positive_records,
            negative=negative_records,
            enabled=use_ai,
            settings=learning_settings,
            groups_path=groups_path,
            groups=groups,
        )
        from intergenic_regions.regional_analysis import regional_outputs

        regional_summary = regional_outputs(
            directory=stage / "regions",
            positive=positive_records,
            negative=negative_records,
            motif_directory=stage / "motifs",
            enabled=use_regions,
            use_ai=use_ai,
            learning_settings=learning_settings,
            groups_path=groups_path,
            parent_groups=dict(
                zip(
                    [*positive_records, *negative_records], groups, strict=True
                )
            ),
            regions=retained,
            genes=genes,
            evidence_settings=evidence_settings,
            both_strands=(enrichment_settings or {}).get("both_strands", True),
            **(regional_settings or {}),
        )
        evidence_by_id = {r["sequence_id"]: r for r in evidence}
        prediction_by_id = {r["sequence_id"]: r for r in predictions}
        candidate_rows: list[dict[str, Any]] = []
        for region in retained:
            overlap_sources = sorted(
                {
                    r["reference_name"]
                    for r in references
                    if r["sequence_id"] == region.sequence_id
                    and r["overlap_bp"] > 0
                }
            )
            candidate_rows.append(
                {
                    **evidence_by_id[region.sequence_id],
                    "label": "positive"
                    if region.sequence_id in positive_records
                    else "negative",
                    "held_out_signature_score": prediction_by_id.get(
                        region.sequence_id, {}
                    ).get("held_out_signature_score"),
                    "reference_support_count": len(overlap_sources),
                    "overlapping_reference_names": json.dumps(
                        obj=overlap_sources
                    ),
                }
            )
        candidate_rows = prioritise_candidates(rows=candidate_rows)
        write_tsv(
            path=stage / "candidate_summary.tsv",
            rows=candidate_rows,
            fields=tuple(candidate_rows[0]),
        )
        write_tsv(
            path=stage / "candidate_priorities.tsv",
            rows=candidate_rows,
            fields=tuple(candidate_rows[0]),
        )
        summary: dict[str, Any] = {
            "positive_regions": len(positive),
            "negative_regions": len(negative),
            "leakage_exclusions": len(audit),
            "matched_background": matching_settings is not None,
            "motifs": native,
            "ai": ai_summary,
            "multiscale": regional_summary,
            "optional_evidence_supplied": any(
                (evidence_settings or {}).get(k)
                for k in (
                    "accessibility_bed",
                    "enhancer_bed",
                    "evidence_tsv",
                    "references",
                )
            ),
        }
        if scan_genome:
            from intergenic_regions.scanning import genome_scan_outputs

            targets, screen_settings = scan_targets(
                enrichment_path=stage / "motifs" / "motif_enrichment.tsv",
                settings=scan_settings or {},
            )
            summary["scan"] = genome_scan_outputs(
                directory=stage / "genome_scan",
                genome_path=genome_path,
                annotation_path=annotation_path,
                annotation_format=annotation_format,
                positive_genes=positive_genes,
                negative_genes=negative_genes,
                motifs=targets,
                selection_settings=serialise_settings(
                    settings=scan_settings or {}
                ),
                **screen_settings,
            )
        write_json(path=stage / "summary.json", data=summary)
        inputs = (
            [genome_path, annotation_path, positive_genes, negative_genes]
            + [path for path in (motif_path, groups_path) if path is not None]
            + evidence_input_paths(settings=evidence_settings or {})
        )
        settings = {
            "annotation_format": annotation_format,
            "extraction": extraction_settings or {},
            "enrichment": enrichment_settings or {},
            "learning": learning_settings or {},
            "use_ai": use_ai,
            "group_by": group_by,
            "matching": matching_settings,
            "evidence": serialise_settings(settings=evidence_settings or {}),
            "scan_genome": scan_genome,
            "scan": scan_settings or {},
            "use_regions": use_regions,
            "regional_search": regional_settings or {},
        }
        write_json(
            path=stage / "manifest.json",
            data=provenance(inputs=inputs, settings=settings),
        )
        analysis_report(
            directory=stage,
            summary=summary,
            candidates=candidate_rows,
            motif_directory=stage / "motifs",
            extra_tables={
                "Leakage exclusions": audit,
                "Reference overlaps": references,
            },
        )
    return summary


def scan_targets(
    *,
    enrichment_path: Path | None = None,
    motif_path: Path | None = None,
    settings: dict[str, Any] | None = None,
) -> tuple[list[Any], dict[str, Any]]:
    """Separate motif selection from genome-screen settings.

    Args:
        enrichment_path: Generated positive-versus-negative enrichment TSV.
        motif_path: Alternatively, a known-motif file.
        settings: Selection and screening keyword arguments.

    Returns:
        Validated consensus targets and remaining screening settings.
    """
    from intergenic_regions.scanning import read_scan_motifs

    remaining = dict(settings or {})
    selection = {
        key: remaining.pop(key)
        for key in ("motif_format", "q_threshold", "max_motifs", "motif_ids")
        if key in remaining
    }
    targets = read_scan_motifs(
        enrichment_path=enrichment_path,
        motif_path=motif_path,
        **selection,
    )
    return targets, remaining


def genome_scan_workflow(
    *,
    genome_path: Path,
    output: Path,
    enrichment_path: Path | None = None,
    motif_path: Path | None = None,
    annotation_path: Path | None = None,
    annotation_format: str = "auto",
    positive_genes: Path | None = None,
    negative_genes: Path | None = None,
    settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Publish an atomic genome-wide consensus screen and offline dashboard.

    Args:
        genome_path: Complete genome FASTA.
        output: New output directory.
        enrichment_path: Motif-enrichment TSV for positive-enriched selection.
        motif_path: Alternatively, an IUPAC, MEME or JASPAR motif file.
        annotation_path: Optional full annotation for distances and overlaps.
        annotation_format: Annotation input format.
        positive_genes: Optional discovery foreground gene identifiers.
        negative_genes: Optional discovery control gene identifiers.
        settings: Motif-selection and genome-screen options.

    Returns:
        Scan summary. Matches are unvalidated sequence hypotheses.
    """
    from intergenic_regions.scanning import genome_scan_outputs

    targets, screen_settings = scan_targets(
        enrichment_path=enrichment_path,
        motif_path=motif_path,
        settings=settings,
    )
    with output_bundle(path=output) as stage:
        summary = genome_scan_outputs(
            directory=stage,
            genome_path=genome_path,
            motifs=targets,
            annotation_path=annotation_path,
            annotation_format=annotation_format,
            positive_genes=positive_genes,
            negative_genes=negative_genes,
            selection_settings=serialise_settings(settings=settings or {}),
            **screen_settings,
        )
        inputs = [genome_path] + [
            p
            for p in (
                enrichment_path,
                motif_path,
                annotation_path,
                positive_genes,
                negative_genes,
            )
            if p is not None
        ]
        write_json(
            path=stage / "manifest.json",
            data=provenance(
                inputs=inputs,
                settings={
                    "annotation_format": annotation_format,
                    **serialise_settings(settings=settings or {}),
                },
            ),
        )
    return summary
