"""Optional RNA-start, mapping, differential-expression and motif workflows."""

import csv
import json
import logging
import math
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

from intergenic_regions.annotation import read_annotation
from intergenic_regions.extraction import GeneIndex
from intergenic_regions.genome import Genome
from intergenic_regions.io import (
    output_bundle,
    provenance,
    read_identifiers,
    write_json,
    write_tsv,
)
from intergenic_regions.models import Gene
from intergenic_regions.reporting import write_report
from intergenic_regions.rna_expression import (
    DE_FIELDS,
    differential_expression,
    read_count_matrix,
    read_de_results,
    select_de_gene_sets,
)
from intergenic_regions.rna_inputs import (
    RnaSample,
    open_rna_coverage,
    read_gene_exons,
    read_rna_samples,
)
from intergenic_regions.rna_reporting import (
    RNA_FIGURE_NOTES,
    RNA_NOTES,
    plot_rna_results,
    rna_motif_profiles,
)
from intergenic_regions.rna_starts import (
    compare_rna_starts,
    consensus_rna_starts,
    detect_coverage_start,
    rna_search_interval,
)
from intergenic_regions.rna_tools import count_rna_genes, map_rna_reads

LOGGER = logging.getLogger(__name__)


def rna_input_paths(*, samples: Sequence[RnaSample]) -> list[Path]:
    """Collect source and alignment-index files for reproducibility hashes.

    Args:
        samples: Biological sample objects.

    Returns:
        Existing distinct paths, without including generated mapping outputs.
    """
    paths: list[Path] = []
    for sample in samples:
        for source in (
            sample.bam,
            sample.coverage,
            sample.coverage_plus,
            sample.coverage_minus,
            sample.read1,
            sample.read2,
        ):
            if source is not None:
                paths.append(source)
        if sample.bam:
            paths += [
                index
                for index in (
                    Path(str(sample.bam) + ".bai"),
                    sample.bam.with_suffix(".bai"),
                    Path(str(sample.bam) + ".crai"),
                    sample.bam.with_suffix(".crai"),
                    Path(str(sample.bam) + ".csi"),
                )
                if index.is_file()
            ]
    return list(dict.fromkeys(paths))


def rna_start_outputs(
    *,
    directory: Path,
    genome_path: Path,
    annotation_path: Path,
    samples: Sequence[RnaSample],
    annotation_format: str = "auto",
    genes: Sequence[Gene] | None = None,
    identifiers: Sequence[str] | None = None,
    upstream: int = 2000,
    downstream: int = 300,
    window: int = 20,
    min_depth: float = 2.0,
    min_fold_rise: float = 4.0,
    min_covered_fraction: float = 0.6,
    pseudocount: float = 0.5,
    min_replicates: int = 2,
    agreement_tolerance: int = 40,
    min_support_fraction: float = 0.5,
    min_mapq: int = 20,
    max_reads: int = 100000,
    exclude_duplicates: bool = True,
    profile_bin_width: int = 20,
    max_genes: int = 10000,
    max_profile_rows: int = 2000000,
    condition: str | None = None,
    reference: str | None = None,
) -> dict[str, Any]:
    """Write bounded replicate searches and condition-specific start audits.

    Args:
        directory: New/staged destination; caller controls atomic publication.
        genome_path: Genomic FASTA.
        annotation_path: Full annotation, including all blockers.
        samples: Indexed BAM/CRAM or transcript-strand/unstranded bedGraphs.
        annotation_format: Original supported annotation format.
        genes: Optional already parsed complete genes.
        identifiers: Target subset, or all genes up to the declared cap.
        upstream: Maximum safe upstream search.
        downstream: Maximum first-exon/gene-span search.
        window: Coverage-edge resolution/sustained-block size.
        min_depth: Raw downstream depth floor.
        min_fold_rise: Required descriptive rise ratio.
        min_covered_fraction: Sustained covered-base fraction.
        pseudocount: Positive depth pseudocount.
        min_replicates: Independent observations needed for support.
        agreement_tolerance: Maximum agreeing coordinate span.
        min_support_fraction: Required fraction of all condition samples.
        min_mapq: Alignment quality threshold.
        max_reads: Hard alignment cap per locus.
        exclude_duplicates: Exclude duplicate-marked reads.
        profile_bin_width: Bounded exported plot resolution.
        max_genes: Hard target-gene cap.
        max_profile_rows: Hard row cap, checked before reading coverage.
        condition: Optional target condition for descriptive shift comparisons.
        reference: Optional reference condition for those comparisons.

    Returns:
        Summary, calls, consensus starts, profiles and shifts.

    Raises:
        ValueError: Inputs are invalid or resource limits are exceeded.
    """
    if (
        not samples
        or len({s.sample_id for s in samples}) != len(samples)
        or len({(s.condition, s.replicate_id) for s in samples})
        != len(samples)
        or isinstance(profile_bin_width, bool)
        or profile_bin_width < 1
        or isinstance(max_genes, bool)
        or max_genes < 1
        or isinstance(max_profile_rows, bool)
        or max_profile_rows < 1
    ):
        raise ValueError("Invalid RNA samples/profile resource limits")
    detect_coverage_start(
        depth=[],
        genomic_start=0,
        strand="+",
        window=window,
        min_depth=min_depth,
        min_fold_rise=min_fold_rise,
        min_covered_fraction=min_covered_fraction,
        pseudocount=pseudocount,
    )
    consensus_rna_starts(
        calls=[],
        min_replicates=min_replicates,
        agreement_tolerance=agreement_tolerance,
        min_support_fraction=min_support_fraction,
    )
    conditions = {sample.condition for sample in samples}
    if (
        condition is not None
        and condition not in conditions
        or reference is not None
        and (
            condition is None
            or reference not in conditions
            or reference == condition
        )
    ):
        raise ValueError(
            "RNA selected/reference conditions are absent or equal"
        )
    full = (
        list(genes)
        if genes is not None
        else read_annotation(
            path=annotation_path, annotation_format=annotation_format
        )
    )
    exons = read_gene_exons(
        path=annotation_path, annotation_format=annotation_format, genes=full
    )
    by_id = {gene.gene_id: gene for gene in full}
    selected = list(identifiers) if identifiers is not None else list(by_id)
    if (
        not selected
        or len(set(selected)) != len(selected)
        or set(selected) - set(by_id)
    ):
        raise ValueError("RNA target genes are empty, repeated or unknown")
    if len(selected) > max_genes:
        raise ValueError(
            "RNA target count exceeds max_genes; use a subset or raise th"
            "e cap explicitly"
        )
    calls: list[dict[str, Any]] = []
    profiles: list[dict[str, Any]] = []
    with Genome(path=genome_path) as genome:
        index = GeneIndex(genes=full, lengths=genome.lengths)
        intervals: dict[str, tuple[int, int, str]] = {}
        for gene_id in selected:
            gene = by_id[gene_id]
            intervals[gene_id] = (
                rna_search_interval(
                    index=index,
                    gene=gene,
                    exons=exons.get(gene_id, []),
                    upstream=upstream,
                    downstream=downstream,
                )
                if gene.strand != "."
                else (gene.start, gene.start, "unknown_strand")
            )
        upper_bound = len(samples) * sum(
            math.ceil((end - start) / profile_bin_width) + 1
            for start, end, _ in intervals.values()
            if end > start
        )
        if upper_bound > max_profile_rows:
            raise ValueError(
                "RNA profile output exceeds max_profile_rows; increase bi"
                "n width or reduce genes"
            )
        for sample in samples:
            LOGGER.info(
                "Searching candidate RNA starts: sample=%s, genes=%d",
                sample.sample_id,
                len(selected),
            )
            with open_rna_coverage(
                sample=sample,
                lengths=genome.lengths,
                genome_path=genome_path,
                min_mapq=min_mapq,
                max_reads=max_reads,
                exclude_duplicates=exclude_duplicates,
            ) as coverage:
                for gene_id in selected:
                    gene = by_id[gene_id]
                    left, right, exon_source = intervals[gene_id]
                    anchor = gene.start if gene.strand != "-" else gene.end - 1
                    depth = (
                        coverage.depth(
                            contig=gene.contig,
                            start=left,
                            end=right,
                            strand=gene.strand,
                        )
                        if right > left
                        else []
                    )
                    call = detect_coverage_start(
                        depth=depth,
                        genomic_start=left,
                        strand=gene.strand if gene.strand != "." else "+",
                        window=window,
                        min_depth=min_depth,
                        min_fold_rise=min_fold_rise,
                        min_covered_fraction=min_covered_fraction,
                        pseudocount=pseudocount,
                    )
                    if not depth:
                        call["status"] = exon_source
                    calls.append(
                        {
                            "gene_id": gene_id,
                            "contig": gene.contig,
                            "strand": gene.strand,
                            "annotation_anchor": anchor,
                            "sample_id": sample.sample_id,
                            "condition": sample.condition,
                            "replicate_id": sample.replicate_id,
                            "strandedness": sample.strandedness,
                            "search_start": left,
                            "search_end": right,
                            "exon_source": exon_source,
                            **call,
                        }
                    )
                    oriented = (
                        depth if gene.strand == "+" else list(reversed(depth))
                    )
                    relative_first = (
                        left - anchor
                        if gene.strand == "+"
                        else anchor - (right - 1)
                    )
                    grouped: dict[int, list[float]] = defaultdict(list)
                    for offset, value in enumerate(oriented):
                        relative = relative_first + offset
                        grouped[
                            math.floor(relative / profile_bin_width)
                            * profile_bin_width
                        ].append(value)
                    for bin_start, values in grouped.items():
                        mean = sum(values) / len(values)
                        profiles.append(
                            {
                                "gene_id": gene_id,
                                "sample_id": sample.sample_id,
                                "condition": sample.condition,
                                "replicate_id": sample.replicate_id,
                                "relative_start": bin_start,
                                "relative_end": bin_start + profile_bin_width,
                                "observed_bases": len(values),
                                "raw_depth": mean,
                                "divided_depth": mean / sample.size_factor,
                                "size_factor": sample.size_factor,
                                "coordinate_system": (
                                    "annotated_5_prime_gene_base"
                                ),
                            }
                        )
    starts = consensus_rna_starts(
        calls=calls,
        min_replicates=min_replicates,
        agreement_tolerance=agreement_tolerance,
        min_support_fraction=min_support_fraction,
    )
    shifts = (
        compare_rna_starts(
            starts=starts, condition=condition, reference=reference
        )
        if reference and condition
        else []
    )
    directory.mkdir(parents=True, exist_ok=True)
    for filename, rows, empty in (
        ("replicate_calls.tsv", calls, ("gene_id", "sample_id", "status")),
        ("candidate_starts.tsv", starts, ("gene_id", "condition", "status")),
        (
            "coverage_profiles.tsv",
            profiles,
            ("gene_id", "sample_id", "raw_depth"),
        ),
        (
            "start_shifts.tsv",
            shifts,
            (
                "gene_id",
                "condition",
                "reference",
                "oriented_shift_bp",
                "status",
            ),
        ),
    ):
        write_tsv(
            path=directory / filename,
            rows=rows,
            fields=tuple(rows[0]) if rows else empty,
        )
    summary = {
        "status": "completed",
        "biological_samples": len(samples),
        "searched_genes": len(selected),
        "supported_gene_conditions": sum(
            row["status"] == "replicate_supported_candidate" for row in starts
        ),
        "candidate_status_counts": dict(
            Counter(row["status"] for row in starts)
        ),
        "profile_rows": len(profiles),
        "profile_bin_width": profile_bin_width,
        "conditions": sorted(conditions),
        "anchor_condition": condition,
        "unstranded_samples": [
            s.sample_id for s in samples if s.strandedness == "unstranded"
        ],
        "genes_without_annotated_exons": sorted(set(selected) - set(exons)),
        "start_method": (
            "descriptive_sustained_coverage_rise_independent_replicates"
        ),
        "tss_status": "unvalidated_coverage_edges",
        "edge_p_values": None,
        "depth_unit": (
            "aligned_fragment_depth_for_BAM_or_supplied_bedGraph_depth"
        ),
        "coverage_divisors": {s.sample_id: s.size_factor for s in samples},
    }
    write_json(path=directory / "summary.json", data=summary)
    return {
        "summary": summary,
        "calls": calls,
        "starts": starts,
        "profiles": profiles,
        "shifts": shifts,
    }


def _write_rna_report(
    *,
    directory: Path,
    data: Mapping[str, Any],
    labels: Mapping[str, str] | None = None,
    de_rows: Sequence[Mapping[str, Any]] = (),
    motif_profiles: Sequence[Mapping[str, Any]] = (),
    counts: Mapping[str, Mapping[str, int]] | None = None,
    max_loci: int = 6,
) -> None:
    summary = data["summary"]
    images = plot_rna_results(
        profiles=data["profiles"],
        starts=data["starts"],
        calls=data["calls"],
        shifts=data["shifts"],
        de_rows=de_rows,
        motif_profiles=motif_profiles,
        labels=labels or {},
        directory=directory / "figures",
        counts=counts,
        size_factors=summary.get("coverage_divisors"),
        max_loci=max_loci,
    )
    write_report(
        path=directory / "report.html",
        title="RNA evidence · candidate starts and upstream motifs",
        summary={"rna": summary},
        tables={
            "Candidate start neighbourhoods": data["starts"],
            "Independent biological replicate calls": data["calls"],
            "Descriptive candidate start shifts": data["shifts"],
            "Gene differential expression": de_rows,
            "Upstream motif distance profiles": motif_profiles,
        },
        images=images,
        notes=[
            *RNA_NOTES,
            *[
                f"DE model warning: {warning}"
                for warning in (summary.get("de") or {}).get("warnings", [])
            ],
        ],
        figure_notes=RNA_FIGURE_NOTES,
        links={
            "Candidate starts TSV": "candidate_starts.tsv",
            "Replicate calls TSV": "replicate_calls.tsv",
            "Coverage profiles TSV": "coverage_profiles.tsv",
            "Start shifts TSV": "start_shifts.tsv",
        },
    )


def rna_starts_workflow(
    *,
    samples_path: Path,
    genome_path: Path,
    annotation_path: Path,
    output: Path,
    identifiers_path: Path | None = None,
    annotation_format: str = "auto",
    settings: dict[str, Any] | None = None,
    max_loci: int = 6,
) -> dict[str, Any]:
    """Run a complete standalone candidate-start workflow with an HTML report.

    Args:
        samples_path: Mapped/coverage biological sample TSV.
        genome_path: Genomic FASTA.
        annotation_path: Full annotation.
        output: New result directory, published atomically.
        identifiers_path: Optional exact target-gene list.
        annotation_format: Supported annotation format.
        settings: rna_start_outputs keyword settings.
        max_loci: Bounded locus plot count.

    Returns:
        Candidate-start summary.
    """
    samples = read_rna_samples(path=samples_path)
    with output_bundle(path=output) as stage:
        data = rna_start_outputs(
            directory=stage,
            genome_path=genome_path,
            annotation_path=annotation_path,
            samples=samples,
            annotation_format=annotation_format,
            identifiers=read_identifiers(path=identifiers_path)
            if identifiers_path
            else None,
            **(settings or {}),
        )
        _write_rna_report(directory=stage, data=data, max_loci=max_loci)
        write_json(
            path=stage / "manifest.json",
            data=provenance(
                inputs=[
                    samples_path,
                    genome_path,
                    annotation_path,
                    *rna_input_paths(samples=samples),
                ]
                + ([identifiers_path] if identifiers_path else []),
                settings={
                    "command": "rna-starts",
                    "annotation_format": annotation_format,
                    **(settings or {}),
                },
            ),
        )
    return data["summary"]


def rna_mapping_workflow(
    *,
    samples_path: Path,
    genome_path: Path,
    output: Path,
    settings: dict[str, Any] | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Control mapping or export a dry-run plan using argument vectors.

    Args:
        samples_path: FASTQ sample/condition TSV.
        genome_path: Genomic reference.
        output: New result directory.
        settings: Controlled STAR/samtools settings.
        dry_run: Plan only, requiring no installed aligner.

    Returns:
        Explicit mapping completion or dry-run status.
    """
    samples = read_rna_samples(path=samples_path)
    with output_bundle(path=output) as stage:
        _, summary = map_rna_reads(
            samples=samples,
            genome_path=genome_path,
            directory=stage,
            dry_run=dry_run,
            **(settings or {}),
        )
        summary["working_directory"] = str(output.resolve())
        write_json(
            path=stage / "mapping_plan.json",
            data={
                key: value for key, value in summary.items() if key != "status"
            },
        )
        write_json(path=stage / "summary.json", data=summary)
        write_json(
            path=stage / "manifest.json",
            data=provenance(
                inputs=[
                    samples_path,
                    genome_path,
                    *rna_input_paths(samples=samples),
                ],
                settings={
                    "command": "rna-map",
                    "dry_run": dry_run,
                    **{
                        name: str(value) if isinstance(value, Path) else value
                        for name, value in (settings or {}).items()
                    },
                },
            ),
        )
    return summary


def rna_de_workflow(
    *,
    samples_path: Path,
    counts_path: Path,
    output: Path,
    condition: str,
    reference: str,
    settings: dict[str, Any] | None = None,
    selection_settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Fit replicate-aware DE from raw counts and export transparent gene sets.

    Args:
        samples_path: Biological observation/design metadata TSV.
        counts_path: Raw integer gene count TSV.
        output: New result directory.
        condition: Numerator/target condition.
        reference: Reference condition.
        settings: differential_expression keyword settings.
        selection_settings: Predeclared DE foreground/control thresholds.

    Returns:
        DE method/QC summary plus selected gene counts.
    """
    samples = read_rna_samples(path=samples_path, require_source=False)
    counts = read_count_matrix(path=counts_path, samples=samples)
    with output_bundle(path=output) as stage:
        rows, summary = differential_expression(
            counts=counts,
            samples=samples,
            condition=condition,
            reference=reference,
            **(settings or {}),
        )
        positive, negative, audit = select_de_gene_sets(
            rows=rows, **(selection_settings or {})
        )
        _write_gene_sets(
            directory=stage, positive=positive, negative=negative, audit=audit
        )
        write_tsv(path=stage / "de_results.tsv", rows=rows, fields=DE_FIELDS)
        summary.update(
            {
                "positive_genes": len(positive),
                "screening_controls": len(negative),
            }
        )
        write_json(path=stage / "summary.json", data=summary)
        images = plot_rna_results(
            profiles=[],
            starts=[],
            calls=[],
            shifts=[],
            de_rows=rows,
            motif_profiles=[],
            labels={
                **dict.fromkeys(positive, "positive"),
                **dict.fromkeys(negative, "negative"),
            },
            directory=stage / "figures",
            counts=counts,
            size_factors=summary["size_factors"],
            max_loci=0,
        )
        write_report(
            path=stage / "report.html",
            title="RNA gene differential expression",
            summary={
                "rna": {"biological_samples": len(samples), "de": summary}
            },
            tables={
                "Gene differential expression": rows,
                "Gene-set selection audit": audit,
            },
            images=images,
            notes=RNA_NOTES,
            figure_notes=RNA_FIGURE_NOTES,
            links={
                "DE results TSV": "de_results.tsv",
                "Gene-set audit TSV": "gene_set_audit.tsv",
            },
        )
        write_json(
            path=stage / "manifest.json",
            data=provenance(
                inputs=[samples_path, counts_path],
                settings={
                    "command": "rna-de",
                    "condition": condition,
                    "reference": reference,
                    "de": settings or {},
                    "selection": selection_settings or {},
                },
            ),
        )
    return summary


def _write_gene_sets(
    *,
    directory: Path,
    positive: Sequence[str],
    negative: Sequence[str],
    audit: Sequence[Mapping[str, Any]],
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "positive_genes.txt").write_text(
        data="".join(f"{gene}\n" for gene in positive), encoding="utf-8"
    )
    (directory / "negative_genes.txt").write_text(
        data="".join(f"{gene}\n" for gene in negative), encoding="utf-8"
    )
    write_tsv(
        path=directory / "gene_set_audit.tsv",
        rows=audit,
        fields=tuple(audit[0]) if audit else ("gene_id", "label", "reason"),
    )


def _read_rows(*, path: Path) -> list[dict[str, Any]]:
    with path.open(mode="r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(f=stream, delimiter="\t"))


def rna_pipeline_workflow(
    *,
    samples_path: Path,
    genome_path: Path,
    annotation_path: Path,
    output: Path,
    condition: str,
    reference: str | None = None,
    positive_genes: Path | None = None,
    negative_genes: Path | None = None,
    de_path: Path | None = None,
    counts_path: Path | None = None,
    run_de: bool = False,
    annotation_format: str = "auto",
    map_reads: bool = False,
    mapping_settings: dict[str, Any] | None = None,
    counting_settings: dict[str, Any] | None = None,
    de_settings: dict[str, Any] | None = None,
    selection_settings: dict[str, Any] | None = None,
    rna_settings: dict[str, Any] | None = None,
    pipeline_settings: dict[str, Any] | None = None,
    fallback: str = "annotation",
    allow_provisional: bool = False,
    max_loci: int = 6,
) -> dict[str, Any]:
    """Run RNA evidence through safe motifs, automatic ML/SHAP and HTML.

    Args:
        samples_path: Biological sample TSV with BAM, tracks or reads.
        genome_path: Full genomic FASTA.
        annotation_path: Full annotation, including all gene blockers.
        output: New directory, published after the complete workflow succeeds.
        condition: One start-anchor condition used for both gene classes.
        reference: Required reference for DE, optional for supplied gene lists.
        positive_genes: Optional exact supplied positive IDs.
        negative_genes: Exact supplied negative IDs; supply both lists.
        de_path: Optional imported DESeq2-style result TSV with known contrast.
        counts_path: Raw count matrix; otherwise count BAM exons if needed.
        run_de: Also count BAM exons and fit DE while retaining curated labels.
        annotation_format: Supported annotation format.
        map_reads: Explicitly control mapping of any FASTQ samples first.
        mapping_settings: Bounded STAR/samtools settings.
        counting_settings: Bounded featureCounts settings.
        de_settings: PyDESeq2 contrast/design settings.
        selection_settings: Predeclared DE foreground/control selection.
        rna_settings: Coverage/start detection and resource limits.
        pipeline_settings: Existing motif/ML/SHAP/evidence/scan settings.
        fallback: Explicit annotation fallback or exclusion for unusable edges.
        allow_provisional: Explicitly permit a singleton unvalidated edge.
        max_loci: Bounded RNA locus plot count.

    Returns:
        Combined RNA, motif, model, multi-scale and optional scan summary.

    Raises:
        ValueError: Inputs, contrasts or gene classes cannot support analysis.
    """
    if bool(positive_genes) != bool(negative_genes) or (
        de_path and (counts_path or run_de)
    ):
        raise ValueError(
            "Supply both gene lists; choose imported DE or raw counts"
        )
    if (
        not positive_genes or de_path or counts_path or run_de
    ) and not reference:
        raise ValueError(
            "DE input/selection requires an explicit reference condition"
        )
    original_samples = read_rna_samples(path=samples_path)
    conditions = {sample.condition for sample in original_samples}
    if (
        condition not in conditions
        or reference
        and (reference not in conditions or reference == condition)
    ):
        raise ValueError("RNA anchor/contrast conditions are absent or equal")
    genes = read_annotation(
        path=annotation_path, annotation_format=annotation_format
    )
    gene_ids = {gene.gene_id for gene in genes}
    exons = read_gene_exons(
        path=annotation_path, annotation_format=annotation_format, genes=genes
    )
    from intergenic_regions.workflows import (
        analysis_report,
        pipeline_workflow,
        read_region_table,
    )

    with output_bundle(path=output) as stage:
        samples = list(original_samples)
        mapping: dict[str, Any] = {"status": "not_required"}
        raw = [sample for sample in samples if sample.read1]
        if raw:
            if not map_reads:
                raise ValueError(
                    "FASTQ samples need --map-reads or a precomputed inde"
                    "xed BAM"
                )
            mapped, mapping = map_rna_reads(
                samples=raw,
                genome_path=genome_path,
                directory=stage / "mapping",
                **(mapping_settings or {}),
            )
            replacements = {sample.sample_id: sample for sample in mapped}
            samples = [
                replacements.get(sample.sample_id, sample)
                for sample in samples
            ]
            mapping["working_directory"] = str((output / "mapping").resolve())
            write_json(
                path=stage / "mapping" / "mapping_plan.json",
                data={
                    key: value
                    for key, value in mapping.items()
                    if key != "status"
                },
            )
        counts = (
            read_count_matrix(path=counts_path, samples=samples)
            if counts_path
            else None
        )
        de_rows: list[dict[str, Any]] = []
        de_summary: dict[str, Any] = {
            "status": "not_supplied",
            "reason": "User-supplied gene lists",
        }
        if de_path:
            de_rows = read_de_results(path=de_path)
            de_summary = {
                "status": "imported",
                "condition": condition,
                "reference": reference,
                "method": "user_supplied_DESeq2_style_results",
                "contrast_verification": "user_declared",
            }
        elif not positive_genes or counts_path or run_de:
            if counts is None:
                counts = count_rna_genes(
                    samples=samples,
                    genes=genes,
                    exons=exons,
                    directory=stage / "expression" / "counting",
                    **(counting_settings or {}),
                )
            assert reference is not None
            de_rows, de_summary = differential_expression(
                counts=counts,
                samples=samples,
                condition=condition,
                reference=reference,
                **(de_settings or {}),
            )
            samples = [
                replace(
                    sample,
                    size_factor=de_summary["size_factors"].get(
                        sample.sample_id, sample.size_factor
                    ),
                )
                for sample in samples
            ]
        if de_rows and {row["gene_id"] for row in de_rows} - gene_ids:
            raise ValueError("DE gene IDs do not match the full annotation")
        if positive_genes and negative_genes:
            positive, negative = (
                read_identifiers(path=positive_genes),
                read_identifiers(path=negative_genes),
            )
            audit = [
                {
                    "gene_id": gene,
                    "label": label,
                    "reason": "user_supplied_gene_set",
                }
                for label, ids in (
                    ("positive", positive),
                    ("negative", negative),
                )
                for gene in ids
            ]
        else:
            positive, negative, audit = select_de_gene_sets(
                rows=de_rows, **(selection_settings or {})
            )
        if (
            not positive
            or not negative
            or set(positive) & set(negative)
            or (set(positive) | set(negative)) - gene_ids
        ):
            raise ValueError(
                "RNA motif analysis needs two non-empty, disjoint annotat"
                "ed gene sets"
            )
        _write_gene_sets(
            directory=stage / "expression",
            positive=positive,
            negative=negative,
            audit=audit,
        )
        if de_rows:
            write_tsv(
                path=stage / "expression" / "de_results.tsv",
                rows=de_rows,
                fields=DE_FIELDS,
            )
        if counts:
            write_tsv(
                path=stage / "expression" / "gene_counts.tsv",
                rows=[
                    {"gene_id": gene, **row} for gene, row in counts.items()
                ],
                fields=("gene_id", *(s.sample_id for s in samples)),
            )
        write_json(path=stage / "expression" / "summary.json", data=de_summary)
        data = rna_start_outputs(
            directory=stage / "rna",
            genome_path=genome_path,
            annotation_path=annotation_path,
            samples=samples,
            genes=genes,
            identifiers=[*positive, *negative],
            annotation_format=annotation_format,
            condition=condition,
            reference=reference,
            **(rna_settings or {}),
        )
        sequence_directory = stage / "sequence_analysis"
        summary = pipeline_workflow(
            genome_path=genome_path,
            annotation_path=annotation_path,
            positive_genes=stage / "expression" / "positive_genes.txt",
            negative_genes=stage / "expression" / "negative_genes.txt",
            output=sequence_directory,
            annotation_format=annotation_format,
            rna_starts=data["starts"],
            rna_condition=condition,
            rna_fallback=fallback,
            allow_provisional_rna=allow_provisional,
            **(pipeline_settings or {}),
        )
        sequence_manifest = json.loads(
            (sequence_directory / "manifest.json").read_text(encoding="utf-8")
        )
        for entry in sequence_manifest["inputs"]:
            source = Path(entry["path"])
            if source.is_relative_to(stage):
                entry["path"] = str(source.relative_to(stage))
        for child in sequence_directory.iterdir():
            child.rename(target=stage / child.name)
        sequence_directory.rmdir()
        anchor_audit = _read_rows(path=stage / "rna_anchor_audit.tsv")
        for row in anchor_audit:
            row["applied_anchor"] = (
                int(row["applied_anchor"]) if row["applied_anchor"] else None
            )
        regions = [
            *read_region_table(
                path=stage / "positive_extraction" / "regions.tsv"
            ),
            *read_region_table(
                path=stage / "negative_extraction" / "regions.tsv"
            ),
        ]
        # The standard audit reader omits sequence bases; recover lengths for
        # opportunity denominators from the explicitly generated class FASTAs.
        from intergenic_regions.io import read_fasta

        positive_records = read_fasta(path=stage / "positive.fasta")
        negative_records = read_fasta(path=stage / "negative.fasta")
        records = {**positive_records, **negative_records}
        regions = [
            replace(region, sequence=records[region.sequence_id])
            for region in regions
            if region.sequence_id in records
        ]
        sequence_labels = {
            **dict.fromkeys(positive_records, "positive"),
            **dict.fromkeys(negative_records, "negative"),
        }
        motif_profiles = rna_motif_profiles(
            regions=regions,
            anchors=anchor_audit,
            sites=_read_rows(path=stage / "motifs" / "motif_sites.tsv"),
            motifs=_read_rows(path=stage / "motifs" / "motif_enrichment.tsv"),
            labels=sequence_labels,
        )
        write_tsv(
            path=stage / "rna" / "motif_distance_profiles.tsv",
            rows=motif_profiles,
            fields=tuple(motif_profiles[0])
            if motif_profiles
            else (
                "motif_id",
                "label",
                "anchor_source",
                "site_count",
                "eligible_motif_centres",
            ),
        )
        candidates = _read_rows(path=stage / "candidate_priorities.tsv")
        de_by_id = {row["gene_id"]: row for row in de_rows}
        for row in candidates:
            expression = de_by_id.get(row["gene_id"], {})
            row.update(
                {
                    f"de_{name}": expression.get(name)
                    for name in ("baseMean", "log2FoldChange", "padj")
                }
            )
        for filename in ("candidate_summary.tsv", "candidate_priorities.tsv"):
            write_tsv(
                path=stage / filename,
                rows=candidates,
                fields=tuple(candidates[0]),
            )
        rna_summary = data["summary"]
        rna_summary.update(
            {
                "de": de_summary,
                "mapping": mapping,
                "applied_rna_anchors": sum(
                    row["anchor_source"] == "rna_coverage_candidate"
                    for row in anchor_audit
                ),
                "anchor_sources": dict(
                    Counter(row["anchor_source"] for row in anchor_audit)
                ),
                "positive_gene_set": len(positive),
                "negative_gene_set": len(negative),
                "gene_sets": "user_supplied"
                if positive_genes
                else "DE_threshold_screening_sets",
                "regulatory_status": (
                    "unvalidated_promoter_or_regulatory_hypotheses"
                ),
                "ML_features": (
                    "sequence_only_RNA_label_defining_features_excluded"
                ),
            }
        )
        summary["rna"] = rna_summary
        write_json(path=stage / "rna" / "summary.json", data=rna_summary)
        _write_rna_report(
            directory=stage / "rna",
            data=data,
            labels={
                **dict.fromkeys(positive, "positive"),
                **dict.fromkeys(negative, "negative"),
            },
            de_rows=de_rows,
            motif_profiles=motif_profiles,
            counts=counts,
            max_loci=max_loci,
        )
        write_json(path=stage / "summary.json", data=summary)
        inputs = [
            samples_path,
            genome_path,
            annotation_path,
            *rna_input_paths(samples=original_samples),
        ] + [
            path
            for path in (positive_genes, negative_genes, de_path, counts_path)
            if path is not None
        ]
        manifest = provenance(
            inputs=list(dict.fromkeys(inputs)),
            settings={
                "command": "rna-pipeline",
                "condition": condition,
                "reference": reference,
                "rna": rna_settings or {},
                "fallback": fallback,
                "allow_provisional": allow_provisional,
                "de": de_settings or {},
                "selection": selection_settings or {},
                "sequence_analysis": sequence_manifest,
            },
        )
        write_json(path=stage / "manifest.json", data=manifest)
        analysis_report(
            directory=stage,
            summary=summary,
            candidates=candidates,
            motif_directory=stage / "motifs",
            extra_tables={
                title: _read_rows(path=stage / filename)
                for title, filename in (
                    ("Leakage exclusions", "background_audit.tsv"),
                    ("Reference overlaps", "reference_overlaps.tsv"),
                )
                if (stage / filename).is_file()
            },
        )
    return summary
