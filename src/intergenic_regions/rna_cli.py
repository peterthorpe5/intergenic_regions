"""Named RNA commands with curated gene sets as primary inputs."""

import argparse
from pathlib import Path
from typing import Any


def add_rna_start_options(*, parser: argparse.ArgumentParser) -> None:
    """Add explicit coverage, replicate-agreement and resource limits.

    Args:
        parser: RNA start/pipeline command parser.
    """
    parser.add_argument("--rna-upstream", type=int, default=2000)
    parser.add_argument("--rna-downstream", type=int, default=300)
    parser.add_argument("--rna-window", type=int, default=20)
    parser.add_argument("--rna-min-depth", type=float, default=2.0)
    parser.add_argument("--rna-min-fold-rise", type=float, default=4.0)
    parser.add_argument("--rna-min-covered-fraction", type=float, default=0.6)
    parser.add_argument("--rna-pseudocount", type=float, default=0.5)
    parser.add_argument("--rna-min-replicates", type=int, default=2)
    parser.add_argument("--rna-agreement-tolerance", type=int, default=40)
    parser.add_argument("--rna-min-support-fraction", type=float, default=0.5)
    parser.add_argument("--rna-min-mapq", type=int, default=20)
    parser.add_argument("--rna-max-reads", type=int, default=100000)
    parser.add_argument("--rna-include-duplicates", action="store_true")
    parser.add_argument("--rna-profile-bin-width", type=int, default=20)
    parser.add_argument("--rna-max-genes", type=int, default=10000)
    parser.add_argument("--rna-max-profile-rows", type=int, default=2000000)
    parser.add_argument("--rna-max-locus-plots", type=int, default=6)


def rna_start_options(*, args: argparse.Namespace) -> dict[str, Any]:
    """Convert named CLI coverage settings to tested API keyword arguments.

    Args:
        args: Parsed RNA command.

    Returns:
        rna_start_outputs keyword settings, excluding plot/output arguments.
    """
    names = (
        "upstream",
        "downstream",
        "window",
        "min_depth",
        "min_fold_rise",
        "min_covered_fraction",
        "pseudocount",
        "min_replicates",
        "agreement_tolerance",
        "min_support_fraction",
        "min_mapq",
        "max_reads",
        "profile_bin_width",
        "max_genes",
        "max_profile_rows",
    )
    return {
        **{name: getattr(args, f"rna_{name}") for name in names},
        "exclude_duplicates": not args.rna_include_duplicates,
    }


def _add_mapping_options(*, parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--star-index", type=Path)
    parser.add_argument("--star-executable", default="STAR")
    parser.add_argument("--samtools-executable", default="samtools")
    parser.add_argument("--sort-ram", type=int, default=2000000000)
    parser.add_argument("--index-ram", type=int, default=8000000000)
    parser.add_argument("--read-length", type=int)
    parser.add_argument("--tool-timeout", type=float, default=86400)


def _add_de_options(*, parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--de-covariates", nargs="+", default=[])
    parser.add_argument("--de-alpha", type=float, default=0.05)
    parser.add_argument("--de-min-total-count", type=int, default=10)
    parser.add_argument(
        "--de-direction", choices=("up", "down", "both"), default="up"
    )
    parser.add_argument(
        "--de-min-abs-log2-fold-change", type=float, default=1.0
    )
    parser.add_argument("--de-min-base-mean", type=float, default=10.0)
    parser.add_argument(
        "--de-control-max-abs-log2-fold-change", type=float, default=0.25
    )
    parser.add_argument("--de-control-min-padj", type=float, default=0.5)


def add_rna_commands(*, subcommands: Any) -> None:
    """Add optional RNA commands without importing BAM or DE libraries.

    Args:
        subcommands: Root argparse subparser collection.
    """
    from intergenic_regions.cli import (
        add_evidence_options,
        add_extraction_options,
        add_learning_options,
        add_motif_options,
        add_region_options,
        add_scan_options,
    )

    for name, description in (
        (
            "rna-pipeline",
            (
                "Curated gene sets or DE sets, RNA starts, motifs and automat"
                "ic ML/SHAP"
            ),
        ),
        (
            "rna-starts",
            "Infer unvalidated coverage-start neighbourhoods from mapped RNA",
        ),
        (
            "rna-map",
            "Control splice-aware STAR mapping, or export a dry-run plan",
        ),
        (
            "rna-de",
            "Fit optional replicate-aware DE from raw integer gene counts",
        ),
    ):
        parser = subcommands.add_parser(name=name, help=description)
        parser.add_argument("--verbose", action="store_true")
        parser.add_argument("--log-file", type=Path)
        parser.add_argument("--output-dir", "-o", type=Path, required=True)
        parser.add_argument(
            "--samples",
            type=Path,
            required=True,
            help="Sample/condition/biological-replicate TSV",
        )
        if name == "rna-pipeline":
            add_extraction_options(parser=parser)
            add_evidence_options(parser=parser)
            add_motif_options(parser=parser)
            add_learning_options(parser=parser)
            add_region_options(parser=parser)
            add_scan_options(parser=parser)
            parser.add_argument(
                "--positive-genes",
                type=Path,
                help="Curated foreground list; no DE analysis required",
            )
            parser.add_argument(
                "--negative-genes",
                type=Path,
                help=(
                    "Curated control list; preserved when DE context is added"
                ),
            )
            source = parser.add_mutually_exclusive_group()
            source.add_argument("--de-tsv", type=Path)
            source.add_argument("--counts-tsv", type=Path)
            parser.add_argument("--map-reads", action="store_true")
            parser.add_argument(
                "--run-de",
                action="store_true",
                help="Count BAM exons and fit DE; keep supplied gene labels",
            )
            parser.add_argument(
                "--featurecounts-executable", default="featureCounts"
            )
            parser.add_argument(
                "--no-ml", dest="ai", action="store_false", default=True
            )
            parser.add_argument(
                "--group-by", choices=("gene", "contig"), default="gene"
            )
            parser.add_argument("--match-background", action="store_true")
            parser.add_argument("--background-ratio", type=int, default=1)
            parser.add_argument("--max-gc-difference", type=float, default=0.1)
            parser.add_argument("--max-length-ratio", type=float, default=1.5)
            parser.add_argument("--scan-genome", action="store_true")
            parser.add_argument(
                "--rna-start-fallback",
                choices=("annotation", "exclude"),
                default="annotation",
            )
            parser.add_argument("--allow-provisional-rna", action="store_true")
        elif name != "rna-de":
            parser.add_argument("--genome", type=Path, required=True)
            if name == "rna-starts":
                parser.add_argument(
                    "--annotation", "--gff", type=Path, required=True
                )
                parser.add_argument(
                    "--annotation-format",
                    choices=("auto", "gff3", "gtf", "tsv"),
                    default="auto",
                )
                parser.add_argument("--genes", type=Path)
        if name in {"rna-pipeline", "rna-starts"}:
            add_rna_start_options(parser=parser)
        if name in {"rna-pipeline", "rna-map"}:
            _add_mapping_options(parser=parser)
        if name in {"rna-pipeline", "rna-de"}:
            _add_de_options(parser=parser)
        if name != "rna-map":
            parser.add_argument("--condition", required=name != "rna-starts")
            parser.add_argument(
                "--reference-condition", required=name == "rna-de"
            )
        if name == "rna-de":
            parser.add_argument("--counts-tsv", type=Path, required=True)
            parser.add_argument("--threads", type=int, default=1)
        if name == "rna-map":
            parser.add_argument("--dry-run", action="store_true")


def _mapping_settings(*, args: argparse.Namespace) -> dict[str, Any]:
    return {
        "star_index": args.star_index,
        "threads": args.threads,
        "sort_ram": args.sort_ram,
        "index_ram": args.index_ram,
        "read_length": args.read_length,
        "star_executable": args.star_executable,
        "samtools_executable": args.samtools_executable,
        "timeout": args.tool_timeout,
    }


def _de_settings(*, args: argparse.Namespace) -> dict[str, Any]:
    return {
        "covariates": args.de_covariates,
        "alpha": args.de_alpha,
        "min_total_count": args.de_min_total_count,
        "threads": args.threads,
    }


def _selection_settings(*, args: argparse.Namespace) -> dict[str, Any]:
    return {
        "direction": args.de_direction,
        "alpha": args.de_alpha,
        "min_abs_log2_fold_change": args.de_min_abs_log2_fold_change,
        "min_base_mean": args.de_min_base_mean,
        (
            "control_max_abs_log2_fold_change"
        ): args.de_control_max_abs_log2_fold_change,
        "control_min_padj": args.de_control_min_padj,
    }


def dispatch_rna(*, args: argparse.Namespace) -> dict[str, Any]:
    """Execute a named RNA command while keeping curated labels unchanged.

    Args:
        args: Parsed RNA command options.

    Returns:
        JSON-compatible workflow summary.

    Raises:
        ValueError: The command or inputs are invalid.
    """
    from intergenic_regions.rna_workflows import (
        rna_de_workflow,
        rna_mapping_workflow,
        rna_pipeline_workflow,
        rna_starts_workflow,
    )

    if args.command == "rna-map":
        return rna_mapping_workflow(
            samples_path=args.samples,
            genome_path=args.genome,
            output=args.output_dir,
            settings=_mapping_settings(args=args),
            dry_run=args.dry_run,
        )
    if args.command == "rna-de":
        return rna_de_workflow(
            samples_path=args.samples,
            counts_path=args.counts_tsv,
            output=args.output_dir,
            condition=args.condition,
            reference=args.reference_condition,
            settings=_de_settings(args=args),
            selection_settings=_selection_settings(args=args),
        )
    if args.command == "rna-starts":
        return rna_starts_workflow(
            samples_path=args.samples,
            genome_path=args.genome,
            annotation_path=args.annotation,
            annotation_format=args.annotation_format,
            identifiers_path=args.genes,
            output=args.output_dir,
            settings={
                **rna_start_options(args=args),
                "condition": args.condition,
                "reference": args.reference_condition,
            },
            max_loci=args.rna_max_locus_plots,
        )
    if args.command == "rna-pipeline":
        from intergenic_regions.cli import (
            evidence_options,
            extraction_options,
            learning_options,
            motif_options,
            region_options,
            scan_options,
        )

        return rna_pipeline_workflow(
            samples_path=args.samples,
            genome_path=args.genome,
            annotation_path=args.annotation,
            annotation_format=args.annotation_format,
            output=args.output_dir,
            condition=args.condition,
            reference=args.reference_condition,
            positive_genes=args.positive_genes,
            negative_genes=args.negative_genes,
            de_path=args.de_tsv,
            counts_path=args.counts_tsv,
            run_de=args.run_de,
            map_reads=args.map_reads,
            mapping_settings=_mapping_settings(args=args),
            counting_settings={
                "executable": args.featurecounts_executable,
                "threads": args.threads,
                "min_mapq": args.rna_min_mapq,
                "exclude_duplicates": not args.rna_include_duplicates,
                "timeout": args.tool_timeout,
            },
            de_settings=_de_settings(args=args),
            selection_settings=_selection_settings(args=args),
            rna_settings=rna_start_options(args=args),
            fallback=args.rna_start_fallback,
            allow_provisional=args.allow_provisional_rna,
            max_loci=args.rna_max_locus_plots,
            pipeline_settings={
                "extraction_settings": extraction_options(args=args),
                "motif_path": args.motifs,
                "enrichment_settings": motif_options(args=args),
                "use_ai": args.ai,
                "learning_settings": learning_options(args=args),
                "groups_path": args.groups_tsv,
                "group_by": args.group_by,
                "use_regions": not args.no_regions,
                "regional_settings": region_options(args=args),
                "evidence_settings": evidence_options(args=args),
                "scan_genome": args.scan_genome,
                "scan_settings": scan_options(args=args),
                "matching_settings": {
                    "ratio": args.background_ratio,
                    "max_gc_difference": args.max_gc_difference,
                    "max_length_ratio": args.max_length_ratio,
                }
                if args.match_background
                else None,
            },
        )
    raise ValueError("Unknown RNA command")
