"""Named-option command-line interface with explicit, non-zero failures."""

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from intergenic_regions import __version__
from intergenic_regions.io import (
    output_bundle,
    provenance,
    read_fasta,
    write_json,
    write_tsv,
)
from intergenic_regions.reporting import write_report
from intergenic_regions.workflows import (
    attach_evidence,
    enrichment_workflow,
    extract_workflow,
    genome_scan_workflow,
    learning_workflow,
    pipeline_workflow,
    read_region_table,
)

LOGGER = logging.getLogger(__name__)


def add_extraction_options(*, parser: argparse.ArgumentParser) -> None:
    """Add genome, annotation and flank settings to a command.

    Args:
        parser: Command parser.
    """
    parser.add_argument("--genome", "-g", type=Path, required=True)
    parser.add_argument("--annotation", "--gff", type=Path, required=True)
    parser.add_argument(
        "--annotation-format",
        choices=("auto", "gff3", "gtf", "tsv"),
        default="auto",
    )
    limits = parser.add_mutually_exclusive_group()
    limits.add_argument(
        "--length",
        "--upstream",
        "-u",
        type=int,
        default=1000,
        help="Maximum flank length (default: 1000 bp)",
    )
    limits.add_argument(
        "--full-gap",
        action="store_true",
        help="Return the complete contiguous intergenic flank",
    )
    parser.add_argument(
        "--direction",
        choices=("upstream", "downstream", "both"),
        default="upstream",
    )
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--min-length", "-m", type=int, default=1)
    parser.add_argument("--max-ambiguous-fraction", type=float, default=1.0)
    parser.add_argument("--mask-lowercase", action="store_true")


def add_evidence_options(*, parser: argparse.ArgumentParser) -> None:
    """Add optional assay, functional and regulatory-reference inputs.

    Args:
        parser: Command parser.
    """
    parser.add_argument("--accessibility-bed", type=Path)
    parser.add_argument("--enhancer-bed", type=Path)
    parser.add_argument("--evidence-tsv", type=Path)
    parser.add_argument(
        "--reference-dir", type=Path, action="append", default=[]
    )
    parser.add_argument(
        "--organism", help="Scientific-name token, e.g. Homo_sapiens"
    )
    parser.add_argument(
        "--assembly", help="Genome build, e.g. GRCh38 or TAIR10"
    )


def add_motif_options(*, parser: argparse.ArgumentParser) -> None:
    """Add known-motif and exhaustive k-mer enrichment settings.

    Args:
        parser: Command parser.
    """
    parser.add_argument("--motifs", type=Path)
    parser.add_argument(
        "--motif-format",
        choices=("auto", "meme", "jaspar", "iupac"),
        default="auto",
    )
    parser.add_argument("--kmer-lengths", type=int, nargs="+", default=[])
    parser.add_argument("--forward-only", action="store_true")
    parser.add_argument("--site-p-value", type=float, default=1e-4)


def add_learning_options(*, parser: argparse.ArgumentParser) -> None:
    """Add reproducible, fixed-hyperparameter AI settings.

    Args:
        parser: Command parser.
    """
    parser.add_argument(
        "--ai-kmer-lengths", type=int, nargs="+", default=[4, 5, 6]
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--permutations", type=int, default=99)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-features", type=int, default=100000)
    parser.add_argument(
        "--regularisation",
        type=float,
        default=1.0,
        help="Fixed inverse L2 strength C; no validation-driven tuning",
    )
    parser.add_argument("--groups-tsv", type=Path)
    parser.add_argument(
        "--no-shap",
        action="store_true",
        help="Disable automatic held-out SHAP explanations and graphics",
    )
    parser.add_argument("--shap-max-sequences", type=int, default=1000)
    parser.add_argument("--shap-max-features", type=int, default=20)


def add_region_options(*, parser: argparse.ArgumentParser) -> None:
    """Add automatic multi-scale regulatory window search settings.

    Args:
        parser: Motif-analysis command parser.
    """
    parser.add_argument(
        "--no-regions",
        action="store_true",
        help="Disable automatic multi-scale regulatory region search",
    )
    parser.add_argument(
        "--region-lengths",
        type=int,
        nargs="+",
        default=[100, 200, 400, 800, 1600],
        help="Candidate window lengths, distinct from short motif lengths",
    )
    parser.add_argument("--region-step", type=int, default=50)
    parser.add_argument("--region-max-windows", type=int, default=50000)
    parser.add_argument("--region-score-threshold", type=float, default=0.75)
    parser.add_argument("--region-max-motifs", type=int, default=20)
    parser.add_argument("--region-motif-q-value", type=float, default=0.05)
    parser.add_argument("--region-position-bin-width", type=int, default=100)
    parser.add_argument("--region-position-limit", type=int, default=5000)
    parser.add_argument("--region-max-locus-plots", type=int, default=6)


def add_scan_options(*, parser: argparse.ArgumentParser) -> None:
    """Add explicit consensus-screen and gene-start profile settings.

    Args:
        parser: Genome-screen or pipeline command parser.
    """
    parser.add_argument(
        "--scan-mismatches",
        "--max-mismatches",
        type=int,
        default=0,
        help="Maximum substitutions, excluding indels",
    )
    parser.add_argument("--scan-q-value", type=float, default=0.05)
    parser.add_argument("--scan-max-motifs", type=int, default=20)
    parser.add_argument("--scan-motif-ids", nargs="+", default=[])
    parser.add_argument("--scan-intergenic-only", action="store_true")
    parser.add_argument("--scan-forward-only", action="store_true")
    parser.add_argument("--scan-mask-lowercase", action="store_true")
    parser.add_argument("--scan-chunk-size", type=int, default=250000)
    parser.add_argument("--scan-max-hits", type=int, default=2000000)
    parser.add_argument("--scan-upstream", type=int, default=2000)
    parser.add_argument("--scan-downstream", type=int, default=2000)
    parser.add_argument("--scan-bin-width", type=int, default=100)


def build_parser() -> argparse.ArgumentParser:
    """Build the complete CLI and its help text.

    Returns:
        The root argument parser.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Strictly intergenic, strand-aware extraction and "
            "regulatory sequence analysis"
        )
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"intergenic-regions {__version__}",
    )
    subcommands = parser.add_subparsers(dest="command", required=True)
    commands: dict[str, argparse.ArgumentParser] = {}
    for name, help_text in (
        ("extract", "Extract strictly intergenic flanks"),
        ("pipeline", "Gene lists to motifs, optional AI and evidence"),
        ("enrich", "Native motif enrichment from FASTA collections"),
        ("ai", "Learn and validate interpretable sequence signatures"),
        ("predict", "Score new candidates using an exported JSON model"),
        (
            "scan-genome",
            "Screen motif consensuses genome-wide with substitutions",
        ),
        (
            "annotate",
            "Attach optional assay/reference evidence to regions.tsv",
        ),
        ("homer", "Run HOMER with an explicit FASTA background"),
        ("import-reference", "Import a regulatory BED with assembly metadata"),
        ("references", "List supported regulatory-reference sources"),
    ):
        command = subcommands.add_parser(name=name, help=help_text)
        command.add_argument("--verbose", action="store_true")
        command.add_argument("--log-file", type=Path)
        if name != "references":
            command.add_argument(
                "--output-dir",
                "-o",
                type=Path,
                required=True,
                help="New directory; existing paths are refused",
            )
        commands[name] = command
    for name in ("extract", "pipeline"):
        add_extraction_options(parser=commands[name])
        add_evidence_options(parser=commands[name])
    commands["extract"].add_argument("--genes", type=Path)
    pipeline = commands["pipeline"]
    pipeline.add_argument("--positive-genes", type=Path, required=True)
    pipeline.add_argument("--negative-genes", type=Path, required=True)
    pipeline.add_argument("--scan-genome", action="store_true")
    add_scan_options(parser=pipeline)
    scan = commands["scan-genome"]
    scan.add_argument("--genome", type=Path, required=True)
    scan.add_argument("--annotation", "--gff", type=Path)
    scan.add_argument(
        "--annotation-format",
        choices=("auto", "gff3", "gtf", "tsv"),
        default="auto",
    )
    sources = scan.add_mutually_exclusive_group(required=True)
    sources.add_argument("--enrichment-tsv", type=Path)
    sources.add_argument("--motifs", type=Path)
    scan.add_argument(
        "--motif-format",
        choices=("auto", "meme", "jaspar", "iupac"),
        default="auto",
    )
    scan.add_argument("--positive-genes", type=Path)
    scan.add_argument("--negative-genes", type=Path)
    add_scan_options(parser=scan)
    for name in ("pipeline", "enrich"):
        controls = commands[name].add_mutually_exclusive_group()
        controls.add_argument(
            "--ai",
            dest="ai",
            action="store_true",
            default=True,
            help="Run ML alongside motifs (enabled by default)",
        )
        controls.add_argument(
            "--no-ml",
            dest="ai",
            action="store_false",
            help="Explicitly disable automatic sequence learning",
        )
    pipeline.add_argument(
        "--group-by", choices=("gene", "contig"), default="gene"
    )
    pipeline.add_argument("--match-background", action="store_true")
    pipeline.add_argument("--background-ratio", type=int, default=1)
    pipeline.add_argument("--max-gc-difference", type=float, default=0.1)
    pipeline.add_argument("--max-length-ratio", type=float, default=1.5)
    for name in ("enrich", "ai", "homer"):
        commands[name].add_argument(
            "--positive-fasta", type=Path, required=True
        )
        commands[name].add_argument(
            "--negative-fasta", type=Path, required=True
        )
    for name in ("enrich", "pipeline"):
        add_motif_options(parser=commands[name])
        add_region_options(parser=commands[name])
    commands["enrich"].add_argument("--mask-lowercase", action="store_true")
    for name in ("ai", "pipeline", "enrich"):
        add_learning_options(parser=commands[name])
    commands["ai"].add_argument("--candidate-fasta", type=Path)
    commands["predict"].add_argument("--model", type=Path, required=True)
    commands["predict"].add_argument("--fasta", type=Path, required=True)
    annotate = commands["annotate"]
    annotate.add_argument("--regions-tsv", type=Path, required=True)
    add_evidence_options(parser=annotate)
    homer = commands["homer"]
    homer.add_argument("--executable", default="findMotifs.pl")
    homer.add_argument(
        "--motif-lengths", type=int, nargs="+", default=[8, 10, 12]
    )
    homer.add_argument("--threads", type=int, default=1)
    homer.add_argument("--known-motifs", type=Path)
    homer.add_argument("--known-only", action="store_true")
    homer.add_argument("--timeout", type=float, default=3600)
    homer.add_argument("--dry-run", action="store_true")
    reference = commands["import-reference"]
    from intergenic_regions.references import CATALOGUE

    reference.add_argument(
        "--source", choices=tuple(CATALOGUE), default="custom"
    )
    source = reference.add_mutually_exclusive_group()
    source.add_argument("--bed", type=Path)
    source.add_argument("--url")
    reference.add_argument("--organism")
    reference.add_argument("--assembly")
    reference.add_argument("--name")
    reference.add_argument("--expected-sha256")
    return parser


def extraction_options(*, args: argparse.Namespace) -> dict[str, Any]:
    """Translate parsed extraction options to explicit API settings.

    Args:
        args: Parsed command options.

    Returns:
        Extraction keyword arguments.
    """
    return {
        "direction": args.direction,
        "length": None if args.full_gap else args.length,
        "offset": args.offset,
        "min_length": args.min_length,
        "max_ambiguous_fraction": args.max_ambiguous_fraction,
        "mask_lowercase": args.mask_lowercase,
    }


def evidence_options(*, args: argparse.Namespace) -> dict[str, Any]:
    """Translate optional evidence flags to API settings.

    Args:
        args: Parsed command options.

    Returns:
        Evidence keyword arguments.
    """
    return {
        "accessibility_bed": args.accessibility_bed,
        "enhancer_bed": args.enhancer_bed,
        "evidence_tsv": args.evidence_tsv,
        "references": args.reference_dir,
        "organism": args.organism,
        "assembly": args.assembly,
    }


def motif_options(*, args: argparse.Namespace) -> dict[str, Any]:
    """Translate motif flags to native analysis settings.

    Args:
        args: Parsed command options.

    Returns:
        Native enrichment keyword arguments.
    """
    return {
        "motif_format": args.motif_format,
        "lengths": args.kmer_lengths,
        "both_strands": not args.forward_only,
        "site_p_value": args.site_p_value,
    }


def learning_options(*, args: argparse.Namespace) -> dict[str, Any]:
    """Translate AI flags to fixed model settings.

    Args:
        args: Parsed command options.

    Returns:
        Sequence model keyword arguments.
    """
    return {
        "lengths": args.ai_kmer_lengths,
        "folds": args.folds,
        "seed": args.seed,
        "permutations": args.permutations,
        "max_features": args.max_features,
        "regularisation": args.regularisation,
        "shap": not args.no_shap,
        "shap_max_sequences": args.shap_max_sequences,
        "shap_max_features": args.shap_max_features,
    }


def scan_options(*, args: argparse.Namespace) -> dict[str, Any]:
    """Translate genome-screen flags to selection and scan settings.

    Args:
        args: Parsed command options.

    Returns:
        Validated workflow keyword settings.
    """
    return {
        "max_mismatches": args.scan_mismatches,
        "q_threshold": args.scan_q_value,
        "max_motifs": args.scan_max_motifs,
        "motif_ids": args.scan_motif_ids,
        "both_strands": not args.scan_forward_only,
        "intergenic_only": args.scan_intergenic_only,
        "mask_lowercase": args.scan_mask_lowercase,
        "chunk_size": args.scan_chunk_size,
        "max_hits": args.scan_max_hits,
        "upstream": args.scan_upstream,
        "downstream": args.scan_downstream,
        "bin_width": args.scan_bin_width,
    }


def region_options(*, args: argparse.Namespace) -> dict[str, Any]:
    """Translate candidate-window flags into multi-scale search settings.

    Args:
        args: Parsed command options.

    Returns:
        Window-generation, triage and reporting keyword arguments.
    """
    return {
        "lengths": args.region_lengths,
        "step": args.region_step,
        "max_windows": args.region_max_windows,
        "score_threshold": args.region_score_threshold,
        "max_motifs": args.region_max_motifs,
        "motif_q_value": args.region_motif_q_value,
        "position_bin_width": args.region_position_bin_width,
        "position_limit": args.region_position_limit,
        "max_locus_plots": args.region_max_locus_plots,
    }


def dispatch(*, args: argparse.Namespace) -> dict[str, Any]:
    """Execute a parsed command.

    Args:
        args: Valid parsed options.

    Returns:
        A JSON-compatible command summary.

    Raises:
        ValueError: Input data or analysis settings are invalid.
    """
    if args.command == "references":
        from intergenic_regions.references import CATALOGUE

        return {
            "sources": CATALOGUE,
            "catalogue_version": "2026-10-01",
            "note": (
                "Preset datasets retain their documented build. Custom and "
                "plant tracks require explicit organism/assembly and a "
                "selected BED."
            ),
        }
    output = args.output_dir
    if args.command == "extract":
        return extract_workflow(
            genome_path=args.genome,
            annotation_path=args.annotation,
            output=output,
            identifiers_path=args.genes,
            annotation_format=args.annotation_format,
            extraction_settings=extraction_options(args=args),
            evidence_settings=evidence_options(args=args),
        )
    if args.command == "pipeline":
        matching = (
            {
                "ratio": args.background_ratio,
                "max_gc_difference": args.max_gc_difference,
                "max_length_ratio": args.max_length_ratio,
            }
            if args.match_background
            else None
        )
        return pipeline_workflow(
            genome_path=args.genome,
            annotation_path=args.annotation,
            positive_genes=args.positive_genes,
            negative_genes=args.negative_genes,
            output=output,
            annotation_format=args.annotation_format,
            extraction_settings=extraction_options(args=args),
            motif_path=args.motifs,
            enrichment_settings=motif_options(args=args),
            use_ai=args.ai,
            learning_settings=learning_options(args=args),
            groups_path=args.groups_tsv,
            group_by=args.group_by,
            matching_settings=matching,
            evidence_settings=evidence_options(args=args),
            scan_genome=args.scan_genome,
            scan_settings=scan_options(args=args),
            use_regions=not args.no_regions,
            regional_settings=region_options(args=args),
        )
    if args.command == "enrich":
        return enrichment_workflow(
            positive_path=args.positive_fasta,
            negative_path=args.negative_fasta,
            output=output,
            motif_path=args.motifs,
            settings=motif_options(args=args),
            mask_lowercase=args.mask_lowercase,
            use_ai=args.ai,
            learning_settings=learning_options(args=args),
            groups_path=args.groups_tsv,
            use_regions=not args.no_regions,
            regional_settings=region_options(args=args),
        )
    if args.command == "ai":
        return learning_workflow(
            positive_path=args.positive_fasta,
            negative_path=args.negative_fasta,
            output=output,
            settings=learning_options(args=args),
            groups_path=args.groups_tsv,
            candidates_path=args.candidate_fasta,
        )
    if args.command == "scan-genome":
        return genome_scan_workflow(
            genome_path=args.genome,
            output=output,
            annotation_path=args.annotation,
            annotation_format=args.annotation_format,
            enrichment_path=args.enrichment_tsv,
            motif_path=args.motifs,
            positive_genes=args.positive_genes,
            negative_genes=args.negative_genes,
            settings={
                "motif_format": args.motif_format,
                **scan_options(args=args),
            },
        )
    if args.command == "predict":
        from intergenic_regions.learning import predict_sequences

        model = json.loads(args.model.read_text(encoding="utf-8"))
        rows = predict_sequences(
            model=model, sequences=read_fasta(path=args.fasta)
        )
        with output_bundle(path=output) as stage:
            write_tsv(
                path=stage / "candidate_scores.tsv",
                rows=rows,
                fields=tuple(rows[0]),
            )
            write_json(
                path=stage / "manifest.json",
                data=provenance(
                    inputs=[args.model, args.fasta],
                    settings={"command": "predict"},
                ),
            )
            write_report(
                path=stage / "report.html",
                title="Candidate sequence signature scores",
                summary={"candidates": len(rows)},
                tables={"Candidates": rows},
                notes=[
                    "Candidate scores are positive-class sequence signatures, "
                    "not enhancer probabilities. Exact training duplicates "
                    "are flagged."
                ],
            )
        return {"candidates": len(rows)}
    if args.command == "annotate":
        with output_bundle(path=output) as stage:
            evidence, references = attach_evidence(
                directory=stage,
                regions=read_region_table(path=args.regions_tsv),
                **evidence_options(args=args),
            )
            write_report(
                path=stage / "report.html",
                title="Optional regulatory evidence",
                summary={
                    "regions": len(evidence),
                    "reference_queries": len(references),
                },
                tables={
                    "Supplied experimental evidence": evidence,
                    "Regulatory reference overlaps": references,
                },
                notes=[
                    "Linked gene-level evidence is distinct from experimental "
                    "evidence on the extracted genomic interval."
                ],
            )
        return {"regions": len(evidence), "reference_queries": len(references)}
    if args.command == "import-reference":
        from intergenic_regions.references import import_reference

        return import_reference(
            output=output,
            source=args.source,
            local_bed=args.bed,
            url=args.url,
            organism=args.organism,
            assembly=args.assembly,
            name=args.name,
            expected_sha256=args.expected_sha256,
        )
    if args.command == "homer":
        from intergenic_regions.homer import run_homer

        return run_homer(
            positive=args.positive_fasta,
            negative=args.negative_fasta,
            output=output,
            executable=args.executable,
            lengths=tuple(args.motif_lengths),
            threads=args.threads,
            known_motifs=args.known_motifs,
            known_only=args.known_only,
            timeout=args.timeout,
            dry_run=args.dry_run,
        )
    raise ValueError(f"Unsupported command: {args.command}")


def main(*, argv: Sequence[str] | None = None) -> int:
    """Run the CLI with reliable exit status and stderr diagnostics.

    Args:
        argv: Optional explicit arguments, excluding the program name.

    Returns:
        Zero on success, two for invalid inputs/failed operations, 130 when
        interrupted. argparse handles help/version and invalid flags itself.
    """
    parser = build_parser()
    args = parser.parse_args(args=argv)
    handlers: list[logging.Handler] = [
        logging.StreamHandler(stream=sys.stderr)
    ]
    if args.log_file:
        try:
            handlers.append(
                logging.FileHandler(filename=args.log_file, encoding="utf-8")
            )
        except OSError as exc:
            print(f"Cannot open log file: {exc}", file=sys.stderr)
            return 2
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=handlers,
        force=True,
    )
    try:
        result = dispatch(args=args)
    except ImportError as exc:
        LOGGER.error(
            "Analysis dependency missing; install "
            "intergenic-regions[analysis]: %s",
            exc,
        )
        return 2
    except (OSError, ValueError, RuntimeError) as exc:
        LOGGER.error("%s", exc, exc_info=args.verbose)
        return 2
    except KeyboardInterrupt:
        LOGGER.warning(
            "Interrupted; no incomplete output bundle was published"
        )
        return 130
    print(json.dumps(obj=result, sort_keys=True, allow_nan=False))
    return 0
