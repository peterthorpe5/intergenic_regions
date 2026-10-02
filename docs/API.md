# Python API

Project-defined functions use keyword-only arguments, type annotations and
Google-style docstrings. Paths are pathlib.Path objects. Data model constructors
also use named arguments. Sequence analysis dependencies are imported only when
the relevant modules/workflows are used.

| Module | Public interface |
| --- | --- |
| Package root | Gene, Region, Genome, GeneIndex, read_annotation, extract_regions, reverse_complement |
| annotation | parse_attributes, parse_feature, resolve_roots, aggregate_gene, read_coordinates, read_annotation |
| extraction | merge_intervals, GeneIndex.free_flank, extract_regions |
| genome | Genome context manager/close/fetch/lengths, normalise_genome, reverse_complement, sequence_composition |
| io | text/FASTA/ID readers, FASTA/TSV/JSON writers, output_bundle, file_fingerprint, provenance |
| motifs | motif readers, exact/PWM scanners, kmer_counts, pooled_background, prepare_pwm, analyse_motifs |
| statistics | enrichment_test, adjust_fdr, kmer_family_size |
| background | canonical_sequence, check_sequence_sets, collapse_overlaps, sanitise_regions, match_background |
| evidence | PeakIndex, read_bed, read_functional_evidence, annotate_evidence |
| learning | group reader, features, validation/model fitting, permutations, fit_sequence_model, predict_sequences |
| explanations | explanation_sample, linear_shap_values, summarise_shap |
| shap_reporting | shap_explanation, plot_shap (official SHAP plots) |
| scanning | ScanMotif, read_scan_motifs, consensus_sites, genome_scan_outputs |
| scan_reporting | plot_scan, scan_report |
| windows | RegulatoryWindow, make_windows, assign_window_groups, merge_window_candidates |
| window_learning | parent_window_weights, fit_window_model |
| regional_analysis | window_motif_counts, window_position_profiles, regional_outputs |
| window_reporting | plot_regional_results |
| prioritisation | prioritise_candidates |
| references | canonical_assembly, download_bed, import_reference, query_references |
| homer | homer_command, run_homer |
| reporting | plot helpers, save_figure, write_report |
| workflows | extract_workflow, enrichment_workflow, learning_workflow, pipeline_workflow and output/report helpers |
| cli | build_parser, option helpers, dispatch, main |

## Complete workflow

```python
from pathlib import Path
from intergenic_regions.workflows import pipeline_workflow

summary = pipeline_workflow(
    genome_path=Path("genome.fasta"),
    annotation_path=Path("genes.gff3"),
    positive_genes=Path("positive.txt"),
    negative_genes=Path("negative.txt"),
    output=Path("results/run_1"),
    extraction_settings={"length": 1000, "min_length": 50},
    enrichment_settings={"lengths": [5, 6]},
    learning_settings={"folds": 5, "permutations": 99, "seed": 42},
    evidence_settings={"accessibility_bed": Path("atac.bed")},
)
```

Omit evidence_settings for sequence-only analysis. ML is enabled by default;
use_ai=False is an explicit opt-out. enrichment_workflow has the same default
and accepts learning_settings/groups_path. fit_sequence_model can be used
directly when arrays/metrics are required rather than files.
SHAP is enabled by default. `learning_settings={"shap": False}` is the opt-out;
`shap_max_sequences` and `shap_max_features` bound explanation output. The
five-value return contract of `fit_sequence_model` is unchanged. Supply a mutable
`explanation_outputs={}` to receive long rows, all-feature importance and a
method summary in that output dictionary.

## Genome screening

```python
from intergenic_regions.workflows import genome_scan_workflow

summary = genome_scan_workflow(
    genome_path=Path("genome.fasta"),
    annotation_path=Path("genes.gff3"),
    enrichment_path=Path("results/run_1/motifs/motif_enrichment.tsv"),
    output=Path("results/genome_screen"),
    settings={
        "max_mismatches": 1,
        "max_motifs": 20,
        "intergenic_only": True,
        "upstream": 2000,
        "downstream": 500,
        "bin_width": 100,
    },
)
```

Alternatively, `pipeline_workflow(scan_genome=True, scan_settings={...})`
screens newly enriched motifs in the same bundle. `motif_path` replaces
`enrichment_path` for supplied patterns/PWMs. `scan_targets` separates selection
from scan settings. `consensus_sites(sequence="ACGT", pattern="ARG",
max_mismatches=1)` is the in-memory overlapping-site API.

`genome_scan_outputs` is a low-level streaming writer: callers must provide
their own `output_bundle` context for atomic publication. Public workflows
already do so. Full site tables stream to disk; preview rows are capped at 500.
Hits are physical windows, merged when both orientations match; the matching
method is consensus substitutions, including for explicitly converted PWMs.

## Data and failure contracts

Gene spans are non-empty; coordinates are zero-based, half-open. Region records
include audited exclusions: sequence is empty for excluded flanks. A region
table reloaded for evidence annotation has coordinates but no bases. sequence_id
escapes reserved gene-ID characters and includes the direction.

Genome must be used inside a context manager. It keeps its temporary index
private and closes handles on failure. Reusing a closed accessor raises an error.
All genes must belong to the indexed assembly and fit their contig lengths.
The target list never controls which genes block extraction.

Malformed data generally raise ValueError, unavailable files raise OSError
subclasses, and failed external operations raise RuntimeError. CLI main turns
these into logged non-zero exits. output_bundle stages a complete new directory
and refuses overwrite; its context manager removes partial output on failure.

No root/API import requires scipy, matplotlib or scikit-learn. Direct motif
imports require [motifs], direct learning imports [ml]. Workflow defaults attempt
ML and report unavailable/not_estimable status explicitly when appropriate.
Numerical explanations require [ml]; official SHAP graphics require [analysis].
The root extraction API and CLI help remain usable without SHAP or NumPy.

## Multi-scale regulatory windows

`pipeline_workflow` and `enrichment_workflow` accept `use_regions=False` to opt
out, otherwise they run regional search automatically. Configure it with
`regional_settings={"lengths": [100, 200, 400, 800, 1600, 3200], "step": 50}`.
Use `extraction_settings={"length": None}` for complete safe gaps, or supply
a larger maximum length. Existing short-word lengths remain separate settings.

```python
from intergenic_regions.windows import make_windows, assign_window_groups
from intergenic_regions.window_learning import fit_window_model

windows, availability = make_windows(
    positive=positive_sequences,
    negative=negative_sequences,
    lengths=[100, 200, 400, 800],
    step=50,
)
groups = assign_window_groups(windows=windows)
predictions, parents, metrics, summary, model, explanations = fit_window_model(
    windows=[w for w in windows if w.width == 200],
    parent_groups=groups,
    lengths=[4, 5, 6],
    folds=5,
)
```

Supply matching `regions` and full `genes` to `make_windows` for exact genomic
coordinates and signed annotated-start distances. FASTA-only windows retain
oriented source offsets, with genomic coordinates/distance None. Too-short
sources receive audit rows and no truncated windows. `assign_window_groups`
unions optional parent-family groups and identical/reverse-complement windows
across every requested scale. Partial homology still needs user grouping.

`fit_window_model` fits exactly one scale, with parent/class-balanced weights
and source-level fold metrics. Its six-value return includes numerical SHAP
separately from the reusable model. It inherits source labels and is not a
model trained on measured enhancer intervals. Its JSON is compatible with
`predict_sequences`, which checks the additional `window_length` field.

`regional_outputs` writes low-level results inside an existing staged output
directory; public workflows provide atomic `output_bundle` publication.
`merge_window_candidates` joins overlapping scored windows as exploratory
unions. Neither their boundaries nor their peak scores have region-level FDR.
