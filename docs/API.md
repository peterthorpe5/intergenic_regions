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
