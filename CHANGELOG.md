# Changelog

## 1.3.0

- RNA-seq interrogation with curated positive and negative gene lists as primary
  inputs; optional fitted/imported DE context never replaces their labels.
- Bounded BAM/CRAM or bedGraph coverage with explicit library conventions,
  spliced-block/paired-fragment depth and independent biological replicate calls.
- Descriptive sustained coverage rises, replicate-supported neighbourhoods and
  condition shifts, preserving uncertainty and the original safe DNA boundaries.
- Optional raw FASTQ mapping via controlled two-pass STAR/samtools and raw exon
  counts via featureCounts; explicit threads/RAM/time limits, logs and versions.
- Optional PyDESeq2 raw-count DE with replicate/design validation, transparent
  screening sets, missing-test handling and recorded model warnings.
- Automatic existing sequence ML/SHAP and multi-scale analysis for RNA-defined
  comparisons, excluding RNA quantities used to define labels from predictors.
- RNA coverage/support/locus/motif-distance and optional DE graphics, offline
  reports with interpretation notes and embedded PNG/vector PDF downloads.
- Separate Conda/pip RNA environment, input schemas, TransStart migration notes,
  reproducible synthetic fixtures and biological/API/external-tool regression tests.
- Exactly stratified splits for singleton independent-gene groups, preventing
  heuristic class loss during permutation validation in the Conda environment.

RNA starts, promoter hypotheses and enhancer candidates remain unvalidated.
Sequence-only workflows and the original strand-aware blocker behaviour remain
compatible with version 1.2.

## 1.2.0

Multi-scale regulatory search runs automatically with `pipeline` and `enrich`.
Separate search-window lengths from motif/AI word lengths; default scales are
100, 200, 400, 800 and 1,600 bp. Safe flank coordinates, negative-strand sequence
orientation, end-aligned coverage and an explicit per-scale availability audit
preserve the original biological contract.

Per-scale weak-label models hold out complete source/group units, union exact
shared windows across scales, balance training mass by parent and report
source-level metrics alongside a composition baseline. Exact held-out SHAP
and official graphics run automatically when available. Variable-length window
unions carry peak/scale provenance, local optional evidence, FASTA/BED/TSV
exports and explicit uncertainty, without region-level significance claims.

HTML adds scale coverage/validation, distance and density heatmaps, position
coverage, candidate span distributions and locus views. Allocation/plot/profile
limits are configurable; `--no-regions` disables the new layer. Existing
whole-flank models, motif tests, genome motif screening and extraction commands
retain their existing interfaces. Window models enforce their training length
when reused through `predict`.

Add a Conda environment file for the compiled analysis dependencies and
regression/property tests covering every new public function.
