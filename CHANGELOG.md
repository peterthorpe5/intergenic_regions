# Changelog

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
