# Verification

Run pytest from the repository root after installing .[dev].
The suite uses small deterministic datasets and no mandatory network calls.

| Tests | Protection |
| --- | --- |
| test_extraction.py | Strand orientation, boundaries, offsets, nested/overlapping genes, contig ends, masks, per-base property oracle |
| test_annotation.py | GFF3/GTF/legacy parsing, unsorted hierarchy, UTRs, shared/missing parents, escapes, discontinuous IDs, invalid data |
| test_demo_legacy_api.py | Original 2020 fixture, public imports and reproducible example generator |
| test_io_models.py | Immutable model/ID and text/FASTA/TSV/JSON contracts, hashes, atomic bundles |
| test_background_evidence.py | Shared-promoter leakage, duplicates, GC/length matching, evidence union/unknown/negative cases |
| test_motifs_statistics.py | Motif formats, known/exact/PWM scans, hypothesis family, Fisher/FDR checks, planted motifs |
| test_learning.py | Fold-only vocabulary/scaling, grouped splits, null labels, permutations, inert model schema and prediction |
| test_shap_explanations.py | Agreement with official LinearExplainer, fold backgrounds, subsampling, additive reconstruction, official PNG/PDF plots and missing plotting dependency |
| test_genome_scanning.py | Brute-force matching property oracle, substitutions/IUPAC/strands, chunk boundaries, intergenic blockers, signed distances, tied starts, denominator heatmaps, hit caps and atomic outputs |
| test_automatic_ml_prioritisation.py | Automatic ML defaults, unsupported data, transparent ranks, uncertainty and HTML interactivity |
| test_workflows_cli_reporting.py | Focused workflows, plots, escaping, option handling, failures and end-to-end subprocess analysis |
| test_references_homer.py | Reference build/checksum validation, controlled downloads, HOMER arguments/errors/timeouts |

Every public function is exercised directly or through a focused workflow test.
The property oracle independently walks occupied bases rather than copying the
binary-search implementation. Legacy off-by-one output is corrected explicitly:
a neighbour's first base must never be part of the returned sequence.

Tests require no live HOMER installation or downloaded enhancer database.
Those integrations use controlled executables/mocks; real external integration
availability depends on the user's environment.

Optional browser verification requires Node.js, Playwright and Chromium:

```bash
node tests/check_report.cjs \
  --report results/demo_analysis/report.html \
  --screenshot results/dashboard.png \
  --result-json results/browser_checks.json
```

It checks offline rendering, embedded images, table filtering/reset, numerical
sorting in both directions, mobile overflow and JavaScript errors. An explicit
--browser-executable or --chromium-module can select an installed alternative.
