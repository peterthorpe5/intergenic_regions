# Output reference

All package-owned tabular files are TSV. Null evidence/undefined values appear
as empty fields (HTML displays NA). JSON is strict: no NaN or Infinity.
HTML embeds PNG/PDF figure downloads, CSS and JavaScript; no external runtime
is required. Data/model links remain relative to the companion folder.

| File | Meaning |
| --- | --- |
| report.html | Unified overview for pipeline/enrich; searchable/sortable previews up to 500 rows per table |
| summary.json | Workflow summary and explicit ML status |
| manifest.json | Package/Python/software versions, settings, SHA-256 input fingerprints |
| regions.fasta | Retained, transcription-oriented sequence |
| regions.tsv | Every requested flank, coordinate/status/boundary/composition audit |
| regions.bed | Retained half-open genomic intervals |
| regions.gff3 | Retained one-based inclusive genomic intervals |
| positive.fasta / negative.fasta | Pipeline sequence representatives actually used in analysis |
| background_audit.tsv | Removed overlaps and exact/reverse-complement duplicates |
| background_matching.tsv | Optional matched controls, GC difference and length ratio |
| evidence.tsv | Accessibility/regulatory overlaps and linked functional evidence |
| reference_overlaps.tsv | Imported-reference source/build/checksum and overlap |
| candidate_priorities.tsv | Ranked candidates, evidence tiers, held-out scores, uncertainty |
| candidate_summary.tsv | Pipeline companion with the same detailed ranked records |
| motifs/motif_enrichment.tsv | Known motifs/words, presence counts, effects, p/q-values |
| motifs/motif_sites.tsv | Known motif occurrences in sequence-relative coordinates |
| motifs/figures/*.png / *.pdf | Motif prevalence, enrichment, logos and composition |
| ai/held_out_predictions.tsv | Scores for each training observation from its held-out fold |
| ai/sequence_signatures.tsv | Full-data model word coefficients and fold stability |
| ai/fold_metrics.tsv | Sample counts, vocabulary size, ROC/AP and baseline per fold |
| ai/model.json | Reusable inert model, weights/scaling/training fingerprints |
| ai/summary.json | Model metrics or an explicit disabled/unavailable/not_estimable status |
| ai/figures/*.png / *.pdf | Held-out ROC/PR, signatures and permutation results |
| ai/shap_values.tsv | Leading fold-specific contributions plus an additive remainder per explained sequence; base value, logit, score and scope |
| ai/shap_importance.tsv | Mean absolute/signed SHAP across all features and the explained sample; fold availability counts |
| ai/figures/shap_*.png / *.pdf | Official SHAP bar, beeswarm, heatmap, highest/lowest-score waterfalls and leading-feature dependence |
| genome_scan/report.html | Offline whole-genome screen with distance, gene burden and mismatch graphics |
| genome_scan/genome_motif_sites.tsv / *.bed | Complete streaming physical hits, substitutions, reference-strand bases, nearest-gene proxy and explicit uncertainty |
| genome_scan/scan_motifs.tsv | Selected consensus patterns, source kind and optional discovery q-value |
| genome_scan/distance_profiles.tsv | Per-motif/cohort signed distance bins, physical counts, eligible windows and density per million |
| genome_scan/scan_burden.tsv | Counts by motif, contig and minimum substitution number |
| genome_scan/gene_motif_counts.tsv | Nearest-gene physical motif burden; tied sites use a stable representative |
| genome_scan/summary.json | Screen method/settings/counts and TSS-proxy interpretation |
| regions/report.html | Offline multi-scale region dashboard and links to scale-specific SHAP reports |
| regions/window_availability.tsv | Complete-window counts and too-short status per source/scale |
| regions/window_scores.tsv | All window offsets, optional genomic coordinates, density, held-out scores and local optional evidence |
| regions/windows.bed | All windows with genomic coordinates; empty for FASTA-only enrich |
| regions/parent_predictions.tsv | One mean/peak held-out score per source/scale, group and fold |
| regions/fold_metrics.tsv | Source counts and source-class validation metrics at each scale |
| regions/candidate_regions.tsv | Ranked window unions, peak/scale support, evidence and explicit uncertainty; region_q_value is empty |
| regions/candidate_regions.bed / *.fasta | Genomic union intervals and transcription-oriented candidate sequences |
| regions/position_profiles.tsv | Window-length/position bins, parent coverage and parent-weighted resemblance/density means |
| regions/scale_N/model.json | Reusable weak-label window model, enforcing its training length |
| regions/scale_N/shap_values.tsv / shap_importance.tsv | Exact held-out window explanations against weighted training-only backgrounds |
| regions/scale_N/report.html / figures/ | Official SHAP plots and interpretation when the scale is estimable |
| regions/evidence.tsv / reference_overlaps.tsv | Evidence re-queried at each window's own genomic coordinates |
| regions/candidate_evidence/ | Evidence re-queried at each union's own genomic coordinates |
| regions/figures/region_*.png / *.pdf | Coverage, source validation, distance/density heatmaps, union spans and individual locus plots |

In standalone enrich output, motif tables/figures are at the root rather than
under motifs/. Prediction-only output uses candidate_scores.tsv. Extraction
bundles include evidence tables even when no evidence was supplied.
Standalone scan-genome writes its scan files at the root rather than under
genome_scan/. Its preview contains at most 500 hits; the TSV/BED remains complete.

The top-level report links full data; table previews are not an alternative
to the TSV when analysing hundreds or thousands of candidates. Image captions
identify plots; PNG and PDF copies remain independently reusable.

## Status interpretation

retained: extracted and passed filters. unknown_strand: no orientation.
no_intergenic_space: the adjacent flank is occupied or touches a contig end.
offset_exceeds_gap: requested skip consumes the available gap.
below_minimum_length / excess_ambiguity: explicit filtering reasons.

stop_reason names the limiting boundary of the entire available gap; a maximum
length may return a shorter piece before that boundary. available_length is
before offset and retained-length filtering.

ML completed means internal held-out modelling finished; it does not indicate
external biological validation. not_estimable/unavailable/disabled report why
no model scores were produced. Existing result paths are never reused.

Multi-scale `completed` means windows were searched; inspect the separate status
for every model scale. `no_complete_windows` means every requested window was
longer than its source. `parent_too_short` is a per-source availability status.
Disabling ML retains density/availability output; candidate unions require
estimable scores. Scores unavailable at a scale are empty/NA, never zero.
Position profiles use annotated 5-prime bases as TSS proxies or oriented FASTA
offsets, with an explicit coordinate-system field. The configured profile limit
clips these summaries only; full window/candidate tables remain complete.
