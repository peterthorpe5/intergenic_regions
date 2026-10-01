# Reproducible example

Run make_demo.py with --output-dir, as shown in the root README. It produces
both positive- and negative-strand target genes and explicit flanking blockers.
Random sequence, planted G-box motifs and all assay/reference metadata are
synthetic. --seed and --sequences-per-class control reproducibility/sample size.

Expected flank length is 240 bases for every target, even with --length 1000,
because the blockers truncate extraction. Positive promoters have planted
CACGTG occurrences. Reverse-complemented minus-strand promoters must match the
supplied oriented positive.fasta/negative.fasta exactly; a unit test verifies it.

The pipeline example in the README runs motifs plus automatic ML and produces
report.html. Results demonstrate functioning software, not biological validity.
All options are named and output paths must be new.
The version 1.1.0 README example also runs automatic held-out SHAP and a genome
screen with one substitution. Open the unified dashboard and the linked scan
report; inspect exact counts/denominators in distance_profiles.tsv. Official
SHAP plots and position/gene heatmaps have PNG/PDF downloads embedded in HTML.
