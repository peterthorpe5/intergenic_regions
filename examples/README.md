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

## Variable-length planted region demonstration

```bash
python examples/make_demo.py \
  --output-dir results/multiscale_inputs \
  --flank-length 2400 --module-lengths 120 360 720 \
  --sequences-per-class 12
python -m intergenic_regions pipeline \
  --genome results/multiscale_inputs/genome.fasta \
  --annotation results/multiscale_inputs/genes.gff3 \
  --positive-genes results/multiscale_inputs/positive_genes.txt \
  --negative-genes results/multiscale_inputs/negative_genes.txt \
  --groups-tsv results/multiscale_inputs/groups.tsv \
  --motifs results/multiscale_inputs/motifs.tsv \
  --accessibility-bed results/multiscale_inputs/accessibility.bed \
  --full-gap --region-lengths 100 200 400 800 1600 --region-step 100 \
  --folds 3 --permutations 19 \
  --output-dir results/multiscale_analysis
```

The fixture contains motif clusters of several lengths at varying source
positions. Negative controls carry a same-composition alternate pattern.
`synthetic_modules.tsv` records planted source offsets; these are software
fixtures, not measured enhancer intervals. Neither model validation nor
agreement with their spans establishes biological enhancer detection.
Open the root and `regions/report.html`, plus individual scale SHAP reports.
