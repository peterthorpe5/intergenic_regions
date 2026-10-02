# Migration from the 2020 script

The starting point is master commit c36293fc8f40117a8401a3e42aba255e172c914f.
Version 1.0.0 is a major interface change. The old modules and nose tests are
replaced; retain an old checkout if a historical pipeline requires that API.

| Old interface | New interface |
| --- | --- |
| `python intergenic_regions.py` | `python -m intergenic_regions extract` or `intergenic-regions extract` |
| `--gff file.gff` | `--annotation file.gff3`; `--gff` remains an alias |
| `-g` / `--genome` | Same named options |
| `-u` / `--upstream` | Aliases for `--length`; direction defaults to upstream |
| `-o output.fasta` | `--output-dir new_directory` / `-o new_directory` |
| `-m` / `--min_len` | `--min-length` / `-m` |
| `-z` / `--user_defined_genic` | Removed: this package guarantees strictly intergenic output |
| Space-separated coordinate tables | Tab-separated five-column input |
| nose | pytest, branch coverage, property tests and subprocess workflows |

The output is a directory: sequences are in `regions.fasta`, with a complete
`regions.tsv` audit, BED/GFF3 coordinates, evidence, input checksums and HTML.
Output directories must not exist. This prevents silently overwriting results.

## Deliberate corrections

GFF3 coordinates are converted from one-based inclusive to zero-based half-open.
No first base of a neighbouring gene is returned, including on the negative
strand. The original simplified fixture's gap between positions 50 and 55
contains positions 51–54, four bases: a five-base request must return four,
not five. A dedicated regression uses this original fixture.

The old minimum length excluded lengths at or below the threshold. The new
threshold is inclusive, default 1. To reproduce old `--min_len 3` filtering,
use `--min-length 4`.

The old parser relied on CDS-derived coordinates and simplified identifiers.
The new parser uses the complete gene/descendant span, including UTRs/introns,
and preserves exact suffixes. This can shorten a flank or shift its anchor when
the original tool started at a coding boundary instead of the gene boundary.
For historical CDS-boundary reproducibility, supply an explicit five-column
table; its coordinates become the declared gene spans. Only those declared
spans can then block extraction, so this is not equivalent to a full annotation.

Input IDs must match the parsed gene IDs exactly. Do not automatically strip
`.1`, `.t1` or similar suffixes. GFF3 Parent hierarchies determine gene roots;
legacy five-column rows preserve their identifiers.

## Updated example

```bash
python -m intergenic_regions extract \
  --gff genes.gff3 --genome genome.fasta \
  --upstream 1000 --min-length 4 --output-dir results/upstream
```

Flanks are strand-oriented. Coordinates still refer to the original genomic
interval; they do not reverse numerical ordering on the minus strand.

For new motif analysis use `pipeline` or `enrich`: ML runs by default.
There is no requirement to run a second command. Existing `--ai` flags are
accepted. A small or ungroupable dataset keeps its motif results and explicitly
reports why the model cannot be estimated.

The source archive is a complete package, not a legacy-file overlay. Prefer
applying the supplied Git bundle/patch so removed scripts really are removed.

## Version 1.0.0 to 1.1.0

Existing extraction/enrichment/model APIs and named options remain compatible.
Reinstall `python -m pip install --editable '.[analysis]'` to obtain SHAP.
Motif workflows now generate exact held-out SHAP automatically when ML is
estimable. `--no-shap` disables this addition; `--no-ml` disables the model.
Old schema-version-1 JSON models still predict normally. New models additionally
record a full-training feature background mean; held-out SHAP never uses it.

`--scan-genome` adds a consensus screen to `pipeline`; `scan-genome` reuses
enrichment TSVs or supplied motifs independently. Genome screening is opt-in.
Substitutions default to zero; `--scan-mismatches 1` allows one. Extraction's
strict boundary behaviour and native discovery statistics are unchanged.

Figures now provide embedded PNG/PDF downloads in HTML. TSV/BED/model links
still require the companion output folder. New SHAP/scan TSVs are additional
outputs; they do not change candidate uncertainty or evidence tiers.

## Version 1.1.0 to 1.2.0

Motif workflows now also run multi-scale region analysis automatically. Existing
commands still work; `--no-regions` restores the earlier output scope. New
results live in `regions/`. The original whole-flank candidate files, motif
q-values and genome-screening interface retain their existing meanings.

Window lengths are independent of `--kmer-lengths` and `--ai-kmer-lengths`.
The default flank cap remains 1,000 bp: larger scales are unavailable unless
you increase `--length` or select `--full-gap`. Neighbouring genes still limit
all flanks. Short gaps are audited, without padding or genic extension.

Activate your working Conda environment and reinstall:
`python -m pip install --only-binary=:all: --editable '.[dev]'`.
For a new environment, `conda env create --file environment.yml` installs
compiled analysis dependencies through Conda. Historical schema-version-1
whole-flank models continue to predict normally. New window models additionally
check that prediction sequences match their recorded training length.
