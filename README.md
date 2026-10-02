# Intergenic regions

**Strand-aware intergenic extraction, motif enrichment, multi-scale regulatory
region search, automatic ML with SHAP, genome motif screening and HTML reports.**

**Project creator, original author and maintainer:**
[Peter Thorpe (@peterthorpe5)](https://github.com/peterthorpe5).

Peter Thorpe developed the original package, its biological concept and its
strand-aware approach to extracting upstream and intergenic sequence while
respecting neighbouring gene boundaries. He continues to lead the scientific
design, development priorities and maintenance of this project.

The version 1.0.0 overhaul was developed with AI assistance from OpenAI Codex
for implementation, testing and documentation, under Peter Thorpe's direction.
The original Git history is retained, preserving attribution for the earlier
package and its development.
The version 1.1.0 SHAP and genome-screening upgrade continues this work under
Peter Thorpe's scientific direction, with AI-assisted implementation and testing.
Version 1.2.0 adds the multi-scale search proposed by Peter Thorpe to distinguish
short motif patterns from larger, variable-length candidate regulatory regions.

This modern Python 3.11+ package replaces the script at master commit
`c36293fc8f40117a8401a3e42aba255e172c914f` (4 September 2020).

The biological contract is unchanged: stop before any neighbouring annotated
gene, on either strand; never jump a blocker; reverse-complement negative-strand
targets into transcriptional orientation. Overlapping genes, unknown strands
and unavailable flanks receive explicit exclusion records.

## Installation


### using conda:

Alternatively, create the same analysis environment from the supplied file:
`conda env create --file environment.yml`. For an existing working environment,
activate it and use the pip install/check commands below.

```bash
mamba create -n intergenic_regions \
  --override-channels -c conda-forge \
  --strict-channel-priority \
  python=3.11 pip \
  numpy=1.26.4 scipy=1.13.1 matplotlib=3.9.2 \
  scikit-learn=1.5.2 shap=0.49.1 \
  numba=0.61.2 llvmlite=0.44.0 \
  -y
```

```bash
conda activate intergenic_regions
```

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install --only-binary=:all: --editable '.[dev]'
```

Then check the installation:

```bash
python -m pip check
python -m intergenic_regions --version
python -m pytest --cov=intergenic_regions --cov-report=term-missing
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python -m build
```


### Other Python environments with binary wheels

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --only-binary=:all: --editable '.[analysis]'
python -m intergenic_regions --help
```

The installed command is `intergenic-regions`. Extraction-only installation:
`python -m pip install --editable .`. The `[analysis]` extra includes motifs,
plots, ML and the official SHAP package; `[motifs]` and `[ml]` remain separate
internally. Extraction does not load ML libraries. This project is not assumed
to be published on PyPI.

## Complete gene-list analysis

Each gene list contains **one exact annotation gene ID per line**, with optional
blank lines and `#` comments. Always provide the **full annotation**, including
genes outside the target/control lists: they must still block extraction.

```bash
python -m intergenic_regions pipeline \
  --genome genome.fasta \
  --annotation genes.gff3 \
  --positive-genes genes_of_interest.txt \
  --negative-genes control_genes.txt \
  --length 1000 \
  --kmer-lengths 5 6 \
  --output-dir results/my_analysis \
  --log-file analysis.log
```

Open `results/my_analysis/report.html`. The unified dashboard contains result
cards, motif and model plots, searchable/sortable tables, ranked candidates,
uncertainty labels and links to full TSV files. Images, styles and scripts are
embedded: the HTML works offline. Figure PNG/PDF downloads are embedded too;
full TSV/BED/model links need the companion folder. PNG and vector PDF plots
are also supplied separately.

**Both `pipeline` and `enrich` automatically run ML alongside motif analysis.**
A separate `ai` command supports reuse and experiments. `--no-ml` explicitly
disables learning; `--ai` remains an accepted, redundant compatibility flag.
If the observations/groups cannot support validation, motif results are kept
and ML is marked `not_estimable`, with a reason. Missing ML dependencies are
marked `unavailable`. Neither status fabricates scores. Invalid settings or
malformed grouping files still fail.
**SHAP runs automatically whenever the model is estimable.** `--no-shap`
explicitly disables explanations without disabling ML. Without the SHAP plotting
dependency, exact numerical explanations are retained and plot status explains
how to install `[analysis]`.

## Enhancer-sized regions: search multiple lengths

A transcription-factor motif and an enhancer are different analysis units.
Short words describe local sequence patterns; larger regions can contain several
sites, their spacing and surrounding context. There is no universal enhancer
length. In a human MPRA study, testing the same candidate elements at **192,
354 and 678 bp** gave substantial activity differences. A plant study found
activity in **169 bp** segments and context-dependent cooperation between their
functional elements. These are experimental examples, not universal boundaries
([Klein et al., 2020](https://doi.org/10.1038/s41592-020-0965-y);
[Plant enhancer study, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11218779/)).

**`pipeline` and `enrich` now search overlapping 100, 200, 400, 800 and 1,600 bp
windows automatically**, stepping by 50 bp. Configure other lengths, including
larger ones, without turning a whole enhancer into an exact long k-mer:

```bash
python -m intergenic_regions pipeline \
  --genome genome.fasta \
  --annotation genes.gff3 \
  --positive-genes positive.txt \
  --negative-genes negative.txt \
  --full-gap \
  --kmer-lengths 4 6 8 \
  --region-lengths 100 200 400 800 1600 3200 \
  --region-step 50 \
  --accessibility-bed tissue_atac_peaks.bed \
  --output-dir results/multiscale
```

Omit `--accessibility-bed` for a complete sequence-only run. `--length 5000`
can replace `--full-gap` to set a maximum search extent. The extraction default
is still 1,000 bp: scales larger than a retained flank are marked unavailable.
Every window stays inside the original free flank, including on the negative
strand. End-aligned windows cover the last bases; short gaps are audited without
padding or silently shortening a requested scale.

| Setting | Meaning |
| --- | --- |
| `--kmer-lengths 4 6 8` | Short exact patterns tested for source-level enrichment |
| `--ai-kmer-lengths 4 5 6` | Short word features used by the models |
| `--region-lengths 100 200 400 800 1600` | Larger candidate search windows |
| `--length 5000` / `--full-gap` | Extent of the strictly intergenic source flank |

Each available scale gets its own window model automatically. Whole source
flanks, supplied families/contigs and exact shared windows stay together during
validation. Training weights give each source equal mass within its class;
reported ROC/AP use one mean held-out score per source. Window labels are
inherited from gene/flank labels, so this evaluates **class resemblance**, not
experimentally validated enhancer detection. Each scale also supplies exact
held-out SHAP tables and, when installed, official SHAP plots.

Above-threshold overlapping windows form variable-length **candidate unions**.
The default `--region-score-threshold 0.75` is an exploratory triage setting,
not an enhancer probability or significance threshold. Scale scores are not
assumed to share calibration. Searching more windows gives more chances for a
high peak; the unions have **no region-level q-value**, and their edges are not
measured enhancer boundaries. All candidates remain `unvalidated_candidate`.

The HTML includes scale availability, source-level validation, window-length
versus gene-start-distance heatmaps, motif density, position coverage, candidate
span distributions and individual locus views with optional accessibility
overlap. Signed gene-start distances use the annotation's 5′ base as a TSS
proxy; FASTA-only `enrich` uses oriented sequence offsets. See
`regions/report.html`, `regions/window_scores.tsv`, `regions/candidate_regions.tsv`,
`regions/candidate_regions.bed` and `regions/candidate_regions.fasta`.

Use `--no-regions` to disable this layer explicitly. `--no-ml` keeps the window
availability and motif-density analysis while leaving model scores empty;
`--no-shap` disables explanations. Limits include `--region-max-windows 50000`,
`--region-max-locus-plots 6` and a configurable profile distance range.
Larger searches that exceed the window limit fail with guidance; increase the
step, restrict the input or intentionally raise the limit. Short-word models
still lack explicit motif-spacing/cooperativity features. Regulatory activity,
enhancer-to-gene linkage and functional boundaries require independent evidence.

The regional search covers extracted intergenic flanks or supplied FASTA.
The genome-screening command continues to locate motif sites across the genome.
Enhancers within genes or beyond an intervening gene are outside a contiguous
strictly intergenic flank, even when `--full-gap` is used.

Accessibility, regulatory annotations and functional evidence are optional:

```bash
python -m intergenic_regions pipeline \
  --genome genome.fasta \
  --annotation genes.gff3 \
  --positive-genes genes_of_interest.txt \
  --negative-genes control_genes.txt \
  --motifs motifs.meme \
  --accessibility-bed tissue_atac_peaks.bed \
  --enhancer-bed regulatory_annotations.bed \
  --evidence-tsv functional_results.tsv \
  --groups-tsv independent_groups.tsv \
  --output-dir results/with_evidence
```

Evidence must use the same assembly and contig names. Missing evidence is NA,
distinct from zero overlap in supplied data. Open chromatin, motif presence
and positive-class resemblance support hypotheses; none proves enhancer
function. Gene-level assays remain distinct from evidence on the actual flank.

## Reproducible end-to-end demo

No downloads are needed. The generator makes a synthetic genome, unsorted GFF3,
positive/negative genes, planted G-box motifs, groups and optional evidence.

```bash
python examples/make_demo.py --output-dir results/demo_inputs
python -m intergenic_regions pipeline \
  --genome results/demo_inputs/genome.fasta \
  --annotation results/demo_inputs/genes.gff3 \
  --positive-genes results/demo_inputs/positive_genes.txt \
  --negative-genes results/demo_inputs/negative_genes.txt \
  --motifs results/demo_inputs/motifs.tsv \
  --kmer-lengths 6 \
  --accessibility-bed results/demo_inputs/accessibility.bed \
  --evidence-tsv results/demo_inputs/functional.tsv \
  --groups-tsv results/demo_inputs/groups.tsv \
  --folds 3 \
  --permutations 19 \
  --scan-genome --scan-mismatches 1 --scan-max-motifs 8 \
  --scan-intergenic-only --scan-upstream 300 --scan-downstream 300 \
  --scan-bin-width 20 \
  --output-dir results/demo_analysis
```

Open `results/demo_analysis/report.html`. Planted demo results illustrate
software behaviour, **not biological discoveries or enhancer validation**.
Default settings are five folds and 99 permutations; the demo uses fewer
replicates for speed. See [example guide](examples/README.md).

## Extraction only

```bash
python -m intergenic_regions extract \
  --genome genome.fasta --annotation genes.gff3 \
  --genes targets.txt --length 1000 --min-length 50 \
  --output-dir results/upstream
```

Omit `--genes` to process all genes.

| Option | Behaviour |
| --- | --- |
| `--direction upstream` | Default transcription-relative upstream flank. |
| `--direction downstream` | Downstream flank with identical boundary protection. |
| `--direction both` | Both flanks, oriented to the target gene. |
| `--full-gap` | Entire contiguous gap; mutually exclusive with `--length`. |
| `--offset 100` | Skip nearest bases within the same gap; never leap a blocker. |
| `--min-length 50` | Inclusive retained-length threshold. |
| `--mask-lowercase` | Replace soft-masked bases with N. |
| `--max-ambiguous-fraction 0.1` | Audit/exclude excessively ambiguous regions. |
| `--annotation-format gff3` | Explicit GFF3, GTF or legacy TSV; default auto. |

The package handles unsorted hierarchies, multiple transcripts, shared children,
percent-escaped identifiers, omitted roots and discontinuous features. Full
gene/descendant spans include introns and UTRs. Invalid coordinates, conflicting
IDs/strands and cyclic hierarchies fail. Annotation completeness determines what
can be recognised as genic; linear contigs are assumed. Circular genomes and
transcript-specific TSS selection are not currently modelled.

Extraction writes FASTA, BED, GFF3, a region audit, evidence tables, summary,
input checksums and HTML. **TSV/BED/API are zero-based, half-open; GFF3 is
one-based, inclusive.** FASTA region IDs escape reserved gene-ID characters
and append `|upstream` or `|downstream`; the original gene ID is in the TSV.

## Existing FASTA: motifs and automatic ML

```bash
python -m intergenic_regions enrich \
  --positive-fasta positive_promoters.fasta \
  --negative-fasta negative_promoters.fasta \
  --motifs motifs.jaspar --kmer-lengths 5 6 \
  --groups-tsv independent_groups.tsv \
  --output-dir results/fasta_analysis
```

Known motifs: MEME DNA probability matrices, JASPAR count matrices, or an IUPAC
TSV with `motif_id`, `pattern` and optional `name`. With neither a motif file
nor lengths, six-base words are tested. Motif reverse complements are combined
unless `--forward-only` is supplied; ML always uses strand-canonical words.

Native enrichment uses a one-sided sequence-presence Fisher exact test and
Benjamini–Hochberg FDR over known motifs **and all possible tested words**,
including unobserved words. Outputs include p/q-values, log-space values,
prevalence, effect sizes and approximate effect-size confidence intervals.
PWM site p-values and sequence-level enrichment p-values are distinct.
Native discovery searches exact words, rather than reproducing HOMER's full
de novo optimisation. See [methods](docs/METHODS.md).

ML uses regularised logistic regression on canonical k-mer frequencies and
GC/length/ambiguity covariates, with a composition-only baseline. Vocabulary
selection, scaling and fitting are repeated inside each validation fold.
Supply homology/family/chromosome groups when needed; pipeline also offers
`--group-by contig`. Grouping is never silently relaxed. Scores measure
positive-class resemblance, **not calibrated enhancer probabilities**.

Pipeline filters overlapping/shared genomic regions and exact or reverse-
complement sequence duplicates, with an audit. Standalone FASTA analysis
rejects duplicates instead of silently modifying inputs. Optional deterministic
greedy GC/length matching: `--match-background --background-ratio 1
--max-gc-difference 0.1 --max-length-ratio 1.5`. All positives must receive
controls; this is not a globally optimal assignment.

## Automatic SHAP explanations and graphics

The official SHAP package produces a global importance bar, beeswarm,
sequence-by-feature heatmap, waterfall examples for the highest/lowest model
scores, and a dependence plot for the leading displayed feature. All appear
in the unified HTML, with PNG/PDF downloads and exact TSV values.

Explanations use the **held-out fold's model and training-only background**.
For this linear model, interventional SHAP is exact. The baseline plus all
contributions reconstructs the model's positive-class **log-odds**, not changes
in probability. Fold backgrounds differ; the heatmap's top trace is labelled
as the sum of SHAP contributions.

By default, a reproducible, approximately class-proportional sample of up to
1,000 sequences is explained. `--shap-max-sequences 2000` changes that limit;
`--shap-max-features 30` changes the leading features retained per sequence
(supported range 1–100). Global importance includes **all** model features.
An explicit additive `other_features` row preserves omitted contributions.
Grey points have no feature value, including that remainder and words absent
from a fold's vocabulary.

Correlated, overlapping words affect attribution. SHAP explains the classifier,
**not causal regulation or validated enhancer function**. It does not explain
the separate contextual-evidence ranking. New-candidate prediction remains
separate from held-out explanations. See [methods](docs/METHODS.md).

## Screen identified motifs across the whole genome

Add `--scan-genome` to `pipeline` to screen identified motifs in the same run.
The screen is opt-in because a genome can produce millions of matches.
Alternatively, reuse completed enrichment results:

```bash
python -m intergenic_regions scan-genome \
  --genome genome.fasta --annotation genes.gff3 \
  --enrichment-tsv results/my_analysis/motifs/motif_enrichment.tsv \
  --positive-genes genes_of_interest.txt --negative-genes control_genes.txt \
  --scan-q-value 0.05 --scan-max-motifs 20 \
  --scan-mismatches 1 --scan-intergenic-only \
  --scan-upstream 2000 --scan-downstream 500 --scan-bin-width 100 \
  --output-dir results/genome_screen
```

Use `--motifs motifs.tsv` instead of `--enrichment-tsv` for supplied IUPAC,
MEME or JASPAR motifs. No annotation is required for matching alone; it is
required for `--scan-intergenic-only`, gene cohorts and positional heatmaps.
Always provide the full annotation to protect all genic boundaries.

| Setting | Interpretation |
| --- | --- |
| `--scan-q-value 0.05` | Select motifs with positive prevalence greater than negative prevalence and discovery q ≤ threshold. This is not a genome-hit p-value. |
| `--scan-max-motifs 20` | Pattern limit; enriched motifs sort by q then ID. Supplied motifs retain file order. Selected targets are written to TSV. |
| `--scan-motif-ids Gbox motif_2` | Restrict exact IDs, subject to discovery selection and the pattern limit. |
| `--scan-mismatches 1` | Allow substitutions only, excluding indels. IUPAC degeneracy does not consume a mismatch. The limit must be less than every selected motif's length. |
| `--scan-forward-only` | Search the reference strand only; otherwise both orientations. |
| `--scan-intergenic-only` | Exclude any overlap with an annotated gene, including introns/UTRs and unknown-strand genes. |
| `--scan-mask-lowercase` | Exclude soft-masked windows. Ambiguous genomic windows are always excluded. |
| `--scan-max-hits 2000000` | Abort atomically if exceeded; never publish a silently truncated scan. |
| `--scan-chunk-size 250000` | Chunk core size; overlapping tails protect cross-boundary matches. |

This is **consensus/IUPAC Hamming screening**, not PWM-score scanning. PWMs
become explicitly labelled maximum-probability `pwm_consensus` patterns;
discovery PWM thresholds are not reproduced. Discovery q-values are provenance,
never genome-hit significance. Overlapping sites are retained. A physical site
matching both orientations appears once with strand `.` and its minimum mismatch
count. `matched_sequence` contains reference-strand bases; coordinates are
zero-based, half-open.

The report adds raw-count and opportunity-normalised distance heatmaps,
positive/control/other cohort profiles, a gene-by-motif burden heatmap, mismatch
counts and contig counts. Signed distance connects the site's midpoint to the
**nearest annotated 5′ gene base**: `start` on `+`, `end - 1` on `-`; upstream
is negative in transcriptional orientation. This is a **TSS proxy**, not a
measured TSS or proof of a motif-to-gene link. Equidistant/coincident starts
are flagged and excluded from positional profiles; unknown-strand genes have
no start proxy. Gene burden retains a stable representative for tied assignments.

Density is sites per million eligible, unambiguous start windows of the same
motif length, under the same scan restrictions. Zero opportunity gives NA.
The profiles reuse discovery motifs and cohorts descriptively; they are
**not independent enrichment tests**. Every hit retains
`enhancer_status=unvalidated_sequence_match`.

## Evidence-aware prioritisation

Candidates are ranked by contextual support tier, then held-out ML score.
All retain `enhancer_status=unvalidated_candidate`.

| Tier | Support |
| --- | --- |
| 3 | Accessibility plus regulatory annotation/reference overlap. |
| 2 | Accessibility or annotation/reference overlap. |
| 1 | Positive linked gene evidence without contradictory gene evidence. |
| 0 | Sequence-only hypothesis or insufficient contextual support. |

This is a transparent triage rule, not a calibrated statistical fusion.
Gene-level positive/negative assays are counted separately. Explicit positive
values: `positive`, `supported`, `1`, `true`; negative: `negative`,
`not_supported`, `0`, `false`. Other values remain uninterpreted.

Functional-evidence TSV example:

```text
gene_id	evidence_type	value	source
GENE001	reporter_assay	positive	tissue-specific experiment
GENE002	CRISPR_perturbation	negative	independent experiment
```

Grouping TSV uses analysed FASTA identifiers:

```text
sequence_id	group
GENE001|upstream	family_1
GENE002|upstream	family_1
```

Reuse the inert JSON model on independent sequences:

```bash
python -m intergenic_regions predict \
  --model results/my_analysis/ai/model.json \
  --fasta independent_candidates.fasta \
  --output-dir results/predictions
```

Training duplicates are flagged. Candidate predictions are separate from
held-out validation scores.

## HOMER and regulatory reference tracks

HOMER is installed separately. The wrapper passes an explicit background,
avoids shell interpretation, and checks errors/timeouts:

```bash
python -m intergenic_regions homer \
  --positive-fasta results/my_analysis/positive.fasta \
  --negative-fasta results/my_analysis/negative.fasta \
  --motif-lengths 8 10 12 --threads 8 --output-dir results/homer
```

Use `--dry-run` for a clearly labelled command plan. External HOMER result
formats remain HOMER's own; package-native tables are TSV. Automatic ML
accompanies native `pipeline`/`enrich`; `homer` delegates to the external tool.

Local reference tracks can be imported and checked against organism/build:

```bash
python -m intergenic_regions import-reference \
  --source custom --bed regulatory_regions.bed \
  --organism Arabidopsis_thaliana --assembly TAIR10 \
  --name tissue_track --output-dir references/tissue_track
python -m intergenic_regions annotate \
  --regions-tsv results/upstream/regions.tsv \
  --reference-dir references/tissue_track \
  --organism Arabidopsis_thaliana --assembly TAIR10 \
  --output-dir results/reference_overlap
```

`references` lists source metadata. Tracks are hashed, different builds and
changed bundles are rejected, and no implicit liftover is performed. External
downloads are optional. Core workflows and local imports work offline.

## Python API

All project functions use named arguments:

```python
from pathlib import Path
from intergenic_regions import (
    Genome,
    GeneIndex,
    extract_regions,
    read_annotation,
)

genes = read_annotation(path=Path("genes.gff3"))
with Genome(path=Path("genome.fasta")) as genome:
    index = GeneIndex(genes=genes, lengths=genome.lengths)
    regions = extract_regions(
        index=index,
        genome=genome,
        identifiers=["GENE001", "GENE002"],
        length=1000,
    )
retained = [r for r in regions if r.status == "retained"]
```

[API guide](docs/API.md), [migration from the 2020 script](docs/MIGRATION.md),
[methods and uncertainty](docs/METHODS.md), [output reference](docs/OUTPUTS.md).

## Tests and development

```bash
python -m pip install --editable '.[dev]'
python -m pytest --cov=intergenic_regions --cov-report=term-missing
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python -m build
```

Unit/focused workflow tests exercise every public function, including failures.
Boundary regressions use exact plus/minus cases and an independent per-base
oracle for random overlapping/nested genes. Subprocess tests cover complete
CLI analyses, errors, refused overwrites and atomic publication.
[Test guide](tests/README.md) · [Contributing](CONTRIBUTING.md).

Genomes use private temporary on-disk indices; flank queries use sorted merged
intervals instead of rescanning genes. Gzip genomes need temporary uncompressed
disk space. On HPC set `TMPDIR` to writable scratch. Annotation/selected
sequences remain in memory; ML/permutation cost scales with observations,
features and folds. No numerical speed-up is claimed without benchmarking.

Each command needs a **new output directory**; failed bundles are not published
and existing results are refused. Progress is logged to stderr and optional
`--log-file`; summaries go to stdout as JSON. Exit codes: 0 success, 2 failed
operation/input, 130 interrupted. A successful motif workflow may carry an
explicitly unavailable model: inspect the model-status card and summary.
