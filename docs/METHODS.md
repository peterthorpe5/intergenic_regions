# Methods, evidence and limitations

## Strictly intergenic extraction

Internally all spans are zero-based and half-open. A gene's span is the union
of its explicit record and recognised genic descendants, including introns.
All annotated genes are blockers, independently of target selection and strand.
Merged touching/overlapping spans form a sorted interval index. The free flank
touching the selected boundary is found by binary search; an occupied boundary
returns no sequence. It never advances to a farther gap.

For upstream on +, the anchor is the gene start and extraction proceeds left.
For upstream on −, the anchor is the gene end and extraction proceeds right.
Downstream swaps these sides. Every minus-strand flank is reverse-complemented.
Offset and length are applied inside the original free gap, then length and
ambiguity filters produce retained/excluded audit rows.

GFF3 Parent lists are split before percent decoding, preserving escaped commas
inside IDs. Records may be unsorted. Shared children contribute to each root.
Missing root records are inferred conservatively from Parent IDs. Discontinuous
consistent records are combined. Conflict/cycle checks prevent silently unsafe
gene models. GTF gene_id and legacy coordinate TSV are supported. Recognised
feature types are documented by GENIC_FEATURES in the annotation module.
Unrecognised custom feature types are not automatically assumed to be genes.

No sequence-only algorithm can detect every missing gene annotation. Circular
coordinates, transcript-specific TSSs and splicing are outside this version.
Annotation errors and assembly mismatches must be addressed upstream.

## Native motif statistics

Each sequence contributes presence/absence for each tested motif, irrespective
of its number of overlapping occurrences. Positive hits versus negative hits
form a 2×2 contingency table. The enrichment test is the one-sided Fisher exact
test, calculated as a hypergeometric upper tail. Its log representation is
retained to protect tiny probabilities from numerical underflow.

Effect estimates include prevalence, fold enrichment and a Haldane 0.5 corrected
odds ratio with approximate Wald 95% limits. These limits are approximate,
especially with small/sparse counts. Pseudocounts never enter the Fisher test.
A zero control prevalence gives an undefined/infinite fold enrichment represented
as null/NA, distinct from a p-value of zero.

Benjamini–Hochberg uses the complete requested hypothesis family: known motifs
plus all possible requested words, with reverse complements combined when
both strands are searched. Unobserved words count as p=1 hypotheses. FDR is
calculated in log space. Overlapping words are dependent; q-values are nominal
BH estimates and do not substitute for experimental replication.

IUPAC motifs use exact degenerate matching, with no statistical site p-value.
MEME/JASPAR matrices are normalised and scored against pooled ACGT background
with pseudocounts. A discretised score distribution yields zero-order per-site
tail probabilities. Site thresholds are per position, not genome-wide
significance. PWM site p-values and between-sequence enrichment q-values are
different quantities. Ambiguous windows are skipped; lowercase can be masked.
Native exact-word discovery does not model spaced/composite/long de novo motifs.

Genomic overlaps and exact/reverse-complement duplicate sequences are removed
before pipeline testing, favouring foreground representatives. Related loci
may still be dependent: use appropriate groups and biological study design.
Standalone FASTA analysis lacks coordinates and cannot detect genomic overlaps.

## Automatic ML

The separate learning module is an interpretable, supervised L2 logistic
classifier, not a pretrained enhancer foundation model. Features are canonical
word frequencies plus GC, log length and ambiguity. A composition-only model
provides a baseline. Fixed regularisation avoids selecting hyperparameters on
reported test folds. Vocabulary selection, scaling and model fitting are
performed inside every training fold. Signature coefficients from the full
model are exploratory; their fold stability is supplied, without per-feature
significance claims.

Validation uses stratified folds or stratified group folds. Groups are never
split or silently relaxed. At least five sequences per class and enough
independent, two-class folds are required. Family/homology/chromosome grouping
must reflect the dependence structure; grouping by gene alone cannot rule out
paralogous leakage. Out-of-fold scores evaluate labelled-sequence resemblance,
not enhancer validity or causal regulation.

Label permutations repeat the entire modelling procedure. Without groups,
labels are shuffled globally. In mixed-class groups, labels are shuffled within
each group. Pure-class groups exchange their labels only with equal-sized
pure groups. This preserves class totals but assumes the specified blocks are
exchangeable; it is not a universal null model. No usable exchanges means no
permutation p-value. The empirical p-value is
(1 + null scores at least observed)/(1 + completed permutations); its minimum
resolution is 1/(B+1). With 99 permutations this is 0.01. Use more replicates
when stronger p-value resolution is required.

ROC AUC, average precision, fold metrics, baseline performance and permutation
results should be inspected together. They are exploratory internal validation;
there is no claim of external validation or calibrated probabilities. The JSON
model contains weights/scaling/schema/checksums, never executable pickle data.
Exact training-sequence duplicates are flagged when predicting new candidates;
near duplicates and homologues still require independent validation.

## Exact held-out SHAP

For a linear logit model with intercept b, weights w and training background
mean μ, the expected logit is `b + w·μ`; feature i contributes
`w_i × (x_i − μ_i)`. This is exact interventional SHAP, checked against
`shap.LinearExplainer` with an independent masker. Each held-out observation
uses its own fold's fitted classifier, vocabulary, scaling and training-only
feature mean. Test observations and the full-data model never define that
background. Permutation models do not generate additional SHAP output.

All features enter mean-absolute global importance. A feature absent from a
fold contributes zero and has an unavailable feature value. Importance averages
over the explained sample; seeded sampling is approximately proportional by
class and can affect feature ranking. Feature blocks bound sparse-to-dense
memory. Leading local rows plus an explicit additive remainder reconstruct
logits and held-out signature scores. Sampling limits, scope and reconstruction
error are recorded. This sample does not imply a population estimate.

Official `shap.plots` generates bar, beeswarm, heatmap, waterfall and dependence
graphics. Heatmap sequences sort by actual held-out logits; its top trace sums
contributions without fold-specific baselines and is labelled accordingly.
Waterfalls show the highest/lowest predicted scores. Dependence-plot fold
variation is not an interaction test. Grey values include the remainder and
features absent from a fold. Units are positive-class log-odds, not probability
changes, enrichment significance or enhancer validation.

Interventional attribution treats perturbations independently despite correlated
or overlapping words. Correlation-dependent attribution, causal effects and
explanations of the contextual-evidence ranking are outside this method.
Chromatin evidence remains a separate support layer. Missing official plotting
dependencies produce an explicit status while preserving numerical SHAP and ML.

## Whole-genome consensus screening

Enrichment-driven selection requires positive prevalence greater than negative
prevalence and discovery q ≤ the configured limit, ordered by q then ID and
capped by a pattern count. Supplied motifs retain file order. Requested IDs
further restrict selection. No eligible motif gives `no_selected_motifs` and
empty site files, rather than invented hits.

IUPAC positions accept their allowed bases without consuming substitutions.
Hamming distance counts other substitutions; indels are excluded. Ambiguous
genomic windows are always skipped; soft masking is optional. PWMs become
explicit maximum-probability consensuses, so this screen does not reproduce
their discovery scoring threshold. Discovery q-values only preserve provenance.
No genome-hit significance is inferred from them or from mismatch counts.

Both orientations are searched by default. A coordinate matching both is one
physical site with strand `.` and minimum mismatch count. Overlapping windows
are retained. Chunk tails cover the longest pattern; core-start ownership
avoids boundary omissions/duplicates. Full TSV/BED sites stream to disk; an
exceeded hit limit fails the atomic workflow. Intergenic-only filtering rejects
any overlap with the complete gene spans, regardless of target membership,
strand or availability of a gene-start proxy.

Association uses the nearest known-strand annotated 5′ gene base: zero-based
`start` on + and `end−1` on −. For midpoint c and proxy t, distance is `c−t`
on + and `t−c` on −. This is a TSS proxy, not a measured TSS or a regulatory
link. Equidistant/coincident starts select the first coordinate/ID representative,
retain a tie flag and are excluded from positional numerators and denominators.
Gene burden retains that representative; unknown-strand genes still block hits.

Density divides hits by eligible scanned start windows of the same pattern
length, nearest-gene cohort and distance bin, scaled by 10^6. Masking, ambiguity,
intergenic restriction, ties and contig edges affect counts and opportunities.
The final bin may be narrower; plots use true edges. Zero opportunity gives NA;
raw counts are also supplied. Contig/gene burden is descriptive, not normalised
enrichment. Discovery-selected motifs/cohorts are reused and cannot establish
independent validation. Every hit remains `unvalidated_sequence_match`.

## Optional evidence and candidate triage

BED evidence measures union overlap, preventing double counting. Missing tracks
are null; supplied tracks with no overlap are zero. Assembly and contig names
must agree. User tracks have no automatic liftover. Imported references carry
organism/build/source/checksum metadata and are checked on query.

Accessibility supports open chromatin in that assay's tissue/condition.
Regulatory-reference overlap is an annotation match. Gene-level assays are
linked records, not validation of an extracted interval. Explicit positive
and negative values are counted separately; arbitrary values are uninterpreted.

Ranking first uses contextual-evidence tiers, then held-out signature scores.
It does not learn a joint evidence model or manufacture a combined probability.
Every candidate remains unvalidated. A high score, motif enrichment, accessibility
or catalogue overlap alone does not establish enhancer activity, enhancer-to-gene
linkage, tissue specificity or mechanism.

## Primary documentation

- [Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
- [SciPy Fisher exact test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)
- [scikit-learn grouped cross-validation](https://scikit-learn.org/stable/modules/cross_validation.html)
- [scikit-learn permutation-test documentation](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.permutation_test_score.html)
- [HOMER FASTA-mode documentation](https://homer.ucsd.edu/homer/microarray/fasta.html)
- [SHAP linear explainer](https://shap.readthedocs.io/en/latest/generated/shap.LinearExplainer.html)
- [SHAP plotting API](https://shap.readthedocs.io/en/latest/api.html#plots)

Our grouped permutation procedure is described above explicitly; it extends
within-group shuffling to equal-sized pure groups and is not presented as
identical to scikit-learn's permutation_test_score.
