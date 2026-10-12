# RNA-seq evidence, candidate starts and curated gene comparisons

Version 1.3 extends Peter Thorpe's TransStart concept. A spatial RNA coverage
rise can help locate a candidate transcript-start neighbourhood and focus the
upstream motif comparison for a biological question. It is separate from a
between-condition DE test. The primary input is a curated positive gene list and
a curated negative list; DE can supply context or generate lists when absent.

## Analysis units and evidence boundaries

| Layer | Observation and method | Meaning and limits |
| --- | --- | --- |
| Start candidate | Sustained aligned-base coverage rise, called per biological sample | A descriptive edge; standard RNA-seq does not precisely validate transcript ends. |
| Replicate support | Coordinate agreement within each gene/condition | Independent observations and a recorded resolution/range; no confidence interval or edge p-value. |
| Gene DE | Raw exonic counts, declared contrast and categorical design, PyDESeq2 | Replicate-aware gene expression test; not a differential promoter-usage or enhancer-activity test. |
| Motif enrichment | Presence in RNA-refined/fallback positive versus negative flanks | Fisher/BH sequence comparison; gene-set and anchor-selection bias still require assessment. |
| Sequence ML/SHAP | Source-grouped held-out DNA signatures and additive log-odds explanations | Exploratory class resemblance. Label-defining RNA quantities are not predictors; SHAP is not causal evidence. |
| Functional evidence | Optional accessibility peaks, catalogues and supplied experiments | Contextual support, kept distinct from sequence/DE/start evidence; no invented joint enhancer probability. |

Curated labels never change when optional DE is fitted or imported. This includes
genes with small, missing or oppositely directed DE estimates. Selected lists and
the reason for every selected label are exported. If lists are absent, explicit
DE thresholds generate positive and expressed screening-control sets. Missing
tests and low-expression genes cannot become controls. Non-significant tests do
not establish equivalent expression. Users can inspect and trim those lists.

## Coverage access and biological conventions

All coordinates inside the API and bedGraph outputs are zero-based half-open.
GFF3/GTF parsing uses the established hierarchy resolver; exon records can be
unsorted, multi-parent and escaped. Each gene's annotated exons are merged for
coverage search limits and exonic counting. This union does not identify a unique
isoform or manufacture splice junctions. Missing exons remain absent. Coordinate
TSVs have no exon information and cannot support automatic exonic counting.

Alignments must declare coordinate sorting, have a usable index and match genome
contig lengths. CRAM uses the supplied genomic reference; counting here requires
BAM. Read blocks, rather than full alignment spans, provide depth. Skipped
introns, deletions and unsequenced mate gaps stay zero. Overlapping mates count
once per fragment, keyed by read name and read group. Reads marked unmapped,
secondary, supplementary or failed QC are excluded. MAPQ, duplicate inclusion
and a hard fetched-read cap are explicit settings; NH-tagged multi-mappers are
excluded. Alignments lacking NH cannot prove uniqueness from that tag, so the
mapper/quality convention must be assessed for externally produced BAMs.

`forward` means read 1 is sense to the transcript; `reverse` means antisense.
Read 2 follows the opposite convention. Unstranded depth can support a spatial
edge but cannot assign transcript origin on its own and remains marked ambiguous.
Transcript-plus/minus bedGraphs are already assigned to transcript strands and
remain in genomic coordinate order; they are never reversed again. A single
track must be labelled unstranded. Omitted bases mean zero; truncated tracks
must not be mistaken for measured zero-expression regions.

Read counts/depth filters are related but distinct. Paired featureCounts uses
paired-fragment counting, both mapped mates and same-contig fragment constraints;
depth considers observed aligned blocks even when a mate is absent. Record and
assess library QC rather than treating the two outputs as identical quantities.

## Candidate-start detection

The search begins in the original contiguous upstream gap, bounded by **all**
annotated genes on both strands. It may extend into the target's first 5-prime
exonic union, capped at the configured downstream limit. A different gene inside
that extension also stops the search. Missing exons allow a separately audited
gene-span search fallback; they never supply substitute exonic counts.

Depth is oriented in the gene's transcriptional direction. With default block
size W=20, every candidate requires one complete block before and two complete
blocks after it. Both downstream blocks must meet the raw mean-depth floor,
covered-base fraction and pseudocount-adjusted rise threshold. The smaller of the
two downstream means/fractions controls the check, rejecting an isolated spike.
The strongest log2 rise multiplied by covered fraction selects the edge; ties
use the first oriented position. Alternative separated edge clusters are counted.
This is a transparent heuristic, not a trained TSS detector or formal test.

Calls are made independently in each biological sample. Consensus chooses the
largest coordinate cluster whose full span is at most the predeclared tolerance
(default 40 bp), using rise strength and annotation proximity for ties. Sliding
cluster windows avoid chaining widely separated observations. The consensus is
the integer median coordinate. Support must meet both the replicate count floor
(default two) and the fraction of **all** condition observations (default 0.5).
Missing calls stay in that denominator. Technical libraries are not replicates.
Singletons remain provisional; discordant or insufficient support stays explicit.

The displayed neighbourhood includes block resolution and observed replicate
range. It is not a statistical confidence interval; bases within a coverage curve
are correlated and receive no invented p-values or TSS probabilities. Standard
RNA-seq can be biased by fragmentation, degradation, 3-prime enrichment,
transcript length, mappability, isoforms, overlapping transcription and
read-through. Specialised 5-prime evidence is needed to validate a start.

Between-condition shifts require supported candidates in both conditions. Signs
are transcription-oriented, positive downstream, and neighbourhood overlap is
recorded. Changes in detectability can reflect depth rather than promoter usage.
There is no promoter-usage p-value or test of alternative first-exon usage here.

## Strict extraction and selection bias

The original annotation is never rewritten. A supported edge can refine upstream
DNA only inside the same original contiguous gap. On the positive strand the
candidate's first transcribed base is the exclusive upper boundary; on the
negative strand extraction starts at candidate+1 and is reverse-complemented
once. Offsets and maximum lengths stay inside the available gap.

An edge inside the target's annotated UTR/exon remains a reported candidate,
but does not authorise returning genic sequence. An edge beyond a blocker cannot
jump that blocker. Both use an explicitly recorded annotation fallback or the
requested exclusion. Exact annotated 5-prime boundary candidates are usable.
Provisional singletons require an explicit opt-in; they remain unvalidated.

Both curated classes use the same condition and detection/extraction rules.
Nevertheless expression differences can affect edge availability, source length
and sequence inclusion. Inspect `rna_anchor_audit.tsv`, gene-set/exclusion audits
and class coverage before interpreting motif enrichment. Match controls for the
biological question, including expression, coverage quality, GC/length,
mappability, family and experimental design where appropriate. Automatic
GC/length matching does not establish expression matching.

## Optional DE and mapping

PyDESeq2 receives raw integer sample-by-gene counts for the declared contrast.
At least two independent biological observations per condition and an estimable
full-rank design with residual replication are required. Explicit batch/subject
terms are categorical. Missing/constant, confounded or saturated designs fail.
At least 20 genes must pass the declared total-count floor, and every selected
library must have counts. The model uses size factors, dispersion estimation,
Cook's refitting/filtering, Wald tests against zero log2 effect and independent
filtering. An additional effect-size threshold selects a screening set; it is
not itself a test against that effect-size threshold. No LFC shrinkage or complex
interaction/continuous-covariate contrast is implemented in this version; import
an externally fitted result table when those methods are required.

All input genes are exported, with missing results for filtered/untested genes.
Warnings and finite adjusted-test counts are recorded; a fit with no estimable
adjusted tests is marked `not_estimable`, retaining NA rather than zero p-values.
Estimated size factors divide displayed coverage; detection depth thresholds
remain explicitly raw. These displays are not TPM/FPKM.

Controlled mapping uses two-pass STAR on a genomic reference, unique alignments,
coordinate-sorted BAM and samtools indexing. A new unannotated index has zero
splice-junction overhang; the read-length overhang is supplied for two-pass
alignment. No junctions are invented by joining a gene's union of isoform exons.
Safe read links protect paths containing spaces/metacharacters; subprocesses use
argument vectors, explicit threads/RAM/time limits and bounded error excerpts.
Successful output records executable versions, commands and logs; raw-read links
are removed. Existing result directories are never overwritten. External indexes
without this package's genome fingerprint remain labelled identity-unverified.
Read QC/adapter trimming and protocol choice remain upstream responsibilities.

featureCounts receives a one-based SAF of merged exonic unions, counts each
sample with its declared layout/strand convention, and excludes ambiguous gene
assignment and multi-mapping under the configured defaults. Missing-exon genes
are audited rather than counted over an invented whole-gene feature.

## HTML and motif-position graphics

Coverage heatmaps use log1p of condition mean divided depth, with unobserved
positions grey. Locus panels show each biological sample, condition means and
candidate resolution intervals alongside the annotated gene boundary. Panels
are bounded and balanced between supplied classes; they are previews.

RNA motif-distance profiles use **actual native detected sites**, not substituted
PWM consensuses, and the separately audited applied anchor. Palindromic opposite
strand records at the same physical site are counted once. Site counts are
divided by eligible motif-centre positions for that motif width in the actual
returned sequence. RNA-derived and fallback anchor classes remain separate;
unavailable positions are not padded. Discovery-enriched motifs may be chosen
from the same complete data, so the heatmap is descriptive, with no positional
significance test. Existing genome-scan and multi-scale plots retain their
annotated-gene coordinate system; the two anchors are explicitly distinguished.

DE volcanoes and row-centred log1p divided-count heatmaps are visual summaries,
separate from fitting. RNA evidence and optional DE columns appear alongside
candidate records without changing evidence tiers or inventing an activity
probability. Embedded PNG/PDF figures, expandable notes and sortable/searchable
tables work offline. Complete TSVs and linked reports need the companion folder.

## Primary sources

- [Peter Thorpe's original TransStart](https://github.com/peterthorpe5/public_scripts/tree/master/TransStart)
- [Comprehensive comparison of 5-prime RNA-sequencing methods](https://pubmed.ncbi.nlm.nih.gov/29867192/)
- [SEASTAR alternative-start-site evaluation using RNA-seq](https://pubmed.ncbi.nlm.nih.gov/29546410/)
- [DESeq2 official vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [PyDESeq2 official pipeline](https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_minimal_pydeseq2_pipeline.html)
- [pysam alignment API](https://pysam.readthedocs.io/en/latest/api.html)
- [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [Subread/featureCounts official documentation](https://subread.sourceforge.net/)
