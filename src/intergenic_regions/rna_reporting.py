"""Inspectable RNA coverage, uncertainty and upstream motif profiles."""

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from intergenic_regions.models import Region
from intergenic_regions.reporting import save_figure

RNA_NOTES = [
    "RNA-seq coverage edges are descriptive candidate transcription-start "
    "neighbourhoods. Fragmentation, degradation, 3-prime bias, isoforms, "
    "read-through, mappability and unstranded libraries can move or obscure "
    "them. CAGE/RAMPAGE, 5-prime RACE or appropriate independent evidence is "
    "needed to validate a transcription start.",
    "Detection uses a full upstream block followed by two sustained covered "
    "downstream blocks, independently in biological replicates. The reported "
    "interval combines window resolution and observed replicate range; it "
    "is not a statistical confidence interval. Coverage bases are correlated "
    "and receive no invented edge p-values.",
    "Only supported candidate anchors within the original contiguous "
    "intergenic gap refine motif extraction by default. Starts inside the "
    "target annotation or beyond any neighbouring gene fall back explicitly "
    "to the annotation, or are excluded if requested. Negative-strand DNA "
    "remains reverse-complemented into transcriptional orientation.",
    "Both gene classes use the same predeclared condition and start-selection "
    "rules. Inspect anchor sources and excluded genes: differing expression "
    "or coverage quality can create selection bias. Original multi-scale "
    "plots retain annotated-gene distances; RNA motif profiles use their "
    "separately audited applied anchors.",
    "Gene-level DE uses raw exonic counts and independent biological "
    "replicates, with a declared contrast and optional categorical nuisance "
    "terms. NA tests and low-expression genes cannot become negatives. "
    "Small observed effects with non-significant tests define screening "
    "controls, not proven equivalent expression.",
    "Automatically fitted sequence ML/SHAP measures held-out resemblance "
    "to the selected gene class. RNA quantities used to define DE labels "
    "are not predictor features. Keep gene families/chromosomes together "
    "with explicit groups when required; cross-validation is exploratory "
    "and does not replace validation in independent experiments.",
    "RNA evidence supports gene expression and candidate promoter location. "
    "It does not validate a promoter sequence, a TF binding event or an "
    "enhancer. Accessibility and functional evidence remain separate, "
    "optional evidence fields; no joint enhancer probability is invented.",
    "Coverage plots divide depths by the manifest size_factor when supplied "
    "or by fitted DE size factors. A unit divisor leaves raw depth. These "
    "are coverage visualisations, not TPM/FPKM, and depth alone is not used "
    "for differential-expression significance.",
    "Motif distance profiles are descriptive. Motifs may be selected using "
    "the same complete discovery dataset, so these heatmaps are not an "
    "independent confirmation of enrichment or a positional significance "
    "test.",
]

RNA_FIGURE_NOTES = {
    "rna_coverage_heatmap": (
        "Each row is one gene in one condition. Colours show log1p mean depth"
        " after the recorded sample divisors. Gaps remain unobserved. Distanc"
        "es are relative to the annotated 5-prime gene base, not a validated "
        "TSS. Compare classes alongside their coverage and anchor-availabilit"
        "y audits."
    ),
    "rna_start_support": (
        "Independent biological replicate agreement is shown explicitly. A mi"
        "ssing edge may reflect weak coverage rather than absent transcriptio"
        "n. Supported edges still require orthogonal TSS validation."
    ),
    "rna_start_shifts": (
        "These are descriptive shifts between replicate-supported candidate c"
        "overage edges. Positive distances are downstream in transcriptional "
        "orientation. No differential promoter-usage p-values are assigned."
    ),
    "rna_de_volcano": (
        "Each point is a tested gene. Adjusted p-values and fold changes come"
        " from the recorded DE method/contrast. Grey genes can include low-po"
        "wer comparisons; non-significance does not prove unchanged expressio"
        "n."
    ),
    "rna_expression_heatmap": (
        "The displayed subset has the largest variance in log1p counts divide"
        "d by recorded size factors. Rows are centred for visual comparison. "
        "This display is separate from model fitting and does not establish a"
        " regulatory mechanism."
    ),
    "rna_motif_distance_heatmap": (
        "Counts use actual detected motif sites and are divided by eligible m"
        "otif-centre positions in the returned sequences. RNA anchors and ann"
        "otation fallbacks remain separate. Discovery-selected motifs and cov"
        "erage-dependent gene selection can bias this display."
    ),
}


def rna_motif_profiles(
    *,
    regions: Sequence[Region],
    anchors: Sequence[Mapping[str, Any]],
    sites: Sequence[Mapping[str, Any]],
    motifs: Sequence[Mapping[str, Any]],
    labels: Mapping[str, str],
    bin_width: int = 50,
    upstream: int = 2000,
    max_motifs: int = 20,
    q_value: float = 0.05,
) -> list[dict[str, Any]]:
    """Count actual motif sites relative to separately audited applied anchors.

    Args:
        regions: Original-safe, strand-oriented upstream sequences.
        anchors: Per-sequence applied anchor and source audit.
        sites: Native motif-site records, not substituted PWM consensuses.
        motifs: Native enrichment records; discovery selection is descriptive.
        labels: Positive/negative class per sequence ID.
        bin_width: Distance bin width.
        upstream: Maximum displayed upstream distance.
        max_motifs: Maximum positively enriched discovery motifs.
        q_value: Predeclared discovery q threshold.

    Returns:
        Actual site counts, covered-parent counts and eligible motif-centre
        opportunities per motif/class/anchor-source/distance bin. Palindromic
        opposite-strand records at the same site are counted once.

    Raises:
        ValueError: Profile bounds or anchor metadata are invalid.
    """
    if (
        isinstance(bin_width, bool)
        or bin_width < 1
        or upstream < bin_width
        or math.ceil(upstream / bin_width) > 500
        or not 1 <= max_motifs <= 100
        or not math.isfinite(q_value)
        or not 0 < q_value <= 1
    ):
        raise ValueError("Invalid RNA motif profile bounds")
    selected = [
        r
        for r in motifs
        if float(r["q_value"]) <= q_value
        and (float(r["positive_fraction"]) > float(r["negative_fraction"]))
    ][:max_motifs]
    by_motif = {r["motif_id"]: r for r in selected}
    by_region = {r.sequence_id: r for r in regions if r.status == "retained"}
    by_anchor = {r["sequence_id"]: r for r in anchors}
    if len(by_anchor) != len(anchors) or len(by_region) != sum(
        r.status == "retained" for r in regions
    ):
        raise ValueError("RNA motif anchors/regions repeat")
    for identifier in by_region:
        if identifier not in labels or labels[identifier] not in {
            "positive",
            "negative",
        }:
            raise ValueError("RNA motif class metadata is missing")
        if (
            identifier not in by_anchor
            or by_anchor[identifier].get("applied_anchor") is None
        ):
            raise ValueError("RNA motif source lacks an applied anchor")
    counts: dict[tuple[str, str, str, int], int] = defaultdict(int)
    seen: set[tuple[str, str, int, int]] = set()
    for site in sites:
        identifier, motif = site["sequence_id"], site["motif_id"]
        if identifier not in by_region or motif not in by_motif:
            continue
        first, last = int(site["start"]), int(site["end"])
        region = by_region[identifier]
        if first < 0 or last <= first or last > len(region.sequence):
            raise ValueError("RNA motif site is outside its source sequence")
        site_key = (identifier, motif, first, last)
        if site_key in seen:
            continue
        seen.add(site_key)
        anchor = by_anchor.get(identifier)
        if anchor is None or anchor.get("applied_anchor") is None:
            raise ValueError(
                "Retained RNA motif region lacks an applied anchor"
            )
        midpoint = (
            region.start + (first + last - 1) / 2
            if region.strand == "+"
            else (region.end - (first + last + 1) / 2)
        )
        distance = (midpoint - int(anchor["applied_anchor"])) * (
            1 if region.strand == "+" else -1
        )
        if -upstream <= distance < 0:
            bin_start = math.floor(distance / bin_width) * bin_width
            counts[
                (motif, labels[identifier], anchor["anchor_source"], bin_start)
            ] += 1
    opportunities: dict[tuple[str, str, str, int], int] = defaultdict(int)
    parents: dict[tuple[str, str, str, int], int] = defaultdict(int)
    bins = range(-math.ceil(upstream / bin_width) * bin_width, 0, bin_width)
    for motif, record in by_motif.items():
        width = len(record["consensus"])
        for identifier, region in by_region.items():
            if identifier not in labels or labels[identifier] not in {
                "positive",
                "negative",
            }:
                raise ValueError("RNA motif class metadata is missing")
            anchor = by_anchor.get(identifier)
            if anchor is None or anchor.get("applied_anchor") is None:
                raise ValueError("RNA motif source lacks an applied anchor")
            if width < 1:
                raise ValueError("RNA motif consensus has no width")
            possible = max(0, len(region.sequence) - width + 1)
            midpoint = (
                region.start + (width - 1) / 2
                if region.strand == "+"
                else region.end - 1 - (width - 1) / 2
            )
            first_distance = (midpoint - int(anchor["applied_anchor"])) * (
                1 if region.strand == "+" else -1
            )
            for left in bins:
                low = max(0, math.ceil(max(left, -upstream) - first_distance))
                high = min(
                    possible, math.ceil(left + bin_width - first_distance)
                )
                key = (
                    motif,
                    labels[identifier],
                    anchor["anchor_source"],
                    left,
                )
                if high > low:
                    opportunities[key] += high - low
                    parents[key] += 1
    return [
        {
            "motif_id": motif,
            "label": label,
            "anchor_source": source,
            "bin_start": max(left, -upstream),
            "bin_end": left + bin_width,
            "site_count": counts[(motif, label, source, left)],
            "covered_parents": parents[(motif, label, source, left)],
            "eligible_motif_centres": opportunity,
            "sites_per_1000_eligible_centres": 1000
            * counts[(motif, label, source, left)]
            / opportunity,
            "position_test_p_value": None,
        }
        for (motif, label, source, left), opportunity in sorted(
            opportunities.items()
        )
    ]


def plot_rna_results(
    *,
    profiles: Sequence[Mapping[str, Any]],
    starts: Sequence[Mapping[str, Any]],
    calls: Sequence[Mapping[str, Any]],
    shifts: Sequence[Mapping[str, Any]],
    de_rows: Sequence[Mapping[str, Any]],
    motif_profiles: Sequence[Mapping[str, Any]],
    labels: Mapping[str, str],
    directory: Path,
    counts: Mapping[str, Mapping[str, int]] | None = None,
    size_factors: Mapping[str, float] | None = None,
    max_loci: int = 6,
) -> list[Path]:
    """Plot RNA evidence without concealing missing coverage or uncertainty.

    Args:
        profiles: Binned per-replicate raw/divided coverage.
        starts: Consensus starts and non-detections.
        calls: Independent replicate edge records.
        shifts: Descriptive between-condition shifts.
        de_rows: Tested DE results, if available.
        motif_profiles: Actual sites per eligible upstream centre.
        labels: Gene ID to selected class.
        directory: Figure destination; PNG and vector PDF are written.
        counts: Optional raw integer counts for a visual expression heatmap.
        size_factors: Supplied/fitted count divisors; default one means raw.
        max_loci: Maximum locus figures, balanced between supplied classes.

    Returns:
        PNG paths with companion PDF files. Empty inputs produce no plots.

    Raises:
        ValueError: Locus bounds or size factors are invalid.
    """
    if (
        isinstance(max_loci, bool)
        or not 0 <= max_loci <= 20
        or any(
            not math.isfinite(value) or value <= 0
            for value in (size_factors or {}).values()
        )
    ):
        raise ValueError("Invalid RNA plot limits or sample factors")
    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    images: list[Path] = []

    def finish(*, figure: Any, name: str) -> None:
        figure.tight_layout()
        images.append(
            save_figure(figure=figure, directory=directory, name=name)
        )

    if profiles:
        genes = sorted(
            {row["gene_id"] for row in profiles},
            key=lambda gene: (labels.get(gene, "unselected"), gene),
        )[:60]
        conditions = sorted({row["condition"] for row in profiles})
        columns = sorted({int(row["relative_start"]) for row in profiles})
        pairs = [
            (gene, condition) for condition in conditions for gene in genes
        ]
        values: dict[tuple[str, str, int], list[float]] = defaultdict(list)
        for row in profiles:
            values[
                (row["gene_id"], row["condition"], int(row["relative_start"]))
            ].append(float(row["divided_depth"]))
        matrix = np.full((len(pairs), len(columns)), np.nan)
        for i, (gene, condition) in enumerate(pairs):
            for j, left in enumerate(columns):
                observed = values.get((gene, condition, left), [])
                if observed:
                    matrix[i, j] = math.log1p(sum(observed) / len(observed))
        figure, axis = plt.subplots(
            figsize=(12, max(4, min(24, len(pairs) * 0.19)))
        )
        heat = axis.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=plt.get_cmap("viridis").with_extremes(bad="#e5e7eb"),
        )
        ticks = np.linspace(
            0, len(columns) - 1, min(9, len(columns)), dtype=int
        )
        axis.set(
            xticks=ticks,
            xticklabels=[columns[j] for j in ticks],
            yticks=np.arange(len(pairs)),
            yticklabels=[
                f"{gene} · {condition} · {labels.get(gene, 'unselected')}"
                for gene, condition in pairs
            ],
            xlabel="Distance from annotated 5-prime gene base (bp)",
            title="RNA coverage · supplied/fitted depth divisors",
        )
        axis.tick_params(axis="y", labelsize=6)
        figure.colorbar(heat, ax=axis, label="log1p mean divided depth")
        finish(figure=figure, name="rna_coverage_heatmap")
    if starts:
        categories: dict[tuple[str, str], int] = defaultdict(int)
        for row in starts:
            categories[(row["condition"], row["status"])] += 1
        conditions = sorted({key[0] for key in categories})
        statuses = sorted({key[1] for key in categories})
        figure, axis = plt.subplots(figsize=(10, 5))
        bottom = np.zeros(len(conditions))
        for status in statuses:
            bar_values = np.array(
                [categories[(condition, status)] for condition in conditions]
            )
            axis.bar(
                conditions,
                bar_values,
                bottom=bottom,
                label=status.replace("_", " "),
            )
            bottom += bar_values
        axis.set(
            ylabel="Gene-condition observations",
            title="Biological replicate support for candidate starts",
        )
        axis.legend(
            frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1, 1)
        )
        finish(figure=figure, name="rna_start_support")
    estimable = [row for row in shifts if row["oriented_shift_bp"] is not None]
    if estimable:
        figure, axis = plt.subplots(figsize=(9, 4))
        axis.hist(
            [row["oriented_shift_bp"] for row in estimable],
            bins=min(30, max(1, len(estimable))),
            color="#187a8a",
        )
        axis.axvline(x=0, color="#c67638", linestyle="--")
        axis.set(
            xlabel=(
                "Condition − reference start shift (bp; downstream positive)"
            ),
            ylabel="Genes with supported candidates in both conditions",
            title="Candidate start shifts · descriptive, no usage p-values",
        )
        finish(figure=figure, name="rna_start_shifts")
    tested = [
        row
        for row in de_rows
        if row.get("padj") is not None
        and row.get("log2FoldChange") is not None
    ]
    if tested:
        figure, axis = plt.subplots(figsize=(9, 5))
        axis.scatter(
            [row["log2FoldChange"] for row in tested],
            [
                min(300, -math.log10(max(float(row["padj"]), 1e-300)))
                for row in tested
            ],
            c=[
                "#187a8a"
                if labels.get(row["gene_id"]) == "positive"
                else "#c67638"
                if labels.get(row["gene_id"]) == "negative"
                else "#a8adb3"
                for row in tested
            ],
            s=18,
            alpha=0.7,
        )
        axis.axvline(x=0, color="#b6b6b6", linewidth=0.7)
        axis.set(
            xlabel="Condition/reference log2 fold change",
            ylabel="−log10 imported/fitted adjusted p-value (cap 300)",
            title="Gene differential expression · recorded contrast",
        )
        finish(figure=figure, name="rna_de_volcano")
    if counts:
        samples = list(next(iter(counts.values())))
        genes = list(counts)
        if samples and genes:
            matrix = np.array(
                [
                    [
                        math.log1p(
                            counts[gene][sample]
                            / (size_factors or {}).get(sample, 1)
                        )
                        for sample in samples
                    ]
                    for gene in genes
                ]
            )
            chosen = np.argsort(np.var(matrix, axis=1))[
                -min(40, len(genes)) :
            ][::-1]
            centred = matrix[chosen] - matrix[chosen].mean(
                axis=1, keepdims=True
            )
            figure, axis = plt.subplots(
                figsize=(
                    max(7, len(samples) * 0.5),
                    max(4, len(chosen) * 0.22),
                )
            )
            bound = max(0.1, float(np.max(np.abs(centred))))
            heat = axis.imshow(
                centred, aspect="auto", cmap="RdBu_r", vmin=-bound, vmax=bound
            )
            axis.set(
                xticks=range(len(samples)),
                xticklabels=samples,
                yticks=range(len(chosen)),
                yticklabels=[genes[int(j)] for j in chosen],
                title="Expression overview · row-centred log1p divided counts",
            )
            axis.tick_params(axis="x", rotation=45, labelsize=7)
            axis.tick_params(axis="y", labelsize=7)
            figure.colorbar(
                heat,
                ax=axis,
                label="Difference from gene mean log1p divided count",
            )
            finish(figure=figure, name="rna_expression_heatmap")
    if motif_profiles:
        rows = sorted(
            {
                (r["motif_id"], r["label"], r["anchor_source"])
                for r in motif_profiles
            }
        )
        columns = sorted({int(r["bin_start"]) for r in motif_profiles})
        lookup = {
            (
                r["motif_id"],
                r["label"],
                r["anchor_source"],
                int(r["bin_start"]),
            ): float(r["sites_per_1000_eligible_centres"])
            for r in motif_profiles
        }
        matrix = np.array(
            [
                [lookup.get((*row, left), np.nan) for left in columns]
                for row in rows
            ]
        )
        figure, axis = plt.subplots(figsize=(12, max(4, len(rows) * 0.35)))
        heat = axis.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=plt.get_cmap("magma").with_extremes(bad="#e5e7eb"),
        )
        ticks = np.linspace(
            0, len(columns) - 1, min(9, len(columns)), dtype=int
        )
        axis.set(
            xticks=ticks,
            xticklabels=[columns[j] for j in ticks],
            yticks=range(len(rows)),
            yticklabels=[" · ".join(row) for row in rows],
            xlabel="Distance from applied candidate/fallback start (bp)",
            title="Upstream motif positions · discovery-selected, descriptive",
        )
        axis.tick_params(axis="y", labelsize=7)
        figure.colorbar(
            heat, ax=axis, label="Sites per 1,000 eligible motif centres"
        )
        finish(figure=figure, name="rna_motif_distance_heatmap")
    loci: list[str] = []
    classes = sorted(set(labels.values())) or ["unselected"]
    available = {row["gene_id"] for row in profiles}
    ordered = {
        label: sorted(
            gene
            for gene in available
            if labels.get(gene, "unselected") == label
        )
        for label in classes
    }
    while len(loci) < max_loci and any(ordered.values()):
        for label in classes:
            if ordered[label] and len(loci) < max_loci:
                loci.append(ordered[label].pop(0))
    for number, gene_id in enumerate(loci, start=1):
        selected = [row for row in profiles if row["gene_id"] == gene_id]
        samples = sorted({row["sample_id"] for row in selected})
        conditions = sorted({row["condition"] for row in selected})
        colours = {
            condition: plt.get_cmap("tab10")(i % 10)
            for i, condition in enumerate(conditions)
        }
        figure, axes = plt.subplots(
            nrows=2,
            figsize=(11, 6),
            sharex=True,
            gridspec_kw={"height_ratios": [3, 2]},
        )
        condition_values: dict[tuple[str, float], list[float]] = defaultdict(
            list
        )
        for sample in samples:
            locus_rows = sorted(
                [row for row in selected if row["sample_id"] == sample],
                key=lambda r: r["relative_start"],
            )
            condition = locus_rows[0]["condition"]
            positions = [
                (row["relative_start"] + row["relative_end"]) / 2
                for row in locus_rows
            ]
            depths = [row["divided_depth"] for row in locus_rows]
            axes[0].plot(
                positions,
                depths,
                color=colours[condition],
                linewidth=0.8,
                alpha=0.35,
            )
            for position, value in zip(positions, depths, strict=True):
                condition_values[(condition, position)].append(value)
        for condition in conditions:
            positions = sorted(
                position
                for level, position in condition_values
                if level == condition
            )
            axes[0].plot(
                positions,
                [
                    sum(condition_values[(condition, position)])
                    / len(condition_values[(condition, position)])
                    for position in positions
                ],
                color=colours[condition],
                label=condition,
                linewidth=2,
            )
        for i, sample in enumerate(samples):
            found = next(
                (
                    row
                    for row in calls
                    if row["gene_id"] == gene_id and row["sample_id"] == sample
                ),
                None,
            )
            if found and found["candidate_start"] is not None:
                sign = 1 if found["strand"] == "+" else -1
                point = (
                    found["candidate_start"] - found["annotation_anchor"]
                ) * sign
                low, high = sorted(
                    [
                        (found[name] - found["annotation_anchor"]) * sign
                        for name in (
                            "neighbourhood_start",
                            "neighbourhood_end",
                        )
                    ]
                )
                axes[1].plot(
                    [low, high],
                    [i, i],
                    color=colours[found["condition"]],
                    linewidth=3,
                )
                axes[1].scatter(
                    [point], [i], color=colours[found["condition"]], s=22
                )
        for axis in axes:
            axis.axvline(x=0, color="#b6b6b6", linestyle="--", linewidth=1)
        axes[0].set(
            ylabel="Depth / recorded sample divisor",
            title=(
                f"{gene_id} · {labels.get(gene_id, 'unselected')}"
                " · candidate start evidence"
            ),
        )
        axes[0].legend(frameon=False)
        axes[1].set(
            yticks=range(len(samples)),
            yticklabels=samples,
            xlabel="Distance from annotated 5-prime gene base (bp)",
            ylabel="Biological sample",
        )
        finish(figure=figure, name=f"rna_locus_{number:02d}")
    return images
