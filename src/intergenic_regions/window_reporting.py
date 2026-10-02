"""Offline graphics for search scales, source validation and locus support."""

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from intergenic_regions.reporting import save_figure


def plot_regional_results(
    *,
    rows: Sequence[Mapping[str, Any]],
    audit: Sequence[Mapping[str, Any]],
    candidates: Sequence[Mapping[str, Any]],
    profiles: Sequence[Mapping[str, Any]],
    scales: Sequence[Mapping[str, Any]],
    directory: Path,
    max_loci: int = 6,
) -> list[Path]:
    """Plot availability, scale comparisons, position maps and candidate loci.

    Args:
        rows: Full window scores and optional local evidence.
        audit: Availability per source and scale.
        candidates: Ranked unions of supported search windows.
        profiles: Parent-weighted position summaries.
        scales: Per-scale source-level validation metrics and fitting status.
        directory: PNG/PDF destination.
        max_loci: Maximum illustrated loci, zero to 20.

    Returns:
        Embedded-report PNG paths, each with a companion vector PDF.

    Raises:
        ValueError: The locus limit is invalid.
    """
    if (
        not isinstance(max_loci, int)
        or isinstance(max_loci, bool)
        or not 0 <= max_loci <= 20
    ):
        raise ValueError("Locus plot limit must be between zero and 20")
    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    paths: list[Path] = []
    widths = sorted({r["window_length"] for r in audit})
    colours = {"positive": "#187a8a", "negative": "#c67638"}
    figure, axis = plt.subplots(figsize=(10, 4))
    for label, shift in (("positive", -0.18), ("negative", 0.18)):
        counts = Counter(
            r["window_length"]
            for r in audit
            if r["label"] == label and r["windows"] > 0
        )
        axis.bar(
            np.arange(len(widths)) + shift,
            [counts[w] for w in widths],
            width=0.35,
            color=colours[label],
            label=label.title(),
        )
    axis.set(
        xticks=np.arange(len(widths)),
        xticklabels=[str(w) for w in widths],
        xlabel="Search-window length (bp)",
        ylabel="Source flanks with complete windows",
        title="Scale availability · short gaps are never padded",
    )
    axis.legend(frameon=False)
    figure.tight_layout()
    paths.append(
        save_figure(
            figure=figure, directory=directory, name="region_scale_coverage"
        )
    )
    fitted = [s for s in scales if s["status"] == "completed"]
    if fitted:
        figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
        for axis, field, baseline, title in (
            (axes[0], "roc_auc", "baseline_roc_auc", "Source-class ROC AUC"),
            (
                axes[1],
                "average_precision",
                "baseline_average_precision",
                "Source-class average precision",
            ),
        ):
            positions = [s["window_length"] for s in fitted]
            axis.plot(
                positions,
                [s[field] for s in fitted],
                marker="o",
                color="#187a8a",
                label="Sequence words + composition",
            )
            axis.plot(
                positions,
                [s[baseline] for s in fitted],
                marker="o",
                color="#c67638",
                label="Composition baseline",
            )
            axis.set(
                xscale="log",
                xticks=positions,
                xticklabels=[str(n) for n in positions],
                xlabel="Window length (bp)",
                ylabel=title,
                ylim=(-0.02, 1.02),
                title=title + " · one mean score per source",
            )
            axis.legend(frameon=False, fontsize=8)
        figure.tight_layout()
        paths.append(
            save_figure(
                figure=figure,
                directory=directory,
                name="region_scale_validation",
            )
        )
    if profiles:
        profile_widths = sorted({r["window_length"] for r in profiles})
        boundaries = sorted(
            {r["bin_start"] for r in profiles}
            | {r["bin_end"] for r in profiles}
        )
        column = {left: i for i, left in enumerate(boundaries[:-1])}
        row_index = {width: i for i, width in enumerate(profile_widths)}
        coordinate = profiles[0]["coordinate_system"]
        xlabel = (
            "Window midpoint from annotated gene 5′ base (bp; TSS proxy)"
            if coordinate == "annotated_gene_5prime_base"
            else "Window midpoint in oriented source sequence (bp)"
        )
        for field, name, title, maximum in (
            (
                "mean_held_out_signature_score",
                "region_score_distance",
                "Mean held-out sequence resemblance",
                1.0,
            ),
            (
                "mean_motif_sites_per_kb",
                "region_motif_distance",
                "Descriptive enriched-motif sites per kb",
                None,
            ),
            (
                "parents",
                "region_position_coverage",
                "Sources contributing to each position bin",
                None,
            ),
        ):
            if not any(r[field] is not None for r in profiles):
                continue
            figure, axes = plt.subplots(
                nrows=2,
                ncols=1,
                figsize=(11, 7),
                sharex=True,
                constrained_layout=True,
            )
            upper = (
                maximum
                if maximum is not None
                else max(1, max(r[field] for r in profiles))
            )
            cmap = plt.get_cmap("viridis").with_extremes(bad="#e6eaed")
            for axis, label in zip(
                axes, ("positive", "negative"), strict=True
            ):
                matrix = np.full(
                    (len(profile_widths), len(boundaries) - 1), np.nan
                )
                for profile in profiles:
                    if (
                        profile["label"] == label
                        and profile[field] is not None
                    ):
                        matrix[
                            row_index[profile["window_length"]],
                            column[profile["bin_start"]],
                        ] = profile[field]
                mesh = axis.pcolormesh(
                    boundaries,
                    np.arange(len(profile_widths) + 1),
                    np.ma.masked_invalid(matrix),
                    cmap=cmap,
                    vmin=0,
                    vmax=upper,
                    shading="flat",
                )
                axis.set(
                    yticks=np.arange(len(profile_widths)) + 0.5,
                    yticklabels=[str(w) for w in profile_widths],
                    ylabel="Window length (bp)",
                    title=label.title() + " source flanks",
                )
                figure.colorbar(mesh, ax=axis, label=title)
            axes[-1].set_xlabel(xlabel)
            figure.suptitle(
                title + " · average within each source, then across sources"
            )
            paths.append(
                save_figure(figure=figure, directory=directory, name=name)
            )
    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
    if candidates:
        axes[0].hist(
            [r["length"] for r in candidates],
            bins=min(15, len(candidates)),
            color="#187a8a",
        )
        counts = Counter(
            width
            for r in candidates
            for width in r["supporting_lengths"].split(",")
        )
        sizes = sorted(counts, key=int)
        axes[1].bar(sizes, [counts[w] for w in sizes], color="#187a8a")
    else:
        axes[0].text(
            0.5,
            0.5,
            "No above-threshold window unions",
            ha="center",
            va="center",
            transform=axes[0].transAxes,
        )
    axes[0].set(
        xlabel="Union length (bp; exploratory boundaries)",
        ylabel="Candidate unions",
        title="Candidate span distribution",
    )
    axes[1].set(
        xlabel="Supporting search scale (bp)",
        ylabel="Candidate unions",
        title="Scale support · not independent evidence",
    )
    figure.tight_layout()
    paths.append(
        save_figure(
            figure=figure, directory=directory, name="region_candidate_lengths"
        )
    )
    by_parent: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        by_parent[row["parent_id"]].append(row)
    order = list(dict.fromkeys([r["parent_id"] for r in candidates]))
    order.extend(p for p in sorted(by_parent) if p not in order)
    for number, parent in enumerate(order[:max_loci], start=1):
        items = by_parent[parent]
        widths = sorted({r["window_length"] for r in items})
        loci = [r for r in candidates if r["parent_id"] == parent]
        accessible = any(
            r.get("accessibility_overlap_fraction") is not None for r in items
        )
        figure, axes = plt.subplots(
            nrows=4 if accessible else 3,
            ncols=1,
            figsize=(11, 8 if accessible else 7),
            sharex=True,
            constrained_layout=True,
            gridspec_kw={
                "height_ratios": [0.7, 2, 1, 1] if accessible else [0.7, 2, 1]
            },
        )
        for candidate in loci:
            axes[0].broken_barh(
                [(candidate["sequence_start"], candidate["length"])],
                (0, 1),
                facecolors="#187a8a",
                alpha=0.55,
            )
        axes[0].set(
            ylim=(0, 1),
            yticks=[],
            title=parent + " · candidate window unions (unvalidated)",
        )
        if not loci:
            axes[0].text(
                0.5,
                0.5,
                "No above-threshold union",
                ha="center",
                va="center",
                transform=axes[0].transAxes,
            )
        cmap = plt.get_cmap("viridis").with_extremes(bad="#e6eaed")
        axes[1].set_facecolor("#e6eaed")
        for index, width in enumerate(widths):
            selected = sorted(
                [r for r in items if r["window_length"] == width],
                key=lambda r: r["sequence_start"],
            )
            centers = [
                (r["sequence_start"] + r["sequence_end"]) / 2 for r in selected
            ]
            if len(centers) == 1:
                boundaries = [centers[0] - 0.5, centers[0] + 0.5]
            else:
                boundaries = [
                    centers[0] - (centers[1] - centers[0]) / 2,
                    *[
                        (a + b) / 2
                        for a, b in zip(centers[:-1], centers[1:], strict=True)
                    ],
                    centers[-1] + (centers[-1] - centers[-2]) / 2,
                ]
            scores = np.asarray(
                [[r["held_out_signature_score"] for r in selected]],
                dtype=float,
            )
            mesh = axes[1].pcolormesh(
                boundaries,
                [index, index + 1],
                np.ma.masked_invalid(scores),
                cmap=cmap,
                vmin=0,
                vmax=1,
            )
        axes[1].set(
            yticks=np.arange(len(widths)) + 0.5,
            yticklabels=[str(w) for w in widths],
            ylabel="Window length (bp)",
            title="Window midpoint resemblance · grey = unavailable",
        )
        figure.colorbar(mesh, ax=axes[1], label="Source-class resemblance")
        for width in widths:
            selected = sorted(
                [r for r in items if r["window_length"] == width],
                key=lambda r: r["sequence_start"],
            )
            x = [
                (r["sequence_start"] + r["sequence_end"]) / 2 for r in selected
            ]
            axes[2].plot(
                x,
                [r["motif_sites_per_kb"] for r in selected],
                label=f"{width} bp",
            )
            if accessible:
                axes[3].plot(
                    x,
                    [
                        r.get("accessibility_overlap_fraction")
                        for r in selected
                    ],
                    label=f"{width} bp",
                )
        axes[2].set(
            ylabel="Motif sites / kb",
            title="Descriptive motif density by search scale",
        )
        axes[2].legend(frameon=False, ncol=min(5, len(widths)), fontsize=8)
        if accessible:
            axes[3].set(
                ylabel="Peak overlap fraction",
                ylim=(-0.02, 1.02),
                title="Supplied accessibility overlap · contextual support",
            )
        axes[-1].set(
            xlabel="Transcription-oriented source offset (bp)",
            xlim=(0, max(r["sequence_end"] for r in items)),
        )
        first = items[0]
        if first["distance_to_gene_start"] is not None:
            offset = (
                first["distance_to_gene_start"]
                - (first["sequence_start"] + first["sequence_end"]) / 2
            )
            secondary = axes[0].secondary_xaxis(
                "top",
                functions=(
                    lambda x, offset=offset: x + offset,
                    lambda x, offset=offset: x - offset,
                ),
            )
            secondary.set_xlabel(
                "Distance from annotated gene 5′ base (bp; TSS proxy)"
            )
        paths.append(
            save_figure(
                figure=figure,
                directory=directory,
                name=f"region_locus_{number:02d}",
            )
        )
    return paths
