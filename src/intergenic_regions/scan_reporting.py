"""Genome scan graphics with explicit opportunity denominators."""

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from intergenic_regions.reporting import save_figure, write_report


def plot_scan(
    *,
    profiles: Sequence[Mapping[str, Any]],
    burden: Sequence[Mapping[str, Any]],
    genes: Sequence[Mapping[str, Any]],
    directory: Path,
) -> list[Path]:
    """Plot positional density, cohort profiles, burden and mismatches.

    Args:
        profiles: Per-motif/cohort distance-bin counts and denominators.
        burden: Complete counts by motif, contig and mismatch number.
        genes: Complete nearest-gene motif counts.
        directory: Figure destination for paired PNG/PDF graphics.

    Returns:
        PNG paths for non-empty summaries; empty scans need no invented plots.
    """
    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt

    images = []
    if profiles and any(int(row["eligible_windows"]) > 0 for row in profiles):
        grouped: dict[tuple[str, float, float], list[float]] = defaultdict(
            lambda: [0.0, 0.0]
        )
        for row in profiles:
            key = (
                row["motif_id"],
                float(row["bin_start"]),
                float(row["bin_end"]),
            )
            grouped[key][0] += float(row["site_count"])
            grouped[key][1] += float(row["eligible_windows"])
        motifs = sorted({key[0] for key in grouped})
        bins = sorted({(key[1], key[2]) for key in grouped})
        counts = np.asarray(
            [[grouped[m, a, b][0] for a, b in bins] for m in motifs]
        )
        denominators = np.asarray(
            [[grouped[m, a, b][1] for a, b in bins] for m in motifs]
        )
        density = np.divide(
            1e6 * counts,
            denominators,
            out=np.full_like(counts, np.nan),
            where=denominators > 0,
        )
        for values, name, title, units in (
            (
                density,
                "genome_distance_density",
                "Motif density around annotated gene-start proxies",
                "Sites per 10⁶ eligible windows",
            ),
            (
                counts,
                "genome_distance_counts",
                "Raw motif counts around annotated gene-start proxies",
                "Physical sites",
            ),
        ):
            figure, axis = plt.subplots(
                figsize=(12, max(3, len(motifs) * 0.38))
            )
            mesh = axis.pcolormesh(
                [bins[0][0], *[b for _, b in bins]],
                np.arange(len(motifs) + 1) - 0.5,
                values,
                cmap="magma",
                shading="flat",
            )
            axis.axvline(x=0, color="#66ddd0", linestyle="--", linewidth=1)
            axis.set(
                yticks=np.arange(len(motifs)),
                yticklabels=motifs,
                xlabel=(
                    "Signed midpoint distance (bp): upstream < 0; "
                    "downstream > 0"
                ),
                title=title,
            )
            figure.colorbar(mappable=mesh, ax=axis, label=units)
            figure.tight_layout()
            images.append(
                save_figure(figure=figure, directory=directory, name=name)
            )
        first_motif = max(
            motifs, key=lambda m: sum(grouped[m, a, b][0] for a, b in bins)
        )
        figure, axis = plt.subplots(figsize=(11, 4))
        for cohort, colour in (
            ("positive", "#187a8a"),
            ("negative", "#c67638"),
            ("other", "#6c7686"),
        ):
            records = sorted(
                (
                    r
                    for r in profiles
                    if r["motif_id"] == first_motif and r["cohort"] == cohort
                ),
                key=lambda r: float(r["bin_start"]),
            )
            if any(int(r["eligible_windows"]) for r in records):
                axis.plot(
                    [
                        (float(r["bin_start"]) + float(r["bin_end"])) / 2
                        for r in records
                    ],
                    [
                        r["hits_per_million_windows"]
                        if r["hits_per_million_windows"] is not None
                        else np.nan
                        for r in records
                    ],
                    label=cohort,
                    color=colour,
                )
        axis.axvline(x=0, color="#acb6bd", linestyle="--")
        axis.set(
            xlabel="Signed distance from annotated gene-start proxy (bp)",
            ylabel="Sites per million eligible start windows",
            title=(
                f"{first_motif}: descriptive profiles; "
                "independent validation required"
            ),
        )
        axis.legend(frameon=False)
        figure.tight_layout()
        images.append(
            save_figure(
                figure=figure,
                directory=directory,
                name="genome_cohort_profiles",
            )
        )
    if genes:
        totals: Counter[str] = Counter()
        for row in genes:
            totals[row["gene_id"]] += int(row["site_count"])
        leading = [
            g
            for g, _ in sorted(
                totals.items(), key=lambda pair: (-pair[1], pair[0])
            )[:30]
        ]
        motifs = sorted({r["motif_id"] for r in genes})
        lookup = {
            (r["gene_id"], r["motif_id"]): int(r["site_count"]) for r in genes
        }
        matrix = np.asarray(
            [[lookup.get((g, m), 0) for m in motifs] for g in leading]
        )
        figure, axis = plt.subplots(
            figsize=(max(8, len(motifs) * 0.6), max(4, len(leading) * 0.26))
        )
        image = axis.imshow(np.log1p(matrix), aspect="auto", cmap="viridis")
        axis.set(
            yticks=np.arange(len(leading)),
            yticklabels=leading,
            xticks=np.arange(len(motifs)),
            xticklabels=motifs,
            title="Nearest-gene motif burden: top 30 genes (descriptive)",
        )
        axis.tick_params(axis="x", rotation=65)
        figure.colorbar(
            mappable=image, ax=axis, label="log(1 + physical site count)"
        )
        figure.tight_layout()
        images.append(
            save_figure(
                figure=figure,
                directory=directory,
                name="genome_gene_motif_heatmap",
            )
        )
    if burden:
        mismatches: Counter[int] = Counter()
        contigs: Counter[str] = Counter()
        for row in burden:
            mismatches[int(row["mismatches"])] += int(row["site_count"])
            contigs[row["contig"]] += int(row["site_count"])
        figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
        axes[0].bar(
            list(sorted(mismatches)),
            [mismatches[d] for d in sorted(mismatches)],
            color="#187a8a",
        )
        axes[0].set(
            xlabel="Substitutions (minimum across matching orientations)",
            ylabel="Physical sites",
            title="Exact and approximate matches",
        )
        top_contigs = sorted(contigs, key=lambda c: (-contigs[c], c))[:15]
        axes[1].barh(
            top_contigs[::-1],
            [contigs[c] for c in top_contigs[::-1]],
            color="#c67638",
        )
        axes[1].set(
            xlabel="Physical sites (not normalised by contig length)",
            title="Site counts by contig: top 15",
        )
        figure.tight_layout()
        images.append(
            save_figure(
                figure=figure, directory=directory, name="genome_scan_burden"
            )
        )
    return images


def scan_report(
    *,
    directory: Path,
    summary: Mapping[str, Any],
    sites: Sequence[Mapping[str, Any]],
    profiles: Sequence[Mapping[str, Any]],
    burden: Sequence[Mapping[str, Any]],
    genes: Sequence[Mapping[str, Any]],
) -> None:
    """Write the offline screen, graphics and complete-result links.

    Args:
        directory: Existing scan result directory.
        summary: Audited complete scan summary.
        sites: At most 500 site preview rows; the full TSV is streamed to disk.
        profiles: Complete signed distance-bin summaries.
        burden: Complete motif/contig/substitution counts.
        genes: Complete nearest-gene motif counts.
    """
    images = plot_scan(
        profiles=profiles,
        burden=burden,
        genes=genes,
        directory=directory / "figures",
    )
    write_report(
        path=directory / "report.html",
        title="Genome-wide motif screening",
        summary={"scan": dict(summary)},
        tables={
            "Genome motif site preview": sites,
            "Distance profiles": profiles,
            "Nearest-gene motif burden": genes,
        },
        images=images,
        links={
            "Complete genome sites TSV": "genome_motif_sites.tsv",
            "Genome browser BED": "genome_motif_sites.bed",
            "Distance counts and denominators TSV": "distance_profiles.tsv",
            "Selected consensus targets TSV": "scan_motifs.tsv",
            "Nearest-gene motif counts TSV": "gene_motif_counts.tsv",
        },
        notes=[
            str(summary["interpretation"]),
            str(summary["distance_anchor"]),
            str(summary["profile_normalisation"]),
            "The site table previews at most 500 records; complete counts "
            "and all sites are in the TSV/BED. Gene-start ties are flagged "
            "and excluded from position profiles. Unknown gene strands "
            "have no start proxy.",
            "Discovery-selected motifs and discovery gene cohorts are "
            "reused descriptively. Genome match counts and position "
            "heatmaps do not constitute an independent enrichment test.",
        ],
    )
