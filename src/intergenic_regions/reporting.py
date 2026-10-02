"""Portable HTML reports and publication-friendly PNG/PDF graphics."""

import base64
import html
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from intergenic_regions.genome import sequence_composition


def save_figure(*, figure: Any, directory: Path, name: str) -> Path:
    """Save a plot in PNG and vector PDF and close its resources.

    Args:
        figure: Matplotlib figure.
        directory: Destination directory.
        name: Output basename without an extension.

    Returns:
        PNG path suitable for the HTML report.
    """
    import matplotlib.pyplot as plt

    directory.mkdir(parents=True, exist_ok=True)
    png = directory / f"{name}.png"
    figure.savefig(fname=png, dpi=180, bbox_inches="tight")
    figure.savefig(fname=directory / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig=figure)
    return png


def plot_composition(
    *,
    positive: Mapping[str, str],
    negative: Mapping[str, str],
    directory: Path,
) -> Path:
    """Plot control/foreground length and GC distributions.

    Args:
        positive: Foreground sequences.
        negative: Control sequences.
        directory: Plot destination.

    Returns:
        PNG plot path.
    """
    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
    for label, records, colour in (
        ("Positive", positive, "#187a8a"),
        ("Negative", negative, "#c67638"),
    ):
        composition = [
            sequence_composition(sequence=s) for s in records.values()
        ]
        lengths = [c["length"] for c in composition]
        gc = [c["gc_fraction"] for c in composition]
        axes[0].scatter(
            lengths, gc, label=label, color=colour, alpha=0.65, s=18
        )
        axes[1].hist(lengths, bins=20, alpha=0.5, color=colour, label=label)
    axes[0].set(
        xlabel="Region length (bp)",
        ylabel="GC fraction (ACGT bases)",
        title="Composition comparison",
    )
    axes[1].set(
        xlabel="Region length (bp)",
        ylabel="Sequences",
        title="Length comparison",
    )
    for axis in axes:
        axis.legend(frameon=False)
    figure.tight_layout()
    return save_figure(
        figure=figure, directory=directory, name="background_composition"
    )


def plot_motif_results(
    *, rows: Sequence[Mapping[str, Any]], directory: Path
) -> list[Path]:
    """Plot sequence prevalence and enrichment/FDR relationship.

    Args:
        rows: Sorted native enrichment table.
        directory: Plot destination.

    Returns:
        PNG paths; matching PDF files are saved alongside them.
    """
    import math

    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if not rows:
        return []
    leading = list(rows[:12])[::-1]
    figure, axis = plt.subplots(figsize=(10, max(4, len(leading) * 0.4)))
    positions = np.arange(len(leading))
    axis.barh(
        positions - 0.18,
        [100 * r["positive_fraction"] for r in leading],
        height=0.35,
        color="#187a8a",
        label="Positive",
    )
    axis.barh(
        positions + 0.18,
        [100 * r["negative_fraction"] for r in leading],
        height=0.35,
        color="#c67638",
        label="Negative",
    )
    axis.set(
        yticks=positions,
        yticklabels=[
            f"{r['motif_id']}  q={r['q_value']:.2g}" for r in leading
        ],
        xlabel="Sequences containing motif (%)",
        title="Leading sequence signatures",
    )
    axis.legend(frameon=False)
    figure.tight_layout()
    first = save_figure(
        figure=figure, directory=directory, name="motif_prevalence"
    )
    figure, axis = plt.subplots(figsize=(8, 5))
    x = [math.log2(r["odds_ratio_haldane"]) for r in rows]
    y = [r["minus_log10_q"] for r in rows]
    axis.scatter(
        x,
        y,
        c=["#187a8a" if r["q_value"] <= 0.05 else "#a8adb3" for r in rows],
        s=14,
        alpha=0.6,
    )
    axis.axhline(
        y=-math.log10(0.05), color="#c67638", linestyle="--", label="q = 0.05"
    )
    axis.axvline(x=0, color="#c6c6c6", linewidth=0.8)
    axis.set(
        xlabel="log2 odds ratio (Haldane effect estimate)",
        ylabel="-log10 BH q-value",
        title="Sequence-level motif enrichment",
    )
    axis.legend(frameon=False)
    figure.tight_layout()
    return [
        first,
        save_figure(
            figure=figure, directory=directory, name="motif_enrichment"
        ),
    ]


def plot_logos(
    *,
    rows: Sequence[Mapping[str, Any]],
    motifs: Sequence[Any],
    directory: Path,
) -> Path | None:
    """Draw DNA information-content logos for the leading motif results.

    Args:
        rows: Sorted enrichment results.
        motifs: Imported known motifs.
        directory: Plot destination.

    Returns:
        PNG path or None when no motifs were tested.
    """
    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D

    from intergenic_regions.motifs import pattern_motif

    if not rows:
        return None
    known = {motif.motif_id: motif for motif in motifs}
    leading = rows[:6]
    figure, axes = plt.subplots(
        nrows=len(leading), figsize=(10, 2 * len(leading)), squeeze=False
    )
    colours = {"A": "#339e65", "C": "#2876b8", "G": "#dc9b32", "T": "#c64b59"}
    letters = {
        base: TextPath(
            xy=(0, 0), s=base, size=1, prop=FontProperties(weight="bold")
        )
        for base in "ACGT"
    }
    for axis, record in zip(axes[:, 0], leading, strict=True):
        motif = known.get(record["motif_id"]) or pattern_motif(
            motif_id=record["motif_id"], pattern=record["consensus"]
        )
        matrix = np.asarray(motif.matrix)
        entropy = -(matrix * np.log2(np.clip(matrix, 1e-300, 1))).sum(axis=1)
        heights = matrix * (2 - entropy)[:, None]
        for position, values in enumerate(heights):
            bottom = 0.0
            for base_index in np.argsort(values):
                height = max(0.0, float(values[base_index]))
                if height < 1e-6:
                    continue
                base = "ACGT"[int(base_index)]
                path = letters[base]
                bounds = path.get_extents()
                transform = (
                    Affine2D()
                    .translate(-bounds.x0, -bounds.y0)
                    .scale(0.85 / bounds.width, height / bounds.height)
                    .translate(position + 0.075, bottom)
                )
                axis.add_patch(
                    PathPatch(
                        path=path,
                        transform=transform + axis.transData,
                        color=colours[base],
                        linewidth=0,
                    )
                )
                bottom += height
        axis.set(
            xlim=(0, len(matrix)),
            ylim=(0, 2.05),
            xticks=np.arange(len(matrix)) + 0.5,
            xticklabels=np.arange(1, len(matrix) + 1),
            ylabel="Bits",
            title=f"{record['motif_id']}  (q={record['q_value']:.2g})",
        )
        axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    return save_figure(figure=figure, directory=directory, name="motif_logos")


def plot_learning(
    *,
    predictions: Sequence[Mapping[str, Any]],
    signatures: Sequence[Mapping[str, Any]],
    summary: Mapping[str, Any],
    directory: Path,
) -> list[Path]:
    """Plot held-out predictive performance and interpretable coefficients.

    Args:
        predictions: Held-out sequence predictions.
        signatures: Ranked k-mer coefficients.
        summary: AI summary including permutation scores.
        directory: Plot destination.

    Returns:
        PNG paths, with accompanying PDF figures.
    """
    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.metrics import precision_recall_curve, roc_curve

    labels = np.asarray([r["label"] for r in predictions])
    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
    for field, label, colour in (
        ("held_out_signature_score", "Sequence model", "#187a8a"),
        (
            "held_out_composition_score",
            "GC/length/ambiguity baseline",
            "#c67638",
        ),
    ):
        scores = np.asarray([r[field] for r in predictions])
        false_positive, true_positive, _ = roc_curve(
            y_true=labels, y_score=scores
        )
        precision, recall, _ = precision_recall_curve(
            y_true=labels, y_score=scores
        )
        axes[0].plot(false_positive, true_positive, label=label, color=colour)
        axes[1].plot(recall, precision, label=label, color=colour)
    axes[0].plot([0, 1], [0, 1], color="#a8adb3", linestyle="--")
    axes[1].axhline(
        y=summary["positive_prevalence"], color="#a8adb3", linestyle="--"
    )
    axes[0].set(
        xlabel="False-positive rate",
        ylabel="True-positive rate",
        title="Held-out ROC",
    )
    axes[1].set(
        xlabel="Recall", ylabel="Precision", title="Held-out precision-recall"
    )
    for axis in axes:
        axis.legend(frameon=False, fontsize=8)
    figure.tight_layout()
    paths = [
        save_figure(figure=figure, directory=directory, name="ai_validation")
    ]
    leading = list(signatures[:10])[::-1]
    figure, axis = plt.subplots(figsize=(10, 5))
    axis.barh(
        [r["kmer"] for r in leading],
        [r["coefficient"] for r in leading],
        color="#187a8a",
    )
    axis.set(
        xlabel="Regularised model coefficient per frequency percentage point",
        title="Positive-class sequence signatures (exploratory)",
    )
    figure.tight_layout()
    paths.append(
        save_figure(figure=figure, directory=directory, name="ai_signatures")
    )
    if summary["permutation_scores"]:
        figure, axis = plt.subplots(figsize=(8, 4))
        axis.hist(summary["permutation_scores"], bins=15, color="#a8adb3")
        axis.axvline(
            x=summary["average_precision"],
            color="#187a8a",
            label="Observed held-out AP",
        )
        axis.set(
            xlabel="Average precision under label permutation",
            ylabel="Replicates",
            title=f"Permutation p = {summary['permutation_p_value']:.3g}",
        )
        axis.legend(frameon=False)
        paths.append(
            save_figure(
                figure=figure, directory=directory, name="ai_permutation"
            )
        )
    return paths


def write_report(
    *,
    path: Path,
    title: str,
    summary: Mapping[str, Any],
    tables: Mapping[str, Sequence[Mapping[str, Any]]] | None = None,
    images: Sequence[Path] = (),
    notes: Sequence[str] = (),
    links: Mapping[str, str] | None = None,
) -> None:
    """Write a self-contained HTML report with escaped user-supplied content.

    Args:
        path: HTML destination.
        title: Report title.
        summary: Analysis summary.
        tables: Named result tables; first 500 rows form searchable previews.
        images: PNG files embedded directly in the report.
        notes: Interpretation and method notes.
        links: Named companion files, using relative paths or HTTPS URLs.

    Raises:
        ValueError: A supplied link uses an unsafe URL scheme.
    """
    escape = html.escape
    motif_summary = summary.get("motifs", summary)
    ai_summary = summary.get("ai") or summary
    shap_summary = ai_summary.get("shap") or {}
    scan_summary = summary.get("scan") or {}
    regional_summary = summary.get("multiscale") or {}
    metrics = [
        (
            "Positive regions",
            motif_summary.get(
                "positive_sequences",
                summary.get(
                    "positive_regions", summary.get("retained_regions", "—")
                ),
            ),
        ),
        (
            "Control regions",
            motif_summary.get(
                "negative_sequences", summary.get("negative_regions", "—")
            ),
        ),
        ("Motifs at q ≤ 0.05", motif_summary.get("significant_q_0_05", "—")),
        (
            "Held-out ROC AUC",
            ai_summary.get("roc_auc", summary.get("roc_auc", "—")),
        ),
        (
            "Held-out average precision",
            ai_summary.get(
                "average_precision", summary.get("average_precision", "—")
            ),
        ),
        (
            "ML status",
            ai_summary.get("status", summary.get("status", "not applicable")),
        ),
        (
            "Composition baseline ROC AUC",
            ai_summary.get(
                "baseline_roc_auc", summary.get("baseline_roc_auc", "—")
            ),
        ),
        (
            "Permutation p-value",
            ai_summary.get(
                "permutation_p_value", summary.get("permutation_p_value", "—")
            ),
        ),
    ]
    if shap_summary:
        metrics.extend(
            [
                ("SHAP status", shap_summary.get("status", "—")),
                (
                    "Sequences explained",
                    shap_summary.get("explained_sequences", "—"),
                ),
            ]
        )
    if regional_summary:
        if set(summary) == {"multiscale"}:
            metrics = []
        metrics.extend(
            [
                ("Multi-scale search", regional_summary.get("status", "—")),
                (
                    "Search lengths (bp)",
                    ", ".join(
                        str(n)
                        for n in regional_summary.get("window_lengths", [])
                    ),
                ),
                (
                    "Regulatory search windows",
                    regional_summary.get("total_windows", "—"),
                ),
                (
                    "Candidate window unions",
                    regional_summary.get("candidate_regions", "—"),
                ),
            ]
        )
    if scan_summary:
        if set(summary) == {"scan"}:
            metrics = []
        metrics.extend(
            [
                ("Genome scan status", scan_summary.get("status", "—")),
                (
                    "Screened motif consensuses",
                    scan_summary.get("selected_motifs", "—"),
                ),
                ("Genome motif sites", scan_summary.get("total_sites", "—")),
                (
                    "Genic overlap sites",
                    scan_summary.get("genic_overlap_sites", "—"),
                ),
                ("Genome bases", scan_summary.get("genome_bases", "—")),
                (
                    "Allowed substitutions",
                    scan_summary.get("max_mismatches", "—"),
                ),
            ]
        )
    body = [
        "<header><p class='eyebrow'>INTERGENIC REGIONS · RESULTS</p>",
        f"<h1>{escape(title)}</h1>",
        "<p class='subtitle'>Strand-aware sequences, regulatory patterns "
        "and evidence you can inspect.</p></header>",
        "<section class='metrics' aria-label='Key results'>",
    ]
    for label, value in metrics:
        if value is None:
            value = "not estimated"
        displayed = f"{value:.3f}" if isinstance(value, float) else str(value)
        body.append(
            f"<article class='metric'><span>{escape(label)}</span>"
            f"<strong>{escape(displayed)}</strong></article>"
        )
    body.append("</section><section id='interpretation' class='callout'>")
    body.append("<h2>What these results support</h2>")
    for note in notes:
        body.append(f"<p>{escape(note)}</p>")
    if ai_summary.get("reason"):
        body.append(
            f"<p><strong>ML: {escape(ai_summary['reason'])}</strong></p>"
        )
    body.append("</section>")
    if links:
        body.append("<nav class='downloads' aria-label='Result files'>")
        for label, href in links.items():
            parsed = urlparse(href)
            if parsed.scheme not in {"", "https"} or (
                parsed.netloc and parsed.scheme != "https"
            ):
                raise ValueError(
                    "Report links require relative paths or HTTPS"
                )
            body.append(
                f"<a href='{escape(href, quote=True)}'>{escape(label)}</a>"
            )
        body.append("</nav>")
    body.append("<section id='figures'>")
    previous_category = ""
    for image in images:
        category = (
            "SHAP · held-out model explanations"
            if image.stem.startswith("shap_")
            else "Multi-scale regions · exploratory regulatory hypotheses"
            if image.stem.startswith("region_")
            else "Genome scan · unvalidated sequence matches"
            if image.stem.startswith("genome_")
            else "Model performance and sequence signatures"
            if image.stem.startswith("ai_")
            else "Motifs and sequence balance"
        )
        if category != previous_category:
            if previous_category:
                body.append("</div>")
            body.append(f"<h2>{escape(category)}</h2><div class='gallery'>")
            previous_category = category
        encoded = base64.b64encode(image.read_bytes()).decode("ascii")
        pdf = image.with_suffix(".pdf")
        figure_links = (
            "<div class='figure-links'>"
            f"<a href='data:image/png;base64,{encoded}' "
            f"download='{escape(image.name, quote=True)}'>Download PNG</a>"
        )
        if pdf.is_file():
            # Embed both formats so a copied HTML remains entirely portable.
            encoded_pdf = base64.b64encode(pdf.read_bytes()).decode("ascii")
            figure_links += (
                f"<a href='data:application/pdf;base64,{encoded_pdf}' "
                f"download='{escape(pdf.name, quote=True)}'>Download PDF</a>"
            )
        figure_links += "</div>"
        body.append(
            f'<figure><img alt="{escape(image.stem)}" '
            f'src="data:image/png;base64,{encoded}"><figcaption>'
            f"{escape(image.stem.replace('_', ' '))}"
            f"</figcaption>{figure_links}</figure>"
        )
    if previous_category:
        body.append("</div>")
    body.append("</section>")
    for table_number, (name, rows) in enumerate((tables or {}).items()):
        body.append(
            f"<section id='table-{table_number}' class='results'>"
            f"<h2>{escape(name)}</h2>"
        )
        if not rows:
            body.append("<p>No records.</p></section>")
            continue
        fields = list(rows[0])
        preview = rows[:500]
        body.append(
            "<label class='search'>Filter this table "
            f"<input type='search' data-table='results-{table_number}' "
            "placeholder='Gene, motif, evidence or score…'></label>"
            f"<p class='table-status'>{len(preview)} preview rows of "
            f"{len(rows)} records. Select a column heading to sort.</p>"
        )
        body.append(
            "<div class='table'>"
            f"<table id='results-{table_number}'><thead><tr>"
            + "".join(
                f"<th scope='col'><button type='button' data-sort='{i}'>"
                f"{escape(field.replace('_', ' '))}</button></th>"
                for i, field in enumerate(fields)
            )
            + "</tr></thead><tbody>"
        )
        for row in preview:
            cells: list[str] = []
            for field in fields:
                value = row.get(field)
                text = (
                    "NA"
                    if value is None
                    else f"{value:.5g}"
                    if isinstance(value, float)
                    else str(value)
                )
                cells.append(f"<td>{escape(text)}</td>")
            body.append("<tr>" + "".join(cells) + "</tr>")
        body.append(
            "</tbody></table></div><p class='muted'>Full results are "
            "available "
            "in the accompanying TSV files.</p></section>"
        )
    body.append(
        "<details id='provenance'><summary>Methods, metrics "
        "and run details</summary><pre>"
        + escape(json.dumps(obj=dict(summary), indent=2, allow_nan=False))
        + "</pre></details>"
    )
    assets = Path(__file__).parent / "assets"
    css = (assets / "report.css").read_text(encoding="utf-8")
    script = (assets / "report.js").read_text(encoding="utf-8")
    document = (
        "<!doctype html><html lang='en-GB'><meta charset='utf-8'>"
        "<meta name='viewport' "
        "content='width=device-width,initial-scale=1'><title>"
        + escape(title)
        + "</title><style>"
        + css
        + "</style><body><main>"
        + "\n".join(body)
        + "</main><footer>Intergenic regions · Reproducible sequence analysis "
        "· Supporting evidence remains explicit</footer>"
        + "<script type='application/javascript'>"
        + script
        + "</script></body></html>"
    )
    path.write_text(data=document, encoding="utf-8")
