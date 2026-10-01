"""Official SHAP package graphics for audited out-of-fold explanations."""

import warnings
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from intergenic_regions.reporting import save_figure


def shap_explanation(*, rows: Sequence[Mapping[str, Any]]) -> Any:
    """Construct a SHAP Explanation from complete, additive long-table rows.

    Args:
        rows: One row per displayed feature and sequence, including remainder.

    Returns:
        A SHAP package Explanation with fold-specific base values.

    Raises:
        ValueError: Features duplicate, samples are incomplete or sums differ.
        ImportError: The SHAP package is not installed.
    """
    import shap

    if not rows:
        raise ValueError("No SHAP rows to plot")
    grouped: dict[str, dict[str, Mapping[str, Any]]] = defaultdict(dict)
    for row in rows:
        identifier, feature = row["sequence_id"], row["feature"]
        if feature in grouped[identifier]:
            raise ValueError("Duplicated SHAP feature/sequence row")
        grouped[identifier][feature] = row
    identifiers = sorted(grouped)
    names = sorted(
        grouped[identifiers[0]],
        key=lambda name: (name == "other_features", name),
    )
    if any(set(records) != set(names) for records in grouped.values()):
        raise ValueError("Incomplete SHAP feature rows")
    values = np.asarray(
        [[grouped[i][n]["shap_value"] for n in names] for i in identifiers],
        dtype=float,
    )
    data = np.asarray(
        [
            [
                grouped[i][n]["feature_value"]
                if grouped[i][n]["feature_value"] is not None
                else np.nan
                for n in names
            ]
            for i in identifiers
        ],
        dtype=float,
    )
    baselines = np.asarray(
        [grouped[i][names[0]]["base_value"] for i in identifiers], dtype=float
    )
    logits = np.asarray(
        [grouped[i][names[0]]["model_logit"] for i in identifiers], dtype=float
    )
    if (
        not np.isfinite(values).all()
        or not np.isfinite(baselines).all()
        or np.isinf(data).any()
        or not np.isfinite(logits).all()
        or not np.allclose(values.sum(axis=1) + baselines, logits, atol=1e-8)
    ):
        raise ValueError("SHAP rows do not reconstruct finite model logits")
    displayed = {
        "gc_fraction": "GC fraction (training-scaled)",
        "log_length": "Log length (training-scaled)",
        "ambiguous_fraction": "Ambiguity (training-scaled)",
        "other_features": "Other features (additive remainder)",
    }
    return shap.Explanation(
        values=values,
        base_values=baselines,
        data=data,
        feature_names=[displayed.get(n, n) for n in names],
        instance_names=identifiers,
    )


def plot_shap(
    *,
    rows: Sequence[Mapping[str, Any]],
    importance: Sequence[Mapping[str, Any]],
    directory: Path,
) -> list[Path]:
    """Produce official SHAP global and local explanation graphics.

    Args:
        rows: Audited long explanations; displayed remainder preserves sums.
        importance: Importance across all features, including unplotted words.
        directory: Figure destination; each PNG has a vector PDF companion.

    Returns:
        Embedded-report PNG paths, using SHAP's own plotting functions.

    Raises:
        ValueError: Explanation rows are invalid or importance is empty.
        ImportError: SHAP is unavailable; callers retain the numerical results.
    """
    import matplotlib

    matplotlib.use(backend="Agg")
    import matplotlib.pyplot as plt
    import shap

    if not importance:
        raise ValueError("SHAP importance is empty")
    explanation = shap_explanation(rows=rows)
    paths = []
    global_values = shap.Explanation(
        values=np.asarray([[r["mean_absolute_shap"] for r in importance]]),
        feature_names=[r["feature"] for r in importance],
    )
    shap.plots.bar(shap_values=global_values, max_display=16, show=False)
    plt.gca().set(
        xlabel="Mean |held-out SHAP| (positive-class log-odds)",
        title="Which features influence held-out model scores?",
    )
    paths.append(
        save_figure(
            figure=plt.gcf(), directory=directory, name="shap_global_bar"
        )
    )
    with warnings.catch_warnings():
        # The additive remainder has no feature value and is shown grey.
        warnings.filterwarnings(
            action="ignore",
            message="All-NaN slice encountered",
            category=RuntimeWarning,
        )
        shap.plots.beeswarm(
            shap_values=explanation,
            max_display=16,
            show=False,
            group_remaining_features=False,
        )
    plt.gca().set(
        xlabel="Held-out SHAP contribution (positive-class log-odds)",
        title="Direction and distribution of feature contributions",
    )
    paths.append(
        save_figure(
            figure=plt.gcf(), directory=directory, name="shap_beeswarm"
        )
    )
    logits = explanation.values.sum(axis=1) + explanation.base_values
    heatmap_axis = shap.plots.heatmap(
        shap_values=explanation,
        max_display=16,
        show=False,
        instance_order=np.argsort(logits),
        plot_width=11,
    )
    # SHAP's top trace sums contributions without adding base values.
    ticks = heatmap_axis.get_yticklabels()
    heatmap_axis.set_yticks(
        heatmap_axis.get_yticks(),
        ["Sum of SHAP", *[tick.get_text() for tick in ticks[1:]]],
    )
    heatmap_axis.set(
        xlabel="Sequences ordered by their held-out model log-odds",
        title="Fold-specific contributions and per-sequence baselines",
    )
    paths.append(
        save_figure(figure=plt.gcf(), directory=directory, name="shap_heatmap")
    )
    for name, position in (
        ("positive", int(np.argmax(logits))),
        ("negative", int(np.argmin(logits))),
    ):
        local_explanation = shap.Explanation(
            values=explanation.values[position],
            base_values=explanation.base_values[position],
            data=np.asarray(
                [
                    float(value) if np.isfinite(value) else "NA"
                    for value in explanation.data[position]
                ],
                dtype=object,
            ),
            feature_names=explanation.feature_names,
        )
        shap.plots.waterfall(
            shap_values=local_explanation, max_display=12, show=False
        )
        plt.gca().set_title(
            f"{explanation.instance_names[position]}: explained model log-odds"
        )
        paths.append(
            save_figure(
                figure=plt.gcf(),
                directory=directory,
                name=f"shap_waterfall_{name}",
            )
        )
    available = np.flatnonzero(np.isfinite(explanation.data).any(axis=0))
    if len(available):
        feature = int(
            available[
                np.argmax(
                    np.abs(explanation.values[:, available]).mean(axis=0)
                )
            ]
        )
        shap.plots.scatter(
            shap_values=explanation[:, feature],
            show=False,
            x_jitter=0,
            ylabel="Held-out SHAP contribution (log-odds)",
            title=(
                "Leading-feature dependence; "
                "fold variation is not an interaction test"
            ),
        )
        paths.append(
            save_figure(
                figure=plt.gcf(), directory=directory, name="shap_dependence"
            )
        )
    return paths
