"""Replicate-aware count modelling and explicit, non-circular DE gene sets."""

import csv
import logging
import math
import re
import warnings as warning_tools
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from intergenic_regions.io import open_text
from intergenic_regions.rna_inputs import RnaSample

LOGGER = logging.getLogger(__name__)
DE_FIELDS = (
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)


def read_count_matrix(
    *, path: Path, samples: Sequence[RnaSample]
) -> dict[str, dict[str, int]]:
    """Read raw integer gene counts, never TPM/FPKM or scaled coverage.

    Args:
        path: TSV with gene_id and one column per exact manifest sample_id.
        samples: Metadata specifying all expected sample columns.

    Returns:
        Exact gene/sample count mapping, preserving gene order.

    Raises:
        ValueError: IDs, columns or non-negative integer counts are invalid.
    """
    expected = [sample.sample_id for sample in samples]
    if not expected or len(set(expected)) != len(expected):
        raise ValueError("Count matrix needs unique sample metadata")
    counts: dict[str, dict[str, int]] = {}
    with open_text(path=path) as stream:
        reader = csv.DictReader(f=stream, delimiter="\t")
        headers = reader.fieldnames or []
        if set(headers) != {"gene_id", *expected} or len(headers) != len(
            set(headers)
        ):
            raise ValueError(
                "Count TSV columns must match sample_id metadata exactly"
            )
        for number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{number}: malformed count row")
            gene = row["gene_id"]
            if not gene or any(c.isspace() for c in gene) or gene in counts:
                raise ValueError(
                    "Count matrix gene IDs are invalid or duplicated"
                )
            values: dict[str, int] = {}
            for sample in expected:
                text = row[sample].strip()
                if not re.fullmatch(r"\d+", text):
                    raise ValueError(
                        "DE requires raw non-negative integer counts"
                    )
                value = int(text)
                if value > 2**53 - 1:
                    raise ValueError(
                        "Gene count exceeds exact numeric precision"
                    )
                values[sample] = value
            counts[gene] = values
    if not counts:
        raise ValueError("Gene count matrix is empty")
    return counts


def read_de_results(*, path: Path) -> list[dict[str, Any]]:
    """Import a standard DESeq2-style TSV with explicit missing tests.

    Args:
        path: gene_id, baseMean, log2FoldChange and padj columns, with optional
            lfcSE/stat/pvalue. Export imported adjusted p-values unchanged.

    Returns:
        Canonical result rows; NA tests stay None and cannot select controls.

    Raises:
        ValueError: Required fields, IDs or finite numeric values are invalid.
    """
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    with open_text(path=path) as stream:
        reader = csv.DictReader(f=stream, delimiter="\t")
        headers = reader.fieldnames or []
        if not {"gene_id", "baseMean", "log2FoldChange", "padj"} <= set(
            headers
        ) or len(headers) != len(set(headers)):
            raise ValueError(
                "DE TSV requires gene_id/baseMean/log2FoldChange/padj"
            )
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed DE TSV row")
            gene = row["gene_id"].strip()
            if not gene or any(c.isspace() for c in gene) or gene in seen:
                raise ValueError(
                    "DE gene identifiers are invalid or duplicated"
                )
            seen.add(gene)
            parsed: dict[str, Any] = {"gene_id": gene}
            for name in DE_FIELDS[1:]:
                text = row.get(name, "").strip()
                value = (
                    None
                    if text.lower() in {"", "na", "nan", "none", "."}
                    else float(text)
                )
                if value is not None and (
                    not math.isfinite(value)
                    or (name in {"baseMean", "lfcSE"} and value < 0)
                    or (name in {"pvalue", "padj"} and not 0 <= value <= 1)
                ):
                    raise ValueError(f"Invalid DE value: {name}")
                parsed[name] = value
            output.append(parsed)
    if not output:
        raise ValueError("DE result table is empty")
    return output


def differential_expression(
    *,
    counts: Mapping[str, Mapping[str, int]],
    samples: Sequence[RnaSample],
    condition: str,
    reference: str,
    covariates: Sequence[str] = (),
    alpha: float = 0.05,
    min_total_count: int = 10,
    threads: int = 1,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run the optional PyDESeq2 negative-binomial workflow on raw gene counts.

    Args:
        counts: Non-negative raw integers for every manifest sample.
        samples: Independent biological observations; conditions may be paired
            using an explicit subject covariate, never technical replication.
        condition: Numerator/target condition.
        reference: Denominator/reference condition.
        covariates: Categorical nuisance columns such as batch/subject.
        alpha: Independent-filtering/FDR threshold used by PyDESeq2.
        min_total_count: Predeclared total-count filter.
        threads: Bounded worker count.

    Returns:
        All input genes, including filtered/untested NA rows, and method/QC
        metadata with estimated library size factors.

    Raises:
        ValueError: Counts, replication, design rank or settings are invalid.
        ImportError: Install the optional de extra for PyDESeq2.
    """
    if (
        condition == reference
        or not math.isfinite(alpha)
        or not 0 < alpha < 1
        or isinstance(min_total_count, bool)
        or min_total_count < 1
        or isinstance(threads, bool)
        or not 1 <= threads <= 256
        or len(covariates) != len(set(covariates))
        or any(
            not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name)
            for name in covariates
        )
        or "condition" in covariates
    ):
        raise ValueError("Invalid DE contrast/design/settings")
    selected = [s for s in samples if s.condition in {condition, reference}]
    ids = [s.sample_id for s in selected]
    if len(set(s.sample_id for s in samples)) != len(samples) or len(
        {(s.condition, s.replicate_id) for s in samples}
    ) != len(samples):
        raise ValueError("DE sample/biological observation identifiers repeat")
    replicates = Counter(s.condition for s in selected)
    if any(replicates[level] < 2 for level in (condition, reference)):
        raise ValueError(
            "DE testing requires at least two biological replicates per c"
            "ondition"
        )
    all_ids = {s.sample_id for s in samples}
    if not counts or any(set(row) != all_ids for row in counts.values()):
        raise ValueError("DE counts must match all manifest samples")
    if any(
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value <= 2**53 - 1
        for row in counts.values()
        for value in row.values()
    ):
        raise ValueError("DE counts must be raw, finite non-negative integers")
    kept = [
        gene
        for gene, row in counts.items()
        if sum(row[sample] for sample in ids) >= min_total_count
    ]
    if len(kept) < 20 or any(
        sum(counts[gene][sample] for gene in kept) == 0 for sample in ids
    ):
        raise ValueError(
            "DE needs at least 20 expressed genes and non-empty sample li"
            "braries"
        )
    for name in covariates:
        values = {s.covariates.get(name, "") for s in selected}
        if "" in values or len(values) < 2:
            raise ValueError(f"DE covariate is missing or constant: {name}")
    import numpy as np
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    frame = pd.DataFrame(
        data=[[counts[gene][sample] for gene in kept] for sample in ids],
        index=ids,
        columns=kept,
        dtype="int64",
    )
    metadata = pd.DataFrame(
        data={
            "condition": pd.Categorical(
                values=[s.condition for s in selected],
                categories=[reference, condition],
            ),
            **{
                name: pd.Categorical(
                    values=[s.covariates[name] for s in selected]
                )
                for name in covariates
            },
        },
        index=ids,
    )
    dummy = pd.get_dummies(data=metadata, drop_first=True, dtype=float)
    design_matrix = np.column_stack((np.ones(len(ids)), dummy.to_numpy()))
    if (
        np.linalg.matrix_rank(design_matrix) != design_matrix.shape[1]
        or len(ids) <= design_matrix.shape[1]
    ):
        raise ValueError(
            "DE design is confounded or has no residual biological replication"
        )
    inference = DefaultInference(n_cpus=threads)
    design = "~ " + " + ".join([*covariates, "condition"])
    LOGGER.info(
        "Fitting PyDESeq2: %d genes, %d samples, design=%s",
        len(kept),
        len(ids),
        design,
    )
    dataset = DeseqDataSet(
        counts=frame,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    with warning_tools.catch_warnings(record=True) as model_warnings:
        warning_tools.simplefilter(action="always")
        dataset.deseq2()
        statistics = DeseqStats(
            dds=dataset,
            contrast=["condition", condition, reference],
            alpha=alpha,
            inference=inference,
            quiet=True,
        )
        statistics.summary()
    results = statistics.results_df
    output: list[dict[str, Any]] = []
    for gene in counts:
        row: dict[str, Any] = {"gene_id": gene}
        for name in DE_FIELDS[1:]:
            value = (
                float(results.loc[gene, name])
                if gene in results.index
                else float("nan")
            )
            row[name] = value if math.isfinite(value) else None
        output.append(row)
    warnings = list(
        dict.fromkeys(str(item.message) for item in model_warnings)
    )
    if min(replicates.values()) < 3:
        warnings.append(
            "Only two biological replicates in a condition; "
            "dispersion and power are limited."
        )
    tested = sum(row["padj"] is not None for row in output)
    if not tested:
        warnings.append(
            "No adjusted DE tests are estimable; inspect count variation, "
            "dispersion fitting and independent-filtering diagnostics."
        )
    for warning in warnings:
        LOGGER.warning("DE model: %s", warning)
    factors = dict(zip(ids, dataset.obs["size_factors"], strict=True))
    if any(
        not math.isfinite(value) or value <= 0 for value in factors.values()
    ):
        raise ValueError(
            "DE model could not estimate valid library size factors"
        )
    return output, {
        "status": "completed" if tested else "not_estimable",
        "tested_genes": tested,
        "method": "PyDESeq2",
        "design": design,
        "condition": condition,
        "reference": reference,
        "biological_replicates": dict(replicates),
        "input_genes": len(counts),
        "modelled_genes": len(kept),
        "min_total_count": min_total_count,
        "alpha": alpha,
        "size_factors": {
            sample: float(value) for sample, value in factors.items()
        },
        "warnings": warnings,
        "test": (
            "Wald test against zero log2 fold change; effect-size filtering i"
            "s separate"
        ),
    }


def select_de_gene_sets(
    *,
    rows: Sequence[Mapping[str, Any]],
    direction: str = "up",
    alpha: float = 0.05,
    min_abs_log2_fold_change: float = 1.0,
    min_base_mean: float = 10.0,
    control_max_abs_log2_fold_change: float = 0.25,
    control_min_padj: float = 0.5,
) -> tuple[list[str], list[str], list[dict[str, Any]]]:
    """Choose DE foreground and expressed screening controls transparently.

    Args:
        rows: Validated DESeq2-style result records.
        direction: up, down or both DE directions for the positive class.
        alpha: Positive adjusted p-value threshold.
        min_abs_log2_fold_change: Minimum positive effect magnitude.
        min_base_mean: Expression floor for both classes.
        control_max_abs_log2_fold_change: Small observed control-effect limit.
        control_min_padj: Minimum tested control adjusted p-value.

    Returns:
        Exact positive/negative gene IDs and an audit for every DE row.
        Controls are screening comparators, not proven equivalent-expression
        genes. Untested/filtered/low-count rows never become controls.

    Raises:
        ValueError: Thresholds, duplicate genes or row values are invalid.
    """
    values = (
        alpha,
        min_abs_log2_fold_change,
        min_base_mean,
        control_max_abs_log2_fold_change,
        control_min_padj,
    )
    if (
        direction not in {"up", "down", "both"}
        or any(not math.isfinite(value) for value in values)
        or not 0 < alpha < control_min_padj <= 1
        or not 0 <= control_max_abs_log2_fold_change < min_abs_log2_fold_change
        or min_base_mean < 0
    ):
        raise ValueError("Invalid DE foreground/control thresholds")
    positive: list[str] = []
    negative: list[str] = []
    audit: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        gene = row["gene_id"]
        if gene in seen or not gene:
            raise ValueError("DE selection repeats an invalid gene ID")
        seen.add(gene)
        mean, effect, adjusted = (
            row.get(name) for name in ("baseMean", "log2FoldChange", "padj")
        )
        if any(
            value is not None and not math.isfinite(value)
            for value in (mean, effect, adjusted)
        ):
            raise ValueError("DE selection received non-finite values")
        if (
            adjusted is not None
            and not 0 <= adjusted <= 1
            or mean is not None
            and mean < 0
        ):
            raise ValueError("Invalid DE adjusted p-value/expression mean")
        label, reason = "excluded", "untested_or_missing_result"
        if mean is not None and effect is not None and adjusted is not None:
            if mean < min_base_mean:
                reason = "below_expression_floor"
            elif (
                adjusted <= alpha
                and abs(effect) >= min_abs_log2_fold_change
                and (
                    direction == "both"
                    or (direction == "up" and effect > 0)
                    or (direction == "down" and effect < 0)
                )
            ):
                label, reason = (
                    "positive",
                    "differential_expression_foreground",
                )
                positive.append(gene)
            elif (
                adjusted >= control_min_padj
                and abs(effect) <= control_max_abs_log2_fold_change
            ):
                label, reason = (
                    "negative",
                    "expressed_small_observed_effect_screening_control",
                )
                negative.append(gene)
            else:
                reason = "outside_predeclared_gene_set_thresholds"
        audit.append(
            {
                "gene_id": gene,
                "label": label,
                "reason": reason,
                "baseMean": mean,
                "log2FoldChange": effect,
                "padj": adjusted,
            }
        )
    LOGGER.info(
        "DE gene sets: %d positives, %d screening controls",
        len(positive),
        len(negative),
    )
    return positive, negative, audit
