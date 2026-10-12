from dataclasses import replace

import pytest

from intergenic_regions.rna_expression import (
    differential_expression,
    read_count_matrix,
    read_de_results,
    select_de_gene_sets,
)
from intergenic_regions.rna_inputs import RnaSample, read_rna_samples


def metadata(*, replicates=3):
    return [
        RnaSample(
            sample_id=f"{condition}_{i}",
            condition=condition,
            replicate_id=str(i),
            strandedness="unstranded",
        )
        for condition in ("ref", "treated")
        for i in range(1, replicates + 1)
    ]


@pytest.mark.parametrize(
    "body",
    [
        "gene_id\twrong\ng\t1\n",
        "gene_id\tref_1\ttreated_1\ng\t1.2\t3\n",
        "gene_id\tref_1\ttreated_1\ng\t-1\t3\n",
        "gene_id\tref_1\ttreated_1\ng\tNaN\t3\n",
        "gene_id\tref_1\ttreated_1\ng\t9007199254740992\t3\n",
        "gene_id\tref_1\ttreated_1\ng\t1\t3\ng\t1\t3\n",
        "gene_id\tref_1\ttreated_1\ng\t1\n",
        "gene_id\tref_1\ttreated_1\n",
    ],
)
def test_count_matrix_rejects_normalised_or_invalid_data(tmp_path, body):
    path = tmp_path / "c.tsv"
    path.write_text(body)
    with pytest.raises(ValueError):
        read_count_matrix(path=path, samples=metadata(replicates=1))


def test_count_matrix_exact_samples_and_order(tmp_path):
    path = tmp_path / "c.tsv"
    path.write_text("gene_id\ttreated_1\tref_1\ng2\t0\t12\ng1\t8\t7\n")
    assert read_count_matrix(path=path, samples=metadata(replicates=1)) == {
        "g2": {"ref_1": 12, "treated_1": 0},
        "g1": {"ref_1": 7, "treated_1": 8},
    }
    with pytest.raises(ValueError):
        read_count_matrix(path=path, samples=[])


def test_imported_de_missing_tests_never_become_controls(tmp_path):
    path = tmp_path / "de.tsv"
    path.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\n"
        "up\t100\t2\t0.001\ndown\t100\t-2\t0.001\n"
        "control\t100\t0.1\t0.8\nna\t100\t0\tNA\n"
        "low\t1\t0.1\t0.8\nuncertain\t100\t0.4\t0.7\n"
    )
    rows = read_de_results(path=path)
    assert rows[3]["padj"] is None
    for direction, expected in [
        ("up", ["up"]),
        ("down", ["down"]),
        ("both", ["up", "down"]),
    ]:
        positive, negative, audit = select_de_gene_sets(
            rows=rows, direction=direction
        )
        assert positive == expected
        assert negative == ["control"]
        assert len(audit) == 6
        assert audit[3]["reason"] == "untested_or_missing_result"
        assert audit[4]["reason"] == "below_expression_floor"
        assert "screening_control" in audit[2]["reason"]
    assert select_de_gene_sets(rows=[]) == ([], [], [])


@pytest.mark.parametrize(
    "body",
    [
        "gene_id\tpadj\ng\t0.1\n",
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\n",
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\ng\t-1\t1\t0.1\n",
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\ng\t10\tinf\t0.1\n",
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\ng\t10\t1\t1.1\n",
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\ng\t10\t1\t0.1\n"
        "g\t10\t1\t0.2\n",
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\ng\t10\t1\n",
    ],
)
def test_de_import_validation(tmp_path, body):
    path = tmp_path / "de.tsv"
    path.write_text(body)
    with pytest.raises(ValueError):
        read_de_results(path=path)


@pytest.mark.parametrize(
    "settings",
    [
        {"direction": "guess"},
        {"alpha": 0},
        {"alpha": float("nan")},
        {"control_min_padj": 0.01},
        {"min_abs_log2_fold_change": 0},
        {"control_max_abs_log2_fold_change": 2},
        {"min_base_mean": -1},
    ],
)
def test_de_selection_invalid_settings(settings):
    with pytest.raises(ValueError):
        select_de_gene_sets(rows=[], **settings)


@pytest.mark.parametrize(
    "rows",
    [
        [{"gene_id": ""}],
        [{"gene_id": "g"}, {"gene_id": "g"}],
        [{"gene_id": "g", "padj": float("inf")}],
        [{"gene_id": "g", "padj": -1}],
        [{"gene_id": "g", "baseMean": -1}],
    ],
)
def test_de_selection_invalid_rows(rows):
    with pytest.raises(ValueError):
        select_de_gene_sets(rows=rows)


def test_actual_pydeseq2_recovers_planted_gene_effects(rna_dataset):
    samples = read_rna_samples(path=rna_dataset / "samples.tsv")
    counts = read_count_matrix(
        path=rna_dataset / "gene_counts.tsv", samples=samples
    )
    counts["filtered"] = {s.sample_id: 0 for s in samples}
    rows, summary = differential_expression(
        counts=counts,
        samples=samples,
        condition="stimulated",
        reference="reference",
    )
    positive, _, _ = select_de_gene_sets(rows=rows)
    truth = (rna_dataset / "positive_genes.txt").read_text().splitlines()
    assert set(positive) == set(truth)
    assert summary["method"] == "PyDESeq2"
    assert set(summary["size_factors"]) == set(counts["filtered"])
    assert all(value > 0 for value in summary["size_factors"].values())
    assert rows[-1]["padj"] is None
    assert summary["modelled_genes"] == 48
    assert summary["input_genes"] == 49


@pytest.mark.parametrize(
    "settings",
    [
        {"condition": "ref"},
        {"alpha": 1},
        {"threads": 0},
        {"threads": True},
        {"min_total_count": 0},
        {"covariates": ["condition"]},
        {"covariates": ["x", "x"]},
        {"covariates": ["a+b"]},
        {"min_total_count": 99999},
    ],
)
def test_de_model_settings_validation(settings):
    samples = metadata()
    counts = {f"g{i}": {s.sample_id: 100 for s in samples} for i in range(25)}
    with pytest.raises(ValueError):
        differential_expression(
            **{
                "counts": counts,
                "samples": samples,
                "condition": "treated",
                "reference": "ref",
                **settings,
            }
        )


def test_de_requires_biological_replication_raw_counts_and_estimable_design():
    samples = metadata()
    counts = {f"g{i}": {s.sample_id: 100 for s in samples} for i in range(25)}
    args = dict(
        counts=counts, samples=samples, condition="treated", reference="ref"
    )
    for changes in (
        {"samples": metadata(replicates=1)},
        {"samples": [*samples, samples[0]]},
        {"counts": {"g": {"ref_1": 1}}},
        {"counts": {"g": {s.sample_id: 1.5 for s in samples}}},
        {"covariates": ["missing"]},
        {
            "samples": [
                replace(s, covariates={"batch": s.condition}) for s in samples
            ],
            "covariates": ["batch"],
        },
        {
            "samples": [
                replace(s, covariates={"subject": s.sample_id})
                for s in samples
            ],
            "covariates": ["subject"],
        },
    ):
        with pytest.raises(ValueError):
            differential_expression(**{**args, **changes})
