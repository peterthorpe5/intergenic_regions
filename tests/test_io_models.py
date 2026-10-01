"""Unit tests for immutable models and validated input/output functions."""

import gzip
import hashlib
import json

import pytest

from intergenic_regions import __version__
from intergenic_regions.io import (
    file_fingerprint,
    open_text,
    output_bundle,
    provenance,
    read_fasta,
    read_identifiers,
    write_fasta,
    write_json,
    write_tsv,
)
from intergenic_regions.models import Gene, Region


@pytest.mark.parametrize(
    "changes",
    [
        {"gene_id": ""},
        {"gene_id": "bad id"},
        {"contig": ""},
        {"contig": "bad id"},
        {"start": -1},
        {"end": 1},
        {"strand": "?"},
    ],
)
def test_gene_validation(changes):
    options = dict(gene_id="g.1", contig="chr1", start=1, end=5, strand="+")
    options.update(changes)
    with pytest.raises(ValueError):
        Gene(**options)


def test_models_are_frozen_and_exact():
    from dataclasses import FrozenInstanceError

    gene = Gene(gene_id="gene.1", contig="x", start=0, end=3, strand="-")
    with pytest.raises(FrozenInstanceError):
        gene.start = 5
    region = Region(
        gene_id="g|x",
        contig="x",
        start=0,
        end=1,
        strand="+",
        direction="upstream",
        sequence="A",
        status="retained",
        stop_reason="gene_boundary",
        available_length=1,
    )
    assert region.sequence_id == "g%7Cx|upstream"


def test_text_and_identifier_readers(tmp_path):
    path = tmp_path / "ids.txt.gz"
    with gzip.open(filename=path, mode="wt") as stream:
        stream.write("# comment\ngene.1\n\ngene.2\n")
    assert read_identifiers(path=path) == ["gene.1", "gene.2"]
    with open_text(path=path) as stream:
        assert stream.readline().startswith("#")


@pytest.mark.parametrize("text", ["", "#comment\n", "g\ng\n", "g other\n"])
def test_identifier_failures(tmp_path, text):
    path = tmp_path / "ids"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_identifiers(path=path)


def test_fasta_roundtrip_and_soft_mask(tmp_path):
    path = tmp_path / "data.fa"
    write_fasta(path=path, records={"g.1": "ACGtN" * 30, "g2": "RY"})
    assert read_fasta(path=path)["g.1"] == "ACGTN" * 30
    assert read_fasta(path=path, mask_lowercase=True)["g.1"] == "ACGNN" * 30
    assert len(path.read_text().splitlines()[1]) == 80
    compressed = tmp_path / "data.fa.gz"
    with gzip.open(filename=compressed, mode="wb") as stream:
        stream.write(path.read_bytes())
    assert read_fasta(path=compressed) == read_fasta(path=path)


@pytest.mark.parametrize(
    "text",
    [
        "",
        "AAA\n",
        ">\nAAA",
        ">g\n",
        ">g\nAAA\n>g\nCCC",
        ">g\nAAA\n>x\nCCC\n>g\nTTT",
        ">g\nAZZ",
        ">g\nAA AA",
    ],
)
def test_fasta_failures(tmp_path, text):
    path = tmp_path / "bad.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_fasta(path=path)


def test_tables_json_hash_and_provenance(tmp_path, monkeypatch):
    path = tmp_path / "out.tsv"
    write_tsv(path=path, rows=[{"a": "x", "b": 2}], fields=["a", "b"])
    assert path.read_text() == "a\tb\nx\t2\n"
    fingerprint = file_fingerprint(path=path)
    assert (
        fingerprint["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    )
    manifest = provenance(inputs=[path], settings={"length": 10})
    assert manifest["package_version"] == __version__
    assert manifest["inputs"][0]["bytes"] == path.stat().st_size
    write_json(path=tmp_path / "out.json", data=manifest)
    assert json.loads((tmp_path / "out.json").read_text())["settings"] == {
        "length": 10
    }
    with pytest.raises(ValueError):
        write_json(path=tmp_path / "bad.json", data={"bad": float("nan")})


def test_atomic_bundle_success_and_failure(tmp_path):
    output = tmp_path / "results"
    with output_bundle(path=output) as stage:
        (stage / "x").write_text("success")
        assert not output.exists()
    assert (output / "x").read_text() == "success"
    with pytest.raises(FileExistsError), output_bundle(path=output):
        pass
    failed = tmp_path / "failed"
    with pytest.raises(RuntimeError), output_bundle(path=failed) as stage:
        (stage / "x").write_text("partial")
        raise RuntimeError("failed")
    assert not failed.exists()
    with pytest.raises(FileExistsError), output_bundle(path=failed):
        failed.mkdir()
    assert list(failed.iterdir()) == []
