"""Reference provenance/build safeguards and shell-free HOMER integration."""

import io
import json
import subprocess
from pathlib import Path

import pytest

from intergenic_regions.homer import homer_command, run_homer
from intergenic_regions.models import Region
from intergenic_regions.references import (
    CATALOGUE,
    canonical_assembly,
    download_bed,
    import_reference,
    query_references,
)


def test_reference_local_import_query_assembly_and_checksum(tmp_path):
    bed = tmp_path / "source.bed"
    bed.write_text("x\t2\t8\telement\n")
    output = tmp_path / "reference"
    metadata = import_reference(
        output=output,
        source="custom",
        local_bed=bed,
        organism="test_species",
        assembly="hg38",
    )
    assert metadata["assembly"] == "GRCh38"
    assert metadata["source_url"] is None
    region = Region(
        gene_id="g",
        contig="x",
        start=0,
        end=5,
        strand="+",
        direction="upstream",
        sequence="AAAAA",
        status="retained",
        stop_reason="contig_boundary",
        available_length=5,
    )
    rows = query_references(
        regions=[region],
        reference_directories=[output],
        organism="test_species",
        assembly="GRCh38",
    )
    assert rows[0]["overlap_bp"] == 3
    assert rows[0]["overlap_fraction"] == 0.6
    for organism, assembly in [
        ("wrong", "GRCh38"),
        ("test_species", "GRCh37"),
    ]:
        with pytest.raises(ValueError):
            query_references(
                regions=[region],
                reference_directories=[output],
                organism=organism,
                assembly=assembly,
            )
    (output / "regions.bed").write_text("x\t1\t3\n")
    with pytest.raises(ValueError, match="checksum"):
        query_references(
            regions=[region],
            reference_directories=[output],
            organism="test_species",
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"source": "unknown"},
        {"source": "custom"},
        {"source": "custom", "organism": "two words", "assembly": "build"},
        {"source": "custom", "organism": "test", "assembly": "build"},
        {"source": "screen-human-enhancers", "assembly": "GRCh37"},
        {"source": "screen-human-enhancers", "organism": "mouse"},
        {
            "source": "custom",
            "url": "https://example.com/data",
            "local_bed": Path("bed"),
        },
    ],
)
def test_reference_invalid_metadata(tmp_path, kwargs):
    with pytest.raises(ValueError):
        import_reference(output=tmp_path / "out", **kwargs)


def test_reference_checksums_and_empty_bed(tmp_path):
    path = tmp_path / "in.bed"
    path.write_text("x\t0\t2\n")
    with pytest.raises(ValueError, match="SHA"):
        import_reference(
            output=tmp_path / "out",
            source="custom",
            local_bed=path,
            organism="test",
            assembly="build",
            expected_sha256="bad",
        )
    assert not (tmp_path / "out").exists()
    path.write_text("")
    with pytest.raises(ValueError, match="no intervals"):
        import_reference(
            output=tmp_path / "out",
            source="custom",
            local_bed=path,
            organism="test",
            assembly="build",
        )
    assert canonical_assembly(assembly="TAIR10") == "TAIR10"
    assert canonical_assembly(assembly="mm10") == "GRCm38"
    for build in ("", "bad build"):
        with pytest.raises(ValueError):
            canonical_assembly(assembly=build)


class FakeResponse(io.BytesIO):
    def geturl(self):
        return "https://example.org/data.bed"


def test_reference_download_limits_and_atomic_failures(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "intergenic_regions.references.urlopen",
        lambda **kwargs: FakeResponse(b"x\t0\t2\n"),
    )
    destination = tmp_path / "download.bed"
    download_bed(url="https://example.org/data.bed", path=destination)
    assert destination.read_text() == "x\t0\t2\n"
    with pytest.raises(ValueError):
        download_bed(
            url="https://example.org/data.bed", path=destination, max_bytes=2
        )
    assert not destination.exists()
    for url in (
        "http://example.org/data",
        "file:///tmp/x",
        "https://u:pass@example.org/data",
    ):
        with pytest.raises(ValueError):
            download_bed(url=url, path=destination)
    monkeypatch.setattr(
        "intergenic_regions.references.urlopen",
        lambda **kwargs: FakeResponse(b""),
    )
    with pytest.raises(ValueError, match="empty"):
        download_bed(url="https://example.org/data", path=destination)
    monkeypatch.setattr(
        FakeResponse, "geturl", lambda self: "http://example.org/data"
    )
    with pytest.raises(ValueError, match="redirected"):
        download_bed(url="https://example.org/data", path=destination)


def test_preset_download_import_uses_verified_url(tmp_path, monkeypatch):
    def fake_download(*, url, path):
        assert url == CATALOGUE["screen-human-enhancers"]["url"]
        path.write_text("chr1\t0\t2\n")

    monkeypatch.setattr(
        "intergenic_regions.references.download_bed", fake_download
    )
    metadata = import_reference(
        output=tmp_path / "out", source="screen-human-enhancers"
    )
    assert metadata["organism"] == "Homo_sapiens"


def test_homer_argv_dry_run_and_missing_executable(
    dataset, tmp_path, monkeypatch
):
    command = homer_command(
        positive=dataset["positive"],
        negative=dataset["negative"],
        output=tmp_path / "a space;literal",
        known_motifs=dataset["motifs"],
        known_only=True,
        threads=2,
    )
    assert command[2] == "fasta" and "-fasta" in command
    assert "-mknown" in command and "-nomotif" in command
    assert "a space;literal" in command[3]
    for kwargs in (
        {"threads": 0},
        {"lengths": ()},
        {"lengths": (1,)},
        {"executable": ""},
    ):
        with pytest.raises(ValueError):
            homer_command(
                positive=dataset["positive"],
                negative=dataset["negative"],
                output=tmp_path / "out",
                **kwargs,
            )
    monkeypatch.setattr(
        "intergenic_regions.homer.shutil.which", lambda **kwargs: None
    )
    summary = run_homer(
        positive=dataset["positive"],
        negative=dataset["negative"],
        output=tmp_path / "dry",
        dry_run=True,
    )
    assert summary["status"] == "planned"
    assert (
        json.loads((tmp_path / "dry" / "command.json").read_text())["status"]
        == "planned"
    )
    with pytest.raises(FileNotFoundError):
        run_homer(
            positive=dataset["positive"],
            negative=dataset["negative"],
            output=tmp_path / "missing",
        )
    with pytest.raises(ValueError):
        run_homer(
            positive=dataset["positive"],
            negative=dataset["negative"],
            output=tmp_path / "bad",
            timeout=0,
        )
    with pytest.raises(FileNotFoundError):
        run_homer(
            positive=dataset["positive"],
            negative=dataset["negative"],
            output=tmp_path / "bad",
            known_motifs=tmp_path / "missing",
        )


def test_homer_success_failure_timeout_cleanup(dataset, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "intergenic_regions.homer.shutil.which",
        lambda **kwargs: "/usr/bin/mock-homer",
    )

    def completed(**kwargs):
        assert kwargs["shell"] is False
        Path(kwargs["args"][3]).mkdir()
        return subprocess.CompletedProcess(args=kwargs["args"], returncode=0)

    monkeypatch.setattr("intergenic_regions.homer.subprocess.run", completed)
    assert (
        run_homer(
            positive=dataset["positive"],
            negative=dataset["negative"],
            output=tmp_path / "complete",
        )["status"]
        == "completed"
    )
    monkeypatch.setattr(
        "intergenic_regions.homer.subprocess.run",
        lambda **kwargs: subprocess.CompletedProcess(args=[], returncode=1),
    )
    with pytest.raises(RuntimeError, match="exited"):
        run_homer(
            positive=dataset["positive"],
            negative=dataset["negative"],
            output=tmp_path / "failed",
        )
    assert not (tmp_path / "failed").exists()
    monkeypatch.setattr(
        "intergenic_regions.homer.subprocess.run",
        lambda **kwargs: subprocess.CompletedProcess(args=[], returncode=0),
    )
    with pytest.raises(RuntimeError, match="without producing"):
        run_homer(
            positive=dataset["positive"],
            negative=dataset["negative"],
            output=tmp_path / "empty",
        )

    def timeout(**kwargs):
        raise subprocess.TimeoutExpired(cmd=kwargs["args"], timeout=1)

    monkeypatch.setattr("intergenic_regions.homer.subprocess.run", timeout)
    with pytest.raises(RuntimeError, match="timeout"):
        run_homer(
            positive=dataset["positive"],
            negative=dataset["negative"],
            output=tmp_path / "timeout",
        )
