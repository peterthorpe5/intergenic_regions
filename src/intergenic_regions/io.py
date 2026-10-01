"""Validated text input, atomic output bundles and reproducibility records."""

import csv
import gzip
import hashlib
import json
import os
import platform
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, TextIO

from intergenic_regions._version import __version__

DNA = frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn")


def open_text(*, path: Path) -> TextIO:
    """Open UTF-8 text, transparently reading gzip files.

    Args:
        path: An existing readable text file.

    Returns:
        An open stream owned by the caller.
    """
    if path.suffix.lower() == ".gz":
        return gzip.open(filename=path, mode="rt", encoding="utf-8")
    return path.open(mode="r", encoding="utf-8")


def read_identifiers(*, path: Path) -> list[str]:
    """Read one exact identifier per line, ignoring comments and blank lines.

    Args:
        path: Identifier list, optionally gzip compressed.

    Returns:
        Identifiers in input order.

    Raises:
        ValueError: The list is empty, duplicated or has multiple columns.
    """
    identifiers: list[str] = []
    seen: set[str] = set()
    with open_text(path=path) as stream:
        for number, raw in enumerate(stream, start=1):
            item = raw.strip()
            if not item or item.startswith("#"):
                continue
            if len(item.split()) != 1 or item in seen:
                raise ValueError(f"{path}:{number}: invalid or duplicate ID")
            identifiers.append(item)
            seen.add(item)
    if not identifiers:
        raise ValueError(f"Empty identifier list: {path}")
    return identifiers


def read_fasta(*, path: Path, mask_lowercase: bool = False) -> dict[str, str]:
    """Read a sequence collection, checking identifiers and DNA symbols.

    Args:
        path: FASTA, optionally gzip compressed.
        mask_lowercase: Replace soft-masked bases with ``N``.

    Returns:
        An insertion-ordered dictionary of identifiers and sequences.

    Raises:
        ValueError: FASTA is malformed, duplicated, empty or non-DNA.
    """
    records: dict[str, str] = {}
    pieces: list[str] = []
    current: str | None = None
    with open_text(path=path) as stream:
        for number, raw in enumerate(stream, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(pieces)
                fields = line[1:].split()
                if not fields or fields[0] in records or fields[0] == current:
                    raise ValueError(f"{path}:{number}: invalid FASTA ID")
                current = fields[0]
                pieces = []
            else:
                if current is None or not set(line) <= DNA:
                    raise ValueError(f"{path}:{number}: invalid DNA/FASTA")
                if mask_lowercase:
                    line = "".join("N" if c.islower() else c for c in line)
                pieces.append(line.upper())
    if current is not None:
        records[current] = "".join(pieces)
    if not records or any(not seq for seq in records.values()):
        raise ValueError(
            f"FASTA contains no sequences or empty records: {path}"
        )
    return records


def write_fasta(*, path: Path, records: Mapping[str, str]) -> None:
    """Write a wrapped FASTA collection.

    Args:
        path: Output file.
        records: Identifier-to-sequence mapping.

    Raises:
        ValueError: An identifier or DNA sequence is invalid; validation
            occurs before the destination is opened.
    """
    for identifier, sequence in records.items():
        if (
            not isinstance(identifier, str)
            or not identifier
            or any(c.isspace() for c in identifier)
            or not isinstance(sequence, str)
            or not sequence
            or not set(sequence) <= DNA
        ):
            raise ValueError("Invalid FASTA identifier or DNA sequence")
    with path.open(mode="w", encoding="utf-8") as stream:
        for identifier, sequence in records.items():
            stream.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                stream.write(sequence[start : start + 80] + "\n")


def write_tsv(
    *, path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]
) -> None:
    """Write a tab-separated table with a header, including empty tables.

    Args:
        path: Output file.
        rows: Records with keys matching the field names.
        fields: Ordered column names.
    """
    with path.open(mode="w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            f=stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rowdicts=rows)


def write_json(*, path: Path, data: Any) -> None:
    """Write strict JSON, preventing non-standard NaN/Infinity values.

    Args:
        path: Output file.
        data: JSON-serialisable content.
    """
    path.write_text(
        data=json.dumps(obj=data, indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )


@contextmanager
def output_bundle(*, path: Path) -> Iterator[Path]:
    """Stage a new output directory and publish it only after success.

    Args:
        path: New output directory; existing paths are refused.

    Yields:
        A staging directory on the same filesystem.

    Raises:
        FileExistsError: The destination already exists.
    """
    if path.exists():
        raise FileExistsError(f"Output already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{path.name}-", dir=path.parent
    ) as temporary:
        stage = Path(temporary) / "bundle"
        stage.mkdir()
        yield stage
        if path.exists():
            raise FileExistsError(f"Output appeared during execution: {path}")
        os.rename(src=stage, dst=path)


def file_fingerprint(*, path: Path) -> dict[str, Any]:
    """Hash an input without loading it into memory.

    Args:
        path: Input file.

    Returns:
        Absolute path, byte size and SHA-256 digest.
    """
    digest = hashlib.sha256()
    with path.open(mode="rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def provenance(
    *, inputs: Sequence[Path], settings: Mapping[str, Any]
) -> dict[str, Any]:
    """Record input hashes, settings and installed software versions.

    Args:
        inputs: Files used for the analysis.
        settings: Explicit analysis parameters.

    Returns:
        A strict JSON-compatible manifest.
    """
    packages: dict[str, str] = {}
    for name in (
        "pyfaidx",
        "numpy",
        "scipy",
        "matplotlib",
        "scikit-learn",
        "shap",
    ):
        try:
            packages[name] = version(distribution_name=name)
        except PackageNotFoundError:
            continue
    return {
        "package_version": __version__,
        "python": platform.python_version(),
        "software": packages,
        "inputs": [file_fingerprint(path=p) for p in inputs],
        "settings": dict(settings),
    }
