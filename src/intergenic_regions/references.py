"""Versioned local imports of regulatory databases with assembly checks."""

import json
import logging
import shutil
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from intergenic_regions.evidence import read_bed
from intergenic_regions.io import file_fingerprint, output_bundle, write_json
from intergenic_regions.models import Region

LOGGER = logging.getLogger(__name__)
CATALOGUE = {
    "screen-human-enhancers": {
        "name": "ENCODE SCREEN V4 human candidate enhancers (pELS and dELS)",
        "organism": "Homo_sapiens",
        "assembly": "GRCh38",
        "url": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed",
        "homepage": "https://screen.wenglab.org/downloads",
        "evidence_type": "candidate_enhancer",
    },
    "fantom5-human-enhancers": {
        "name": "FANTOM5 human transcribed enhancers",
        "organism": "Homo_sapiens",
        "assembly": "GRCh38",
        "url": "https://dbarchive.biosciencedbc.jp/data/fantom5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz",
        "homepage": "https://fantom.gsc.riken.go.jp/5/",
        "evidence_type": "transcribed_enhancer",
    },
    "plantregmap": {
        "name": "PlantRegMap user-selected regulatory element track",
        "homepage": "https://plantregmap.gao-lab.org/download.php",
        "evidence_type": "regulatory_element",
    },
    "custom": {
        "name": "User-supplied regulatory BED",
        "evidence_type": "user_annotation",
    },
}


def canonical_assembly(*, assembly: str) -> str:
    """Normalise a small explicit set of genome-build aliases.

    Args:
        assembly: Genome build identifier.

    Returns:
        Canonical identifier, retaining unknown builds exactly.

    Raises:
        ValueError: The identifier is empty or contains whitespace.
    """
    if not assembly or any(c.isspace() for c in assembly):
        raise ValueError("Assembly must be a non-empty token")
    aliases = {
        "hg38": "GRCh38",
        "hg19": "GRCh37",
        "mm10": "GRCm38",
        "mm39": "GRCm39",
    }
    return aliases.get(assembly, assembly)


def download_bed(
    *,
    url: str,
    path: Path,
    timeout: float = 60,
    max_bytes: int = 2_000_000_000,
) -> None:
    """Stream an explicitly selected public BED download to disk.

    Args:
        url: HTTPS dataset URL.
        path: Local destination.
        timeout: Per-operation network timeout in seconds.
        max_bytes: Maximum compressed/download size.

    Raises:
        ValueError: URL, response or size limit is invalid.
        OSError: Download fails. Partial output is removed.
    """
    parsed = urlparse(url=url)
    if (
        parsed.scheme != "https"
        or not parsed.netloc
        or parsed.username
        or parsed.password
        or timeout <= 0
        or max_bytes < 1
    ):
        raise ValueError(
            "A public HTTPS URL and positive download limits are required"
        )
    request = Request(
        url=url, headers={"User-Agent": "intergenic-regions/1.0"}
    )
    try:
        with urlopen(url=request, timeout=timeout) as response:
            if urlparse(response.geturl()).scheme != "https":
                raise ValueError(
                    "Reference download redirected away from HTTPS"
                )
            total = 0
            with path.open(mode="wb") as stream:
                while block := response.read(1024 * 1024):
                    total += len(block)
                    if total > max_bytes:
                        raise ValueError(
                            "Reference download exceeds max_bytes"
                        )
                    stream.write(block)
            if not total:
                raise ValueError("Reference download is empty")
    except Exception:
        path.unlink(missing_ok=True)
        raise


def import_reference(
    *,
    output: Path,
    source: str,
    local_bed: Path | None = None,
    url: str | None = None,
    organism: str | None = None,
    assembly: str | None = None,
    name: str | None = None,
    expected_sha256: str | None = None,
) -> dict[str, Any]:
    """Create a portable, hashed regulatory-reference bundle.

    Args:
        output: New reference directory.
        source: Catalogue entry or ``custom``.
        local_bed: Optional local BED/BED.gz instead of downloading.
        url: Optional explicit HTTPS BED URL.
        organism: Scientific-name token; required for custom/plant tracks.
        assembly: Build identifier; required for custom/plant tracks.
        name: Optional descriptive dataset name.
        expected_sha256: Optional expected digest for reproducible imports.

    Returns:
        Validated dataset metadata.

    Raises:
        ValueError: Source, assembly, checksum or BED content is invalid.
    """
    if source not in CATALOGUE or (local_bed is not None and url is not None):
        raise ValueError("Choose a known source and only one BED input method")
    defaults = CATALOGUE[source]
    selected_assembly = canonical_assembly(
        assembly=assembly or defaults.get("assembly", "")
    )
    selected_organism = organism or defaults.get("organism", "")
    if not selected_organism or any(c.isspace() for c in selected_organism):
        raise ValueError("Supply an organism token, e.g. Arabidopsis_thaliana")
    selected_url = url or defaults.get("url")
    if local_bed is None and not selected_url:
        raise ValueError("This source needs a selected BED URL or --bed")
    if (
        local_bed is None
        and url is None
        and selected_assembly != defaults.get("assembly")
    ):
        raise ValueError(
            "Preset URL assembly does not match requested assembly"
        )
    if (
        local_bed is None
        and url is None
        and selected_organism != defaults.get("organism")
    ):
        raise ValueError(
            "Preset URL organism does not match requested organism"
        )
    with output_bundle(path=output) as stage:
        compressed = (
            local_bed.suffix.lower() == ".gz"
            if local_bed
            else urlparse(url=selected_url or "").path.endswith(".gz")
        )
        bed = stage / ("regions.bed.gz" if compressed else "regions.bed")
        if local_bed is not None:
            shutil.copyfile(src=local_bed, dst=bed)
        else:
            download_bed(url=selected_url or "", path=bed)
        fingerprint = file_fingerprint(path=bed)
        if (
            expected_sha256
            and fingerprint["sha256"] != expected_sha256.lower()
        ):
            raise ValueError(
                "Reference SHA-256 does not match expected digest"
            )
        index = read_bed(path=bed)
        if not index.intervals:
            raise ValueError("Reference BED contains no intervals")
        metadata = {
            "schema_version": 1,
            "source": source,
            "name": name or defaults["name"],
            "organism": selected_organism,
            "assembly": selected_assembly,
            "evidence_type": defaults["evidence_type"],
            "source_url": selected_url if local_bed is None else None,
            "local_source": str(local_bed.resolve()) if local_bed else None,
            "bed_file": bed.name,
            "sha256": fingerprint["sha256"],
            "bytes": fingerprint["bytes"],
            "imported_utc": datetime.now(tz=UTC).isoformat(),
            "coordinate_system": "zero-based half-open",
            "merged_intervals": sum(
                len(spans) for spans in index.intervals.values()
            ),
            "homepage": defaults.get("homepage"),
        }
        write_json(path=stage / "reference.json", data=metadata)
    LOGGER.info(
        "Imported %s reference for %s/%s",
        source,
        selected_organism,
        selected_assembly,
    )
    return metadata


def query_references(
    *,
    regions: Sequence[Region],
    reference_directories: Sequence[Path],
    organism: str,
    assembly: str,
) -> list[dict[str, Any]]:
    """Query local references only after checking organism/build and checksums.

    Args:
        regions: Retained extracted regions.
        reference_directories: Previously imported reference bundles.
        organism: Query genome's scientific-name token.
        assembly: Query genome assembly.

    Returns:
        Per-region, per-database overlap and coverage records.

    Raises:
        ValueError: Metadata, checksum, organism or assembly does not match.
    """
    canonical = canonical_assembly(assembly=assembly)
    rows: list[dict[str, Any]] = []
    for directory in reference_directories:
        metadata = json.loads(
            (directory / "reference.json").read_text(encoding="utf-8")
        )
        if (
            metadata.get("schema_version") != 1
            or metadata.get("organism") != organism
            or canonical_assembly(assembly=metadata.get("assembly", ""))
            != canonical
        ):
            raise ValueError(
                f"Reference organism/assembly mismatch: {directory}"
            )
        filename = metadata.get("bed_file", "")
        if filename not in {"regions.bed", "regions.bed.gz"}:
            raise ValueError("Invalid reference bundle BED filename")
        bed = directory / filename
        if file_fingerprint(path=bed)["sha256"] != metadata.get("sha256"):
            raise ValueError(f"Reference checksum changed: {directory}")
        index = read_bed(path=bed)
        contigs = {
            region.contig for region in regions if region.status == "retained"
        }
        if contigs and not contigs & index.intervals.keys():
            raise ValueError(
                "Reference and query have no shared contig IDs; "
                "provide matching assembly naming"
            )
        for region in regions:
            if region.status != "retained":
                continue
            bases = index.overlap_bases(
                contig=region.contig, start=region.start, end=region.end
            )
            rows.append(
                {
                    "sequence_id": region.sequence_id,
                    "gene_id": region.gene_id,
                    "reference_source": metadata["source"],
                    "reference_name": metadata["name"],
                    "reference_assembly": canonical,
                    "reference_sha256": metadata["sha256"],
                    "evidence_type": metadata["evidence_type"],
                    "overlap_bp": bases,
                    "overlap_fraction": bases / (region.end - region.start),
                }
            )
    return rows
