"""Chunked consensus screening with substitution limits and gene context."""

import csv
import logging
import math
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from intergenic_regions.annotation import read_annotation
from intergenic_regions.extraction import GeneIndex
from intergenic_regions.genome import Genome, reverse_complement
from intergenic_regions.io import (
    open_text,
    read_identifiers,
    write_json,
    write_tsv,
)
from intergenic_regions.motifs import IUPAC, pattern_motif, read_motifs

LOGGER = logging.getLogger(__name__)
MATCH_STATUS = "unvalidated_sequence_match"
SITE_FIELDS = (
    "hit_id",
    "motif_id",
    "pattern",
    "contig",
    "start",
    "end",
    "site_strand",
    "mismatches",
    "matched_sequence",
    "context",
    "nearest_gene_id",
    "gene_strand",
    "gene_start_proxy",
    "distance_to_gene_start",
    "nearest_gene_tied",
    "gene_cohort",
    "source_enrichment_q_value",
    "enhancer_status",
)
PROFILE_FIELDS = (
    "motif_id",
    "cohort",
    "bin_start",
    "bin_end",
    "site_count",
    "eligible_windows",
    "hits_per_million_windows",
)


@dataclass(frozen=True, slots=True, kw_only=True)
class ScanMotif:
    """An explicit consensus screening target with discovery provenance.

    Attributes:
        motif_id: Stable discovery or supplied-motif identifier.
        pattern: DNA IUPAC consensus, at most 100 bases.
        source_kind: Original k-mer, IUPAC or consensus-converted PWM source.
        source_q_value: Discovery q-value; never a genome-hit p-value.
    """

    motif_id: str
    pattern: str
    source_kind: str = "iupac"
    source_q_value: float | None = None

    def __post_init__(self) -> None:
        """Validate identifiers, patterns and discovery significance.

        Raises:
            ValueError: A field is malformed or outside its supported range.
        """
        pattern_motif(motif_id=self.motif_id, pattern=self.pattern)
        if self.pattern != self.pattern.upper():
            raise ValueError("Scan motifs require uppercase consensus")
        if self.source_q_value is not None and (
            not math.isfinite(self.source_q_value)
            or not 0 <= self.source_q_value <= 1
        ):
            raise ValueError("Invalid source enrichment q-value")


def read_scan_motifs(
    *,
    enrichment_path: Path | None = None,
    motif_path: Path | None = None,
    motif_format: str = "auto",
    q_threshold: float = 0.05,
    max_motifs: int = 20,
    motif_ids: Sequence[str] = (),
) -> list[ScanMotif]:
    """Select enriched consensuses or import user-supplied screening motifs.

    Args:
        enrichment_path: Generated motif_enrichment.tsv; retains positive
            enrichment meeting the threshold, ordered by q-value and ID.
        motif_path: Alternative IUPAC, MEME or JASPAR input. PWMs are converted
            explicitly to maximum-probability consensuses, not PWM scoring.
        motif_format: Input motif format.
        q_threshold: Discovery q-value limit; unrelated to scan significance.
        max_motifs: Maximum selected patterns; the summary records the limit.
        motif_ids: Optional exact IDs to retain, subject to the threshold.

    Returns:
        Deterministically selected consensus targets; possibly empty.

    Raises:
        ValueError: Sources, thresholds, fields or requested IDs are invalid.
    """
    if (
        (enrichment_path is None) == (motif_path is None)
        or not math.isfinite(q_threshold)
        or not 0 <= q_threshold <= 1
        or not isinstance(max_motifs, int)
        or isinstance(max_motifs, bool)
        or max_motifs < 1
        or len(set(motif_ids)) != len(motif_ids)
    ):
        raise ValueError(
            "Supply one scan motif source and valid selection limits"
        )
    selected: list[ScanMotif] = []
    seen: set[str] = set()
    if motif_path is not None:
        for motif in read_motifs(path=motif_path, motif_format=motif_format):
            seen.add(motif.motif_id)
            selected.append(
                ScanMotif(
                    motif_id=motif.motif_id,
                    pattern=motif.pattern
                    or "".join(
                        "ACGT"[int(np.argmax(row))] for row in motif.matrix
                    ),
                    source_kind="iupac" if motif.pattern else "pwm_consensus",
                )
            )
    else:
        assert enrichment_path is not None
        with open_text(path=enrichment_path) as stream:
            reader = csv.DictReader(f=stream, delimiter="\t")
            required = {
                "motif_id",
                "consensus",
                "kind",
                "q_value",
                "positive_fraction",
                "negative_fraction",
            }
            if not required <= set(reader.fieldnames or []):
                raise ValueError("Enrichment TSV lacks scan selection columns")
            for row in reader:
                try:
                    target = ScanMotif(
                        motif_id=row["motif_id"],
                        pattern=row["consensus"],
                        source_kind="pwm_consensus"
                        if row["kind"] == "pwm"
                        else row["kind"],
                        source_q_value=float(row["q_value"]),
                    )
                    fractions = [
                        float(row[field])
                        for field in ("positive_fraction", "negative_fraction")
                    ]
                    if target.motif_id in seen or any(
                        not math.isfinite(x) or not 0 <= x <= 1
                        for x in fractions
                    ):
                        raise ValueError(
                            "Duplicate motif or invalid prevalence"
                        )
                except (TypeError, KeyError) as exc:
                    raise ValueError("Malformed enrichment scan row") from exc
                seen.add(target.motif_id)
                if (
                    target.source_q_value is not None
                    and target.source_q_value <= q_threshold
                    and fractions[0] > fractions[1]
                ):
                    selected.append(target)
        selected.sort(key=lambda item: (item.source_q_value, item.motif_id))
    missing = set(motif_ids) - seen
    if missing:
        raise ValueError(
            f"Unknown requested scan motif IDs: {sorted(missing)}"
        )
    return [m for m in selected if not motif_ids or m.motif_id in motif_ids][
        :max_motifs
    ]


def _site_arrays(
    *, sequence: str, pattern: str, max_mismatches: int, both_strands: bool
) -> tuple[NDArray[np.int64], NDArray[np.uint16], NDArray[np.str_]]:
    target = pattern_motif(motif_id="query", pattern=pattern)
    width = len(target.matrix)
    if (
        not isinstance(max_mismatches, int)
        or isinstance(max_mismatches, bool)
        or not 0 <= max_mismatches < width
    ):
        raise ValueError("Mismatches must be an integer below motif length")
    upper = sequence.upper()
    if not set(upper) <= IUPAC.keys():
        raise ValueError("Scan sequence contains invalid DNA symbols")
    codes = np.frombuffer(upper.encode("ascii"), dtype=np.uint8)
    count = max(0, len(sequence) - width + 1)
    valid = np.isin(codes, np.frombuffer(b"ACGT", dtype=np.uint8))
    prefix = np.concatenate(([0], np.cumsum(~valid)))
    eligible = (
        (prefix[width:] - prefix[:-width]) == 0
        if count
        else np.zeros(0, dtype=bool)
    )
    distances: list[NDArray[np.uint16]] = []
    orientations = [pattern.upper()]
    reverse = reverse_complement(sequence=pattern.upper())
    if both_strands:
        orientations.append(reverse)
    for consensus in orientations:
        mismatches = np.zeros(count, dtype=np.uint16)
        for position, letter in enumerate(consensus):
            accepted = np.frombuffer(
                IUPAC[letter].encode("ascii"), dtype=np.uint8
            )
            mismatches += ~np.isin(
                codes[position : position + count], accepted
            )
        distances.append(mismatches)
    forward = eligible & (distances[0] <= max_mismatches)
    reverse_hit = (
        eligible & (distances[1] <= max_mismatches)
        if both_strands
        else np.zeros(count, dtype=bool)
    )
    positions = np.flatnonzero(forward | reverse_hit).astype(np.int64)
    minimum = (
        np.minimum(distances[0], distances[1])
        if both_strands
        else distances[0]
    )
    strands = np.where(
        forward[positions] & reverse_hit[positions],
        ".",
        np.where(forward[positions], "+", "-"),
    )
    return positions, minimum[positions], strands


def consensus_sites(
    *,
    sequence: str,
    pattern: str,
    max_mismatches: int = 0,
    both_strands: bool = True,
) -> list[dict[str, Any]]:
    """Find overlapping physical sites using IUPAC-aware Hamming distance.

    Args:
        sequence: DNA; windows containing ambiguous genome bases are excluded.
        pattern: IUPAC target; degeneracy does not consume a mismatch.
        max_mismatches: Maximum substitutions, excluding indels.
        both_strands: Search reverse complements too. A site matching both
            orientations appears once with strand ``.`` and minimum distance.

    Returns:
        Half-open positions, strand and substitution count for each site.

    Raises:
        ValueError: DNA, pattern or mismatch limit is invalid.
    """
    positions, differences, strands = _site_arrays(
        sequence=sequence,
        pattern=pattern,
        max_mismatches=max_mismatches,
        both_strands=both_strands,
    )
    return [
        {
            "start": int(start),
            "end": int(start) + len(pattern),
            "site_strand": str(strand),
            "mismatches": int(difference),
        }
        for start, difference, strand in zip(
            positions, differences, strands, strict=True
        )
    ]


class _GeneStarts:
    def __init__(self, *, index: GeneIndex, cohorts: dict[str, str]) -> None:
        self.index = index
        self.cohorts = cohorts
        self.genes: dict[str, list[Any]] = defaultdict(list)
        self.positions: dict[str, NDArray[np.int64]] = {}
        self.strands: dict[str, NDArray[np.int64]] = {}
        self.labels: dict[str, NDArray[np.str_]] = {}
        for gene in index.genes.values():
            if gene.strand != ".":
                self.genes[gene.contig].append(gene)
        for contig, genes in self.genes.items():
            genes.sort(
                key=lambda g: (
                    g.start if g.strand == "+" else g.end - 1,
                    g.gene_id,
                )
            )
            self.positions[contig] = np.asarray(
                [g.start if g.strand == "+" else g.end - 1 for g in genes],
                dtype=np.int64,
            )
            self.strands[contig] = np.asarray(
                [1 if g.strand == "+" else -1 for g in genes],
                dtype=np.int64,
            )
            self.labels[contig] = np.asarray(
                [cohorts.get(g.gene_id, "other") for g in genes],
            )

    def locate(
        self, *, contig: str, centres: NDArray[np.float64]
    ) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.bool_]]:
        positions = self.positions.get(contig)
        if positions is None:
            return (
                np.full(len(centres), -1, dtype=np.int64),
                np.full(len(centres), np.nan),
                np.zeros(len(centres), dtype=bool),
            )
        right = np.searchsorted(positions, centres, side="left")
        right = np.minimum(right, len(positions) - 1)
        left = np.maximum(right - 1, 0)
        left_distance = np.abs(centres - positions[left])
        right_distance = np.abs(centres - positions[right])
        nearest = np.where(left_distance <= right_distance, left, right)
        # Identical starts choose the first stable ID and retain ambiguity.
        nearest = np.searchsorted(positions, positions[nearest], side="left")
        equal_start = (
            np.searchsorted(positions, positions[nearest], side="right")
            - np.searchsorted(positions, positions[nearest], side="left")
            > 1
        )
        tied = equal_start | (
            (left != right) & (left_distance == right_distance)
        )
        signed = (centres - positions[nearest]) * self.strands[contig][nearest]
        return nearest.astype(np.int64), signed, tied

    def overlaps(
        self, *, contig: str, starts: NDArray[np.int64], width: int
    ) -> NDArray[np.bool_]:
        spans = self.index.spans.get(contig, [])
        if not spans:
            return np.zeros(len(starts), dtype=bool)
        beginnings = np.asarray([a for a, _ in spans])
        ends = np.asarray([b for _, b in spans])
        previous = np.searchsorted(beginnings, starts + width, side="left") - 1
        return (previous >= 0) & (ends[np.maximum(previous, 0)] > starts)


def genome_scan_outputs(
    *,
    directory: Path,
    genome_path: Path,
    motifs: Sequence[ScanMotif],
    annotation_path: Path | None = None,
    annotation_format: str = "auto",
    positive_genes: Path | None = None,
    negative_genes: Path | None = None,
    max_mismatches: int = 0,
    both_strands: bool = True,
    intergenic_only: bool = False,
    mask_lowercase: bool = False,
    chunk_size: int = 250000,
    max_hits: int = 2000000,
    upstream: int = 2000,
    downstream: int = 2000,
    bin_width: int = 100,
    selection_settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Stream a complete genome scan, audited sites and positional summaries.

    Args:
        directory: Output directory inside an atomic workflow bundle.
        genome_path: Indexed FASTA, including ordinary gzip input.
        motifs: Explicit consensus targets with optional discovery metadata.
        annotation_path: Full annotation for overlap and gene-start proxies.
        annotation_format: GFF3, GTF, legacy TSV or automatic inference.
        positive_genes: Optional discovery foreground IDs for profile cohorts.
        negative_genes: Optional control IDs, disjoint from foreground.
        max_mismatches: Substitution limit applied to every target.
        both_strands: Include reverse complements, merging physical duplicates.
        intergenic_only: Exclude any site overlapping any annotated gene.
        mask_lowercase: Exclude soft-masked genome windows.
        chunk_size: Core bases per chunk; overlapping tails protect boundaries.
        max_hits: Abort when exceeded, without publishing a truncated scan.
        upstream: Profile range upstream of the annotation-derived start proxy.
        downstream: Profile range downstream of that proxy, exclusive.
        bin_width: Positional histogram bin width in bases.
        selection_settings: Optional selection provenance for the report.

    Returns:
        Scan status, complete counts and coordinate/method interpretation.

    Raises:
        ValueError: Inputs/settings are invalid or the hit cap is exceeded.
    """
    integers = (chunk_size, max_hits, upstream, downstream, bin_width)
    if (
        any(not isinstance(x, int) or isinstance(x, bool) for x in integers)
        or chunk_size < 1
        or max_hits < 1
        or min(upstream, downstream) < 0
        or upstream + downstream < 1
        or bin_width < 1
        or math.ceil((upstream + downstream) / bin_width) > 2000
        or len({m.motif_id for m in motifs}) != len(motifs)
        or not isinstance(max_mismatches, int)
        or isinstance(max_mismatches, bool)
        or max_mismatches < 0
        or any(max_mismatches >= len(m.pattern) for m in motifs)
        or (
            annotation_path is None
            and (intergenic_only or positive_genes or negative_genes)
        )
    ):
        raise ValueError(
            "Invalid genome-scan limits, motifs or annotation settings"
        )
    directory.mkdir(parents=True, exist_ok=True)
    positives = read_identifiers(path=positive_genes) if positive_genes else []
    negatives = read_identifiers(path=negative_genes) if negative_genes else []
    if set(positives) & set(negatives):
        raise ValueError("Genome-scan foreground and control genes overlap")
    cohorts = {
        **dict.fromkeys(positives, "positive"),
        **dict.fromkeys(negatives, "negative"),
    }
    genes = (
        read_annotation(
            path=annotation_path, annotation_format=annotation_format
        )
        if annotation_path
        else []
    )
    edges = np.arange(-upstream, downstream, bin_width, dtype=float)
    edges = np.append(edges, downstream)
    bins = len(edges) - 1
    profiles: dict[tuple[str, str], NDArray[np.int64]] = defaultdict(
        lambda: np.zeros(bins, dtype=np.int64)
    )
    opportunities: dict[tuple[int, str], NDArray[np.int64]] = defaultdict(
        lambda: np.zeros(bins, dtype=np.int64)
    )
    burden: Counter[tuple[str, str, int]] = Counter()
    gene_burden: Counter[tuple[str, str, str]] = Counter()
    preview: list[dict[str, Any]] = []
    hit_count = 0
    ambiguous_ties = 0
    genic_count = 0
    with Genome(path=genome_path) as genome:
        index = GeneIndex(genes=genes, lengths=genome.lengths)
        if set(cohorts) - index.genes.keys():
            raise ValueError(
                "Genome-scan cohort contains unknown annotation genes"
            )
        starts_index = _GeneStarts(index=index, cohorts=cohorts)
        total_bases = sum(genome.lengths.values())
        maximum_width = max((len(m.pattern) for m in motifs), default=1)
        with (directory / "genome_motif_sites.tsv").open(
            mode="w", encoding="utf-8", newline=""
        ) as stream:
            writer = csv.DictWriter(
                f=stream,
                fieldnames=SITE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            with (directory / "genome_motif_sites.bed").open(
                mode="w", encoding="utf-8"
            ) as bed:
                for contig, length in genome.lengths.items():
                    LOGGER.info(
                        "Scanning %s (%d bp), %d motifs",
                        contig,
                        length,
                        len(motifs),
                    )
                    for first in (
                        range(0, length, chunk_size) if motifs else ()
                    ):
                        core = min(chunk_size, length - first)
                        sequence = genome.fetch(
                            contig=contig,
                            start=first,
                            end=min(length, first + core + maximum_width - 1),
                        )
                        sequence = (
                            "".join(
                                "N" if c.islower() else c for c in sequence
                            )
                            if mask_lowercase
                            else sequence
                        )
                        sequence = sequence.upper()
                        codes = np.frombuffer(
                            sequence.encode("ascii"), dtype=np.uint8
                        )
                        valid = np.isin(
                            codes, np.frombuffer(b"ACGT", dtype=np.uint8)
                        )
                        prefix = np.concatenate(([0], np.cumsum(~valid)))
                        eligible_by_width: dict[int, NDArray[np.bool_]] = {}
                        for width in sorted({len(m.pattern) for m in motifs}):
                            count = min(
                                core, max(0, len(sequence) - width + 1)
                            )
                            coords = np.arange(
                                first, first + count, dtype=np.int64
                            )
                            eligible = (
                                prefix[width : width + count] - prefix[:count]
                            ) == 0
                            overlaps = starts_index.overlaps(
                                contig=contig, starts=coords, width=width
                            )
                            if intergenic_only:
                                eligible &= ~overlaps
                            eligible_by_width[width] = eligible
                            indices, distances, tied = starts_index.locate(
                                contig=contig,
                                centres=coords.astype(float) + (width - 1) / 2,
                            )
                            local_genes = starts_index.genes.get(contig, [])
                            for cohort in ("positive", "negative", "other"):
                                mask = (
                                    eligible
                                    & ~tied
                                    & (indices >= 0)
                                    & (distances >= -upstream)
                                    & (distances < downstream)
                                )
                                if local_genes:
                                    labels = starts_index.labels[contig]
                                    mask &= (
                                        labels[np.maximum(indices, 0)]
                                        == cohort
                                    )
                                opportunities[width, cohort] += np.histogram(
                                    distances[mask], bins=edges
                                )[0]
                        for motif in motifs:
                            positions, differences, strands = _site_arrays(
                                sequence=sequence,
                                pattern=motif.pattern,
                                max_mismatches=max_mismatches,
                                both_strands=both_strands,
                            )
                            within = positions < core
                            positions, differences, strands = (
                                positions[within],
                                differences[within],
                                strands[within],
                            )
                            if len(positions):
                                keep = eligible_by_width[len(motif.pattern)][
                                    positions
                                ]
                                positions, differences, strands = (
                                    positions[keep],
                                    differences[keep],
                                    strands[keep],
                                )
                            if hit_count + len(positions) > max_hits:
                                raise ValueError(
                                    "Genome scan exceeded max_hits; use "
                                    "stricter motifs/mismatches or raise "
                                    "the explicit limit"
                                )
                            coords = positions + first
                            indices, distances, tied = starts_index.locate(
                                contig=contig,
                                centres=coords.astype(float)
                                + (len(motif.pattern) - 1) / 2,
                            )
                            genic = starts_index.overlaps(
                                contig=contig,
                                starts=coords,
                                width=len(motif.pattern),
                            )
                            local_genes = starts_index.genes.get(contig, [])
                            source_q = motif.source_q_value
                            if local_genes:
                                labels = starts_index.labels[contig][
                                    np.maximum(indices, 0)
                                ]
                                for cohort in (
                                    "positive",
                                    "negative",
                                    "other",
                                ):
                                    mask = (
                                        ~tied
                                        & (indices >= 0)
                                        & (distances >= -upstream)
                                        & (distances < downstream)
                                        & (labels == cohort)
                                    )
                                    profiles[motif.motif_id, cohort] += (
                                        np.histogram(
                                            distances[mask],
                                            bins=edges,
                                        )[0]
                                    )
                            for local, coordinate in enumerate(coords):
                                hit_count += 1
                                gene = (
                                    local_genes[int(indices[local])]
                                    if indices[local] >= 0
                                    else None
                                )
                                cohort = (
                                    cohorts.get(gene.gene_id, "other")
                                    if gene
                                    else "unassigned"
                                )
                                distance = (
                                    float(distances[local]) if gene else None
                                )
                                ambiguous_ties += bool(tied[local])
                                genic_count += bool(genic[local])
                                row = {
                                    "hit_id": f"site_{hit_count:09d}",
                                    "motif_id": motif.motif_id,
                                    "pattern": motif.pattern,
                                    "contig": contig,
                                    "start": int(coordinate),
                                    "end": int(coordinate)
                                    + len(motif.pattern),
                                    "site_strand": str(strands[local]),
                                    "mismatches": int(differences[local]),
                                    "matched_sequence": sequence[
                                        int(positions[local]) : int(
                                            positions[local]
                                        )
                                        + len(motif.pattern)
                                    ],
                                    "context": (
                                        "genic_overlap"
                                        if genic[local]
                                        else "intergenic"
                                    )
                                    if annotation_path
                                    else "annotation_not_supplied",
                                    "nearest_gene_id": gene.gene_id
                                    if gene
                                    else None,
                                    "gene_strand": gene.strand
                                    if gene
                                    else None,
                                    "gene_start_proxy": (
                                        gene.start
                                        if gene.strand == "+"
                                        else gene.end - 1
                                    )
                                    if gene
                                    else None,
                                    "distance_to_gene_start": distance,
                                    "nearest_gene_tied": bool(tied[local]),
                                    "gene_cohort": cohort,
                                    "source_enrichment_q_value": source_q,
                                    "enhancer_status": MATCH_STATUS,
                                }
                                writer.writerow(rowdict=row)
                                bed.write(
                                    f"{contig}\t{row['start']}\t{row['end']}\t{motif.motif_id}\t0\t{row['site_strand']}\n"
                                )
                                if len(preview) < 500:
                                    preview.append(row)
                                burden[
                                    motif.motif_id,
                                    contig,
                                    int(differences[local]),
                                ] += 1
                                if gene:
                                    gene_burden[
                                        gene.gene_id, motif.motif_id, cohort
                                    ] += 1
    profile_rows = []
    for motif in motifs:
        for cohort in ("positive", "negative", "other"):
            for position in range(bins):
                count = int(profiles[motif.motif_id, cohort][position])
                denominator = int(
                    opportunities[len(motif.pattern), cohort][position]
                )
                profile_rows.append(
                    {
                        "motif_id": motif.motif_id,
                        "cohort": cohort,
                        "bin_start": float(edges[position]),
                        "bin_end": float(edges[position + 1]),
                        "site_count": count,
                        "eligible_windows": denominator,
                        "hits_per_million_windows": 1000000
                        * count
                        / denominator
                        if denominator
                        else None,
                    }
                )
    burden_rows = [
        {"motif_id": m, "contig": c, "mismatches": d, "site_count": n}
        for (m, c, d), n in sorted(burden.items())
    ]
    gene_rows = [
        {"gene_id": g, "motif_id": m, "cohort": c, "site_count": n}
        for (g, m, c), n in sorted(gene_burden.items())
    ]
    write_tsv(
        path=directory / "distance_profiles.tsv",
        rows=profile_rows,
        fields=PROFILE_FIELDS,
    )
    write_tsv(
        path=directory / "scan_burden.tsv",
        rows=burden_rows,
        fields=("motif_id", "contig", "mismatches", "site_count"),
    )
    write_tsv(
        path=directory / "gene_motif_counts.tsv",
        rows=gene_rows,
        fields=("gene_id", "motif_id", "cohort", "site_count"),
    )
    write_tsv(
        path=directory / "scan_motifs.tsv",
        rows=[
            {
                "motif_id": m.motif_id,
                "pattern": m.pattern,
                "source_kind": m.source_kind,
                "source_enrichment_q_value": m.source_q_value,
            }
            for m in motifs
        ],
        fields=(
            "motif_id",
            "pattern",
            "source_kind",
            "source_enrichment_q_value",
        ),
    )
    summary = {
        "status": "completed" if motifs else "no_selected_motifs",
        "scan_method": "IUPAC consensus Hamming screening; substitutions only",
        "selected_motifs": len(motifs),
        "selection": selection_settings or {},
        "genome_bases": total_bases,
        "total_sites": hit_count,
        "genic_overlap_sites": genic_count,
        "tied_nearest_gene_sites": ambiguous_ties,
        "max_mismatches": max_mismatches,
        "both_strands": both_strands,
        "intergenic_only": intergenic_only,
        "mask_lowercase": mask_lowercase,
        "chunk_size": chunk_size,
        "max_hits": max_hits,
        "gene_annotation_supplied": annotation_path is not None,
        "profile_range": [-upstream, downstream],
        "bin_width": bin_width,
        "coordinates": "Zero-based half-open; distance uses the site midpoint",
        "distance_anchor": (
            "Annotated 5-prime gene base: start on +, end-1 on -; TSS proxy"
        ),
        "distance_sign": (
            "Negative upstream; positive downstream in gene orientation"
        ),
        "profile_normalisation": (
            "Physical sites per million eligible, unambiguous, scanned "
            "start windows of the same motif length; gene-start ties excluded"
        ),
        "interpretation": (
            "Matches are sequence hypotheses, not validated enhancers. "
            "Discovery q-values are provenance, not genome-hit significance. "
            "PWM consensus scans do not reproduce PWM thresholds."
        ),
    }
    write_json(path=directory / "summary.json", data=summary)
    from intergenic_regions.scan_reporting import scan_report

    scan_report(
        directory=directory,
        summary=summary,
        sites=preview,
        profiles=profile_rows,
        burden=burden_rows,
        genes=gene_rows,
    )
    return summary
