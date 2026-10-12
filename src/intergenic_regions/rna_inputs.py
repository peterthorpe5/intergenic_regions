"""Explicit RNA metadata, exon hierarchies and bounded coverage reads."""

import csv
import logging
import math
import re
from bisect import bisect_right
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from itertools import accumulate
from pathlib import Path
from typing import Any

from intergenic_regions.annotation import (
    Feature,
    parse_feature,
    read_annotation,
    resolve_roots,
)
from intergenic_regions.extraction import merge_intervals
from intergenic_regions.io import open_text
from intergenic_regions.models import Gene

LOGGER = logging.getLogger(__name__)
STRAND_MODES = {"unstranded", "forward", "reverse"}


@dataclass(frozen=True, slots=True, kw_only=True)
class RnaSample:
    """One biological observation; technical libraries must be merged first.

    Attributes:
        sample_id: Unique safe token used for output names.
        condition: Experimental condition token.
        replicate_id: Biological replicate identifier within a condition.
        strandedness: Library convention, relative to read 1 for paired reads.
        layout: single/paired for BAM/FASTQ, coverage for tracks.
        bam: Optional coordinate-sorted, indexed BAM or CRAM.
        coverage: Optional unstranded non-negative bedGraph.
        coverage_plus: Optional transcript-plus-strand bedGraph.
        coverage_minus: Optional transcript-minus-strand bedGraph.
        read1: Optional raw FASTQ, plain or gzip compressed.
        read2: Optional paired FASTQ.
        size_factor: Positive supplied depth divisor, otherwise raw depth.
        covariates: Additional categorical design metadata, e.g. batch/subject.
    """

    sample_id: str
    condition: str
    replicate_id: str
    strandedness: str
    layout: str = "coverage"
    bam: Path | None = None
    coverage: Path | None = None
    coverage_plus: Path | None = None
    coverage_minus: Path | None = None
    read1: Path | None = None
    read2: Path | None = None
    size_factor: float = 1.0
    covariates: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Reject ambiguous sources, unsafe names and invalid metadata."""
        for value in (self.sample_id, self.condition, self.replicate_id):
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
                raise ValueError(
                    "RNA identifiers require safe non-empty tokens"
                )
        if self.strandedness not in STRAND_MODES:
            raise ValueError("Strandedness must be unstranded/forward/reverse")
        if self.layout not in {"single", "paired", "coverage"}:
            raise ValueError("RNA layout must be single/paired/coverage")
        if (
            isinstance(self.size_factor, bool)
            or not math.isfinite(self.size_factor)
            or self.size_factor <= 0
        ):
            raise ValueError("RNA size_factor must be positive and finite")
        tracks = bool(
            self.coverage or self.coverage_plus or self.coverage_minus
        )
        if sum((bool(self.bam), tracks, bool(self.read1))) > 1:
            raise ValueError("Choose exactly one RNA source per sample")
        if self.read2 and not self.read1:
            raise ValueError("read2 requires read1")
        if tracks:
            if self.layout != "coverage":
                raise ValueError("bedGraph samples require layout=coverage")
            if self.coverage and (self.coverage_plus or self.coverage_minus):
                raise ValueError(
                    "Do not mix unstranded and strand-specific tracks"
                )
            if bool(self.coverage_plus) != bool(self.coverage_minus):
                raise ValueError("Supply both transcript-strand bedGraphs")
            if self.coverage and self.strandedness != "unstranded":
                raise ValueError("A single bedGraph must be marked unstranded")
            if self.coverage_plus and self.strandedness == "unstranded":
                raise ValueError(
                    "Transcript-strand tracks need a stranded label"
                )
        if self.bam or self.read1:
            if self.layout == "coverage":
                raise ValueError(
                    "BAM/FASTQ samples require single/paired layout"
                )
        if self.read1 and bool(self.read2) != (self.layout == "paired"):
            raise ValueError("FASTQ layout does not match read1/read2")


def read_rna_samples(
    *, path: Path, require_source: bool = True
) -> list[RnaSample]:
    """Read a TSV with exact sample/condition/biological-replicate metadata.

    Args:
        path: Manifest; relative input paths resolve against this file.
        require_source: False permits metadata-only samples for count DE.

    Returns:
        Samples in manifest order, with existing absolute input paths.

    Raises:
        ValueError: Rows, sources or biological observation IDs are invalid.
        FileNotFoundError: A named source does not exist.
    """
    samples: list[RnaSample] = []
    identifiers: set[str] = set()
    observations: set[tuple[str, str]] = set()
    path_fields = (
        "bam",
        "coverage",
        "coverage_plus",
        "coverage_minus",
        "read1",
        "read2",
    )
    required = {"sample_id", "condition", "replicate_id", "strandedness"}
    reserved = required | set(path_fields) | {"layout", "size_factor"}
    with open_text(path=path) as stream:
        reader = csv.DictReader(f=stream, delimiter="\t")
        headers = reader.fieldnames or []
        if not required <= set(headers) or len(headers) != len(set(headers)):
            raise ValueError(
                "RNA sample TSV needs unique sample metadata columns"
            )
        for number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{number}: malformed RNA sample row")
            paths: dict[str, Path | None] = {}
            for name in path_fields:
                item = row.get(name, "").strip()
                candidate = Path(item).expanduser() if item else None
                if candidate is not None:
                    if not candidate.is_absolute():
                        candidate = path.resolve().parent / candidate
                    candidate = candidate.resolve()
                    if not candidate.is_file():
                        raise FileNotFoundError(
                            f"RNA input does not exist: {candidate}"
                        )
                paths[name] = candidate
            sample = RnaSample(
                sample_id=row["sample_id"].strip(),
                condition=row["condition"].strip(),
                replicate_id=row["replicate_id"].strip(),
                strandedness=row["strandedness"].strip(),
                layout=row.get("layout", "coverage").strip() or "coverage",
                size_factor=float(row.get("size_factor", "1").strip() or "1"),
                covariates={
                    name: row[name].strip()
                    for name in headers
                    if name not in reserved
                },
                **paths,
            )
            if require_source and not any(
                (
                    sample.bam,
                    sample.coverage,
                    sample.coverage_plus,
                    sample.read1,
                )
            ):
                raise ValueError(
                    f"RNA sample has no source: {sample.sample_id}"
                )
            observation = (sample.condition, sample.replicate_id)
            if sample.sample_id in identifiers or observation in observations:
                raise ValueError(
                    "RNA samples/biological observations repeat; merge te"
                    "chnical "
                    "libraries before analysis"
                )
            identifiers.add(sample.sample_id)
            observations.add(observation)
            samples.append(sample)
    if not samples:
        raise ValueError("RNA sample manifest is empty")
    LOGGER.info("Read %d biological RNA observations", len(samples))
    return samples


def read_gene_exons(
    *,
    path: Path,
    annotation_format: str = "auto",
    genes: Sequence[Gene] | None = None,
) -> dict[str, list[tuple[int, int]]]:
    """Resolve exon-to-gene ancestry independently of annotation record order.

    Args:
        path: Full GFF3/GTF; coordinate TSVs have no exon information.
        annotation_format: Standard annotation format, or auto detection.
        genes: Already validated full gene spans, avoiding duplicate parsing.

    Returns:
        Merged exonic intervals per exact gene ID; missing exons stay absent.

    Raises:
        ValueError: Exons conflict with their gene, or ancestry is invalid.
    """
    full_genes = (
        list(genes)
        if genes is not None
        else read_annotation(path=path, annotation_format=annotation_format)
    )
    by_id = {gene.gene_id: gene for gene in full_genes}
    nodes: dict[str, Feature] = {}
    exons: list[Feature] = []
    detected = annotation_format
    with open_text(path=path) as stream:
        for line in stream:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            if detected == "auto":
                if len(line.rstrip().split("\t")) == 5:
                    return {}
                detected = "gtf" if re.search(r'gene_id\s+"', line) else "gff3"
            if detected == "tsv":
                return {}
            feature = parse_feature(line=line, annotation_format=detected)
            if identifier := feature.attributes.get("ID"):
                nodes[identifier] = feature
            if feature.kind in {"exon", "pseudogenic_exon"}:
                exons.append(feature)
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    cache: dict[str, tuple[str, ...]] = {}
    for exon in exons:
        direct = exon.attributes.get("gene_id")
        roots = (
            (direct,)
            if direct
            else tuple(
                sorted(
                    {
                        root
                        for parent in exon.parents
                        for root in resolve_roots(
                            identifier=parent, nodes=nodes, cache=cache
                        )
                    }
                )
            )
        )
        if not roots:
            raise ValueError("RNA exon lacks a gene_id or resolvable Parent")
        for root in roots:
            gene = by_id.get(root)
            if (
                gene is None
                or exon.contig != gene.contig
                or (
                    exon.strand not in {gene.strand, "."}
                    or exon.start < gene.start
                    or exon.end > gene.end
                )
            ):
                raise ValueError(f"Exon conflicts with gene: {root}")
            grouped[root].append((exon.start, exon.end))
    return {
        gene: merge_intervals(intervals=spans)
        for gene, spans in grouped.items()
    }


def read_bedgraph(
    *, path: Path, lengths: Mapping[str, int]
) -> dict[str, list[tuple[int, int, float]]]:
    """Read non-overlapping, per-contig sorted bedGraph in genomic coordinates.

    Args:
        path: Four-column track; omitted bases represent zero coverage.
        lengths: Exact genome contig names and lengths.

    Returns:
        Non-negative finite coverage spans, without allocating whole contigs.

    Raises:
        ValueError: Coordinates, contigs, ordering or depths are invalid.
    """
    tracks: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    with open_text(path=path) as stream:
        for number, line in enumerate(stream, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            columns = line.split()
            if len(columns) != 4:
                raise ValueError(
                    f"{path}:{number}: bedGraph requires four columns"
                )
            contig, left, right, value = columns
            start, end, depth = int(left), int(right), float(value)
            if (
                contig not in lengths
                or start < 0
                or end <= start
                or end > lengths[contig]
                or not math.isfinite(depth)
                or depth < 0
                or (tracks[contig] and start < tracks[contig][-1][1])
            ):
                raise ValueError(
                    f"{path}:{number}: invalid/overlapping bedGraph span"
                )
            tracks[contig].append((start, end, depth))
    LOGGER.info(
        "Loaded %s with %d coverage spans",
        path,
        sum(map(len, tracks.values())),
    )
    return dict(tracks)


def alignment_transcript_strand(*, read: Any, strandedness: str) -> str:
    """Translate alignment orientation using the declared library convention.

    Args:
        read: A pysam alignment; paired read 2 has the opposite convention.
        strandedness: unstranded, forward (R1 sense) or reverse (antisense).

    Returns:
        Transcript strand, or ``.`` for an unstranded library.

    Raises:
        ValueError: The convention or paired read flags are invalid.
    """
    if strandedness not in STRAND_MODES:
        raise ValueError("Invalid RNA strandedness")
    if strandedness == "unstranded":
        return "."
    if read.is_paired and read.is_read1 == read.is_read2:
        raise ValueError("Paired RNA alignment lacks a unique R1/R2 flag")
    negative = bool(read.is_reverse)
    if strandedness == "reverse":
        negative = not negative
    if read.is_paired and read.is_read2:
        negative = not negative
    return "-" if negative else "+"


class CoverageReader:
    """Reusable local coverage access; spliced gaps are never read coverage."""

    def __init__(
        self,
        *,
        sample: RnaSample,
        lengths: Mapping[str, int],
        genome_path: Path,
        min_mapq: int = 20,
        max_reads: int = 100000,
        exclude_duplicates: bool = True,
    ) -> None:
        """Open one validated source and enforce bounded per-locus work.

        Args:
            sample: A mapped alignment or coverage-track sample.
            lengths: Reference lengths.
            genome_path: Reference for CRAM decoding.
            min_mapq: Minimum mapping quality; NH > 1 is excluded when tagged.
            max_reads: Hard cap on fetched alignments per locus.
            exclude_duplicates: Exclude aligner-marked duplicate reads.

        Raises:
            ValueError: Source, reference or settings are invalid.
            ImportError: BAM access needs the optional pysam dependency.
        """
        if (
            isinstance(min_mapq, bool)
            or not 0 <= min_mapq <= 255
            or isinstance(max_reads, bool)
            or max_reads < 1
        ):
            raise ValueError("Invalid RNA mapping-quality/read cap")
        self.sample = sample
        self.lengths = dict(lengths)
        self.min_mapq = min_mapq
        self.max_reads = max_reads
        self.exclude_duplicates = exclude_duplicates
        self.handle: Any = None
        self.tracks: dict[str, dict[str, list[tuple[int, int, float]]]] = {}
        self.track_starts: dict[str, dict[str, list[int]]] = {}
        if sample.bam:
            import pysam

            self.handle = pysam.AlignmentFile(
                filename=str(sample.bam),
                mode="rc" if sample.bam.suffix.lower() == ".cram" else "rb",
                reference_filename=str(genome_path),
            )
            try:
                if not self.handle.has_index():
                    raise ValueError("RNA BAM/CRAM must be indexed")
                if (
                    self.handle.header.to_dict().get("HD", {}).get("SO")
                    != "coordinate"
                ):
                    raise ValueError(
                        "RNA BAM/CRAM must declare coordinate sorting"
                    )
                reference = dict(
                    zip(
                        self.handle.references,
                        self.handle.lengths,
                        strict=True,
                    )
                )
                if any(
                    reference.get(name) != size
                    for name, size in lengths.items()
                ):
                    raise ValueError(
                        "RNA alignment reference differs from the genome"
                    )
            except BaseException:
                self.close()
                raise
        else:
            sources = (
                {".": sample.coverage}
                if sample.coverage
                else {"+": sample.coverage_plus, "-": sample.coverage_minus}
            )
            if not all(sources.values()):
                raise ValueError("Map raw reads before requesting coverage")
            for strand, path in sources.items():
                assert path is not None
                self.tracks[strand] = read_bedgraph(path=path, lengths=lengths)
                self.track_starts[strand] = {
                    contig: [start for start, _, _ in spans]
                    for contig, spans in self.tracks[strand].items()
                }

    def close(self) -> None:
        """Close any alignment handle; repeated closing is harmless."""
        if self.handle is not None:
            self.handle.close()
            self.handle = None

    def depth(
        self, *, contig: str, start: int, end: int, strand: str
    ) -> list[float]:
        """Fetch genomic-order aligned-fragment depth within one bounded locus.

        Args:
            contig: Exact genome identifier.
            start: Zero-based inclusive start.
            end: Zero-based exclusive end.
            strand: Target transcript strand.

        Returns:
            Raw depths; overlapping paired alignments count once per fragment.
            Deletions, skipped introns and unsequenced mate gaps remain zero.

        Raises:
            ValueError: Invalid coordinates/layout, or the read cap is hit.
        """
        if (
            contig not in self.lengths
            or start < 0
            or end <= start
            or end > self.lengths[contig]
            or strand not in {"+", "-"}
        ):
            raise ValueError("Invalid RNA coverage query")
        differences = [0.0] * (end - start + 1)
        if self.tracks:
            key = "." if "." in self.tracks else strand
            spans = self.tracks[key].get(contig, [])
            starts = self.track_starts[key].get(contig, [])
            first = max(0, bisect_right(a=starts, x=start) - 1)
            for number in range(first, len(spans)):
                left, right, value = spans[number]
                if left >= end:
                    break
                if right <= start:
                    continue
                differences[max(left, start) - start] += value
                differences[min(right, end) - start] -= value
        elif self.handle is not None:
            fragments: dict[tuple[str, str], list[tuple[int, int]]] = (
                defaultdict(list)
            )
            for number, read in enumerate(
                self.handle.fetch(contig=contig, start=start, stop=end),
                start=1,
            ):
                if number > self.max_reads:
                    raise ValueError(
                        "RNA locus exceeds max_reads; narrow the search"
                    )
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_qcfail
                    or read.mapping_quality < self.min_mapq
                    or (self.exclude_duplicates and read.is_duplicate)
                    or (read.has_tag(tag="NH") and read.get_tag(tag="NH") > 1)
                ):
                    continue
                if read.is_paired != (self.sample.layout == "paired"):
                    raise ValueError(
                        "RNA BAM layout differs from sample metadata"
                    )
                transcript_strand = alignment_transcript_strand(
                    read=read, strandedness=self.sample.strandedness
                )
                if transcript_strand not in {".", strand}:
                    continue
                read_group = (
                    str(read.get_tag(tag="RG"))
                    if read.has_tag(tag="RG")
                    else ""
                )
                fragments[(read.query_name, read_group)].extend(
                    read.get_blocks()
                )
            for blocks in fragments.values():
                for left, right in merge_intervals(intervals=blocks):
                    left, right = max(start, left), min(end, right)
                    if right > left:
                        differences[left - start] += 1
                        differences[right - start] -= 1
        else:
            raise ValueError("RNA coverage reader is closed")
        return list(accumulate(differences))[:-1]


@contextmanager
def open_rna_coverage(
    *,
    sample: RnaSample,
    lengths: Mapping[str, int],
    genome_path: Path,
    min_mapq: int = 20,
    max_reads: int = 100000,
    exclude_duplicates: bool = True,
) -> Iterator[CoverageReader]:
    """Own and close a reusable bounded RNA coverage reader.

    Args:
        sample: Mapped or track sample.
        lengths: Genome reference lengths.
        genome_path: Genome path for CRAM.
        min_mapq: Mapping-quality filter.
        max_reads: Hard per-locus alignment cap.
        exclude_duplicates: Exclude duplicate-marked alignments.

    Yields:
        Reader for repeated local depth queries.
    """
    reader = CoverageReader(
        sample=sample,
        lengths=lengths,
        genome_path=genome_path,
        min_mapq=min_mapq,
        max_reads=max_reads,
        exclude_duplicates=exclude_duplicates,
    )
    try:
        yield reader
    finally:
        reader.close()
