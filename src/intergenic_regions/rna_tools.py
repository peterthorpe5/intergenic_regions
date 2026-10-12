"""Auditable STAR/samtools/featureCounts adapters; never invoke a shell."""

import csv
import logging
import math
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

from intergenic_regions.genome import Genome
from intergenic_regions.io import (
    DNA,
    file_fingerprint,
    open_text,
    write_json,
    write_tsv,
)
from intergenic_regions.models import Gene
from intergenic_regions.rna_inputs import RnaSample

LOGGER = logging.getLogger(__name__)


def inspect_fastq(*, path: Path, max_records: int = 100) -> dict[str, Any]:
    """Validate a bounded four-line FASTQ sample before launching mapping.

    Args:
        path: Plain/gzip FASTQ, already adapter/QC processed as appropriate.
        max_records: Number of complete records to inspect, not a full QC run.

    Returns:
        Sampled record count and maximum observed read length.

    Raises:
        ValueError: The inspected records are empty, truncated or malformed.
    """
    if isinstance(max_records, bool) or not 1 <= max_records <= 10000:
        raise ValueError("FASTQ inspection limit must be 1..10000")
    count, maximum = 0, 0
    with open_text(path=path) as stream:
        for _ in range(max_records):
            header = stream.readline()
            if not header:
                break
            sequence, plus, quality = (
                stream.readline().rstrip("\r\n") for _ in range(3)
            )
            if (
                not header.startswith("@")
                or not plus.startswith("+")
                or not sequence
                or not set(sequence) <= DNA
                or len(quality) != len(sequence)
                or any(not 33 <= ord(symbol) <= 126 for symbol in quality)
            ):
                raise ValueError(f"Malformed/truncated FASTQ record: {path}")
            count += 1
            maximum = max(maximum, len(sequence))
    if not count:
        raise ValueError(f"Empty FASTQ: {path}")
    return {
        "sampled_records": count,
        "maximum_sampled_read_length": maximum,
        "quality_control_status": "bounded_format_check_not_complete_read_qc",
    }


def plan_rna_mapping(
    *,
    samples: Sequence[RnaSample],
    genome_path: Path,
    directory: Path,
    star_index: Path | None = None,
    threads: int = 4,
    sort_ram: int = 2000000000,
    index_ram: int = 8000000000,
    star_executable: str = "STAR",
    samtools_executable: str = "samtools",
    read_length: int | None = None,
    junctions: bool = False,
) -> dict[str, Any]:
    """Prepare controlled splice-aware mapping and indexing argument vectors.

    Args:
        samples: Raw FASTQ biological samples only.
        genome_path: Full genomic FASTA, not a transcriptome.
        directory: Mapping working directory, also used as subprocess cwd.
        star_index: Existing STAR index, otherwise build a new index.
        threads: Explicit thread bound.
        sort_ram: Explicit STAR BAM sorting bytes.
        index_ram: Explicit index generation bytes.
        star_executable: STAR executable name or path.
        samtools_executable: samtools executable name or path.
        read_length: Index read length, or a bounded FASTQ estimate.
        junctions: The caller supplies junctions.gtf generated from full exons.

    Returns:
        Strict JSON plan. Read links use safe relative names even if original
        FASTQ paths contain spaces or shell metacharacters.

    Raises:
        ValueError: Settings or raw sample layout are invalid.
    """
    if (
        not samples
        or any(not sample.read1 for sample in samples)
        or len({s.sample_id for s in samples}) != len(samples)
        or isinstance(threads, bool)
        or not 1 <= threads <= 256
        or isinstance(sort_ram, bool)
        or isinstance(index_ram, bool)
        or sort_ram < 10000000
        or index_ram < 10000000
        or (
            read_length is not None
            and (
                isinstance(read_length, bool) or not 2 <= read_length <= 10000
            )
        )
    ):
        raise ValueError(
            "Invalid STAR samples, thread/memory limits or read length"
        )
    sampled = [
        inspect_fastq(path=path)
        for sample in samples
        for path in (sample.read1, sample.read2)
        if path is not None
    ]
    maximum = max(row["maximum_sampled_read_length"] for row in sampled)
    if read_length is not None and read_length < maximum:
        raise ValueError("STAR read_length is shorter than an inspected read")
    overhang = (read_length or maximum) - 1
    with Genome(path=genome_path) as genome:
        total, contigs = sum(genome.lengths.values()), len(genome.lengths)
    commands: list[dict[str, Any]] = []
    index = str(star_index.resolve()) if star_index else "star_index"
    if star_index is not None:
        if not all(
            (star_index / name).is_file()
            for name in ("Genome", "SA", "SAindex", "genomeParameters.txt")
        ):
            raise ValueError("Provided STAR directory is not a complete index")
        manifest = star_index / "intergenic_reference.json"
        if manifest.is_file():
            import json

            saved = json.loads(manifest.read_text(encoding="utf-8"))
            if (
                saved.get("genome_sha256")
                != file_fingerprint(path=genome_path)["sha256"]
            ):
                raise ValueError(
                    "Provided STAR index has a different genome fingerprint"
                )
    else:
        argv = [
            star_executable,
            "--runMode",
            "genomeGenerate",
            "--runThreadN",
            str(threads),
            "--genomeDir",
            index,
            "--genomeFastaFiles",
            str(genome_path.resolve()),
            "--limitGenomeGenerateRAM",
            str(index_ram),
            "--genomeSAindexNbases",
            str(min(14, max(1, int(math.log2(total) / 2 - 1)))),
            "--genomeChrBinNbits",
            str(min(18, max(3, int(math.log2(max(total / contigs, 1)))))),
            "--sjdbOverhang",
            str(overhang if junctions else 0),
        ]
        if junctions:
            argv += ["--sjdbGTFfile", "junctions.gtf"]
        commands.append({"step": "build_index", "argv": argv})
    linked_reads: list[dict[str, str]] = []
    for sample in samples:
        assert sample.read1 is not None
        compressed = sample.read1.suffix.lower() == ".gz"
        if (
            sample.read2
            and (sample.read2.suffix.lower() == ".gz") != compressed
        ):
            raise ValueError(
                "Both reads of a FASTQ pair must use the same compression"
            )
        paths: list[str] = []
        for mate, path in enumerate((sample.read1, sample.read2), start=1):
            if path is not None:
                name = f"reads/{sample.sample_id}_R{mate}.fastq" + (
                    ".gz" if compressed else ""
                )
                linked_reads.append(
                    {"source": str(path.resolve()), "link": name}
                )
                paths.append(name)
        argv = [
            star_executable,
            "--runThreadN",
            str(threads),
            "--genomeDir",
            index,
            "--readFilesIn",
            *paths,
            "--twopassMode",
            "Basic",
            "--sjdbOverhang",
            str(overhang),
            "--outFilterMultimapNmax",
            "1",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--outSAMattributes",
            "NH",
            "HI",
            "AS",
            "nM",
            "--limitBAMsortRAM",
            str(sort_ram),
            "--outFileNamePrefix",
            f"{sample.sample_id}/",
        ]
        if compressed:
            argv += ["--readFilesCommand", "gzip", "-cd"]
        commands += [
            {"step": f"map_{sample.sample_id}", "argv": argv},
            {
                "step": f"index_{sample.sample_id}",
                "argv": [
                    samtools_executable,
                    "index",
                    "-@",
                    str(threads),
                    f"{sample.sample_id}/Aligned.sortedByCoord.out.bam",
                ],
            },
        ]
    return {
        "working_directory": str(directory.resolve()),
        "commands": commands,
        "read_links": linked_reads,
        "read_length_for_index": overhang + 1,
        "threads": threads,
        "sort_ram_bytes": sort_ram,
        "index_ram_bytes": index_ram,
        "reference_verification": "fingerprinted_package_index"
        if star_index and (star_index / "intergenic_reference.json").is_file()
        else "supplied_index_identity_unverified"
        if star_index
        else "new_reference_index",
        "notes": [
            "Mapping uses a genomic reference, splice-aware two-pass STAR"
            " and unique alignments.",
            "Read QC/adapter trimming is an upstream responsibility; boun"
            "ded format checks are recorded.",
            "Dry-run plans are not mapping results. Supplied external ind"
            "ex identities need independent verification.",
        ],
    }


def _run_tool(
    *, argv: Sequence[str], directory: Path, name: str, timeout: float
) -> None:
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError(
            "RNA external-tool timeout must be positive and finite"
        )
    if shutil.which(argv[0]) is None:
        raise FileNotFoundError(
            f"RNA executable unavailable: {argv[0]}; "
            "install the RNA tools environment"
        )
    LOGGER.info("Running RNA tool step %s", name)
    with (directory / f"{name}.stdout.log").open(
        mode="w", encoding="utf-8"
    ) as stdout:
        with (directory / f"{name}.stderr.log").open(
            mode="w", encoding="utf-8"
        ) as stderr:
            try:
                result = subprocess.run(
                    args=list(argv),
                    cwd=directory,
                    stdout=stdout,
                    stderr=stderr,
                    check=False,
                    timeout=timeout,
                )
            except subprocess.TimeoutExpired as exc:
                raise TimeoutError(f"RNA tool step timed out: {name}") from exc
    if result.returncode:
        error = (directory / f"{name}.stderr.log").read_text(encoding="utf-8")[
            -4000:
        ]
        raise RuntimeError(
            f"RNA tool failed ({name}, exit {result.returncode}): {error}"
        )


def _tool_version(*, executable: str, arguments: Sequence[str]) -> str:
    try:
        result = subprocess.run(
            args=[executable, *arguments],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=10,
        )
        lines = result.stdout.decode(errors="replace").strip().splitlines()
        return lines[0][:200] if lines else "version_not_reported"
    except (OSError, subprocess.TimeoutExpired):
        return "version_unavailable"


def map_rna_reads(
    *,
    samples: Sequence[RnaSample],
    genome_path: Path,
    directory: Path,
    timeout: float = 86400,
    dry_run: bool = False,
    **settings: Any,
) -> tuple[list[RnaSample], dict[str, Any]]:
    """Run or export an auditable mapping plan in a staged directory.

    Args:
        samples: FASTQ biological samples.
        genome_path: Genomic FASTA.
        directory: Existing staged destination; caller owns atomic publication.
        timeout: Per-command wall-clock bound in seconds.
        dry_run: Write a plan without invoking aligners or claiming results.
        **settings: Keyword-only settings for plan_rna_mapping.

    Returns:
        Mapped samples and status metadata, or original samples for a dry run.

    Raises:
        RuntimeError: Mapping or indexing fails, with a bounded error excerpt.
    """
    directory.mkdir(parents=True, exist_ok=True)
    plan = plan_rna_mapping(
        samples=samples,
        genome_path=genome_path,
        directory=directory,
        junctions=False,
        **settings,
    )
    write_json(path=directory / "mapping_plan.json", data=plan)
    if dry_run:
        return list(samples), {"status": "dry_run", **plan}
    (directory / "reads").mkdir()
    if settings.get("star_index") is None:
        (directory / "star_index").mkdir()
    for link in plan["read_links"]:
        (directory / link["link"]).symlink_to(target=link["source"])
    for sample in samples:
        (directory / sample.sample_id).mkdir()
    for command in plan["commands"]:
        _run_tool(
            argv=command["argv"],
            directory=directory,
            name=command["step"],
            timeout=timeout,
        )
    if settings.get("star_index") is None:
        write_json(
            path=directory / "star_index" / "intergenic_reference.json",
            data={
                "genome_sha256": file_fingerprint(path=genome_path)["sha256"],
                "junction_source": "two_pass_read_inferred_splice_junctions",
            },
        )
    mapped: list[RnaSample] = []
    for sample in samples:
        bam = directory / sample.sample_id / "Aligned.sortedByCoord.out.bam"
        if not bam.is_file() or not Path(str(bam) + ".bai").is_file():
            raise RuntimeError(
                "Mapping returned success without a complete indexed BAM"
            )
        mapped.append(replace(sample, bam=bam, read1=None, read2=None))
    for link in plan["read_links"]:
        (directory / link["link"]).unlink()
    (directory / "reads").rmdir()
    versions = {
        "STAR": _tool_version(
            executable=settings.get("star_executable", "STAR"),
            arguments=["--version"],
        ),
        "samtools": _tool_version(
            executable=settings.get("samtools_executable", "samtools"),
            arguments=["--version"],
        ),
    }
    plan["tool_versions"] = versions
    write_json(path=directory / "mapping_plan.json", data=plan)
    fields = (
        "sample_id",
        "condition",
        "replicate_id",
        "strandedness",
        "layout",
        "size_factor",
        "bam",
    )
    extra = sorted({name for sample in samples for name in sample.covariates})
    write_tsv(
        path=directory / "mapped_samples.tsv",
        rows=[
            {
                **{
                    name: getattr(sample, name)
                    for name in fields
                    if name != "bam"
                },
                "bam": f"{sample.sample_id}/Aligned.sortedByCoord.out.bam",
                **sample.covariates,
            }
            for sample in mapped
        ],
        fields=(*fields, *extra),
    )
    return mapped, {"status": "completed", **plan}


def count_rna_genes(
    *,
    samples: Sequence[RnaSample],
    genes: Sequence[Gene],
    exons: Mapping[str, Sequence[tuple[int, int]]],
    directory: Path,
    executable: str = "featureCounts",
    threads: int = 4,
    min_mapq: int = 20,
    exclude_duplicates: bool = True,
    timeout: float = 86400,
) -> dict[str, dict[str, int]]:
    """Count annotated exonic fragments using featureCounts.

    Args:
        samples: Indexed BAM samples with explicit layout/strand convention.
        genes: Full annotation; gene spans never substitute for exons.
        exons: Merged exonic unions per exact gene ID.
        directory: Existing staged output directory.
        executable: featureCounts executable.
        threads: Explicit thread bound.
        min_mapq: Minimum mapping quality.
        exclude_duplicates: Exclude duplicates consistently with depth.
        timeout: Per-sample wall-clock bound.

    Returns:
        Raw counts for annotated exonic genes; ambiguity is excluded.

    Raises:
        ValueError: Alignment inputs, exon metadata or settings are invalid.
        RuntimeError: featureCounts fails or does not produce valid counts.
    """
    if (
        not samples
        or any(not s.bam or s.bam.suffix.lower() != ".bam" for s in samples)
        or not exons
        or isinstance(threads, bool)
        or not 1 <= threads <= 256
        or isinstance(min_mapq, bool)
        or not 0 <= min_mapq <= 255
        or len({s.sample_id for s in samples}) != len(samples)
    ):
        raise ValueError(
            "Gene counting needs BAM samples, annotated exons and valid limits"
        )
    by_id = {gene.gene_id: gene for gene in genes}
    if any(identifier not in by_id for identifier in exons):
        raise ValueError("Count exon gene is absent from full annotation")
    if len(by_id) != len(genes) or any(
        not spans
        or any(
            left < by_id[identifier].start
            or right > by_id[identifier].end
            or right <= left
            for left, right in spans
        )
        for identifier, spans in exons.items()
    ):
        raise ValueError("Counting exons conflict with their target gene")
    if any(
        not (
            Path(str(s.bam) + ".bai").is_file()
            or s.bam.with_suffix(".bai").is_file()
        )
        for s in samples
        if s.bam is not None
    ):
        raise ValueError("Gene counting requires indexed BAM inputs")
    directory.mkdir(parents=True, exist_ok=True)
    annotation = directory / "exonic_unions.saf"
    rows: list[dict[str, Any]] = [
        {
            "GeneID": gene_id,
            "Chr": by_id[gene_id].contig,
            "Start": left + 1,
            "End": right,
            "Strand": by_id[gene_id].strand,
        }
        for gene_id, spans in exons.items()
        for left, right in spans
    ]
    if any(row["Start"] < 1 or row["End"] < row["Start"] for row in rows):
        raise ValueError("Invalid counting exon interval")
    write_tsv(
        path=annotation,
        rows=rows,
        fields=("GeneID", "Chr", "Start", "End", "Strand"),
    )
    counts: dict[str, dict[str, int]] = {gene: {} for gene in exons}
    plans: list[dict[str, Any]] = []
    for sample in samples:
        output = directory / f"{sample.sample_id}.featurecounts.tsv"
        argv = [
            executable,
            "-T",
            str(threads),
            "-F",
            "SAF",
            "-a",
            str(annotation.resolve()),
            "-o",
            str(output.resolve()),
            "-s",
            str(
                {"unstranded": 0, "forward": 1, "reverse": 2}[
                    sample.strandedness
                ]
            ),
            "-Q",
            str(min_mapq),
            "--primary",
        ]
        if sample.layout == "paired":
            argv += ["-p", "--countReadPairs", "-B", "-C"]
        if exclude_duplicates:
            argv += ["--ignoreDup"]
        assert sample.bam is not None
        argv += [str(sample.bam.resolve())]
        plans.append({"sample_id": sample.sample_id, "argv": argv})
        _run_tool(
            argv=argv,
            directory=directory,
            name=f"count_{sample.sample_id}",
            timeout=timeout,
        )
        if not output.is_file():
            raise RuntimeError("featureCounts did not produce an output table")
        observed: set[str] = set()
        with output.open(mode="r", encoding="utf-8") as stream:
            reader = csv.reader(
                (line for line in stream if not line.startswith("#")),
                delimiter="\t",
            )
            headers = next(reader, [])
            if len(headers) != 7 or headers[0] != "Geneid":
                raise RuntimeError("Unexpected featureCounts output format")
            for columns in reader:
                if (
                    len(columns) != 7
                    or columns[0] not in counts
                    or columns[0] in observed
                ):
                    raise RuntimeError(
                        "Invalid/duplicate featureCounts gene row"
                    )
                value = int(columns[-1])
                if value < 0:
                    raise RuntimeError("Negative featureCounts count")
                counts[columns[0]][sample.sample_id] = value
                observed.add(columns[0])
        if observed != set(counts):
            raise RuntimeError("featureCounts omitted expected exonic genes")
    write_tsv(
        path=directory / "gene_counts.tsv",
        rows=[{"gene_id": gene, **values} for gene, values in counts.items()],
        fields=("gene_id", *(s.sample_id for s in samples)),
    )
    write_json(
        path=directory / "counting_plan.json",
        data={
            "commands": plans,
            "tool_versions": {
                "featureCounts": _tool_version(
                    executable=executable, arguments=["-v"]
                )
            },
            "unit": "exonic_reads_or_paired_fragments",
            "ambiguous_gene_assignment": "excluded",
            "multi_mapping": "excluded",
            "genes_without_exons": sorted(set(by_id) - set(exons)),
        },
    )
    return counts
