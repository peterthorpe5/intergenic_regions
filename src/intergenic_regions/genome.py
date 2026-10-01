"""Indexed genome access without modifying input directories."""

import gzip
import logging
import shutil
import tempfile
from pathlib import Path
from types import TracebackType

from pyfaidx import Fasta, FastaIndexingError

from intergenic_regions.io import DNA, open_text

LOGGER = logging.getLogger(__name__)

COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn"
)


def normalise_genome(*, source: Path, destination: Path) -> None:
    """Stream irregularly wrapped FASTA into a private indexable copy.

    Args:
        source: Original plain or gzip FASTA.
        destination: Temporary normalised FASTA path.

    Raises:
        ValueError: IDs repeat, records are empty, DNA is malformed or
            source and destination are the same path.
    """
    if source.resolve() == destination.resolve():
        raise ValueError("Normalisation source and destination must differ")
    seen: set[str] = set()
    current: str | None = None
    buffer = ""
    record_bases = 0
    with open_text(path=source) as stream:
        with destination.open(mode="w", encoding="utf-8") as output:
            for raw in stream:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current is not None and not record_bases:
                        raise ValueError("Genome contains an empty record")
                    if buffer:
                        output.write(buffer + "\n")
                    fields = line[1:].split()
                    if not fields or fields[0] in seen:
                        raise ValueError("Invalid or duplicate genome ID")
                    current = fields[0]
                    seen.add(current)
                    buffer, record_bases = "", 0
                    output.write(f">{current}\n")
                else:
                    if current is None or not set(line) <= DNA:
                        raise ValueError("Invalid DNA or genome FASTA")
                    record_bases += len(line)
                    buffer += line
                    stop = len(buffer) - len(buffer) % 80
                    for start in range(0, stop, 80):
                        output.write(buffer[start : start + 80] + "\n")
                    buffer = buffer[stop:]
            if current is None or not record_bases:
                raise ValueError(
                    "Genome contains no sequence or an empty record"
                )
            if buffer:
                output.write(buffer + "\n")


def reverse_complement(*, sequence: str) -> str:
    """Reverse-complement DNA, preserving IUPAC ambiguity and case.

    Args:
        sequence: A DNA sequence.

    Returns:
        The reverse complement.

    Raises:
        ValueError: Non-DNA symbols are present.
    """
    if not set(sequence) <= DNA:
        raise ValueError("Sequence contains non-DNA symbols")
    return sequence.translate(COMPLEMENT)[::-1]


def sequence_composition(*, sequence: str) -> dict[str, float | int]:
    """Calculate GC among ACGT bases and ambiguity among all bases.

    Args:
        sequence: DNA, possibly empty.

    Returns:
        Length, GC fraction and ambiguous fraction; empty GC is zero.
    """
    upper = sequence.upper()
    valid = sum(upper.count(c) for c in "ACGT")
    return {
        "length": len(sequence),
        "gc_fraction": (upper.count("G") + upper.count("C")) / valid
        if valid
        else 0.0,
        "ambiguous_fraction": 1.0 - valid / len(sequence) if sequence else 0.0,
    }


class Genome:
    """A context-managed pyfaidx genome with a private temporary index.

    Ordinary gzip genomes are streamed to temporary disk before indexing.
    Set TMPDIR to a writable scratch directory on HPC systems.
    """

    def __init__(self, *, path: Path) -> None:
        """Prepare a genome accessor; indexing occurs on context entry.

        Args:
            path: Genome FASTA, optionally gzip compressed.
        """
        self.path = path
        self._temporary: tempfile.TemporaryDirectory[str] | None = None
        self._fasta: Fasta | None = None

    def __enter__(self) -> "Genome":
        """Open a genome and build its private index.

        Returns:
            This accessor.

        Raises:
            ValueError: The FASTA cannot be indexed or is empty.
        """
        self._temporary = tempfile.TemporaryDirectory(
            prefix="intergenic-genome-"
        )
        temporary = Path(self._temporary.name)
        path = self.path
        try:
            if path.suffix.lower() == ".gz":
                path = temporary / "genome.fasta"
                with gzip.open(filename=self.path, mode="rb") as source:
                    with path.open(mode="wb") as destination:
                        shutil.copyfileobj(fsrc=source, fdst=destination)
            settings = {
                "indexname": str(temporary / "genome.fai"),
                "as_raw": True,
                "strict_bounds": True,
                "duplicate_action": "stop",
                "sequence_always_upper": False,
            }
            try:
                self._fasta = Fasta(filename=str(path), **settings)
            except FastaIndexingError as exc:
                if "Line length" not in str(exc):
                    raise
                LOGGER.info(
                    "Normalising irregular FASTA wrapping: %s", self.path
                )
                normalised = temporary / "normalised.fasta"
                normalise_genome(source=self.path, destination=normalised)
                self._fasta = Fasta(filename=str(normalised), **settings)
            if not self._fasta.keys():
                raise ValueError("Genome has no FASTA records")
            if any(len(self._fasta[key]) == 0 for key in self._fasta.keys()):
                raise ValueError("Genome contains an empty FASTA record")
        except FastaIndexingError as exc:
            self.close()
            raise ValueError(f"Cannot index genome FASTA: {exc}") from exc
        except Exception:
            self.close()
            raise
        return self

    def close(self) -> None:
        """Close file handles and remove private temporary files."""
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None
        if self._temporary is not None:
            self._temporary.cleanup()
            self._temporary = None

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Release the indexed genome, including after failures."""
        self.close()

    @property
    def lengths(self) -> dict[str, int]:
        """Return contig lengths for an open genome.

        Raises:
            RuntimeError: The accessor is closed.
        """
        if self._fasta is None:
            raise RuntimeError("Genome must be opened with a context manager")
        return {name: len(self._fasta[name]) for name in self._fasta.keys()}

    def fetch(self, *, contig: str, start: int, end: int) -> str:
        """Fetch exactly the requested half-open interval.

        Args:
            contig: FASTA identifier.
            start: Inclusive zero-based start.
            end: Exclusive zero-based end.

        Returns:
            DNA retaining its original case.

        Raises:
            ValueError: Contig, bounds or DNA symbols are invalid.
            RuntimeError: The accessor is closed.
        """
        if self._fasta is None:
            raise RuntimeError("Genome is closed")
        if contig not in self._fasta:
            raise ValueError(f"Annotation contig absent from genome: {contig}")
        if start < 0 or end < start or end > len(self._fasta[contig]):
            raise ValueError(
                f"Out-of-bounds genome interval: {contig}:{start}-{end}"
            )
        sequence = str(self._fasta[contig][start:end])
        if len(sequence) != end - start or not set(sequence) <= DNA:
            raise ValueError(
                f"Invalid genome sequence: {contig}:{start}-{end}"
            )
        return sequence
