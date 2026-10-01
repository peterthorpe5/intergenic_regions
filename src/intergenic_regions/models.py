"""Immutable data models using zero-based, half-open coordinates."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True, kw_only=True)
class Gene:
    """A complete gene span, including introns and annotated UTRs.

    Attributes:
        gene_id: Exact annotation identifier; suffixes are preserved.
        contig: FASTA record identifier.
        start: Zero-based inclusive start.
        end: Zero-based exclusive end.
        strand: ``+``, ``-`` or ``.`` (unknown).
    """

    gene_id: str
    contig: str
    start: int
    end: int
    strand: str

    def __post_init__(self) -> None:
        """Reject invalid spans, identifiers and strands.

        Raises:
            ValueError: A field is invalid.
        """
        if not self.gene_id or any(c.isspace() for c in self.gene_id):
            raise ValueError("Gene identifiers must be non-empty tokens.")
        if not self.contig or any(c.isspace() for c in self.contig):
            raise ValueError("Contig identifiers must be non-empty tokens.")
        if self.start < 0 or self.end <= self.start:
            raise ValueError(f"Invalid gene span: {self.gene_id}")
        if self.strand not in {"+", "-", "."}:
            raise ValueError(f"Invalid strand: {self.strand}")


@dataclass(frozen=True, slots=True, kw_only=True)
class Region:
    """An extracted flank or an audited exclusion.

    Attributes:
        gene_id: Target identifier.
        contig: FASTA record identifier.
        start: Zero-based inclusive start.
        end: Zero-based exclusive end.
        strand: Gene strand.
        direction: Upstream or downstream relative to transcription.
        sequence: Strand-oriented sequence; empty for exclusions.
        status: ``retained`` or an exclusion reason.
        stop_reason: Boundary limiting the available intergenic flank.
        available_length: Length before the offset and length filters.
    """

    gene_id: str
    contig: str
    start: int
    end: int
    strand: str
    direction: str
    sequence: str
    status: str
    stop_reason: str
    available_length: int

    @property
    def sequence_id(self) -> str:
        """Return an unambiguous, reversible FASTA identifier."""
        from urllib.parse import quote

        return f"{quote(self.gene_id, safe='._-')}|{self.direction}"
