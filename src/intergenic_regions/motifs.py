"""Known-motif scanning and strand-aware exhaustive exact k-mer enrichment."""

import csv
import logging
import re
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from intergenic_regions.background import check_sequence_sets
from intergenic_regions.genome import reverse_complement
from intergenic_regions.io import open_text
from intergenic_regions.statistics import (
    adjust_fdr,
    enrichment_test,
    kmer_family_size,
)

LOGGER = logging.getLogger(__name__)
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


@dataclass(frozen=True, slots=True, kw_only=True)
class Motif:
    """A DNA probability matrix with stable identifier and display name."""

    motif_id: str
    name: str
    matrix: tuple[tuple[float, ...], ...]
    pattern: str | None = None

    def __post_init__(self) -> None:
        """Validate motif identifiers and probability matrices.

        Raises:
            ValueError: Identifier, shape or probabilities are invalid.
        """
        matrix = np.asarray(self.matrix, dtype=float)
        if not self.motif_id or any(c.isspace() for c in self.motif_id):
            raise ValueError("Motif IDs must be non-empty tokens")
        if (
            matrix.ndim != 2
            or matrix.shape[1] != 4
            or not 1 <= matrix.shape[0] <= 100
        ):
            raise ValueError(
                "DNA motif matrices need 1-100 rows and four columns"
            )
        if (
            not np.isfinite(matrix).all()
            or (matrix < 0).any()
            or not np.allclose(matrix.sum(axis=1), 1, atol=1e-5)
        ):
            raise ValueError(
                "Motif rows must contain finite probabilities summing to one"
            )
        if self.pattern is not None and (
            len(self.pattern) != len(self.matrix)
            or not set(self.pattern) <= IUPAC.keys()
        ):
            raise ValueError("Invalid IUPAC motif pattern")


def pattern_motif(*, motif_id: str, pattern: str, name: str = "") -> Motif:
    """Convert an IUPAC pattern to a probability matrix for visualisation.

    Args:
        motif_id: Unique identifier.
        pattern: DNA IUPAC consensus.
        name: Optional display name.

    Returns:
        An exact-pattern motif; matching uses the pattern, not its PWM.

    Raises:
        ValueError: The pattern is invalid.
    """
    upper = pattern.upper()
    if not upper or not set(upper) <= IUPAC.keys():
        raise ValueError(f"Invalid IUPAC motif: {pattern}")
    matrix = tuple(
        tuple(
            1 / len(IUPAC[c]) if base in IUPAC[c] else 0.0 for base in "ACGT"
        )
        for c in upper
    )
    return Motif(
        motif_id=motif_id, name=name or motif_id, matrix=matrix, pattern=upper
    )


def read_motifs(*, path: Path, motif_format: str = "auto") -> list[Motif]:
    """Read minimal MEME DNA matrices, JASPAR PFMs or IUPAC consensus TSV.

    Args:
        path: Local motif file, optionally gzip compressed.
        motif_format: ``auto``, ``meme``, ``jaspar`` or ``iupac``.

    Returns:
        Motifs in file order.

    Raises:
        ValueError: Format, matrices or IDs are invalid or duplicated.
    """
    with open_text(path=path) as stream:
        lines = [
            line.strip()
            for line in stream
            if line.strip() and not line.startswith("#")
        ]
    if not lines:
        raise ValueError("Empty motif file")
    detected = motif_format
    if detected == "auto":
        detected = (
            "meme"
            if lines[0].startswith("MEME version")
            else "jaspar"
            if lines[0].startswith(">")
            else "iupac"
        )
    motifs: list[Motif] = []
    if detected == "iupac":
        reader = csv.DictReader(f=lines, delimiter="\t")
        if not {"motif_id", "pattern"} <= set(reader.fieldnames or []):
            raise ValueError("IUPAC TSV requires motif_id and pattern columns")
        for row in reader:
            motifs.append(
                pattern_motif(
                    motif_id=row["motif_id"],
                    pattern=row["pattern"],
                    name=row.get("name", ""),
                )
            )
    elif detected == "jaspar":
        position = 0
        while position < len(lines):
            header = lines[position]
            if not header.startswith(">"):
                raise ValueError("JASPAR motif requires a >identifier header")
            fields = header[1:].split(maxsplit=1)
            rows: dict[str, list[float]] = {}
            for line in lines[position + 1 : position + 5]:
                match = re.fullmatch(
                    r"([ACGT])\s+\[?\s*([0-9.eE+\s-]+)\s*\]?", line
                )
                if not match or match.group(1) in rows:
                    raise ValueError("JASPAR requires unique A/C/G/T rows")
                rows[match.group(1)] = [
                    float(x) for x in match.group(2).split()
                ]
            if (
                set(rows) != set("ACGT")
                or len({len(row) for row in rows.values()}) != 1
            ):
                raise ValueError("Malformed JASPAR matrix")
            matrix = np.asarray([rows[base] for base in "ACGT"], dtype=float).T
            if (matrix < 0).any() or (matrix.sum(axis=1) <= 0).any():
                raise ValueError(
                    "JASPAR counts must be non-negative with positive row sums"
                )
            matrix /= matrix.sum(axis=1, keepdims=True)
            motifs.append(
                Motif(
                    motif_id=fields[0],
                    name=fields[-1],
                    matrix=tuple(tuple(row) for row in matrix),
                )
            )
            position += 5
    elif detected == "meme":
        if not lines[0].startswith("MEME version"):
            raise ValueError("MEME motif file lacks its version header")
        alphabet = next(
            (line for line in lines if line.startswith("ALPHABET")),
            "ALPHABET= ACGT",
        )
        if alphabet.replace(" ", "") != "ALPHABET=ACGT":
            raise ValueError(
                "Only the standard ACGT DNA alphabet is supported"
            )
        position = 0
        while position < len(lines):
            if not lines[position].startswith("MOTIF "):
                position += 1
                continue
            fields = lines[position].split(maxsplit=2)
            position += 1
            if position >= len(lines) or not lines[position].startswith(
                "letter-probability matrix:"
            ):
                raise ValueError(
                    "Minimal MEME motif requires a probability matrix"
                )
            width = re.search(r"\bw\s*=\s*(\d+)", lines[position])
            if not width:
                raise ValueError("Minimal MEME matrices must declare w=width")
            size = int(width.group(1))
            probability_rows = tuple(
                tuple(float(x) for x in row.split())
                for row in lines[position + 1 : position + 1 + size]
            )
            if len(probability_rows) != size:
                raise ValueError("Truncated MEME probability matrix")
            motifs.append(
                Motif(
                    motif_id=fields[1],
                    name=fields[-1],
                    matrix=probability_rows,
                )
            )
            position += size + 1
    else:
        raise ValueError(f"Unknown motif format: {detected}")
    if not motifs or len({m.motif_id for m in motifs}) != len(motifs):
        raise ValueError("No motifs or duplicated motif IDs")
    return motifs


def kmer_counts(
    *, sequence: str, lengths: Sequence[int], both_strands: bool = True
) -> Counter[str]:
    """Count valid exact DNA words, skipping ambiguous windows.

    Args:
        sequence: DNA sequence.
        lengths: K-mer lengths.
        both_strands: Canonicalise each word with its reverse complement.

    Returns:
        Counts keyed by exact or canonical word.
    """
    counts: Counter[str] = Counter()
    upper = sequence.upper()
    complement = upper.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for width in lengths:
        for start in range(len(upper) - width + 1):
            word = upper[start : start + width]
            if set(word) <= set("ACGT"):
                if both_strands:
                    word = min(
                        word,
                        complement[
                            len(upper) - start - width : len(upper) - start
                        ],
                    )
                counts[word] += 1
    return counts


def pooled_background(*, sequences: Sequence[str]) -> NDArray[np.float64]:
    """Estimate a strand-symmetric, zero-order DNA background without labels.

    Args:
        sequences: All analysed sequences.

    Returns:
        Positive A/C/G/T frequencies, with one pseudocount per base.
    """
    counts = np.ones(4, dtype=float)
    for sequence in sequences:
        counts += np.asarray([sequence.upper().count(b) for b in "ACGT"])
    counts = (counts + counts[::-1]) / 2
    return counts / counts.sum()


@dataclass(frozen=True, slots=True, kw_only=True)
class PWMScanner:
    """Discretised log-odds scores and their exact zero-order null tail."""

    weights: NDArray[np.int64]
    threshold: int
    minimum: int
    tail: NDArray[np.float64]
    palindrome: bool


def prepare_pwm(
    *,
    motif: Motif,
    background: NDArray[np.float64],
    site_p_value: float = 1e-4,
    resolution: int = 100,
    pseudocount: float = 0.001,
) -> PWMScanner:
    """Build a dynamic-programming score distribution for PWM scanning.

    Args:
        motif: Known DNA PWM.
        background: Strand-symmetric A/C/G/T null probabilities.
        site_p_value: Maximum per-position null tail probability.
        resolution: Integer score bins per log2 unit.
        pseudocount: Total uniform probability mass added per motif row.

    Returns:
        A scanner using the same discretised scores for calibration and hits.

    Raises:
        ValueError: Background or threshold settings are invalid.
    """
    if (
        not 0 < site_p_value < 1
        or not 1 <= resolution <= 1000
        or not 0 < pseudocount <= 1
    ):
        raise ValueError("Invalid PWM scanning parameters")
    if (
        background.shape != (4,)
        or not np.isfinite(background).all()
        or (background <= 0).any()
        or not np.isclose(background.sum(), 1)
        or not np.allclose(background, background[::-1])
    ):
        raise ValueError(
            "PWM background must be positive, normalised and strand symmetric"
        )
    matrix = (np.asarray(motif.matrix) + pseudocount / 4) / (1 + pseudocount)
    weights = np.rint(resolution * np.log2(matrix / background)).astype(
        np.int64
    )
    distribution = np.ones(1, dtype=float)
    minimum = 0
    for row in weights:
        row_minimum = int(row.min())
        shifted = row - row_minimum
        updated = np.zeros(len(distribution) + int(shifted.max()), dtype=float)
        for base in range(4):
            shift = int(shifted[base])
            updated[shift : shift + len(distribution)] += (
                background[base] * distribution
            )
        distribution = updated
        minimum += row_minimum
    tail = np.clip(np.cumsum(distribution[::-1])[::-1], 0, 1)
    passing = np.flatnonzero(tail <= site_p_value)
    threshold = (
        minimum + int(passing[0]) if len(passing) else minimum + len(tail)
    )
    return PWMScanner(
        weights=weights,
        threshold=threshold,
        minimum=minimum,
        tail=tail,
        palindrome=bool(np.array_equal(weights, weights[::-1, ::-1])),
    )


def scan_pwm(
    *, sequence: str, scanner: PWMScanner, both_strands: bool = True
) -> list[dict[str, Any]]:
    """Scan a DNA sequence with bounded-memory vectorised score accumulation.

    Args:
        sequence: Strand-oriented region sequence.
        scanner: Calibrated score distribution.
        both_strands: Scan reverse-complement motif sites as well.

    Returns:
        Sites with zero-based sequence coordinates, site strand, binned score
        and unadjusted per-site p-value. Palindromic sites are reported once.
    """
    width = len(scanner.weights)
    upper = sequence.upper()
    if len(upper) < width:
        return []
    lookup = np.full(256, 4, dtype=np.int64)
    for index, base in enumerate("ACGT"):
        lookup[ord(base)] = index
    codes = lookup[np.frombuffer(upper.encode("ascii"), dtype=np.uint8)]
    cumulative = np.concatenate(([0], np.cumsum(codes == 4)))
    valid = cumulative[width:] == cumulative[:-width]
    weights_list = [("+", scanner.weights)]
    if both_strands and not scanner.palindrome:
        weights_list.append(("-", scanner.weights[::-1, ::-1]))
    hits: list[dict[str, Any]] = []
    for strand, weights in weights_list:
        scores = np.zeros(len(upper) - width + 1, dtype=np.int64)
        for position, row in enumerate(weights):
            padded = np.append(row, 0)
            scores += padded[codes[position : position + len(scores)]]
        for start in np.flatnonzero(valid & (scores >= scanner.threshold)):
            score = int(scores[start])
            hits.append(
                {
                    "start": int(start),
                    "end": int(start) + width,
                    "site_strand": strand,
                    "binned_score": score,
                    "site_p_value": float(
                        scanner.tail[score - scanner.minimum]
                    ),
                }
            )
    return hits


def scan_pattern(
    *, sequence: str, pattern: str, both_strands: bool = True
) -> list[dict[str, Any]]:
    """Find overlapping exact IUPAC sites without matching ambiguous DNA.

    Args:
        sequence: DNA sequence.
        pattern: Valid IUPAC consensus.
        both_strands: Include reverse-complement consensus sites.

    Returns:
        Half-open sequence coordinates and site strands.

    Raises:
        ValueError: Consensus is invalid.
    """
    motif = pattern_motif(motif_id="pattern", pattern=pattern)
    patterns = [("+", motif.pattern or "")]
    reverse = reverse_complement(sequence=pattern.upper())
    if both_strands and reverse != motif.pattern:
        patterns.append(("-", reverse))
    hits: list[dict[str, Any]] = []
    for strand, consensus in patterns:
        expression = (
            "(?=(" + "".join(f"[{IUPAC[c]}]" for c in consensus) + "))"
        )
        for match in re.finditer(expression, sequence.upper()):
            hits.append(
                {
                    "start": match.start(),
                    "end": match.start() + len(consensus),
                    "site_strand": strand,
                    "binned_score": None,
                    "site_p_value": None,
                }
            )
    return hits


def analyse_motifs(
    *,
    positive: dict[str, str],
    negative: dict[str, str],
    motifs: Sequence[Motif] = (),
    kmer_lengths: Sequence[int] = (),
    both_strands: bool = True,
    site_p_value: float = 1e-4,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Run native enrichment with one joint, fully specified testing family.

    Args:
        positive: Foreground sequences.
        negative: Control sequences.
        motifs: Known PWMs or exact IUPAC patterns.
        kmer_lengths: Exact word lengths to screen exhaustively.
        both_strands: Combine reverse-complement k-mers and scan both strands.
        site_p_value: Label-independent site threshold for known PWMs.

    Returns:
        Enrichment records, known-motif sites and analysis summary.

    Raises:
        ValueError: Sequence sets, motif IDs or analysis settings are invalid.
    """
    check_sequence_sets(positive=positive, negative=negative)
    if not motifs and not kmer_lengths:
        raise ValueError("Supply known motifs or k-mer lengths")
    if len({m.motif_id for m in motifs}) != len(motifs):
        raise ValueError("Duplicated motif IDs")
    family = kmer_family_size(
        lengths=kmer_lengths, both_strands=both_strands
    ) + len(motifs)
    records: list[dict[str, Any]] = []
    sites: list[dict[str, Any]] = []
    counts: list[Counter[str]] = []
    for collection in (positive, negative):
        presence: Counter[str] = Counter()
        for sequence in collection.values():
            presence.update(
                kmer_counts(
                    sequence=sequence,
                    lengths=kmer_lengths,
                    both_strands=both_strands,
                ).keys()
            )
        counts.append(presence)
    for word in sorted(counts[0].keys() | counts[1].keys()):
        records.append(
            {
                "motif_id": f"kmer:{word}",
                "name": word,
                "kind": "kmer",
                "consensus": word,
                **enrichment_test(
                    positive_hits=counts[0][word],
                    negative_hits=counts[1][word],
                    positive_total=len(positive),
                    negative_total=len(negative),
                ),
            }
        )
    background = pooled_background(
        sequences=[*positive.values(), *negative.values()]
    )
    for motif in motifs:
        scanner = (
            prepare_pwm(
                motif=motif, background=background, site_p_value=site_p_value
            )
            if motif.pattern is None
            else None
        )
        hits_by_group: list[int] = []
        for label, collection in (
            ("positive", positive),
            ("negative", negative),
        ):
            hits_count = 0
            for identifier, sequence in collection.items():
                hits = (
                    scan_pwm(
                        sequence=sequence,
                        scanner=scanner,
                        both_strands=both_strands,
                    )
                    if scanner
                    else scan_pattern(
                        sequence=sequence,
                        pattern=motif.pattern or "",
                        both_strands=both_strands,
                    )
                )
                hits_count += bool(hits)
                sites.extend(
                    {
                        "motif_id": motif.motif_id,
                        "sequence_id": identifier,
                        "label": label,
                        **hit,
                    }
                    for hit in hits
                )
            hits_by_group.append(hits_count)
        records.append(
            {
                "motif_id": motif.motif_id,
                "name": motif.name,
                "kind": "iupac" if motif.pattern else "pwm",
                "consensus": motif.pattern
                or "".join(
                    "ACGT"[int(np.argmax(row))] for row in motif.matrix
                ),
                **enrichment_test(
                    positive_hits=hits_by_group[0],
                    negative_hits=hits_by_group[1],
                    positive_total=len(positive),
                    negative_total=len(negative),
                ),
            }
        )
    adjusted = adjust_fdr(
        log_p_values=[r["log_p_value"] for r in records], family_size=family
    )
    for record, correction in zip(records, adjusted, strict=True):
        record.update(correction)
    records.sort(
        key=lambda r: (r["log_q_value"], r["log_p_value"], r["motif_id"])
    )
    summary = {
        "positive_sequences": len(positive),
        "negative_sequences": len(negative),
        "reported_hypotheses": len(records),
        "testing_family_size": family,
        "background_acgt": background.tolist(),
        "site_p_value": site_p_value,
        "both_strands": both_strands,
        "significant_q_0_05": int(
            sum(r["log_q_value"] <= np.log(0.05) for r in records)
        ),
        "test": "one-sided sequence-presence Fisher exact",
        "fdr": "Benjamini-Hochberg across known motifs and all possible k-mers",
    }
    LOGGER.info(
        "Tested %d hypotheses; %d pass q <= 0.05",
        family,
        summary["significant_q_0_05"],
    )
    return records, sites, summary
