"""Parse unsorted GFF3, GTF and legacy five-column coordinate tables."""

import logging
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from intergenic_regions.io import open_text
from intergenic_regions.models import Gene

LOGGER = logging.getLogger(__name__)
GENIC_FEATURES = {
    "gene",
    "pseudogene",
    "transcript",
    "mrna",
    "exon",
    "cds",
    "utr",
    "five_prime_utr",
    "three_prime_utr",
    "ncrna",
    "trna",
    "rrna",
    "snrna",
    "snorna",
    "mirna",
    "lncrna",
    "primary_transcript",
    "pseudogenic_transcript",
    "pseudogenic_exon",
    "pseudogenic_cds",
}


@dataclass(frozen=True, slots=True, kw_only=True)
class Feature:
    """A parsed annotation record used to resolve the parent hierarchy."""

    contig: str
    kind: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str]
    parents: tuple[str, ...]


def parse_attributes(*, text: str, annotation_format: str) -> dict[str, str]:
    """Parse GFF3 key/value pairs or quoted GTF attributes.

    Args:
        text: Ninth annotation column.
        annotation_format: ``gff3`` or ``gtf``.

    Returns:
        Attributes with GFF3 percent escapes decoded.

    Raises:
        ValueError: Attributes are malformed or duplicated.
    """
    attributes: dict[str, str] = {}
    if text == ".":
        return attributes
    if annotation_format == "gtf":
        matches = list(re.finditer(r'([^\s;]+)\s+"([^"\n]*)"\s*;?', text))
        remainder = re.sub(r'([^\s;]+)\s+"([^"\n]*)"\s*;?', "", text)
        if remainder.strip():
            raise ValueError(
                "Malformed GTF attributes; quoted values required"
            )
        pairs = [(m.group(1), m.group(2)) for m in matches]
    elif annotation_format == "gff3":
        pairs = []
        for item in text.rstrip(";").split(";"):
            if "=" not in item:
                raise ValueError(
                    "Malformed GFF3 attribute: expected key=value"
                )
            key, value = item.split("=", maxsplit=1)
            pairs.append((key, unquote(value)))
    else:
        raise ValueError(f"Unsupported attribute format: {annotation_format}")
    for key, value in pairs:
        if not key or key in attributes:
            raise ValueError(f"Empty or repeated attribute: {key}")
        attributes[key] = value
    return attributes


def parse_feature(*, line: str, annotation_format: str) -> Feature:
    """Convert a nine-column annotation line to internal coordinates.

    Args:
        line: A GFF3 or GTF data row.
        annotation_format: ``gff3`` or ``gtf``.

    Returns:
        A validated feature.

    Raises:
        ValueError: Columns, coordinates or strand are invalid.
    """
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != 9:
        raise ValueError("Annotation rows require nine tab-separated columns")
    contig, _, kind, start, end, score, strand, phase, text = fields
    attributes = parse_attributes(
        text=text, annotation_format=annotation_format
    )
    first, last = int(start), int(end)
    if (
        not contig
        or not kind
        or any(c.isspace() for c in contig)
        or phase not in {".", "0", "1", "2"}
        or (score != "." and not math.isfinite(float(score)))
    ):
        raise ValueError("Invalid annotation contig, feature, score or phase")
    if first < 1 or last < first or strand not in {"+", "-", ".", "?"}:
        raise ValueError("Invalid annotation coordinates or strand")
    # Split Parent before decoding so escaped commas remain part of an ID.
    parent_text = (
        next(
            (
                s.split("=", 1)[1]
                for s in text.split(";")
                if s.startswith("Parent=")
            ),
            "",
        )
        if annotation_format == "gff3"
        else ""
    )
    parents = tuple(unquote(p) for p in parent_text.split(",") if p)
    return Feature(
        contig=unquote(contig),
        kind=kind.lower(),
        start=first - 1,
        end=last,
        strand="." if strand == "?" else strand,
        attributes=attributes,
        parents=parents,
    )


def resolve_roots(
    *,
    identifier: str,
    nodes: dict[str, Feature],
    cache: dict[str, tuple[str, ...]],
    visiting: frozenset[str] = frozenset(),
) -> tuple[str, ...]:
    """Resolve gene ancestors, allowing omitted root gene records.

    Args:
        identifier: Feature ID or Parent ID.
        nodes: Annotation features indexed by their exact IDs.
        cache: Memoised ancestor results.
        visiting: Current traversal path for cycle detection.

    Returns:
        Root gene identifiers; a missing root is inferred from its Parent ID.

    Raises:
        ValueError: The hierarchy contains a cycle or exceeds safe depth.
    """
    if identifier in cache:
        return cache[identifier]
    if identifier in visiting or len(visiting) > 100:
        raise ValueError(
            f"Cyclic or excessive annotation hierarchy: {identifier}"
        )
    roots: tuple[str, ...]
    node = nodes.get(identifier)
    if node is None or node.kind.endswith("gene") or not node.parents:
        roots = (identifier,)
    else:
        roots = tuple(
            sorted(
                {
                    root
                    for parent in node.parents
                    for root in resolve_roots(
                        identifier=parent,
                        nodes=nodes,
                        cache=cache,
                        visiting=visiting | {identifier},
                    )
                }
            )
        )
    cache[identifier] = roots
    return roots


def aggregate_gene(*, gene_id: str, features: list[Feature]) -> Gene:
    """Combine all child and gene spans conservatively, including introns.

    Args:
        gene_id: Exact root identifier.
        features: Records belonging to this gene.

    Returns:
        The union span of the gene and its annotated children.

    Raises:
        ValueError: A gene crosses contigs or has contradictory strands.
    """
    contigs = {f.contig for f in features}
    strands = {f.strand for f in features if f.strand != "."}
    if len(contigs) != 1 or len(strands) > 1:
        raise ValueError(f"Conflicting contig/strand for gene {gene_id}")
    return Gene(
        gene_id=gene_id,
        contig=next(iter(contigs)),
        start=min(f.start for f in features),
        end=max(f.end for f in features),
        strand=next(iter(strands), "."),
    )


def read_coordinates(*, path: Path) -> list[Gene]:
    """Read legacy contig/start/end/strand/ID tables with 1-based coordinates.

    Args:
        path: Five-column TSV, with an optional column header.

    Returns:
        Unique genes in input order.

    Raises:
        ValueError: Rows are malformed or identifiers repeat.
    """
    genes: list[Gene] = []
    seen: set[str] = set()
    with open_text(path=path) as stream:
        for number, line in enumerate(stream, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields == ["contig", "start", "end", "strand", "gene_id"]:
                continue
            if len(fields) != 5:
                raise ValueError(f"{path}:{number}: expected five TSV columns")
            contig, start, end, strand, identifier = fields
            identifier = identifier.removeprefix("ID=").split(";", 1)[0]
            if identifier in seen:
                raise ValueError(
                    f"{path}:{number}: duplicate gene {identifier}"
                )
            genes.append(
                Gene(
                    gene_id=identifier,
                    contig=contig,
                    start=int(start) - 1,
                    end=int(end),
                    strand=strand,
                )
            )
            seen.add(identifier)
    if not genes:
        raise ValueError(f"No genes in {path}")
    return genes


def read_annotation(
    *, path: Path, annotation_format: str = "auto"
) -> list[Gene]:
    """Read a complete annotation, irrespective of record order.

    Args:
        path: GFF3, GTF or legacy coordinate TSV, optionally gzip compressed.
        annotation_format: ``auto``, ``gff3``, ``gtf`` or ``tsv``.

    Returns:
        Sorted complete gene spans with exact, unique identifiers.

    Raises:
        ValueError: Annotation content or its parent hierarchy is invalid.
    """
    if annotation_format == "tsv":
        return read_coordinates(path=path)
    if annotation_format not in {"auto", "gff3", "gtf"}:
        raise ValueError(f"Unknown annotation format: {annotation_format}")
    features: list[Feature] = []
    nodes: dict[str, Feature] = {}
    seen_spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    detected = annotation_format
    with open_text(path=path) as stream:
        for number, line in enumerate(stream, start=1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            if detected == "auto":
                columns = line.rstrip().split("\t")
                if len(columns) == 5:
                    return read_coordinates(path=path)
                detected = "gtf" if re.search(r'gene_id\s+"', line) else "gff3"
            try:
                feature = parse_feature(line=line, annotation_format=detected)
            except ValueError as exc:
                raise ValueError(f"{path}:{number}: {exc}") from exc
            identifier = feature.attributes.get("ID")
            if identifier:
                existing = nodes.get(identifier)
                if existing and (
                    existing.kind != feature.kind
                    or existing.contig != feature.contig
                    or existing.strand != feature.strand
                    or existing.parents != feature.parents
                    or any(
                        feature.start < end and feature.end > start
                        for start, end in seen_spans[identifier]
                    )
                ):
                    raise ValueError(
                        f"Conflicting/duplicate feature ID: {identifier}"
                    )
                nodes[identifier] = feature
                seen_spans[identifier].append((feature.start, feature.end))
            if feature.kind in GENIC_FEATURES or feature.kind.endswith("gene"):
                features.append(feature)
    grouped: dict[str, list[Feature]] = defaultdict(list)
    cache: dict[str, tuple[str, ...]] = {}
    for feature in features:
        roots: tuple[str, ...]
        direct = feature.attributes.get("gene_id")
        identifier = feature.attributes.get("ID")
        if direct:
            roots = (direct,)
        elif identifier and feature.kind.endswith("gene"):
            roots = (identifier,)
        elif feature.parents:
            roots = tuple(
                sorted(
                    {
                        root
                        for parent in feature.parents
                        for root in resolve_roots(
                            identifier=parent, nodes=nodes, cache=cache
                        )
                    }
                )
            )
        elif identifier:
            roots = (identifier,)
        else:
            raise ValueError("Genic feature lacks ID, gene_id or Parent")
        for root in roots:
            grouped[root].append(feature)
    if not grouped:
        raise ValueError(f"No recognised genic features in {path}")
    genes = [
        aggregate_gene(gene_id=g, features=fs) for g, fs in grouped.items()
    ]
    LOGGER.info("Parsed %d gene spans from %s", len(genes), path)
    return sorted(genes, key=lambda g: (g.contig, g.start, g.end, g.gene_id))
