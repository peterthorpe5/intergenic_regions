"""Unsorted hierarchy, identifier, UTR and malformed annotation tests."""

import gzip

import pytest

from intergenic_regions.annotation import (
    Feature,
    aggregate_gene,
    parse_attributes,
    parse_feature,
    read_annotation,
    read_coordinates,
    resolve_roots,
)


def test_attribute_formats():
    assert parse_attributes(
        text="ID=gene.1;Name=a%3Bb;", annotation_format="gff3"
    ) == {"ID": "gene.1", "Name": "a;b"}
    assert (
        parse_attributes(
            text='gene_id "g.1"; transcript_id "t1";', annotation_format="gtf"
        )["gene_id"]
        == "g.1"
    )
    assert parse_attributes(text=".", annotation_format="gff3") == {}


def test_discontinuous_gene_and_escaped_contig(tmp_path):
    path = tmp_path / "split.gff3"
    path.write_text(
        "x%3A1\ts\tgene\t1\t5\t.\t+\t.\tID=g\n"
        "x%3A1\ts\tgene\t10\t20\t.\t+\t.\tID=g\n"
    )
    genes = read_annotation(path=path)
    assert len(genes) == 1
    assert (genes[0].contig, genes[0].start, genes[0].end) == ("x:1", 0, 20)


@pytest.mark.parametrize(
    "row",
    [
        "x\ts\tgene\t1\t5\tNaN\t+\t.\tID=g",
        "x\ts\tgene\t1\t5\tbad\t+\t.\tID=g",
        "x\ts\tCDS\t1\t5\t.\t+\t4\tID=g",
        "\ts\tgene\t1\t5\t.\t+\t.\tID=g",
        "x\ts\t\t1\t5\t.\t+\t.\tID=g",
    ],
)
def test_invalid_score_phase_and_identifiers(row):
    with pytest.raises(ValueError):
        parse_feature(line=row, annotation_format="gff3")


@pytest.mark.parametrize(
    "text,format_name",
    [
        ("ID=g;ID=h", "gff3"),
        ("=g", "gff3"),
        ("bad", "gff3"),
        ('gene_id "g"; bad', "gtf"),
        ("ID=g", "bad"),
    ],
)
def test_bad_attributes(text, format_name):
    with pytest.raises(ValueError):
        parse_attributes(text=text, annotation_format=format_name)


def test_feature_coordinates_and_escaped_parent():
    feature = parse_feature(
        line="x\ts\tCDS\t1\t9\t.\t?\t0\tParent=g%2Cx,t2",
        annotation_format="gff3",
    )
    assert (feature.start, feature.end, feature.strand) == (0, 9, ".")
    assert feature.parents == ("g,x", "t2")


@pytest.mark.parametrize(
    "row",
    [
        "wrong",
        "x\ts\tgene\t0\t9\t.\t+\t.\tID=g",
        "x\ts\tgene\t3\t2\t.\t+\t.\tID=g",
        "x\ts\tgene\t1\t9\t.\tx\t.\tID=g",
    ],
)
def test_feature_failures(row):
    with pytest.raises(ValueError):
        parse_feature(line=row, annotation_format="gff3")


def test_hierarchy_and_cycle():
    feature = Feature(
        contig="x",
        kind="mrna",
        start=1,
        end=4,
        strand="+",
        attributes={"ID": "t"},
        parents=("g",),
    )
    cache = {}
    assert resolve_roots(
        identifier="t", nodes={"t": feature}, cache=cache
    ) == ("g",)
    assert resolve_roots(identifier="t", nodes={}, cache=cache) == ("g",)
    cyclic = Feature(
        contig="x",
        kind="mrna",
        start=1,
        end=4,
        strand="+",
        attributes={},
        parents=("t",),
    )
    with pytest.raises(ValueError):
        resolve_roots(identifier="t", nodes={"t": cyclic}, cache={})
    with pytest.raises(ValueError):
        resolve_roots(
            identifier="t",
            nodes={},
            cache={},
            visiting=frozenset(str(i) for i in range(102)),
        )


def test_unsorted_gff_preserves_ids_and_utrs(tmp_path):
    path = tmp_path / "annotation.gff3"
    path.write_text(
        "##gff-version 3\n"
        "x\ts\tCDS\t20\t40\t.\t+\t0\tID=cds;Parent=t.1\n"
        "x\ts\texon\t10\t50\t.\t+\t.\tParent=t.1\n"
        "x\ts\tmRNA\t10\t50\t.\t+\t.\tID=t.1;Parent=g.1\n"
        "x\ts\tgene\t10\t50\t.\t+\t.\tID=g.1\n##FASTA\n>x\nAAAA\n"
    )
    genes = read_annotation(path=path)
    assert [(g.gene_id, g.start, g.end) for g in genes] == [("g.1", 9, 50)]


def test_gtf_without_gene_records_and_multi_parent(tmp_path):
    path = tmp_path / "genes.gtf.gz"
    with gzip.open(filename=path, mode="wt") as stream:
        stream.write(
            'x\ts\tCDS\t20\t30\t.\t-\t0\tgene_id "g.1"; '
            'transcript_id "t.1";\n'
            'x\ts\texon\t10\t40\t.\t-\t.\tgene_id "g.1"; '
            'transcript_id "t.2";\n'
        )
    gene = read_annotation(path=path)[0]
    assert (gene.gene_id, gene.start, gene.end, gene.strand) == (
        "g.1",
        9,
        40,
        "-",
    )
    path = tmp_path / "multi.gff"
    path.write_text("x\ts\texon\t10\t20\t.\t+\t.\tParent=g1,g2\n")
    assert {g.gene_id for g in read_annotation(path=path)} == {"g1", "g2"}


def test_coordinate_table_and_no_sort_requirement(tmp_path):
    path = tmp_path / "coords.tsv"
    path.write_text(
        "#header\ncontig\tstart\tend\tstrand\tgene_id\nx\t20\t30\t-\tID=g.2;Name=hello\nx\t1\t5\t+\tg.1\n"
    )
    assert read_coordinates(path=path)[0].gene_id == "g.2"
    assert len(read_annotation(path=path, annotation_format="tsv")) == 2
    assert len(read_annotation(path=path)) == 2


@pytest.mark.parametrize(
    "text,format_name",
    [
        ("", "auto"),
        ("x\ts\tregion\t1\t5\t.\t+\t.\tID=r", "auto"),
        ("x\ts\texon\t1\t5\t.\t+\t.\t.", "gff3"),
        (
            "x\ts\tgene\t1\t5\t.\t+\t.\tID=g\nx\ts\tgene\t1\t5\t.\t+\t.\tID=g",
            "gff3",
        ),
        (
            "x\ts\tCDS\t1\t5\t.\t+\t.\tID=c;Parent=t\ny\ts\tCDS\t1\t5\t.\t+\t.\tID=c;Parent=t",
            "gff3",
        ),
        ("bad", "auto"),
        ("", "bad"),
    ],
)
def test_annotation_failures(tmp_path, text, format_name):
    path = tmp_path / "bad"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_annotation(path=path, annotation_format=format_name)


@pytest.mark.parametrize("text", ["", "bad", "x\t1\t5\t+\tg\nx\t2\t6\t+\tg\n"])
def test_coordinate_failures(tmp_path, text):
    path = tmp_path / "bad.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_coordinates(path=path)


def test_gene_aggregation_conflicts():
    features = [
        Feature(
            contig="x",
            kind="exon",
            start=1,
            end=4,
            strand=".",
            attributes={},
            parents=(),
        ),
        Feature(
            contig="x",
            kind="cds",
            start=2,
            end=5,
            strand="+",
            attributes={},
            parents=(),
        ),
    ]
    gene = aggregate_gene(gene_id="g", features=features)
    assert (gene.start, gene.end, gene.strand) == (1, 5, "+")
    with pytest.raises(ValueError):
        aggregate_gene(
            gene_id="g",
            features=[
                features[1],
                Feature(
                    contig="y",
                    kind="cds",
                    start=2,
                    end=5,
                    strand="-",
                    attributes={},
                    parents=(),
                ),
            ],
        )
