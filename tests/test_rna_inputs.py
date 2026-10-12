import gzip
from pathlib import Path

import pysam
import pytest

from intergenic_regions.models import Gene
from intergenic_regions.rna_inputs import (
    CoverageReader,
    RnaSample,
    alignment_transcript_strand,
    open_rna_coverage,
    read_bedgraph,
    read_gene_exons,
    read_rna_samples,
)


def sample(**kwargs):
    return RnaSample(
        **{
            "sample_id": "a",
            "condition": "treated",
            "replicate_id": "1",
            "strandedness": "unstranded",
            **kwargs,
        }
    )


@pytest.mark.parametrize(
    "changes",
    [
        {"sample_id": "../bad"},
        {"condition": "a b"},
        {"replicate_id": ""},
        {"strandedness": "guess"},
        {"layout": "guess"},
        {"size_factor": 0},
        {"size_factor": True},
        {"size_factor": float("nan")},
        {"read2": Path("x")},
        {"bam": Path("x")},
        {"read1": Path("x")},
        {"bam": Path("x"), "layout": "single", "coverage": Path("y")},
        {"coverage_plus": Path("x")},
        {"coverage": Path("x"), "strandedness": "forward"},
        {"coverage_plus": Path("x"), "coverage_minus": Path("y")},
        {"coverage": Path("x"), "coverage_plus": Path("y")},
        {"coverage": Path("x"), "layout": "single"},
        {"read1": Path("x"), "layout": "paired"},
    ],
)
def test_sample_validation(changes):
    with pytest.raises(ValueError):
        sample(**changes)


def test_manifest_relative_paths_and_covariates(tmp_path):
    track = tmp_path / "depth.bg"
    track.write_text("chr1\t10\t20\t3\n")
    path = tmp_path / "samples.tsv"
    path.write_text(
        "sample_id\tcondition\treplicate_id\tstrandedness\tcoverage\t"
        "size_factor\tbatch\n"
        "r1\tref\t1\tunstranded\tdepth.bg\t2\tA\n"
    )
    result = read_rna_samples(path=path)
    assert result[0].coverage == track.resolve()
    assert result[0].covariates == {"batch": "A"}
    assert result[0].size_factor == 2
    with open_rna_coverage(
        sample=result[0], lengths={"chr1": 30}, genome_path=tmp_path / "g.fa"
    ) as reader:
        assert reader.depth(contig="chr1", start=8, end=22, strand="-") == (
            [0] * 2 + [3] * 10 + [0] * 2
        )
        with pytest.raises(ValueError):
            reader.depth(contig="chr1", start=20, end=31, strand="+")


@pytest.mark.parametrize(
    "body,exception",
    [
        ("", ValueError),
        ("sample_id\tcondition\n", ValueError),
        ("sample_id\tcondition\treplicate_id\tstrandedness\n", ValueError),
        (
            "sample_id\tcondition\treplicate_id\tstrandedness\n"
            "a\tc\t1\tunstranded\textra\n",
            ValueError,
        ),
        (
            "sample_id\tcondition\treplicate_id\tstrandedness\na\tc\t1\n",
            ValueError,
        ),
        (
            "sample_id\tcondition\treplicate_id\tstrandedness\n"
            "a\tc\t1\tunstranded\n",
            ValueError,
        ),
        (
            "sample_id\tcondition\treplicate_id\tstrandedness\tcoverage\n"
            "a\tc\t1\tunstranded\tabsent.bg\n",
            FileNotFoundError,
        ),
    ],
)
def test_manifest_rejects_invalid_inputs(tmp_path, body, exception):
    path = tmp_path / "samples.tsv"
    path.write_text(body)
    with pytest.raises(exception):
        read_rna_samples(path=path)


def test_metadata_only_and_technical_replicates(tmp_path):
    path = tmp_path / "s.tsv"
    header = "sample_id\tcondition\treplicate_id\tstrandedness\n"
    path.write_text(header + "a\tc\t1\tunstranded\n")
    assert len(read_rna_samples(path=path, require_source=False)) == 1
    path.write_text(header + "a\tc\t1\tunstranded\nb\tc\t1\tunstranded\n")
    with pytest.raises(ValueError, match="biological observations"):
        read_rna_samples(path=path, require_source=False)


def test_unordered_multiparent_exons_and_gtf(tmp_path):
    path = tmp_path / "genes.gff3"
    path.write_text(
        "chr1\tx\texon\t21\t40\t.\t+\t.\tParent=t1,t2\n"
        "chr1\tx\texon\t31\t50\t.\t+\t.\tParent=t1\n"
        "chr1\tx\tmRNA\t11\t80\t.\t+\t.\tID=t2;Parent=g%31\n"
        "chr1\tx\tmRNA\t11\t80\t.\t+\t.\tID=t1;Parent=g%31\n"
        "chr1\tx\tgene\t11\t80\t.\t+\t.\tID=g%31\n"
        "##FASTA\n>chr1\nAAAA\n"
    )
    assert read_gene_exons(path=path) == {"g1": [(20, 50)]}
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tx\texon\t21\t40\t.\t-\t.\tgene_id "g"; transcript_id "t";\n'
    )
    assert read_gene_exons(path=gtf) == {"g": [(20, 40)]}
    tsv = tmp_path / "genes.tsv"
    tsv.write_text("contig\tstart\tend\tstrand\tgene_id\nchr1\t1\t50\t+\tg\n")
    assert read_gene_exons(path=tsv) == {}


@pytest.mark.parametrize(
    "exon",
    [
        "chr1\tx\texon\t1\t9\t.\t+\t.\tParent=g\n",
        "chr1\tx\texon\t21\t30\t.\t-\t.\tParent=g\n",
        "chr2\tx\texon\t21\t30\t.\t+\t.\tParent=g\n",
        "chr1\tx\texon\t21\t30\t.\t+\t.\tID=e\n",
    ],
)
def test_exon_ancestry_conflicts(tmp_path, exon):
    path = tmp_path / "g.gff3"
    path.write_text("chr1\tx\tgene\t11\t80\t.\t+\t.\tID=g\n" + exon)
    genes = [Gene(gene_id="g", contig="chr1", start=10, end=80, strand="+")]
    with pytest.raises(ValueError):
        read_gene_exons(path=path, genes=genes)


@pytest.mark.parametrize(
    "body",
    [
        "chr1\t1\t5\n",
        "chr1\t-1\t5\t1\n",
        "chr1\t1\t1\t1\n",
        "chr2\t1\t5\t1\n",
        "chr1\t1\t101\t1\n",
        "chr1\t1\t5\t-1\n",
        "chr1\t1\t5\tnan\n",
        "chr1\t1\t5\tinf\n",
        "chr1\t1\t5\t1\nchr1\t4\t8\t1\n",
    ],
)
def test_bedgraph_invalid_spans(tmp_path, body):
    path = tmp_path / "x.bg"
    path.write_text(body)
    with pytest.raises(ValueError):
        read_bedgraph(path=path, lengths={"chr1": 100})


def test_strand_tracks_are_already_transcript_oriented(tmp_path):
    plus, minus = tmp_path / "plus.bg", tmp_path / "minus.bg.gz"
    plus.write_text("track type=bedGraph\nchr1\t2\t6\t3\n")
    with gzip.open(filename=minus, mode="wt") as stream:
        stream.write("browser position chr1\nchr1\t4\t8\t7\n")
    reader = CoverageReader(
        sample=sample(
            strandedness="reverse", coverage_plus=plus, coverage_minus=minus
        ),
        lengths={"chr1": 10},
        genome_path=tmp_path / "g.fa",
    )
    assert reader.depth(contig="chr1", start=0, end=10, strand="+") == (
        [0] * 2 + [3] * 4 + [0] * 4
    )
    assert reader.depth(contig="chr1", start=0, end=10, strand="-") == (
        [0] * 4 + [7] * 4 + [0] * 2
    )
    reader.close()
    reader.close()
    for settings in ({"min_mapq": -1}, {"max_reads": 0}):
        with pytest.raises(ValueError):
            CoverageReader(
                sample=sample(),
                lengths={"chr1": 10},
                genome_path=tmp_path / "g.fa",
                **settings,
            )
    with pytest.raises(ValueError, match="Map raw reads"):
        CoverageReader(
            sample=sample(read1=Path("reads"), layout="single"),
            lengths={"chr1": 10},
            genome_path=tmp_path / "g.fa",
        )


@pytest.mark.parametrize(
    "flag,mode,expected",
    [
        (0, "forward", "+"),
        (16, "forward", "-"),
        (0, "reverse", "-"),
        (16, "reverse", "+"),
        (99, "forward", "+"),
        (147, "forward", "+"),
        (83, "reverse", "+"),
        (163, "reverse", "+"),
        (0, "unstranded", "."),
    ],
)
def test_alignment_strand(flag, mode, expected):
    read = pysam.AlignedSegment()
    read.flag = flag
    assert (
        alignment_transcript_strand(read=read, strandedness=mode) == expected
    )


def test_alignment_strand_invalid_flags():
    read = pysam.AlignedSegment()
    read.flag = 1
    with pytest.raises(ValueError):
        alignment_transcript_strand(read=read, strandedness="forward")
    with pytest.raises(ValueError):
        alignment_transcript_strand(read=read, strandedness="guessed")


def write_bam(
    *, path, records, length=200, indexed=True, sorting="coordinate"
):
    header = {
        "HD": {"VN": "1.6", "SO": sorting},
        "SQ": [{"SN": "chr1", "LN": length}],
    }
    with pysam.AlignmentFile(
        filename=str(path), mode="wb", header=header
    ) as f:
        for row in sorted(records, key=lambda item: item["start"]):
            read = pysam.AlignedSegment()
            read.query_name = row["name"]
            read.flag = row.get("flag", 0)
            read.reference_id = 0
            read.reference_start = row["start"]
            read.mapping_quality = row.get("mapq", 60)
            read.cigarstring = row.get("cigar", "10M")
            size = sum(
                n for op, n in read.cigartuples if op in {0, 1, 4, 7, 8}
            )
            read.query_sequence = "A" * size
            read.query_qualities = pysam.qualitystring_to_array("I" * size)
            for tag, value in row.get("tags", [("NH", 1)]):
                read.set_tag(tag=tag, value=value)
            f.write(read)
    if indexed:
        pysam.index(str(path))
    return path


def test_spliced_bam_filters_deletions_and_read_groups(tmp_path):
    path = write_bam(
        path=tmp_path / "x.bam",
        records=[
            {"name": "good", "start": 10, "cigar": "10M30N10M"},
            {"name": "delete", "start": 80, "cigar": "5M5D5M"},
            {"name": "opposite", "start": 10, "flag": 16},
            {"name": "low", "start": 10, "mapq": 1},
            {"name": "multi", "start": 10, "tags": [("NH", 2)]},
            {"name": "dup", "start": 10, "flag": 1024},
            {"name": "sec", "start": 10, "flag": 256},
            {"name": "supp", "start": 10, "flag": 2048},
            {"name": "fail", "start": 10, "flag": 512},
            {"name": "rg", "start": 110, "tags": [("RG", "a")]},
            {"name": "rg", "start": 110, "tags": [("RG", "b")]},
        ],
    )
    mapped = sample(bam=path, layout="single", strandedness="forward")
    with open_rna_coverage(
        sample=mapped, lengths={"chr1": 200}, genome_path=tmp_path / "g.fa"
    ) as reader:
        depth = reader.depth(contig="chr1", start=0, end=130, strand="+")
        assert depth[10:20] == [1] * 10
        assert depth[20:50] == [0] * 30
        assert depth[50:60] == [1] * 10
        assert depth[80:95] == [1] * 5 + [0] * 5 + [1] * 5
        assert depth[110:120] == [2] * 10
    assert reader.handle is None
    with pytest.raises(ValueError, match="closed"):
        reader.depth(contig="chr1", start=0, end=10, strand="+")
    with open_rna_coverage(
        sample=mapped,
        lengths={"chr1": 200},
        genome_path=tmp_path / "g.fa",
        exclude_duplicates=False,
        max_reads=2,
    ) as reader:
        with pytest.raises(ValueError, match="max_reads"):
            reader.depth(contig="chr1", start=0, end=100, strand="+")


def test_paired_mate_overlap_counted_once_and_layout_checked(tmp_path):
    path = write_bam(
        path=tmp_path / "x.bam",
        records=[
            {"name": "pair", "start": 10, "flag": 99, "cigar": "20M"},
            {"name": "pair", "start": 20, "flag": 147, "cigar": "20M"},
        ],
    )
    with open_rna_coverage(
        sample=sample(bam=path, layout="paired", strandedness="forward"),
        lengths={"chr1": 200},
        genome_path=tmp_path / "g.fa",
    ) as reader:
        assert reader.depth(contig="chr1", start=10, end=40, strand="+") == (
            [1] * 30
        )
    with open_rna_coverage(
        sample=sample(bam=path, layout="single"),
        lengths={"chr1": 200},
        genome_path=tmp_path / "g.fa",
    ) as reader:
        with pytest.raises(ValueError, match="layout"):
            reader.depth(contig="chr1", start=10, end=40, strand="+")


@pytest.mark.parametrize(
    "indexed,sorting,lengths",
    [
        (False, "coordinate", {"chr1": 200}),
        (True, "unknown", {"chr1": 200}),
        (True, "coordinate", {"chr1": 201}),
    ],
)
def test_bam_reference_and_index_requirements(
    tmp_path, indexed, sorting, lengths
):
    path = write_bam(
        path=tmp_path / "x.bam", records=[], indexed=indexed, sorting=sorting
    )
    with pytest.raises(ValueError):
        CoverageReader(
            sample=sample(bam=path, layout="single"),
            lengths=lengths,
            genome_path=tmp_path / "g.fa",
        )
