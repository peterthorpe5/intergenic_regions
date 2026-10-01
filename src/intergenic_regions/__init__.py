"""Extract strictly intergenic regions and analyse sequence signatures."""

from intergenic_regions._version import __version__
from intergenic_regions.annotation import read_annotation
from intergenic_regions.extraction import GeneIndex, extract_regions
from intergenic_regions.genome import Genome, reverse_complement
from intergenic_regions.models import Gene, Region

__all__ = [
    "Gene",
    "GeneIndex",
    "Genome",
    "Region",
    "__version__",
    "extract_regions",
    "read_annotation",
    "reverse_complement",
]
