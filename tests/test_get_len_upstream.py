#!/usr/bin/env python
# author: Peter Thorpe October 2019.
# U of St Andrew
# script to test the indexing of the coordinates file.
# pip install nose!!!

import os
from nose.tools import nottest, assert_equal
import unittest
from intergenic_regions_modules.tools import get_len_upstream
from Bio.Seq import Seq
from Bio import SeqIO

# INPUT DATA LOCATION
INPUT = os.path.join("tests", "inputs", "genome.fasta")

Genome_sequence = SeqIO.index(INPUT, "fasta")
Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_1"]


def test_intergenic_region():
    "test intergenic_region 5 "
    intergenic_region = get_len_upstream(Genome_seq_record.seq,
                                         5)
    assert_equal(intergenic_region,
                 "TTTTT")

def test_intergenic_region():
    "test intergenic_region 4 "
    intergenic_region = get_len_upstream(Genome_seq_record.seq,
                                         4)
    assert_equal(intergenic_region,
                 "TTTT")

def test_intergenic_region():
    "test intergenic_region 3 "
    intergenic_region = get_len_upstream(Genome_seq_record.seq,
                                         3)
    assert_equal(intergenic_region,
                 "TTT")

