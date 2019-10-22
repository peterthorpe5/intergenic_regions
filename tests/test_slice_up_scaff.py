#!/usr/bin/env python
# author: Peter Thorpe October 2019.
# U of St Andrew
# script to test the indexing of the coordinates file.
# pip install nose!!!

import os
from nose.tools import nottest, assert_equal
import unittest
from intergenic_regions_modules.tools import slice_up_scaff
from Bio.Seq import Seq
from Bio import SeqIO

# INPUT DATA LOCATION
INPUT = os.path.join("tests", "inputs", "genome.fasta")

Genome_sequence = SeqIO.index(INPUT, "fasta")
Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_1"]


def test_scaffold_seq():
    "test scof is correct"
    assert_equal(Genome_seq_record.seq,
                 "TTTTTATGATGTTTTTTTTTTATGATGATGAAAAAAAAAAATGATGATGATTTTTCCGTGGCTGCGTGCCGTCTCCTTGCGGGGTGGTGCGTGGACTAGTGGTGAG")

def test_slice_up_scaff():
    """run the function slice_up_scaff"""
    intergenic_region = slice_up_scaff(Genome_seq_record.seq,
                                       "1",
                                       "5",
                                       "+")
    assert_equal(intergenic_region, "TTTT")


def test_slice_up_scaff():
    """run the function slice_up_scaff reverse_complment"""
    intergenic_region = slice_up_scaff(Genome_seq_record.seq,
                                       "1",
                                       "5",
                                       "-")
    assert_equal(intergenic_region, "AAAA")


def test_slice_up_scaff():
    """run the function slice_up_scaff NA start"""
    intergenic_region = slice_up_scaff(Genome_seq_record.seq,
                                       "NA",
                                       "5",
                                       "+")
    assert_equal(intergenic_region, "TTTT")


def test_slice_up_scaff():
    """run the function slice_up_scaff NA stop"""
    intergenic_region = slice_up_scaff(Genome_seq_record.seq,
                                       "100",
                                       "NA",
                                       "+")
    # note that the 100 is actually 101, and the end is len(seq) -1
    assert_equal(intergenic_region, "GGTGA")


def test_slice_up_scaff():
    """run the function slice_up_scaff NA stop revrse complement"""
    intergenic_region = slice_up_scaff(Genome_seq_record.seq,
                                       "100",
                                       "NA",
                                       "-")
    # note that the 100 is actually 101, and the end is len(seq) -1
    assert_equal(intergenic_region, "TCACC")

