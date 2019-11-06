#!/usr/bin/env python3
# title: GFF to fasta
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe. The James Hutton Institute, Uk
# imports
import sys
import os
from optparse import OptionParser
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def index_genome_file(genome):
    """index the genome file
    return a seq_record object that can be accessed
    in a dictionary like manner"""
    genome_database = SeqIO.index(genome, "fasta")
    return genome_database


def check_line(line):
    """checks if line is ok
    Starts with a comment of blank ine = not ok. """
    if line.startswith("#"):
        return False
    if not line.strip():
        return False
    return line


def split_line(line):
    """split the gff line
    Takes in a gff line. returns the elements, as str or int"""
    warning_list = ["gene", "exon", "intron"]
    assert len(line.split("\t")) == 9 , "GFF fields wrong length should be 9"
    scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = line.split("\t")
    gene_info = gene_info.rstrip()
    start = int(start)
    stop = int(stop)
    return scaff, source, feature, start, stop, score, \
            direction, frame, gene_info


def gff_to_fasta(gff, genome_database, min_length,
                 outfile, upstream, user_defined_genic,
                 NNN=False):
    """take in gff file. Gets the seq defined by the gff coords.
    If negative direction coding, the reverse complement is generated.
    A min length of seq to return and max len is applied to remove seq
    less than, for example 3 which cant be real and less that e.g., 25k
    which will be flase positives and not informative in downstream analysis
    """
    min_length = int(min_length)
    f_out = open(outfile, "w")
    with open(gff, "r") as f_handle:
        for line in f_handle:
            line = check_line(line)
            if not line:
                continue
            scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = split_line(line)
            seq_record = genome_database[scaff]
            if direction == "+":
                seq_with_genic = seq_record.seq[start:(stop + user_defined_genic)]
            if direction == "-":
                seq_with_genic = reverse_complement(seq_record.seq
                                                    [(start - user_defined_genic):stop])
            info = "Scaff: %s | start-stop: %s:%s %s coding direction  %d bp genic" %(scaff,
                                                                         final_start,
                                                                         final_stop,
                                                                         direction,
                                                                        user_defined_genic)
            
            record = SeqRecord(Seq(seq_with_genic),
                   id=gene_info, name="",
                   description=info)
            SeqIO.write(record, f_out, "fasta")
                                                

