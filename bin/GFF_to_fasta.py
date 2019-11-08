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
    scaff = scaff.rstrip()
    return scaff, source, feature, start, stop, score, \
            direction, frame, gene_info


def return_slice_of_interest(seq_record,
                             start,
                             stop,
                             direction,
                             user_defined_genic,
                             info):
    """func to blah blah blah"""
    if direction == "+":
        stop = stop + user_defined_genic
        seq_with_genic = seq_record.seq[start:stop]
        stop = stop + user_defined_genic
    if direction == "-":
        start = start - user_defined_genic
        seq_with_genic = reverse_complement(seq_record.seq[start:stop])
        info = info + " negative strand has been reverse complmented"
    info = info + " start-stop: " + str(start) + " " + str(stop)
    return seq_with_genic, info


def gff_to_fasta(gff, Genome_sequence, min_length,
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
    # this is already indexed
    # Genome_sequence = SeqIO.index(genome_sequence, "fasta")

    with open(gff, "r") as f_handle:
        for line in f_handle:
            line = check_line(line)
            if not line:
                continue
            scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = split_line(line)
            seq_record = Genome_sequence[scaff]
            text = str(seq_record.seq).index("ATGGCTGAAGTCGAACAATTCAAACAAAACCAGAACGT")
            print("gene2 should be 12274, found at", text)
            part = "Scaff: %s | start-stop: %s:%s " % (scaff, start, stop)
            info =  "%s %s coding direction  %d bp genic" %(part,
                                                            direction,
                                                            user_defined_genic)
            seq_with_genic, info = return_slice_of_interest(seq_record,
                                                            start, stop,
                                                            direction,
                                                            user_defined_genic,
                                                            info)

            ROI_len = len(seq_with_genic)
            # the negtive stand has already been reverse complemtned,
            # so we can treat both strand the same here.
            # this is some logic to make sure we dont return tooo much seq
            # only a chunk of size of interest.
            if ROI_len > upstream + user_defined_genic:
               seq_with_genic = seq_with_genic[(ROI_len -
                                                (upstream + user_defined_genic)):]
            print(seq_with_genic, direction)

            record = SeqRecord(Seq(str(seq_with_genic)),
                   id=gene_info, name="",
                   description=info)
            SeqIO.write(record, f_out, "fasta")                                              


usage = """ stand alone script to blindly extract

regions between start and stop"""

parser = OptionParser(usage=usage)

parser.add_option("--gff",
                  dest="coordinate_file",
                  default="format_for_py_script.out",
                  help="gff file")

parser.add_option("-g", "--genome",
                  dest="genome_sequence",
                  default=None,
                  help="genome_sequence.fasta - this has to be the file " +
                  "used to generate the gene models/GFF file")

parser.add_option("-u", "--upstream",
                  dest="upstream",
                  default=False,
                  help="the amount of nucleotide upstream of the gene " +
                  "start, taking into account gene directions, to " +
                  "return in the outfile by default this will not return " +
                  "sequences of min_lenbp or less.")

parser.add_option("-m", "--min_len",
                  dest="min_len",
                  default="3",
                  help="the min length of seq to return. " +
                  "Any fragments less than this are not returned " +
                  "Default = 3")

parser.add_option("-z", "--user_defined_genic",
                  dest="user_defined_genic",
                  default=0,
                  help="the number of nucleotides from within the " +
                  "gene to return, default is 0")

parser.add_option("-o", "--output",
                  dest="out_file",
                  default="upstream_of_genes.fasta",
                  help="Output filename (fasta file)",
                  metavar="FILE")



# get the user options. TODO. arg parser instead
(options, args) = parser.parse_args()
coordinate_file = options.coordinate_file
genome_sequence = options.genome_sequence
upstream = int(options.upstream)
#genes_file = options.genes_file
min_len = int(options.min_len) + 1
user_defined_genic = int(options.user_defined_genic)
description = "YES"
logfile = options.out_file.split(".fa")[0] + ".WARNINGS.log"

outfile = options.out_file



# Run as script
if __name__ == '__main__':
    Genome_sequence = SeqIO.index(genome_sequence, "fasta")
    gff_to_fasta(coordinate_file,
                 Genome_sequence,
                 min_len,
                 outfile,
                 upstream,
                 user_defined_genic,
                 NNN=False)
