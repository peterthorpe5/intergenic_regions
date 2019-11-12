#!/usr/bin/env python
# Title:
# script to reformt GFF.
# author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

import os
import sys
from collections import OrderedDict

# function

def split_gene_name(gene_info):
    """fucntion to remove all the crap around the gene
    name and return just the name
    e.g. ID+Gene01:go123,pfam,qifjf ...
    return
    Gene01
    """
    gene_info = gene_info.replace("ID=", "").split()[0]
    gene_info = gene_info.split(";")[0]
    gene_info = gene_info.replace("CDS:", "")
    gene_info = gene_info.split("Note=")[0]
    gene_info = gene_info.split(".")[0]
    if "gene_id " in gene_info:
        # this is to deal with gtf files
        gene_info = gene_info.split("gene_id ")[1]
        gene_info = gene_info.replace('"', '')
        gene_info = gene_info.replace(';', '')
    return gene_info.rstrip()


def parse_gff(line):
    """func to return the gff line as elements"""
    assert len(line.split("\t")) ==9 ,"GFF... wrong len should be 9"
    scaf, source, feature, start, stop, score, \
        direction, frame, gene_info = line.split("\t")
    gene = split_gene_name(gene_info)
    return scaf, feature, start, stop, direction, gene.rstrip()


def write_out_gff_data(gene_start_stop_dict, gene_start, gene_stop, out):
    """ func to iterate through the dict and wrtie it out"""
    f_out = open(out, "w")
    for gene, gff_info in gene_start_stop_dict.items():
        scaf, start, stop, direction, gene = gff_info.split()
        actual_start = gene_start[gene]
        actual_stop = gene_stop[gene]
        outdata = "\t".join([scaf,
                             actual_start,
                             actual_stop,
                             direction, gene])

        f_out.write(outdata + "\n")
    f_out.close()


def get_starts_stops(gff, out):
    """function to obtain the true
    starts and stops"""
    f_in = open(gff, "r")
    gene_start_stop_dict = OrderedDict() #  ordered to input
    gene_start = OrderedDict()
    gene_stop = OrderedDict()

    stop = 0 #  temp set this value, will get replaced.
    gene = "not a real gene"  #  temp set this value, will get replaced.
    for line in f_in:
        if line.startswith("#"):
            continue
        last_stop = stop #  the previous iteration will go in this variable
        last_gene = gene
        # call the func to parse line
        scaf, feature, start, stop, direction, gene = parse_gff(line)

        if feature.upper() != "CDS":
            continue # the cds start stop seem to take into account
                     # the UTR region, so it better to use.
        if gene not in gene_start:
            gene_start[gene] = start
            gene_stop[gene] = stop
            outdata = "\t".join([scaf,
                                 "start", "end",
                                 direction, gene])
            gene_start_stop_dict[last_gene] = outdata
        current_stop = gene_stop[gene]
        if int(stop) > int(current_stop):
            gene_stop[gene] = stop

    write_out_gff_data(gene_start_stop_dict, gene_start, gene_stop, out)
    f_in.close()
    

#################################################################################################


#run the program
# Run as script
if __name__ == '__main__':
    if not os.path.isfile(gff):
        print("sorry cannot find you %s file" % gff)
        os._exit(0)
    # print ("warning: GFF files vary! make sure you check the file before using")
    get_starts_stops(gff, out)

