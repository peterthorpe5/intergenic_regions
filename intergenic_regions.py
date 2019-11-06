#!/usr/bin/env python
# TITLE: script to get the upstream regions of genes of interest
# script will return reverse complemnet of negative strand
# coded genes upstream regions.
# author: Peter Thorpe September 2018.
# University of St Andrews

"""
script to return intergenic regions between genes.

"""

# imports
from Bio.Seq import Seq
from Bio import SeqIO
import time
import os
import errno
import logging
import logging.handlers
import sys
from optparse import OptionParser
from intergenic_regions_modules.tools import sys_exit,\
      check_gff, index_gene_scaffold_coordinates, \
      assign_vals_to_list, slice_up_scaff, \
      get_len_upstream, get_coordinate_of_interest, \
      write_out_to_file, write_out_gff




##########################

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.5")
    sys.exit(0)


usage = """Use as follows:

$ python get_upstream_regions.py --coordinates
        coordinate_file.fasta -g genome_sequence.fasta
        -upstream <int> number of nucleotides upstream of strat of gene to return
        e.g.  -u 1000

        -m min length of seq to return
        -o outfile_name

Requirements:
python and biopyton.

This will return (--upstream number) of nucleotides to the start of your genes(s) of
interest (-g) gene_file using data from (-c). Gene file can either be space, tab or \n
separated..

The coordinate file can be generated using a GFF3 file and a linux command line using
or using to other python script in this folder. See the readme:

grep "gene" name.gff3 | grep -v "#" | cut -f1,4,5,7,9 > format_for_py_script.out

This also needs to be sorted.


yeilding this reulting file:

scaffold	start	stop	strand(+/-)	ID=gene_number
GROS_00001	2195	3076	-	ID=GROS_g00002
GROS_00001	8583	10515	+	ID=GROS_g00005.....

"""

parser = OptionParser(usage=usage)

parser.add_option("-c", "--coordinates",
                  dest="coordinate_file",
                  default="format_for_py_script.out",
                  help="NOTE: coordinate_file can generate using " +
                  "linux command line of "
                  "GFF file:  grep 'gene' name.gff3 | grep -v '#' | " +
                  " cut -f1,4,5,7,9 > format_for_py_script.out ."
                  "Default = format_for_py_script.out")

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
description = "YES"
logfile = options.out_file.split(".fa")[0] + "WARNINGS.log"

outfile = open(options.out_file, "w")


# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('intergenic_regions.py: %s'
                               % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     logfile)
        sys.exit(1)
    file_list = [coordinate_file, genome_sequence] #genes_file]
    for user_file in file_list:
        if not os.path.isfile(user_file):
           print("file not found: %s" % user_file)
           os._exit(0)
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    # index the genome with biopython
    logger.info("indexing genome")
    Genome_sequence = SeqIO.index(genome_sequence, "fasta")
    Genome_sequence_time=time.time()
    out = 'import genome file took, %.3f' % (Genome_sequence_time
                                             - start_time)
    # populate the dictionaries
    gene_to_next_gene, gene_to_previous_gene, coordinate_dict, \
            gene_list = index_gene_scaffold_coordinates(coordinate_file)

    logger.info(out)
    # open the gff outfile
    name_gff = options.out_file.split(".")[0] + "_upstream_" + str(upstream)+ "bp.gff"
    gff_outfile = open(name_gff, "w")

    for gene in gene_list:
        # now we need to extract the genic regions, not going into other genes
        final_start, final_stop, \
          direction = get_coordinate_of_interest(gene,
                                                 gene_to_next_gene,
                                                 gene_to_previous_gene,
                                                 coordinate_dict)

        # continue # do we want to exclude these?
        gene_coordinates = coordinate_dict[gene]
        # yes calling this again ..
        scaff, start, stop, direction, \
               gene = assign_vals_to_list(gene_coordinates)
        if final_start == "NA":
            out = "upstream for gene %s falls off start of scaff %s" % (gene,
                                                                        scaff)
            logger.warning(out)
            # continue # do we want to exclude these?
        if final_stop == "NA":
            out = "upstream for gene %s falls off END of scaff %s" % (gene,
                                                                      scaff)
            logger.warning(out)
        Genome_seq_record = Genome_sequence[scaff]
        length_of_contig = len(Genome_seq_record.seq)
        # slice up the scaffold.
        # reverse complement for negative
        # this function gets the entire intergenic regions, regardless of
        # the desired length.
        intergenic_region = slice_up_scaff(Genome_seq_record.seq,
                                           final_start,
                                           final_stop,
                                           direction)
        # this function then get the desired upstream chunk
        upstream_ROI = get_len_upstream(intergenic_region,
                                        upstream,
                                        direction)
        if description.upper() == "YES":
            info = "Scaff: %s | start-stop: %s:%s %s coding direction" %(scaff,
                                                                         final_start,
                                                                         final_stop,
                                                                         direction)
        if len(upstream_ROI) >= min_len:
            write_out_to_file(outfile,
                              gene,
                              upstream,
                              upstream_ROI,
                              info)
        else:
            logger.warn

        write_out_gff(gff_outfile, scaff, final_start,
                      final_stop, direction, gene)
    outfile.close()
    gff_outfile.close()
