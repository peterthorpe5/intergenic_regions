# intergenic_regions

READme info 
==============================================

basic usage:
============

``intergenic_regions.py`` -h 

``Requirements:``
python 3.X and biopyton.  


script to get the upstream regions of genes of interest. script will return up to the gene if the full length falls within that gene. 
Also, script will return reverse complemnet of negative strand coded genes.

Usage: Use as follows:

``python intergenic_regions.py --gff file.gff -g genome_sequence.fasta --upstream <int> number of nucleotides upstream of strat of gene to return 
(e.g. -u 1000) -o outfile_name``

This will return (--upstream number) of nucleotides to the start of your genes(s) of interest (-g) gene_file using data from (-c). Gene file can either be space, tab or  separated.


Options:
  -h, --help            show this help message and exit
  --gff, --coordinates=COORDINATE_FILE
                        gff file of the beastie

  -g GENOME_SEQUENCE, --genome=GENOME_SEQUENCE
                        genome_sequence.fasta - this has to be the file used
                        to generate the gene models/GFF file
  -f GENES_FILE, --gene_names=GENES_FILE
                        a file with a list of gene names to get the upstream
                        regions for
  -u UPSTREAM, --upstream=UPSTREAM
                        the amount of nucleotide upstream of the gene start,
                        taking into account gene directions, to return in the
                        outfile by default this will not return sequences of
                        min_lenbp or less
  -z USER_DEFINED_GENIC, --user_defined_genic=USER_DEFINED_GENIC
                        the number of nucleotides from within the gene to
                        return, default is 0
  -m MIN_LEN, --min_len=MIN_LEN
                        the min length of seq to return. Any fragments less
                        than this are not returned Default = 30
  -o FILE, --output=FILE
                        Output filename (fasta file)

						

# note now this script takes in the full gff file and parses it for you

more usage to come. 



