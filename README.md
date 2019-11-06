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

``python intergenic_regions.py --coordinates coordinate_file.fasta -g genome_sequence.fasta --upstream <int> number of nucleotides upstream of strat of gene to return 
(e.g. -u 1000) -o outfile_name``

This will return (--upstream number) of nucleotides to the start of your genes(s) of interest (-g) gene_file using data from (-c). Gene file can either be space, tab or  separated.


Options:
  -h, --help            show this help message and exit
  -c COORDINATE_FILE, --coordinates=COORDINATE_FILE
                        NOTE: coordinate_file can generate using linux command
                        line of GFF file:  grep 'gene' name.gff3 | grep -v '#'
                        |  cut -f1,4,5,7,9 > format_for_py_script.out .Default
                        = format_for_py_script.out
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

						

The coordinate file using grep or (better) a python script
===========================================================
can be generated using a GFF3 file and a linux command line using:
WARNING: the grep method can return lots of unwanted lines. Use with caution

``grep "gene" name.gff3 | grep -v "#" | cut -f1,4,5,7,9 > format_for_py_script.out``

yeilding this resulting file:
	scaffold        start   stop    strand(+/-)     ID=gene_number
	GROS_00001      2195    3076    -       ID=GROS_g00002
	GROS_00001      8583    10515   +       ID=GROS_g00005.....

The main script will check that the start is always less than the end. GFF file should have starts < stop irrespective of the coding direction. Some GFF files are badly formatted

or use the python script:
``./re_format_GFF_Mcscanx.py`` -h 

Options:
  -h, --help            show this help message and exit
  --gff=FILE            hintsfile
  -m FILE, --mcscan=FILE
                        specific formatting for Mcscan. Default is false,
                        change to True if you need this.
  -s FILE, --species=FILE
                        species prefix to add into column 1 two letters
  -o OUT, --out=OUT     output filenames
  
  The coordinate file must be sorted



yeilding this reulting file:
scaffold	start	stop	strand(+/-)	ID=gene_number
GROS_00001	2195	3076	-	ID=GROS_g00002
GROS_00001	8583	10515	+	ID=GROS_g00005.....




