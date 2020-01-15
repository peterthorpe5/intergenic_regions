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
  
  --gff=GFF_FILE        this is the full gff file  greo for gene coord does
                        not work This will get the starts and stops based on
                        the CDS field.
                        
  -g GENOME_SEQUENCE, --genome=GENOME_SEQUENCE
                        genome_sequence.fasta - this has to be the file used
                        to generate the gene models/GFF file
                        
  -u UPSTREAM, --upstream=UPSTREAM
                        the amount of nucleotide upstream of the gene start,
                        taking into account gene directions, to return in the
                        outfile by default this will not return sequences of
                        min_lenbp or less.
                        
  -m MIN_LEN, --min_len=MIN_LEN
                        the min length of seq to return. Any fragments less
                        than this are not returned Default = 3
                        
  -z USER_DEFINED_GENIC, --user_defined_genic=USER_DEFINED_GENIC
                        the number of nucleotides from within the gene to
                        return, default is 0
                        
  -o FILE, --output=FILE
                        Output filename (fasta file)

						
if -z is called another ouput file is generated. 

# note now this script takes in the full gff file and parses it for you

This is a re write on genomeic_upstream_regions: https://github.com/peterthorpe5/public_scripts/tree/master/genomic_upstream_regions 

Why? A bug was found .. better to re write with tests. The old tool was old complicated code, needed testing.  

``unit tests``
look in the test folder.
run with nose tests. (pip install nose)





