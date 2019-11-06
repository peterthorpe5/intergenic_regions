from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO

def sys_exit(msg, error_level=1):
   """Print error message to stdout and quit with
    given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)


def check_gff(gff_info, gene_list):
      error = "coord file fields wrong length, should be 5"
      assert len(gff_info.split("\t")) == 5 , "%s" % error
      gene = gff_info.split("\t")[4]
      if ";" in gene:
         gene = gene.split(";")[0]
        # check each gene only represented once
      assert gene not in gene_list, ("duplicate found %s. Reformat -C file." % gene)
      gene_list.append(gene)
      return gene, gene_list


def index_gene_scaffold_coordinates(coordinate_file):
    """function to return dictionary genes and coordinates
    without directions
    gene = scaffold_cordinates
    scaffold_cordinates = gff_info.split("\t")[:]
    coordinate_dict[gene] = scaffold_cordinates
    gene_to_next_gene[gene] = next_gene_in_gff
    gene_to_previous_gene[gene] = previous_gene_in_gff

            """
    coordinate_dict = dict()
    gene_to_next_gene = defaultdict(str)
    gene_to_previous_gene = defaultdict(str)
    data = open(coordinate_file, "r")
    genes_coordinate = data.readlines()
    genes_coordinate = [line.replace("ID=", "").rstrip() \
                        for line in genes_coordinate
                        if line.rstrip() != "" if not line.startswith("#")]
    gene_list = []
    data.close()
    current_gene = ""
    last_gene = "NA"

    for gff_info in genes_coordinate:
        gene, gene_list = check_gff(gff_info, gene_list)
        if gene in coordinate_dict.values():
            print ("repeated line in gff sub file")
            continue
        else:
            scaffold_cordinates = gff_info.split("\t")[:]
            coordinate_dict[gene.rstrip()] = scaffold_cordinates
        gene = gene.rstrip()
        gene_to_previous_gene[gene] = last_gene
        gene_to_next_gene[last_gene] = gene
        last_gene = gene
    return gene_to_next_gene, gene_to_previous_gene, coordinate_dict, gene_list

def assign_vals_to_list(gene_coordinates):
   """fun to assign to list. better ways to do this.. but you know.
   This is easy to read"""
   if gene_coordinates ==  "":
      return "end_of_file", "NA", "NA", "NA", "NA"
   scaff = gene_coordinates[0]
   start = gene_coordinates[1]
   stop = gene_coordinates[2]
   direction = gene_coordinates[3]
   gene = gene_coordinates[4]
   return scaff, start, stop, direction, gene.rstrip()


def get_coordinate_of_interest(gene, gene_to_next_gene,
                               gene_to_previous_gene,
                               coordinate_dict):
   """fucn to obtain the full length region between the start of its
   and the next upstream gene, coding direction aware
   input is a gene name:
   gene_to_next_gene[gene] = gene_after_this_gene,
    gene_to_previous_gene[gene] = gene_before_this_gene,
    coordinate_dict[gene] = pathogens_Gpal_scaffold_1	5	10	+	GPLIN_000000100
   returns
   start stop"""
   next_gene = gene_to_next_gene[gene]
   previous_gene = gene_to_previous_gene[gene]
   gene_coordinates = coordinate_dict[gene]
   scaff, start, stop, direction,\
          gene = assign_vals_to_list(gene_coordinates)
   if direction == "+":
      # our gene is positive strand encoing.
      # we need to stop coord in the previous_gene
      # and the start cood from our gene
      if previous_gene != "NA":
         previous_gene_coordinates = coordinate_dict[previous_gene]
         prevscaff, prevstart, prevstop, prevdirection, \
                    prevgene = assign_vals_to_list(previous_gene_coordinates)
         if prevscaff == scaff:
            return prevstop, start, direction
      return "NA", start, direction
   if direction == "-":
      # our gene is negative strand encoing.
      # we need to start coord in the next
      # and the end cood from our gene
      if next_gene != "NA":
         next_gene_coordinates = coordinate_dict[next_gene]
         nextscaff, nextstart, nextstop, nextdirection, \
                    nextgene = assign_vals_to_list(next_gene_coordinates)
         if nextscaff == scaff:
            return stop, nextstart, direction
         if nextscaff == "end_of_file":
            return stop, "NA" , "+"
      return stop, "NA", direction


def slice_up_scaff(Genome_seq_record_seq,
                   final_start,
                   final_stop,
                   direction):
   """func to slice up the scaff region.
   seq{start:stop]
   if negative:
   reverse complement #
   seq{start:stop]"""
   if final_start == "NA": # this is start of the scaff:
      final_start = 1
   if final_stop == "NA": # this is end of the scaff:
      final_stop = len(Genome_seq_record_seq) -1
   intergenic_region = Genome_seq_record_seq[int(final_start):int(final_stop)]
   if direction == "-":
      intergenic_region = intergenic_region.reverse_complement()
   return intergenic_region, final_start, final_stop


def get_len_upstream(intergenic_region, upstream, direction):
   """function to obtain the len of the intergenic region of interets"""
   ROI_len = len(intergenic_region)
   # the negtive stand has already been reverse complemtned,
   # so we can treat both strand the same here.
   if ROI_len > upstream:
      return intergenic_region[(ROI_len - upstream):]
   else:
      return intergenic_region



def write_out_to_file(outfile, gene, upstream,
                      upstream_ROI, info):
   """func to write out to file"""
   data = ">%s\t%dbp upstream\t%s\n%s\n" % (gene, upstream,
                                            info, upstream_ROI)
   outfile.write(data)



def write_out_gff(outfile, scaff, final_start,
                  final_stop, direction, gene,
                  upstream):
   """func to write out the coordinate to a gff file. This is so
   another script can be used to return the genic region if required.
   I didnt want to alter this script to do so as I would have to
   rewrite testsscripts"""
   if int(final_start) > int(final_stop):
      start = final_stop
      final_stop = final_start
      final_start = start
   outfmt  = "\t".join([scaff,
                        "intergenic_regions",
                        "intergenic",
                        final_start,
                        final_stop,
                        "%sbp_upstream_of_gene" % str(upstream),
                        direction,
                        ".",
                        gene + "\n"])
   outfile.write(outfmt)

