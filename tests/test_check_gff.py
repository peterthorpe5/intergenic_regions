#!/usr/bin/env python
# author: Peter Thorpe October 2019.
# U of St Andrew
# script to test the indexing of the coordinates file.
# pip install nose!!!

import os
from nose.tools import nottest, assert_equal, assert_not_equal
import unittest
from intergenic_regions_modules.tools import check_gff

# INPUT DATA LOCATION
INPUT_ok = "pathogens_Gpal_scaffold_1	5	10	+	GPLIN_000000100"
INPUT_too_many_coloumns = "pathogens_Gpal_scaffold_1	5	10	+	GPLIN_000000100	+"
INPUT_gene_name = "pathogens_Gpal_scaffold_1	5	10	+	GPLIN_000000100:GO0000123"


def test_check_gff():
    """test_check_gff return correct names"""
    gene_list = []
    gene, gene_list= check_gff(INPUT_ok, gene_list)
    # print(gene, gene_list)
    assert_equal(gene, "GPLIN_000000100")
    assert_equal(gene_list, ["GPLIN_000000100"])
test_check_gff()
    
##def test_gene_indexing_too_many_coloumns():
##    """test_gene_indexing_too_many_coloumns"""
##    gene_list = []
##    gene, gene_list = check_gff(INPUT_too_many_coloumns, gene_list)
##    assert_not_equal(gene, "GPLIN_000000100")
##    assert_not_equal(gene_list, ["GPLIN_000000100"])


def test_check_alter_names():
    """test_check_alter_names to correct names"""
    gene_list = []
    gene, gene_list= check_gff(INPUT_ok, gene_list)
    # print(gene, gene_list)
    assert_equal(gene, "GPLIN_000000100")
    assert_equal(gene_list, ["GPLIN_000000100"])



