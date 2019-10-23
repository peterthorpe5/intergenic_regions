#!/usr/bin/env python
# author: Peter Thorpe October 2019.
# U of St Andrew
# script to test the indexing of the coordinates file.
# pip install nose!!!

import os
from nose.tools import nottest, assert_equal, assert_not_equal
import unittest
from intergenic_regions_modules.tools import get_coordinate_of_interest

# INPUT DATA LOCATION
INPUT_ok = "pathogens_Gpal_scaffold_1	5	10	+	GPLIN_000000100"
INPUT_too_many_coloumns = "pathogens_Gpal_scaffold_1	5	10	+	GPLIN_000000100	+"
INPUT_gene_name = "pathogens_Gpal_scaffold_1	5	10	+	GPLIN_000000100:GO0000123"

gene_to_next_gene = {'NA': 'GPLIN_000000100', 'GPLIN_000000100': 'GPLIN_000000200', 'GPLIN_000000200': 'GPLIN_000000300', 'GPLIN_000000300': 'GPLIN_000000400', 'GPLIN_000000400': 'GPLIN_000000500', 'GPLIN_000000500': 'GPLIN_000020100', 'GPLIN_000020100': 'GPLIN_000020200', 'GPLIN_000020200': 'GPLIN_000030300', 'GPLIN_000030300': 'GPLIN_000030400', 'GPLIN_000030400': 'GPLIN_000030500'}
gene_to_previous_gene = {'GPLIN_000000100': 'NA', 'GPLIN_000000200': 'GPLIN_000000100', 'GPLIN_000000300': 'GPLIN_000000200', 'GPLIN_000000400': 'GPLIN_000000300', 'GPLIN_000000500': 'GPLIN_000000400', 'GPLIN_000020100': 'GPLIN_000000500', 'GPLIN_000020200': 'GPLIN_000020100', 'GPLIN_000030300': 'GPLIN_000020200', 'GPLIN_000030400': 'GPLIN_000030300', 'GPLIN_000030500': 'GPLIN_000030400'}
coordinate_dict = {'GPLIN_000000100': ['pathogens_Gpal_scaffold_1', '5', '10', '+', 'GPLIN_000000100'], 'GPLIN_000000200': ['pathogens_Gpal_scaffold_1', '20', '30', '+', 'GPLIN_000000200'], 'GPLIN_000000300': ['pathogens_Gpal_scaffold_1', '40', '50', '-', 'GPLIN_000000300'], 'GPLIN_000000400': ['pathogens_Gpal_scaffold_1', '55', '70', '-', 'GPLIN_000000400'], 'GPLIN_000000500': ['pathogens_Gpal_scaffold_1', '100', '110', '+', 'GPLIN_000000500'], 'GPLIN_000020100': ['pathogens_Gpal_scaffold_2', '5', '10', '+', 'GPLIN_000020100'], 'GPLIN_000020200': ['pathogens_Gpal_scaffold_2', '20', '30', '+', 'GPLIN_000020200'], 'GPLIN_000030300': ['pathogens_Gpal_scaffold_3', '40', '50', '-', 'GPLIN_000030300'], 'GPLIN_000030400': ['pathogens_Gpal_scaffold_3', '55', '70', '-', 'GPLIN_000030400'], 'GPLIN_000030500': ['pathogens_Gpal_scaffold_3', '100', '110', '+', 'GPLIN_000030500']}
gene = 'GPLIN_000000200'


def test_get_coordinate_of_interest_positive():
    """test_get_coordinate_of_interest + coding"""
    # gene = 'GPLIN_000000200' + coding

    start, stop, direction = get_coordinate_of_interest(gene,
                                             gene_to_next_gene,
                                             gene_to_previous_gene,
                                             coordinate_dict)
    assert_equal(start, "10")
    assert_equal(stop, "20")
    assert_equal(direction, "+")


def test_get_coordinate_of_interest_negative():
    """test_get_coordinate_of_interest start of scaffold"""
    gene = 'GPLIN_000000100' #- coding
    start, stop, direction = get_coordinate_of_interest(gene,
                                             gene_to_next_gene,
                                             gene_to_previous_gene,
                                             coordinate_dict)
    assert_equal(start, "NA")
    assert_equal(stop, "5")
    assert_equal(direction, "+")


def test_get_coordinate_of_interest_negative2():
    """test_get_coordin_of_interest start of scaffold"""
    gene = 'GPLIN_000000300' #- coding
    start, stop, direction = get_coordinate_of_interest(gene,
                                             gene_to_next_gene,
                                             gene_to_previous_gene,
                                             coordinate_dict)
    assert_equal(start, "50")
    assert_equal(stop, "55")
    assert_equal(direction, "-")



def test_get_coordinate_of_interest_different_scaffold():
    """test_get_coordin_of_interest start of scaffold"""
    gene = 'GPLIN_000020100' #+ coding
    start, stop, direction = get_coordinate_of_interest(gene,
                                             gene_to_next_gene,
                                             gene_to_previous_gene,
                                             coordinate_dict)
    assert_equal(start, "NA")
    assert_equal(stop, "5")
    assert_equal(direction, "+")





