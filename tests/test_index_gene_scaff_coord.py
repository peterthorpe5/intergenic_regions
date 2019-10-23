#!/usr/bin/env python
# author: Peter Thorpe October 2019.
# U of St Andrew
# script to test the indexing of the coordinates file.
# pip install nose!!!

import os
from nose.tools import nottest, assert_equal
import unittest
from intergenic_regions_modules.tools import index_gene_scaffold_coordinates

# INPUT DATA LOCATION
INPUT = os.path.join("tests", "inputs", "tests_gene_indexing.txt")
gene_to_next_gene_test = {'NA': 'GPLIN_000000100', 'GPLIN_000000100': 'GPLIN_000000200', 'GPLIN_000000200': 'GPLIN_000000300', 'GPLIN_000000300': 'GPLIN_000000400', 'GPLIN_000000400': 'GPLIN_000000500', 'GPLIN_000000500': 'GPLIN_000020100', 'GPLIN_000020100': 'GPLIN_000020200', 'GPLIN_000020200': 'GPLIN_000030300', 'GPLIN_000030300': 'GPLIN_000030400', 'GPLIN_000030400': 'GPLIN_000030500'}
gene_to_previous_gene_test = {'GPLIN_000000100': 'NA', 'GPLIN_000000200': 'GPLIN_000000100', 'GPLIN_000000300': 'GPLIN_000000200', 'GPLIN_000000400': 'GPLIN_000000300', 'GPLIN_000000500': 'GPLIN_000000400', 'GPLIN_000020100': 'GPLIN_000000500', 'GPLIN_000020200': 'GPLIN_000020100', 'GPLIN_000030300': 'GPLIN_000020200', 'GPLIN_000030400': 'GPLIN_000030300', 'GPLIN_000030500': 'GPLIN_000030400'}
coordinate_dict_test = {'GPLIN_000000100': ['pathogens_Gpal_scaffold_1', '5', '10', '+', 'GPLIN_000000100'], 'GPLIN_000000200': ['pathogens_Gpal_scaffold_1', '20', '30', '+', 'GPLIN_000000200'], 'GPLIN_000000300': ['pathogens_Gpal_scaffold_1', '40', '50', '-', 'GPLIN_000000300'], 'GPLIN_000000400': ['pathogens_Gpal_scaffold_1', '55', '70', '-', 'GPLIN_000000400'], 'GPLIN_000000500': ['pathogens_Gpal_scaffold_1', '100', '110', '+', 'GPLIN_000000500'], 'GPLIN_000020100': ['pathogens_Gpal_scaffold_2', '5', '10', '+', 'GPLIN_000020100'], 'GPLIN_000020200': ['pathogens_Gpal_scaffold_2', '20', '30', '+', 'GPLIN_000020200'], 'GPLIN_000030300': ['pathogens_Gpal_scaffold_3', '40', '50', '-', 'GPLIN_000030300'], 'GPLIN_000030400': ['pathogens_Gpal_scaffold_3', '55', '70', '-', 'GPLIN_000030400'], 'GPLIN_000030500': ['pathogens_Gpal_scaffold_3', '100', '110', '+', 'GPLIN_000030500']}


def test_gene_to_next_gene_indexing():
    """run the function gene_to_next_
    gene index_gene_scaffold_coordinates
    to check the gene indexing is working"""
    gene_to_next_gene, gene_to_previous_gene, \
    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
    assert_equal(gene_to_next_gene, gene_to_next_gene_test)


def test_next_is_as_expected1():
    """run the function gene_to_next_
    gene index_gene_scaffold_coordinates
    to check the gene indexing is working"""
    gene_to_next_gene, gene_to_previous_gene, \
    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
    result = gene_to_next_gene["GPLIN_000000500"]
    # Alhough this is a differenct scaff's gene. Is it true ..
    assert_equal(result, "GPLIN_000020100")


##def test_next_is_as_expected2():
##    """run the function gene_to_next_
##    gene index_gene_scaffold_coordinates
##    to check the gene indexing is working"""
##    gene_to_next_gene, gene_to_previous_gene, \
##    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
##    result = gene_to_next_gene["GPLIN_000030500"]
##    # Alhough this is a differenct scaff's gene. Is it true ..
##    assert_equal(result, "NA")


def test_next_gene_is_as_expected_end_of_scaff():
    """run the function gene_to_next_
    gene index_gene_scaffold_coordinates
    to check the gene indexing is working"""
    gene_to_next_gene, gene_to_previous_gene, \
    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
    result = gene_to_next_gene["GPLIN_000000100"]
    assert_equal(result, "GPLIN_000000200")


def test_gene_to_prev_gene_indexing():
    """run the function gene_to_previous_gene
    index_gene_scaffold_coordinates
    to check the gene indexing is working"""
    gene_to_next_gene, gene_to_previous_gene, \
    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
    assert_equal(gene_to_previous_gene, gene_to_previous_gene_test)


def test_prev_gene2_is_as_expected():
    """run the function gene_to_next_
    gene index_gene_scaffold_coordinates
    to check the gene indexing is working"""
    gene_to_next_gene, gene_to_previous_gene, \
    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
    result = gene_to_previous_gene["GPLIN_000000200"]
    assert_equal(result, "GPLIN_000000100")


def test_previous3_gene_is_as_expected_scaf_start():
    """run the function gene_to_next_
    gene index_gene_scaffold_coordinates
    to check the gene indexing is working"""
    gene_to_next_gene, gene_to_previous_gene, \
    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
    result = gene_to_previous_gene["GPLIN_000000100"]
    assert_equal(result, "NA")

def test_coordinate_dict_indexing():
    """run the function coordinate_dict
    index_gene_scaffold_coordinates
    to check the gene indexing is working"""
    gene_to_next_gene, gene_to_previous_gene, \
    coordinate_dict, gene_list = index_gene_scaffold_coordinates(INPUT)
    assert_equal(coordinate_dict, coordinate_dict_test)

