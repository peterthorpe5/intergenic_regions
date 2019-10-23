#!/usr/bin/env python

"""Tests of intergenic_regions.py script."""

import os
import shutil
import subprocess

from argparse import Namespace
from nose.tools import nottest, assert_equal, with_setup


DATADIR = os.path.join("tests", "inputs")
OUTDIR = os.path.join("tests", "test_out")
TARGETDIR = os.path.join("tests", "test_targets")
GENOME = os.path.join("tests", "inputs", "genome_simplified.fasta")
COORD = os.path.join("tests", "inputs",
                     "tests_gene_indexing_simplified.txt")


def setup_outdir():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)

setup_outdir()
def test_result_u_5bp():
    """result 5 bp upstream"""
    command = " ".join(["python",
                        "intergenic_regions.py",
                        "-g",
                        GENOME,
                        "-u",
                        "5",
                        "-c",
                        COORD,
                        "-o",
                        os.path.join(OUTDIR,
                                     "test2_5bp.fasta")])
    print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def test_result_u_50bp():
    """result 50 bp upstream"""
    command = " ".join(["python",
                        "intergenic_regions.py",
                        "-g",
                        GENOME,
                        "-u",
                        "50",
                        "-c",
                        COORD,
                        "-o",
                        os.path.join(OUTDIR,
                                     "test2_50bp.fasta")])
    print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

#test_complete_script_notravis()

##def test_logger_creation():
##    """intergenic log creation"""
##    logger = intergenic.construct_logger(NAMESPACE, header=False)



