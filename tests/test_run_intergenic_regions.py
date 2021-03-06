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
GENOME = os.path.join("tests", "inputs", "genome.fasta")
COORD = os.path.join("tests", "inputs",
                     "tests_gene_indexing.txt")


def setup_outdir():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)

setup_outdir()

def test_run_u_5bp():
    """complete script runs 5 bp upstream"""
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
                                     "test_5bp.fasta")])
    print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def test_complete_run_u_50bp():
    """complete script runs 50 bp upstream"""
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
                                     "test_50bp.fasta")])
    print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def test_complete_run_u_500bp():
    """complete script runs 500 bp upstream"""
    command = " ".join(["python",
                        "intergenic_regions.py",
                        "-g",
                        GENOME,
                        "-u",
                        "500",
                        "-c",
                        COORD,
                        "-o",
                        os.path.join(OUTDIR,
                                     "test_50bp.fasta")])
    print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def test_complete_script_notravis_u_50bp_minlen_29():
    """complete script runs 50 bp upstream min len 29"""
    command = " ".join(["python",
                        "intergenic_regions.py",
                        "-g",
                        GENOME,
                        "-u",
                        "50",
                        "-c",
                        COORD,
                        "-m",
                        "29",
                        "-o",
                        os.path.join(OUTDIR,
                                     "test_50bp_min_len29.fasta")])
    print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def test_complete_script_notravis_u_50bp_minlen_10():
    """complete script runs 50 bp upstream min len 10"""
    command = " ".join(["python",
                        "intergenic_regions.py",
                        "-g",
                        GENOME,
                        "-u",
                        "50",
                        "-c",
                        COORD,
                        "-m",
                        "10",
                        "-o",
                        os.path.join(OUTDIR,
                                     "test_50bp_min_len10.fasta")])
    print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)  

#test_complete_script_notravis()

##def test_logger_creation():
##    """intergenic log creation"""
##    logger = intergenic.construct_logger(NAMESPACE, header=False)



