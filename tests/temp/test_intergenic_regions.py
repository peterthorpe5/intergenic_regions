#!/usr/bin/env python3

"""Tests of intergenic_regions.py script."""

import os
import shutil

from argparse import Namespace
from nose.tools import nottest, assert_equal, with_setup

from intergenic_regions import intergenic_regions

NAMESPACE = Namespace(upstream='20',
                      genome='input/genome.fasta')
DATADIR = os.path.join("tests", "inputs")
OUTDIR = os.path.join("tests", "test_out")
TARGETDIR = os.path.join("tests", "test_targets")


def setup_outdir():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_complete_script_notravis():
    """complete script runs"""
    if os.path.isdir(NAMESPACE.outdirname):
        shutil.rmtree(NAMESPACE.outdirname)
    santi_pipeline.run_pycits_main(NAMESPACE)


def test_logger_creation():
    """intergenic log creation"""
    logger = intergenic.construct_logger(NAMESPACE, header=False)



