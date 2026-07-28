"""This module contains unit tests for :mod:`~prody.proteins`."""

from numpy.testing import *
from prody.utilities import importDec
dec = importDec()

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

import os

LOGGER.verbosity = 'none'

PATH_1UBI = pathDatafile('1ubi')
SPLIT_PATH = os.path.split(PATH_1UBI)
EXPECT_PATH = os.path.join(SPLIT_PATH[0], 'addH_' + SPLIT_PATH[1])

class TestFixer(unittest.TestCase):

    def setUp(self):
        """Set PDB file data and parse the PDB file."""

        self.ubi_orig = DATA_FILES['1ubi']
        self.ubi_addH = DATA_FILES['1ubi_addH']

    @dec.slow
    def testPDBFixer(self):

        if prody.PY3K:
            self.filename = addMissingAtoms(PATH_1UBI, method='pdbfixer')
            self.assertEqual(self.filename, EXPECT_PATH,
                            'addMissing atoms wrote the file in the wrong place')
            
            atoms = parsePDB(self.filename)
            self.assertEqual(len(atoms.hydrogen), DATA_FILES['1ubi_addH']['n_h'],
                            'addMissing atoms did not add the right number of hydrogens')
