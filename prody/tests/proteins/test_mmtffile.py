"""This module contains unit tests for :mod:`~prody.proteins`."""

import os

import numpy as np
from numpy.testing import *


from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

class TestParseMMTF(unittest.TestCase):

    def setUp(self):
        """Set PDB file data and parse the PDB file."""

        self.ubi_mmtf = DATA_FILES['1ubi_mmtf']
        self.ubi_pdb = DATA_FILES['1ubi']
        self.nmr_mmtf = DATA_FILES['2k39_mmtf']
        self.biomol_pdb = DATA_FILES['3enl_pdb']
        self.biomol_mmtf = DATA_FILES['3enl_mmtf']
        

    def testUsualCase(self):
        """Test the outcome of a simple parsing scenario."""

        #simple file
        ag = parseDatafile(self.ubi_pdb['file'])
        mag = parseDatafile(self.ubi_mmtf['file'])

        self.assertIsInstance(mag, prody.AtomGroup,
            'parseMMTF failed to return an AtomGroup instance')

        self.assertEqual(mag.numAtoms(), self.ubi_mmtf['n_atoms'],
            'parseMMTF failed to parse correct number of atoms')

        self.assertEqual(mag.numCoordsets(), self.ubi_mmtf['models'],
            'parseMMTF failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(prody.calcRMSD(ag,mag),0,'parseMMTF has non-zero RMSD to parsePDB')
        
        # bioassembly creation
        ag = parseDatafile(self.biomol_pdb['file'], biomol=True)
        mag = parseDatafile(self.biomol_mmtf['file'], biomol=True)

        self.assertIsInstance(mag, prody.AtomGroup,
            'parseMMTF failed to return an AtomGroup instance')

        self.assertEqual(mag.numAtoms(), self.biomol_mmtf['n_atoms'],
            'parseMMTF failed to parse correct number of atoms')

        self.assertEqual(mag.numCoordsets(), self.biomol_mmtf['models'],
            'parseMMTF failed to parse correct number of coordinate sets '
            '(models)')        
        self.assertEqual(prody.calcRMSD(ag,mag),0,'parseMMTF has non-zero RMSD to parsePDB')        
        
        
        # multi-model file (currently fails)
        mag = parseDatafile(self.nmr_mmtf['file'])

        self.assertIsInstance(mag, prody.AtomGroup,
            'parseMMTF failed to return an AtomGroup instance')

        self.assertEqual(mag.numAtoms(), self.nmr_mmtf['n_atoms'],
            'parseMMTF failed to parse correct number of atoms')

        self.assertEqual(mag.numCoordsets(), self.nmr_mmtf['models'],
            'parseMMTF failed to parse correct number of coordinate sets '
            '(models)')        

 
