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
        self.altloc_pdb = DATA_FILES['1pwc_pdb']
        self.altloc_mmtf = DATA_FILES['1pwc_mmtf']
        

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
        
        # altlocs
        ag = parseDatafile(self.altloc_pdb['file'],bonds=True)
        mag = parseDatafile(self.altloc_mmtf['file'],bonds=True)

        self.assertIsInstance(mag, prody.AtomGroup,
            'parseMMTF failed to return an AtomGroup instance')

        self.assertEqual(mag.numAtoms(), self.altloc_mmtf['n_atoms'],
            'parseMMTF failed to parse correct number of atoms')

        self.assertEqual(mag.numCoordsets(), self.altloc_mmtf['models'],
            'parseMMTF failed to parse correct number of coordinate sets '
            '(models)')        
        self.assertEqual(prody.calcRMSD(ag,mag),0,'parseMMTF has non-zero RMSD to parsePDB')       
        self.assertEqual(mag.numBonds(), 2, "parseMMTF has wrong number of bonds")

        ag = parseDatafile(self.altloc_pdb['file'],altloc='all',bonds=True)
        mag = parseDatafile(self.altloc_mmtf['file'],altloc='all',bonds=True)
        
        self.assertEqual(mag.numAtoms(), 3191, "parseMMTF failed to parse correct number of atoms with altloc='any'")
        # mmtf does not interleave altlocs, so can't do all atom rmsd (they don't match up)
        self.assertEqual(prody.calcRMSD(ag.name_C,mag.name_C),0,'parseMMTF has non-zero RMSD to parsePDB')       
        self.assertEqual(mag.numBonds(), 2, "parseMMTF has wrong number of bonds")
        
        ag = parseDatafile(self.altloc_pdb['file'],altloc='B',bonds=True)
        mag = parseDatafile(self.altloc_mmtf['file'],altloc='B',bonds=True)
        
        #the whole ligand is altloc A and is covalently bound
        self.assertEqual(mag.numAtoms(), 3109, "parseMMTF failed to parse correct number of atoms with altloc='any'")
        self.assertEqual(prody.calcRMSD(ag.name_C,mag.name_C),0,'parseMMTF has non-zero RMSD to parsePDB')       
        self.assertEqual(mag.numBonds(), 1, "parseMMTF has wrong number of bonds with altloc='B'")        
        
        # multi-model file (currently fails)
        mag = parseDatafile(self.nmr_mmtf['file'])

        self.assertIsInstance(mag, prody.AtomGroup,
            'parseMMTF failed to return an AtomGroup instance')

        self.assertEqual(mag.numAtoms(), self.nmr_mmtf['n_atoms'],
            'parseMMTF failed to parse correct number of atoms')

        self.assertEqual(mag.numCoordsets(), self.nmr_mmtf['models'],
            'parseMMTF failed to parse correct number of coordinate sets '
            '(models)')        

 
