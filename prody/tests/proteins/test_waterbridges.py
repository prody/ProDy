"""This module contains unit tests for :mod:`~prody.proteins`."""

import numpy as np
from numpy.testing import *
try:
    import numpy.testing.decorators as dec
except ImportError:
    from numpy.testing import dec

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

import os

LOGGER.verbosity = 'none'

PATH_1UBI = pathDatafile('1ubi')
SPLIT_PATH = os.path.split(PATH_1UBI)
TARGET_PATH = os.path.join(SPLIT_PATH[0], 'wb_cluster_' + SPLIT_PATH[1])

class TestWaterBridges(unittest.TestCase):

    def setUp(self):
        """Parse the PDB file and run the core calculations."""
        self.atoms = parseDatafile('1ubi_addH')
        self.waterBridgesChain = calcWaterBridges(self.atoms)
        self.waterBridgesClust = calcWaterBridges(self.atoms, method='cluster')

    @dec.slow
    def testSingleWaterBridgesChain(self):
        if prody.PY3K:
            self.assertIsInstance(self.waterBridgesChain, list,
                                  'calcWaterBridges output should be a list')
            
            self.assertEqual(len(self.waterBridgesChain), 7,
                             'calcWaterBridges chain list should have 7 items')
            
            chain0 = self.waterBridgesChain[1]
            self.assertIsInstance(chain0, prody.proteins.waterbridges.AtomicOutput,
                                  'calcWaterBridges output should contain AtomicOutput')
            
            self.assertEqual(len(chain0.proteins), 2,
                             'water bridges chains should have 2 protein items')
            self.assertIsInstance(chain0.proteins[0], Atom,
                                  'protein items in water bridges chains should be Atom')
            
            self.assertEqual(len(chain0.waters), 1,
                             'water bridges chain should have 1 waters item')
            self.assertIsInstance(chain0.waters[0], Atom,
                                  'waters items in water bridges chains should be Atom')
            
    @dec.slow
    def testSingleWaterBridgesCluster(self):
        if prody.PY3K:
            self.assertIsInstance(self.waterBridgesClust, list,
                                  'calcWaterBridges output should be a list')
            
            self.assertEqual(len(self.waterBridgesClust), 7,
                             'calcWaterBridges clust list should have 7 items')
            
            clust1 = self.waterBridgesClust[1]
            self.assertIsInstance(clust1, prody.proteins.waterbridges.AtomicOutput,
                                  'calcWaterBridges output should contain AtomicOutput')
            
            self.assertEqual(len(clust1.proteins), 3,
                             'some water bridges clusters should have 3 protein items')
            self.assertIsInstance(clust1.proteins[0], Atom,
                                  'protein items in water bridges clusters should be Atom')
            
            self.assertEqual(len(clust1.waters), 1,
                             'water bridges cluster should have 1 waters item')
            self.assertIsInstance(clust1.waters[0], Atom,
                                  'waters items in water bridges clusters should be Atom')

    @dec.slow
    def testSaveWaterBridges(self):
        if prody.PY3K:
            savePDBWaterBridges(self.waterBridgesClust, self.atoms, TARGET_PATH)
            self.saved = parsePDB(TARGET_PATH)

            self.assertFalse(np.allclose(self.atoms.protein.getBetas(), self.saved.protein.getBetas()),
                             'savePDBWaterBridges did not change protein b-factors')
            
            self.assertFalse(np.allclose(self.atoms.protein.getOccupancies(), self.saved.protein.getOccupancies()),
                             'savePDBWaterBridges did not change protein occupancies')
            
            self.assertFalse(np.equal(self.atoms.water.numAtoms(), self.saved.water.numAtoms()),
                             'savePDBWaterBridges did not change the number of water atoms')
            
