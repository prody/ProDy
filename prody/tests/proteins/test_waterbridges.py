"""This module contains unit tests for :mod:`~prody.proteins`."""

import numpy as np
from numpy.testing import *
try:
    import numpy.testing.decorators as dec
except ModuleNotFoundError:
    from numpy.testing import dec

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

import os

LOGGER.verbosity = 'none'

PATH_PDB = pathDatafile('3enl_addH')
SPLIT_PATH = os.path.split(PATH_PDB)
TARGET_PATH = os.path.join(SPLIT_PATH[0], 'wb_cluster_' + SPLIT_PATH[1])

class TestWaterBridges(unittest.TestCase):

    def setUp(self):
        """Parse the PDB file and run the core calculations."""
        self.atoms = parsePDB(PATH_PDB)

        if prody.PY3K:
            self.waterBridgesChain = calcWaterBridges(self.atoms)
            self.chainMultiWater = self.waterBridgesChain[37]

            self.waterBridgesClust = calcWaterBridges(self.atoms, method='cluster')
            self.clustMultiWater = self.waterBridgesClust[36]

            self.waterBridgesClustInfo = calcWaterBridges(self.atoms, method='cluster', 
                                                        output='info')
            self.clustInfoMultiWater = self.waterBridgesClustInfo[36]

            self.pro0 = self.clustMultiWater.proteins[0]
            self.pro1 = self.clustMultiWater.proteins[1]
            self.w0 = self.clustMultiWater.waters[0]
            self.w1 = self.clustMultiWater.waters[1]
        
    @dec.slow
    def testSingleWaterBridgesChain(self):
        if prody.PY3K:
            self.assertIsInstance(self.waterBridgesChain, list,
                                  'calcWaterBridges chain output should be a list')
            
            self.assertEqual(len(self.waterBridgesChain), 77,
                             'calcWaterBridges chain list should have 77 items')
            
            
            self.assertIsInstance(self.chainMultiWater, prody.proteins.waterbridges.AtomicOutput,
                                  'calcWaterBridges output should contain AtomicOutput')
            
            self.assertEqual(len(self.chainMultiWater.proteins), 2,
                             'water bridges chains should have 2 protein items')
            self.assertIsInstance(self.chainMultiWater.proteins[0], Atom,
                                  'protein items in water bridges chains should be Atom')
            
            self.assertEqual(len(self.chainMultiWater.waters), 2,
                             'water bridges chain should have 2 water items')
            self.assertIsInstance(self.chainMultiWater.waters[0], Atom,
                                  'waters items in water bridges chains should be Atom')
            
    @dec.slow
    def testSingleWaterBridgesCluster(self):
        if prody.PY3K:
            self.assertIsInstance(self.waterBridgesClust, list,
                                  'calcWaterBridges clust output should be a list')
            
            self.assertEqual(len(self.waterBridgesClust), 74,
                             'calcWaterBridges clust list should have 74 items')
            
            self.assertIsInstance(self.clustMultiWater, prody.proteins.waterbridges.AtomicOutput,
                                  'calcWaterBridges output should contain AtomicOutput')
            
            self.assertEqual(len(self.clustMultiWater.proteins), 2,
                             'tested water bridges multi-water cluster should have 3 protein items')
            self.assertIsInstance(self.clustMultiWater.proteins[0], Atom,
                                  'protein items in water bridges clusters should be Atom')
            
            self.assertEqual(len(self.clustMultiWater.waters), 2,
                             'tested water bridges multi-water cluster should have 2 waters items')
            self.assertIsInstance(self.clustMultiWater.waters[0], Atom,
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
            
    @dec.slow
    def testSingleWaterBridgesOutputInfo(self):
        if prody.PY3K:
            self.assertIsInstance(self.waterBridgesClustInfo, list,
                                  'calcWaterBridges clust info output should be a list')
            
            self.assertEqual(len(self.waterBridgesClustInfo), 74,
                             'calcWaterBridges clust info list should have 74 items')
            
            self.assertIsInstance(self.clustInfoMultiWater, list,
                                  'tested water bridges multi-water cluster info output should contain lists')
            
            self.assertEqual(len(self.clustInfoMultiWater), 9,
                             'tested water bridges multi-water cluster info should have 9 items')
            self.assertIsInstance(self.clustInfoMultiWater[0], str,
                                  'info items in water bridges clusters should be str')


            self.assertEqual(self.clustInfoMultiWater[0], self.pro0.getResname() + str(self.pro0.getResnum()),
                             'item 0 from selected info should be resname + resnum of atom 0')
            self.assertEqual(self.clustInfoMultiWater[1], self.pro0.getName() + '_' + str(self.pro0.getIndex()),
                             'item 1 from selected info should be name + _ + index of atom 0')
            self.assertEqual(self.clustInfoMultiWater[2], self.pro0.getChid(),
                             'item 2 from selected info should be chid of atom 0')

            self.assertEqual(self.clustInfoMultiWater[3], self.pro1.getResname() + str(self.pro1.getResnum()),
                             'item 3 from selected info should be resname + resnum of atom 1')
            self.assertEqual(self.clustInfoMultiWater[4], self.pro1.getName() + '_' + str(self.pro1.getIndex()),
                             'item 4 from selected info should be name + _ + index of atom 1')
            self.assertEqual(self.clustInfoMultiWater[5], self.pro1.getChid(),
                             'item 5 from selected info should be chid of atom 1')

            self.assertEqual(self.clustInfoMultiWater[6], calcDistance(self.pro0, self.pro1),
                             'item 6 from selected info should be distance from protein atom 0 to protein atom 1')
            self.assertEqual(self.clustInfoMultiWater[7], len(self.clustMultiWater.waters),
                             'item 7 from selected info should be number of waters in the cluster')
            self.assertEqual(self.clustInfoMultiWater[8], 
                             list(map(lambda w: w.getChid() + "_" + str(w.getIndex()), self.clustMultiWater.waters)),
                             'item 8 from selected info should be chid and index for 2 water atoms')
           
