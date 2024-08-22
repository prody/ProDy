"""This module contains unit tests for :mod:`~prody.proteins.interactions`."""

import numpy as np
from numpy import arange
from prody import *
from prody.tests import unittest
from prody.tests.datafiles import *
from numpy.testing import assert_equal

from prody.proteins.interactions import calcHydrogenBondsTrajectory, calcSaltBridgesTrajectory
from prody.proteins.interactions import calcRepulsiveIonicBondingTrajectory, calcPiStackingTrajectory
from prody.proteins.interactions import calcPiCationTrajectory, calcHydrophobicTrajectory
from prody.proteins.interactions import calcDisulfideBondsTrajectory, calcProteinInteractions

import sys

class TestInteractions(unittest.TestCase):

    def setUp(self):
        """Generating new data to compare it with the existing one"""
        
        if prody.PY3K:
            self.ATOMS = parseDatafile('2k39_insty')
            self.ALL_INTERACTIONS = parseDatafile('2k39_all')
            self.ALL_INTERACTIONS2 = parseDatafile('2k39_all2')
            self.HBS_INTERACTIONS = parseDatafile('2k39_hbs')
            self.SBS_INTERACTIONS = parseDatafile('2k39_sbs')
            self.RIB_INTERACTIONS = parseDatafile('2k39_rib')
            self.PISTACK_INTERACTIONS = parseDatafile('2k39_PiStack')
            self.PICAT_INTERACTIONS = parseDatafile('2k39_PiCat')
            self.HPH_INTERACTIONS = parseDatafile('2k39_hph')
            self.HPH_INTERACTIONS2 = parseDatafile('2k39_hph2')
            self.DISU_INTERACTIONS = parseDatafile('2k39_disu')
            
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS)
            np.save('test_2k39_all.npy', self.data_all, allow_pickle=True)

            self.data_hbs = calcHydrogenBondsTrajectory(self.ATOMS)
            np.save('test_2k39_hbs.npy', np.array(self.data_hbs, dtype=object), allow_pickle=True)

            self.data_sbs = calcSaltBridgesTrajectory(self.ATOMS)
            np.save('test_2k39_sbs.npy', np.array(self.data_sbs, dtype=object), allow_pickle=True)

            self.data_rib = calcRepulsiveIonicBondingTrajectory(self.ATOMS)
            np.save('test_2k39_rib.npy', np.array(self.data_rib, dtype=object), allow_pickle=True)

            self.data_PiStack = calcPiStackingTrajectory(self.ATOMS)
            np.save('test_2k39_PiStack.npy', np.array(self.data_PiStack, dtype=object), allow_pickle=True)

            self.data_PiCat = calcPiCationTrajectory(self.ATOMS)
            np.save('test_2k39_PiCat.npy', np.array(self.data_PiCat, dtype=object), allow_pickle=True)

            self.data_hph = calcHydrophobicTrajectory(self.ATOMS)
            np.save('test_2k39_hph.npy', np.array(self.data_hph, dtype=object), allow_pickle=True)

            self.data_disu = calcDisulfideBondsTrajectory(self.ATOMS)
            np.save('test_2k39_disu.npy', np.array(self.data_disu, dtype=object), allow_pickle=True)

    def testAllInteractions(self):
        """Test for all types of interactions."""

        if prody.PY3K:        
            data_test = np.load('test_2k39_all.npy', allow_pickle=True)

            try:
                assert_equal(data_test, self.ALL_INTERACTIONS2,
                         'failed to get correct interactions without hpb.so')
            except:
                assert_equal(data_test, self.ALL_INTERACTIONS,
                         'failed to get correct interactions with hpb.so')
    
    def testHydrogenBonds(self):
        """Test for hydrogen bonds.
        Last column is compared becasue pairs of residues can be reversed and
        order can be also different in the interactions"""

        if prody.PY3K:                
            data_test = np.load('test_2k39_hbs.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HBS_INTERACTIONS]),
                         'failed to get correct hydrogen bonds')        
                     
    def testSaltBridges(self):
        """Test for salt bridges."""

        if prody.PY3K:                
            data_test = np.load('test_2k39_sbs.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.SBS_INTERACTIONS]),
                         'failed to get correct salt bridges')                             

    def testRepulsiveIonicBonding(self):
        """Test for repulsive ionic bonding."""

        if prody.PY3K:                
            data_test = np.load('test_2k39_rib.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.RIB_INTERACTIONS if i]),
                         'failed to get correct repulsive ionic bonding')                             

    def testPiStacking(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = np.load('test_2k39_PiStack.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.PISTACK_INTERACTIONS if i]),
                         'failed to get correct pi-stacking interactions')                             
                     
    def testPiCation(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = np.load('test_2k39_PiCat.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.PICAT_INTERACTIONS if i]),
                         'failed to get correct pi-cation interactions')


    def testHydrophobicInteractions(self):
        """Test for hydrophobic interactions."""

        if prody.PY3K:        
            data_test = np.load('test_2k39_hph.npy', allow_pickle=True)                                                        
            try:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS2]),
                         'failed to get correct hydrophobic interactions without hpb.so')
            except:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS]),
                         'failed to get correct hydrophobic interactions with hpb.so')
        

    def testDisulfideBonds(self):
        """Test for disulfide bonds interactions."""

        if prody.PY3K:               
             data_test = np.load('test_2k39_disu.npy', allow_pickle=True)
             assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if i]),
                          'failed to get correct disulfide bonds')
        
    def testImportHpb(self):

        try:
            import prody.proteins.hpb as hpb
            imported_hpb = True
        except ImportError:
            try:
                import hpb
                imported_hpb = True
            except ImportError:
                imported_hpb = False
            
        if sys.version_info[1] < 11:
            self.assertTrue(imported_hpb)
        else:
            self.assertFalse(imported_hpb)

    @classmethod
    def tearDownClass(cls):
        if prody.PY3K:
            import os
            for filename in ['test_2k39_all.npy', 'test_2k39_hbs.npy', 'test_2k39_sbs.npy', 
                            'test_2k39_rib.npy', 'test_2k39_PiStack.npy', 'test_2k39_PiCat.npy',
                            'test_2k39_hph.npy', 'test_2k39_disu.npy']:
                os.remove(filename)
