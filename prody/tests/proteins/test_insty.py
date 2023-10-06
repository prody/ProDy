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

ATOMS = parseDatafile('2k39_insty')
ALL_INTERACTIONS = parseDatafile('2k39_all')
ALL_INTERACTIONS2 = parseDatafile('2k39_all2')
HBS_INTERACTIONS = parseDatafile('2k39_hbs')
SBS_INTERACTIONS = parseDatafile('2k39_sbs')
RIB_INTERACTIONS = parseDatafile('2k39_rib')
PISTACK_INTERACTIONS = parseDatafile('2k39_PiStack')
PICAT_INTERACTIONS = parseDatafile('2k39_PiCat')
HPH_INTERACTIONS = parseDatafile('2k39_hph')
HPH_INTERACTIONS2 = parseDatafile('2k39_hph2')
DISU_INTERACTIONS = parseDatafile('2k39_disu')

# Generating new data to compare it with the existing one:
INTERACTIONS_ALL = InteractionsTrajectory()
data_all = INTERACTIONS_ALL.calcProteinInteractionsTrajectory(ATOMS)
np.save('test_2k39_all.npy', data_all, allow_pickle=True)

data_hbs = calcHydrogenBondsTrajectory(ATOMS)
np.save('test_2k39_hbs.npy', data_hbs, allow_pickle=True)

data_sbs = calcSaltBridgesTrajectory(ATOMS)
np.save('test_2k39_sbs.npy', data_sbs, allow_pickle=True)

data_rib = calcRepulsiveIonicBondingTrajectory(ATOMS)
np.save('test_2k39_rib.npy', data_rib, allow_pickle=True)

data_PiStack = calcPiStackingTrajectory(ATOMS)
np.save('test_2k39_PiStack.npy', data_PiStack, allow_pickle=True)

data_PiCat = calcPiCationTrajectory(ATOMS)
np.save('test_2k39_PiCat.npy', data_PiCat, allow_pickle=True)

data_hph = calcHydrophobicTrajectory(ATOMS)
np.save('test_2k39_hph.npy', data_hph, allow_pickle=True)

data_disu = calcDisulfideBondsTrajectory(ATOMS)
np.save('test_2k39_disu.npy', data_disu, allow_pickle=True)

class TestInteractions(unittest.TestCase):
    
    def testAllInsteractions(self):
        """Test for all types of interactions."""
        
        data_test = np.load('test_2k39_all.npy', allow_pickle=True)
        
        try:
            assert_equal(data_test, ALL_INTERACTIONS2,
                     'failed to get correct interactions without hpb.so')
        except:
            assert_equal(data_test, ALL_INTERACTIONS,
                     'failed to get correct interactions with hpb.so')                
    
    def testHydrogenBonds(self):
        """Test for hydrogen bonds.
        Last column is compared becasue pairs of residues can be reversed and
        order can be also different in the interactions"""
        
        data_test = np.load('test_2k39_hbs.npy', allow_pickle=True)
        assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in HBS_INTERACTIONS]),
                     'failed to get correct hydrogen bonds')        
                     
    def testSaltBridges(self):
        """Test for salt bridges."""
        
        data_test = np.load('test_2k39_sbs.npy', allow_pickle=True)
        assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in SBS_INTERACTIONS]),
                     'failed to get correct salt bridges')                             

    def testRepulsiveIonicBonding(self):
        """Test for repulsive ionic bonding."""
        
        data_test = np.load('test_2k39_rib.npy', allow_pickle=True)
        assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in RIB_INTERACTIONS if i]),
                     'failed to get correct repulsive ionic bonding')                             

    def testPiStacking(self):
        """Test for pi-stacking interactions."""
        
        data_test = np.load('test_2k39_PiStack.npy', allow_pickle=True)
        assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in PISTACK_INTERACTIONS if i]),
                     'failed to get correct pi-stacking interactions')                             
                     
    def testPiCation(self):
        """Test for pi-stacking interactions."""
        
        data_test = np.load('test_2k39_PiCat.npy', allow_pickle=True)
        assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in PICAT_INTERACTIONS if i]),
                     'failed to get correct pi-cation interactions')                             

    def testHydrophobicInteractions(self):
        """Test for hydrophobic interactions."""
        
        data_test = np.load('test_2k39_hph.npy', allow_pickle=True)

        try:
            assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in HPH_INTERACTIONS2]),
                     'failed to get correct hydrophobic interactions without hpb.so')
        except:
            assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in HPH_INTERACTIONS]),
                     'failed to get correct hydrophobic interactions with hpb.so')                             

    def testDisulfideBonds(self):
        """Test for disulfide bonds interactions."""
        
        data_test = np.load('test_2k39_disu.npy', allow_pickle=True)
        assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in DISU_INTERACTIONS if i]),
                     'failed to get correct disulfide bonds')                                                  
                     
