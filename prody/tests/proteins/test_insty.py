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

    def testAllInteractionsCalc(self):
        """Test for calculating all types of interactions."""

        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS)

            try:
                assert_equal(self.data_all, self.ALL_INTERACTIONS2,
                             'failed to get correct interactions without hpb.so from calculation')
            except AssertionError:
                assert_equal(self.data_all, self.ALL_INTERACTIONS,
                             'failed to get correct interactions with hpb.so from calculation')

    def testAllInteractionsSave(self):
        """Test for saving and loading all types of interactions."""
        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS)

            np.save('test_2k39_all.npy', np.array(self.data_all, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_all.npy', allow_pickle=True)

            try:
                assert_equal(data_test, self.ALL_INTERACTIONS2,
                             'failed to get correct interactions without hpb.so from saving and loading')
            except AssertionError:
                assert_equal(data_test, self.ALL_INTERACTIONS,
                             'failed to get correct interactions with hpb.so from saving and loading')
    
    def testHydrogenBonds(self):
        """Test for hydrogen bonds.
        Last column is compared becasue pairs of residues can be reversed and
        order can be also different in the interactions"""

        if prody.PY3K:                
            data_test = calcHydrogenBondsTrajectory(self.ATOMS)
            assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HBS_INTERACTIONS]),
                         'failed to get correct hydrogen bonds')        
                     
    def testSaltBridgesCalc(self):
        """Test for salt bridges without saving and loading."""

        if prody.PY3K:                
            self.data_sbs = calcSaltBridgesTrajectory(self.ATOMS)
            assert_equal(sorted([i[-1][-1] for i in self.data_sbs]), sorted([i[-1][-1] for i in self.SBS_INTERACTIONS]),
                         'failed to get correct salt bridges')

        
    def testSaltBridgesSave(self):
        """Test for salt bridges with saving and loading (one type with results)."""

        if prody.PY3K:                
            self.data_sbs = calcSaltBridgesTrajectory(self.ATOMS)
            
            np.save('test_2k39_sbs.npy', np.array(self.data_sbs, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_sbs.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.SBS_INTERACTIONS if i]),
                         'failed to get correct salt bridges from saving and loading')


    def testRepulsiveIonicBonding(self):
        """Test for repulsive ionic bonding."""

        if prody.PY3K:                
            data_test = calcRepulsiveIonicBondingTrajectory(self.ATOMS)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.RIB_INTERACTIONS if i]),
                         'failed to get correct repulsive ionic bonding')                             

    def testPiStacking(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = calcPiStackingTrajectory(self.ATOMS)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.PISTACK_INTERACTIONS if i]),
                         'failed to get correct pi-stacking interactions')                             
                     
    def testPiCation(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = calcPiCationTrajectory(self.ATOMS)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.PICAT_INTERACTIONS if i]),
                         'failed to get correct pi-cation interactions')


    def testHydrophobicInteractions(self):
        """Test for hydrophobic interactions."""

        if prody.PY3K:        
            data_test = calcHydrophobicTrajectory(self.ATOMS)
            try:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS2]),
                         'failed to get correct hydrophobic interactions without hpb.so')
            except AssertionError:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS]),
                         'failed to get correct hydrophobic interactions with hpb.so')
        

    def testDisulfideBondsCalc(self):
        """Test for disulfide bonds interactions without saving and loading."""
        if prody.PY3K:
            data_test = calcDisulfideBondsTrajectory(self.ATOMS)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if i]),
                          'failed to get correct disulfide bonds from calculation')
        
    def testDisulfideBondsSave(self):
        """Test for disulfide bonds interactions with saving and loading (one type of interactions with 0)."""
        if prody.PY3K:
            data_test = calcDisulfideBondsTrajectory(self.ATOMS)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if i]),
                          'failed to get correct disulfide bonds from calculation')
            
            np.save('test_2k39_disu.npy', np.array(data_test, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_disu.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if i]), sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if i]),
                         'failed to get correct disulfide bonds from saving and loading')
