"""This module contains unit tests for :mod:`~prody.proteins.interactions`."""

import numpy as np
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
            self.ATOMS = parseDatafile('2k39_insty') # no disulfides
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

            self.ATOMS_FIRST = parseDatafile('2k39_insty_first')
            self.DCD = Trajectory(pathDatafile('2k39_insty_dcd'))
            self.DCD.link(self.ATOMS_FIRST)
            self.DCD.setCoords(self.ATOMS_FIRST)

            self.ATOMS_3O21 = parseDatafile('3o21') # has disulfides & not traj
            self.DISU_INTERACTIONS_3O21 = parseDatafile('3o21_disu')

    def testAllInteractionsCalc(self):
        """Test for calculating all types of interactions."""

        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = np.array(self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS,
                                                                                             stop_frame=13))

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
            self.data_all = np.array(self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS,
                                                                                             stop_frame=13))

            np.save('test_2k39_all.npy', np.array(self.data_all, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_all.npy', allow_pickle=True)

            try:
                assert_equal(data_test, self.ALL_INTERACTIONS2,
                             'failed to get correct interactions without hpb.so from saving and loading')
            except AssertionError:
                assert_equal(data_test, self.ALL_INTERACTIONS,
                             'failed to get correct interactions with hpb.so from saving and loading')
    
    def testAllInteractionsCalcWithTraj(self):
        """Test for calculating all types of interactions."""

        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = np.array(self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS_FIRST,
                                                                                             trajectory=self.DCD,
                                                                                             stop_frame=13))

            try:
                assert_equal(self.data_all, self.ALL_INTERACTIONS2,
                             'failed to get correct interactions without hpb.so from calculation')
            except AssertionError:
                assert_equal(self.data_all, self.ALL_INTERACTIONS,
                             'failed to get correct interactions with hpb.so from calculation')

    def testHydrogenBonds(self):
        """Test for hydrogen bonds.
        Last column is compared becasue pairs of residues can be reversed and
        order can be also different in the interactions"""

        if prody.PY3K:                
            data_test = calcHydrogenBondsTrajectory(self.ATOMS, stop_frame=13)
            assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HBS_INTERACTIONS]),
                         'failed to get correct hydrogen bonds')        
                     
    def testSaltBridgesCalc(self):
        """Test for salt bridges without saving and loading."""

        if prody.PY3K:                
            self.data_sbs = calcSaltBridgesTrajectory(self.ATOMS, stop_frame=13)
            assert_equal(sorted([i[-1][-1] for i in self.data_sbs]), sorted([i[-1][-1] for i in self.SBS_INTERACTIONS]),
                         'failed to get correct salt bridges')

        
    def testSaltBridgesSave(self):
        """Test for salt bridges with saving and loading (one type with results)."""

        if prody.PY3K:                
            self.data_sbs = calcSaltBridgesTrajectory(self.ATOMS, stop_frame=13)
            
            np.save('test_2k39_sbs.npy', np.array(self.data_sbs, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_sbs.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.SBS_INTERACTIONS if len(i) > 0]),
                         'failed to get correct salt bridges from saving and loading')


    def testRepulsiveIonicBonding(self):
        """Test for repulsive ionic bonding."""

        if prody.PY3K:                
            data_test = calcRepulsiveIonicBondingTrajectory(self.ATOMS, stop_frame=13)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.RIB_INTERACTIONS if len(i) > 0]),
                         'failed to get correct repulsive ionic bonding')                             

    def testPiStacking(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = calcPiStackingTrajectory(self.ATOMS, stop_frame=13)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.PISTACK_INTERACTIONS if len(i) > 0]),
                         'failed to get correct pi-stacking interactions')                             
                     
    def testPiCation(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = calcPiCationTrajectory(self.ATOMS, stop_frame=13)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.PICAT_INTERACTIONS if len(i) > 0]),
                         'failed to get correct pi-cation interactions')

    def testHydrophobicInteractions(self):
        """Test for hydrophobic interactions."""

        if prody.PY3K:        
            data_test = calcHydrophobicTrajectory(self.ATOMS, stop_frame=13)
            try:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS2]),
                         'failed to get correct hydrophobic interactions without hpb.so')
            except AssertionError:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS]),
                         'failed to get correct hydrophobic interactions with hpb.so')
        

    def testDisulfideBondsCalcNone(self):
        """Test for disulfide bonds interactions without saving and loading."""
        if prody.PY3K:
            data_test = calcDisulfideBondsTrajectory(self.ATOMS, stop_frame=13)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if len(i) > 0]),
                         'failed to get correct disulfide bonds from 2k39 (None) from calculation')

    def testDisulfideBondsSaveNone(self):
        """Test for disulfide bonds interactions with saving and loading (one type of interactions with 0)."""
        if prody.PY3K:
            data_test = calcDisulfideBondsTrajectory(self.ATOMS, stop_frame=13)
            np.save('test_2k39_disu.npy', np.array(data_test, dtype=object), 
                    allow_pickle=True)

            data_test = np.load('test_2k39_disu.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if len(i) > 0]),
                         'failed to get correct disulfide bonds from 2k39 (None) from saving and loading')

    def testDisulfideBondsCalcSomeNotTraj(self):
        """Test for disulfide bonds interactions without saving and loading."""
        if prody.PY3K:
            data_test = calcDisulfideBonds(self.ATOMS_3O21)
            assert_equal(sorted([i[-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1] for i in self.DISU_INTERACTIONS_3O21 if len(i) > 0]),
                         'failed to get correct disulfide bonds from 3o21 from calculation')

    def testDisulfideBondsSaveSomeNotTraj(self):
        """Test for disulfide bonds interactions with saving and loading (one type of interactions with 0)."""
        if prody.PY3K:
            data_test = calcDisulfideBonds(self.ATOMS_3O21)
            np.save('test_3o21_disu.npy', np.array(data_test, dtype=object), 
                    allow_pickle=True)

            data_test = np.load('test_3o21_disu.npy', allow_pickle=True)
            assert_equal(sorted([i[-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1] for i in self.DISU_INTERACTIONS_3O21 if len(i) > 0]),
                         'failed to get correct disulfide bonds from 3o21 from saving and loading')

    def testPiCationTrajArg(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:
            data_test = calcPiCationTrajectory(self.ATOMS, trajectory=self.ATOMS, stop_frame=13)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.PICAT_INTERACTIONS if len(i) > 0]),
                         'failed to get correct pi-cation interactions')

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
            
        if sys.version_info[1] < 13:
            self.assertTrue(imported_hpb)
        else:
            self.assertFalse(imported_hpb)

    @classmethod
    def tearDownClass(cls):
        if prody.PY3K:
            import os
            for filename in ['test_2k39_all.npy', 'test_2k39_sbs.npy', 'test_2k39_disu.npy']:
                os.remove(filename)
