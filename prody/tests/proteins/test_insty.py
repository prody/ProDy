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
from prody.proteins.interactions import InteractionsTrajectory, Interactions

import sys

class TestInteractions(unittest.TestCase):

    @staticmethod
    def _initialize_data(target):
        """
        Helper method to load data and pre-calculate trajectories.
        Target can be 'cls' (for setUpClass) or 'self' (for setUp).
        """
        # 1. Load Data Files
        target.ATOMS = parseDatafile('2k39_insty') # no disulfides
        target.ALL_INTERACTIONS = parseDatafile('2k39_all')
        target.ALL_INTERACTIONS2 = parseDatafile('2k39_all2')
        target.HBS_INTERACTIONS = parseDatafile('2k39_hbs')
        target.SBS_INTERACTIONS = parseDatafile('2k39_sbs')
        target.RIB_INTERACTIONS = parseDatafile('2k39_rib')
        target.PISTACK_INTERACTIONS = parseDatafile('2k39_PiStack')
        target.PICAT_INTERACTIONS = parseDatafile('2k39_PiCat')
        target.HPH_INTERACTIONS = parseDatafile('2k39_hph')
        target.HPH_INTERACTIONS2 = parseDatafile('2k39_hph2')
        target.DISU_INTERACTIONS = parseDatafile('2k39_disu')

        target.ATOMS_FIRST = parseDatafile('2k39_insty_first')
        target.DCD = Trajectory(pathDatafile('2k39_insty_dcd'))
        target.DCD.link(target.ATOMS_FIRST)
        target.DCD.setCoords(target.ATOMS_FIRST)

        target.ATOMS_3O21 = parseDatafile('3o21') # has disulfides & not traj
        target.DISU_INTERACTIONS_3O21 = parseDatafile('3o21_disu')

        # 2. Pre-calculate expensive trajectories
        target.calc_all_13 = np.array(InteractionsTrajectory().calcProteinInteractionsTrajectory(target.ATOMS, stop_frame=13))
        target.calc_all_traj_13 = np.array(InteractionsTrajectory().calcProteinInteractionsTrajectory(target.ATOMS_FIRST, trajectory=target.DCD, stop_frame=13))
        target.calc_hbs_13 = calcHydrogenBondsTrajectory(target.ATOMS, stop_frame=13)
        target.calc_sbs_13 = calcSaltBridgesTrajectory(target.ATOMS, stop_frame=13)
        target.calc_rib_13 = calcRepulsiveIonicBondingTrajectory(target.ATOMS, stop_frame=13)
        target.calc_pistack_13 = calcPiStackingTrajectory(target.ATOMS, stop_frame=13)
        target.calc_picat_13 = calcPiCationTrajectory(target.ATOMS, stop_frame=13)
        target.calc_picat_traj_arg_13 = calcPiCationTrajectory(target.ATOMS, trajectory=target.ATOMS, stop_frame=13)
        target.calc_hph_13 = calcHydrophobicTrajectory(target.ATOMS, stop_frame=13)
        target.calc_disu_13 = calcDisulfideBondsTrajectory(target.ATOMS, stop_frame=13)
        target.calc_disu_3o21 = calcDisulfideBonds(target.ATOMS_3O21)

    @classmethod
    def setUpClass(cls):
        """Run setup ONCE for Python 3 for performance."""
        if prody.PY3K:
            cls._initialize_data(cls)

    def setUp(self):
        """Run setup PER TEST for Python 2 for compatibility."""
        if not prody.PY3K:
            self._initialize_data(self)

    def testAllInteractionsCalc(self):
        """Test for calculating all types of interactions."""

        if prody.PY3K:
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            self.data_all = self.calc_all_13

            try:
                assert_equal(self.data_all, self.ALL_INTERACTIONS2,
                             'failed to get correct interactions without hpb.so from parallel calculation')
            except AssertionError:
                assert_equal(self.data_all, self.ALL_INTERACTIONS,
                             'failed to get correct interactions with hpb.so from parallel calculation')

    def testAllInteractionsCalcSerial(self):
        """Test for calculating all types of interactions without parallelisation."""

        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            # Note: This is kept as a live calculation because it tests specific flags (max_proc=1) 
            # and uses a shorter stop_frame (3) so it is not a major bottleneck.
            self.data_all = np.array(self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS,
                                                                                             stop_frame=3,
                                                                                             max_proc=1))

            try:
                assert_equal(self.data_all, self.ALL_INTERACTIONS2[:, :4],
                             'failed to get correct interactions without hpb.so from serial calculation')
            except AssertionError:
                assert_equal(self.data_all, self.ALL_INTERACTIONS[:, :4],
                             'failed to get correct interactions with hpb.so from serial calculation')

    def testAllInteractionsSave(self):
        """Test for saving and loading all types of interactions."""
        if prody.PY3K:
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            self.data_all = self.calc_all_13

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
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            self.data_all = self.calc_all_traj_13

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
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_hbs_13
            assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HBS_INTERACTIONS]),
                         'failed to get correct hydrogen bonds')        
                     
    def testSaltBridgesCalc(self):
        """Test for salt bridges without saving and loading."""

        if prody.PY3K:                
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            self.data_sbs = self.calc_sbs_13
            assert_equal(sorted([i[-1][-1] for i in self.data_sbs]), sorted([i[-1][-1] for i in self.SBS_INTERACTIONS]),
                         'failed to get correct salt bridges')

        
    def testSaltBridgesSave(self):
        """Test for salt bridges with saving and loading (one type with results)."""

        if prody.PY3K:                
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            self.data_sbs = self.calc_sbs_13
            
            np.save('test_2k39_sbs.npy', np.array(self.data_sbs, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_sbs.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.SBS_INTERACTIONS if len(i) > 0]),
                         'failed to get correct salt bridges from saving and loading')


    def testRepulsiveIonicBonding(self):
        """Test for repulsive ionic bonding."""

        if prody.PY3K:                
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_rib_13
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.RIB_INTERACTIONS if len(i) > 0]),
                         'failed to get correct repulsive ionic bonding')                             

    def testPiStacking(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_pistack_13
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.PISTACK_INTERACTIONS if len(i) > 0]),
                         'failed to get correct pi-stacking interactions')                             
                     
    def testPiCation(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_picat_13
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), sorted([i[-1][-1] for i in self.PICAT_INTERACTIONS if len(i) > 0]),
                         'failed to get correct pi-cation interactions')

    def testHydrophobicInteractions(self):
        """Test for hydrophobic interactions."""

        if prody.PY3K:        
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_hph_13
            try:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS2]),
                         'failed to get correct hydrophobic interactions without hpb.so')
            except AssertionError:
                assert_equal(sorted([i[-1][-1] for i in data_test]), sorted([i[-1][-1] for i in self.HPH_INTERACTIONS]),
                         'failed to get correct hydrophobic interactions with hpb.so')
        

    def testDisulfideBondsCalcNone(self):
        """Test for disulfide bonds interactions without saving and loading."""
        if prody.PY3K:
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_disu_13
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if len(i) > 0]),
                         'failed to get correct disulfide bonds from 2k39 (None) from calculation')

    def testDisulfideBondsSaveNone(self):
        """Test for disulfide bonds interactions with saving and loading (one type of interactions with 0)."""
        if prody.PY3K:
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_disu_13
            np.save('test_2k39_disu.npy', np.array(data_test, dtype=object), 
                    allow_pickle=True)

            data_test = np.load('test_2k39_disu.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1][-1] for i in self.DISU_INTERACTIONS if len(i) > 0]),
                         'failed to get correct disulfide bonds from 2k39 (None) from saving and loading')

    def testDisulfideBondsCalcSomeNotTraj(self):
        """Test for disulfide bonds interactions without saving and loading."""
        if prody.PY3K:
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_disu_3o21
            assert_equal(sorted([i[-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1] for i in self.DISU_INTERACTIONS_3O21 if len(i) > 0]),
                         'failed to get correct disulfide bonds from 3o21 from calculation')

    def testDisulfideBondsSaveSomeNotTraj(self):
        """Test for disulfide bonds interactions with saving and loading (one type of interactions with 0)."""
        if prody.PY3K:
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_disu_3o21
            np.save('test_3o21_disu.npy', np.array(data_test, dtype=object), 
                    allow_pickle=True)

            data_test = np.load('test_3o21_disu.npy', allow_pickle=True)
            assert_equal(sorted([i[-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1] for i in self.DISU_INTERACTIONS_3O21 if len(i) > 0]),
                         'failed to get correct disulfide bonds from 3o21 from saving and loading')

    def testPiCationTrajArg(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:
            # OPTIMIZATION: Use pre-calculated data from setUpClass
            data_test = self.calc_picat_traj_arg_13
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
            for filename in ['test_2k39_all.npy', 'test_2k39_sbs.npy', 'test_2k39_disu.npy', 'test_3o21_disu.npy']:
                if os.path.exists(filename):
                    os.remove(filename)
