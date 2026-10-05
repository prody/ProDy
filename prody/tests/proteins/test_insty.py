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

    def setUp(self):
        """Generating new data to compare it with the existing one"""
        # Reduce the number of frames to speed up tests while keeping behavior
        self.N_FRAMES = 4
        self.STOP_FRAME = self.N_FRAMES - 1
        
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
                                                                                             stop_frame=self.STOP_FRAME))

            try:
                assert_equal(self.data_all, self.ALL_INTERACTIONS2[:, :self.N_FRAMES],
                             'failed to get correct interactions without hpb.so from parallel calculation')
            except AssertionError:
                assert_equal(self.data_all, self.ALL_INTERACTIONS[:, :self.N_FRAMES],
                             'failed to get correct interactions with hpb.so from parallel calculation')

    def testAllInteractionsCalcSerial(self):
        """Test for calculating all types of interactions without parallelisation."""

        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = np.array(self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS,
                                                                                             stop_frame=self.STOP_FRAME,
                                                                                             max_proc=1))

            try:
                assert_equal(self.data_all, self.ALL_INTERACTIONS2[:, :self.N_FRAMES],
                             'failed to get correct interactions without hpb.so from serial calculation')
            except AssertionError:
                assert_equal(self.data_all, self.ALL_INTERACTIONS[:, :self.N_FRAMES],
                             'failed to get correct interactions with hpb.so from serial calculation')

    def testAllInteractionsSave(self):
        """Test for saving and loading all types of interactions."""
        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = np.array(self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS,
                                                                                             stop_frame=self.STOP_FRAME))

            np.save('test_2k39_all.npy', np.array(self.data_all, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_all.npy', allow_pickle=True)

            try:
                assert_equal(data_test, self.ALL_INTERACTIONS2[:, :self.N_FRAMES],
                             'failed to get correct interactions without hpb.so from saving and loading')
            except AssertionError:
                assert_equal(data_test, self.ALL_INTERACTIONS[:, :self.N_FRAMES],
                             'failed to get correct interactions with hpb.so from saving and loading')
    
    def testAllInteractionsCalcWithTraj(self):
        """Test for calculating all types of interactions."""

        if prody.PY3K:
            self.INTERACTIONS_ALL = InteractionsTrajectory()
            self.data_all = np.array(self.INTERACTIONS_ALL.calcProteinInteractionsTrajectory(self.ATOMS_FIRST,
                                                                                             trajectory=self.DCD,
                                                                                             stop_frame=self.STOP_FRAME))

            try:
                assert_equal(self.data_all, self.ALL_INTERACTIONS2[:, :self.N_FRAMES],
                             'failed to get correct interactions without hpb.so from calculation')
            except AssertionError:
                assert_equal(self.data_all, self.ALL_INTERACTIONS[:, :self.N_FRAMES],
                             'failed to get correct interactions with hpb.so from calculation')

    def testHydrogenBonds(self):
        """Test for hydrogen bonds.
        Last column is compared becasue pairs of residues can be reversed and
        order can be also different in the interactions"""

        if prody.PY3K:                
            data_test = calcHydrogenBondsTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            assert_equal(sorted([i[-1][-1] for i in data_test]),
                         sorted([i[-1][-1] for i in self.HBS_INTERACTIONS[:self.N_FRAMES]]),
                         'failed to get correct hydrogen bonds')        
                     
    def testSaltBridgesCalc(self):
        """Test for salt bridges without saving and loading."""

        if prody.PY3K:                
            self.data_sbs = calcSaltBridgesTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            assert_equal(sorted([i[-1][-1] for i in self.data_sbs]),
                         sorted([i[-1][-1] for i in self.SBS_INTERACTIONS[:self.N_FRAMES]]),
                         'failed to get correct salt bridges')

        
    def testSaltBridgesSave(self):
        """Test for salt bridges with saving and loading (one type with results)."""

        if prody.PY3K:                
            self.data_sbs = calcSaltBridgesTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            
            np.save('test_2k39_sbs.npy', np.array(self.data_sbs, dtype=object), allow_pickle=True)

            data_test = np.load('test_2k39_sbs.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]),
                         sorted([i[-1][-1] for i in self.SBS_INTERACTIONS[:self.N_FRAMES] if len(i) > 0]),
                         'failed to get correct salt bridges from saving and loading')


    def testRepulsiveIonicBonding(self):
        """Test for repulsive ionic bonding."""

        if prody.PY3K:                
            data_test = calcRepulsiveIonicBondingTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]),
                         sorted([i[-1][-1] for i in self.RIB_INTERACTIONS[:self.N_FRAMES] if len(i) > 0]),
                         'failed to get correct repulsive ionic bonding')                             

    def testPiStacking(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = calcPiStackingTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]),
                         sorted([i[-1][-1] for i in self.PISTACK_INTERACTIONS[:self.N_FRAMES] if len(i) > 0]),
                         'failed to get correct pi-stacking interactions')                             
                     
    def testPiCation(self):
        """Test for pi-stacking interactions."""

        if prody.PY3K:                
            data_test = calcPiCationTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]),
                         sorted([i[-1][-1] for i in self.PICAT_INTERACTIONS[:self.N_FRAMES] if len(i) > 0]),
                         'failed to get correct pi-cation interactions')

    def testHydrophobicInteractions(self):
        """Test for hydrophobic interactions.

        This test uses the full trajectory length because results are sensitive
        to the total number of processed frames. Keeping this one full preserves
        correctness while other tests remain shortened for speed.
        """

        if prody.PY3K:        
            data_test = calcHydrophobicTrajectory(self.ATOMS, stop_frame=13)
            try:
                assert_equal(sorted([i[-1][-1] for i in data_test]),
                             sorted([i[-1][-1] for i in self.HPH_INTERACTIONS2]),
                             'failed to get correct hydrophobic interactions without hpb.so')
            except AssertionError:
                assert_equal(sorted([i[-1][-1] for i in data_test]),
                             sorted([i[-1][-1] for i in self.HPH_INTERACTIONS]),
                             'failed to get correct hydrophobic interactions with hpb.so')
        

    def testDisulfideBondsCalcNone(self):
        """Test for disulfide bonds interactions without saving and loading."""
        if prody.PY3K:
            data_test = calcDisulfideBondsTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1][-1] for i in self.DISU_INTERACTIONS[:self.N_FRAMES] if len(i) > 0]),
                         'failed to get correct disulfide bonds from 2k39 (None) from calculation')

    def testDisulfideBondsSaveNone(self):
        """Test for disulfide bonds interactions with saving and loading (one type of interactions with 0)."""
        if prody.PY3K:
            data_test = calcDisulfideBondsTrajectory(self.ATOMS, stop_frame=self.STOP_FRAME)
            np.save('test_2k39_disu.npy', np.array(data_test, dtype=object), 
                    allow_pickle=True)

            data_test = np.load('test_2k39_disu.npy', allow_pickle=True)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]), 
                         sorted([i[-1][-1] for i in self.DISU_INTERACTIONS[:self.N_FRAMES] if len(i) > 0]),
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
            data_test = calcPiCationTrajectory(self.ATOMS, trajectory=self.ATOMS, stop_frame=self.STOP_FRAME)
            assert_equal(sorted([i[-1][-1] for i in data_test if len(i) > 0]),
                         sorted([i[-1][-1] for i in self.PICAT_INTERACTIONS[:self.N_FRAMES] if len(i) > 0]),
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
            for filename in ['test_2k39_all.npy', 'test_2k39_sbs.npy',
                             'test_2k39_disu.npy', 'test_3o21_disu.npy']:
                # only remove what a test actually wrote: a test that failed before
                # saving leaves the file absent, and tearing it down then raised
                # FileNotFoundError, turning one failure into a failure plus an error
                # and hiding which test was the real one
                if os.path.isfile(filename):
                    os.remove(filename)


class _CountingContext(object):
    """A stand-in for the multiprocessing module that counts Process constructions,
    so a test can assert the parallel path really was taken."""

    def __init__(self, wrapped):
        self._wrapped = wrapped
        self.processes = 0

    def __getattr__(self, name):
        return getattr(self._wrapped, name)

    def Process(self, *args, **kwargs):
        self.processes += 1
        return self._wrapped.Process(*args, **kwargs)


class TestParallelInteractions(unittest.TestCase):
    """The serial and parallel paths are tested separately and on purpose.

    Neither can be left to the machine: max_proc defaults to mp.cpu_count()//2,
    which is 1 on a two- or three-core runner and 2 on a four-core one, so a test
    that does not pass max_proc exercises whichever path the hardware happens to
    pick -- which is how a bug that made the parallel path unusable sat in a CI
    matrix passing on some runners and failing on others."""

    def setUp(self):

        self.STOP_FRAME = 5
        if prody.PY3K:
            self.ATOMS = parseDatafile('2k39_insty')

    def counts(self, max_proc, stop_frame=None):
        """Per-frame hydrogen-bond counts, which fingerprint both the result and the
        frame it belongs to."""

        frames = calcHydrogenBondsTrajectory(
            self.ATOMS,
            stop_frame=self.STOP_FRAME if stop_frame is None else stop_frame,
            max_proc=max_proc)
        return [len(frame) for frame in frames]

    def testSerialPath(self):
        """max_proc=1 on purpose.  This is the reference the parallel runs are
        checked against, and on a two- or three-core machine it is the only path
        the default would ever take."""

        if not prody.PY3K:
            return

        counts = self.counts(1)
        self.assertEqual(len(counts), self.STOP_FRAME + 1)
        self.assertTrue(all(count > 0 for count in counts),
                        'no hydrogen bonds found, so the test proves nothing')

    def testParallelPathIsTaken(self):
        """Passing max_proc explicitly is what guarantees the parallel path runs.
        The Process count is asserted too, so the test cannot quietly degrade into
        a second serial run and still pass."""

        if not prody.PY3K:
            return

        from prody.proteins import interactions

        real = interactions.mp
        counting = _CountingContext(real)
        interactions.mp = counting
        try:
            counts = self.counts(2)
        finally:
            interactions.mp = real

        self.assertEqual(counting.processes, self.STOP_FRAME + 1,
                         'the parallel path did not start one process per frame')
        assert_equal(counts, self.counts(1))

    def testParallelKeepsFrameOrder(self):
        """Results must not depend on how many processes computed them, and must come
        back in frame order rather than in the order the processes finished."""

        if not prody.PY3K:
            return

        serial = self.counts(1)
        for max_proc in (2, 3, self.STOP_FRAME + 1):
            assert_equal(self.counts(max_proc), serial,
                         'max_proc={0} disagrees with the serial result'.format(max_proc))

    def testParallelPathUnderSpawn(self):
        """Run the parallel path under the "spawn" start method, which is the default
        on macOS since Python 3.8 and on Windows always, whatever this platform's
        default is.  Spawn pickles the Process target, so this is what catches a
        nested worker on Linux, where "fork" would hide it."""

        if not prody.PY3K:
            return

        import multiprocessing
        from prody.proteins import interactions

        real = interactions.mp
        interactions.mp = multiprocessing.get_context('spawn')
        try:
            # fewer frames: every spawned child re-imports prody
            counts = self.counts(2, stop_frame=2)
        finally:
            interactions.mp = real

        assert_equal(counts, self.counts(1, stop_frame=2))

    def testWorkersArePicklable(self):
        """A Process target is pickled under "spawn", so the workers cannot be nested
        functions.  This needs no multiprocessing at all, so it is the cheap guard
        that fails the moment one is moved back inside its caller."""

        import pickle
        from prody.proteins import interactions, waterbridges

        for worker in (interactions._analyseFrame, interactions._analyseModel,
                       waterbridges._analyseWaterBridgeFrame,
                       waterbridges._analyseWaterBridgeModel,
                       waterbridges._saveBridgesFrame):
            self.assertIs(pickle.loads(pickle.dumps(worker)), worker,
                          '{0} is not picklable'.format(worker.__name__))

    def testNoWaterGivesAnEmptyResult(self):
        """No water means no water bridges, which is an answer and not an error, and
        it has to be the same answer on both paths.  Raising aborted a whole batch on
        one waterless structure; the parallel path meanwhile returned empties, so the
        two paths disagreed."""

        if not prody.PY3K:
            return

        from prody.proteins.waterbridges import calcWaterBridgesTrajectory

        for max_proc in (1, 2, 3):
            frames = calcWaterBridgesTrajectory(self.ATOMS, None, stop_frame=2,
                                                max_proc=max_proc)
            # stop_frame is the last frame, inclusive, so 0..2 is three of them
            self.assertEqual([len(frame) for frame in frames], [0, 0, 0],
                             'max_proc={0} did not give one empty result per frame'
                             .format(max_proc))


    def testStopFrameIsTheLastFrame(self):
        """stop_frame is the index of the last frame to read, inclusive, on both the
        trajectory and the multi-model path and in both modules.  waterbridges sliced
        [start_frame:stop_frame] on its multi-model path, dropping the last frame, and
        did not resolve the -1 default there at all."""

        if not prody.PY3K:
            return

        from prody.proteins.waterbridges import calcWaterBridgesTrajectory

        for max_proc in (1, 2):
            for stop_frame, expected in ((0, 1), (2, 3), (4, 5)):
                frames = calcWaterBridgesTrajectory(self.ATOMS, None,
                                                    stop_frame=stop_frame,
                                                    max_proc=max_proc)
                self.assertEqual(len(frames), expected,
                                 'stop_frame={0} at max_proc={1} gave {2} frames'
                                 .format(stop_frame, max_proc, len(frames)))

            # -1 means all of them, and must not become an empty slice
            frames = calcWaterBridgesTrajectory(self.ATOMS, None, stop_frame=-1,
                                                max_proc=max_proc)
            self.assertEqual(len(frames), self.ATOMS.numCoordsets())

        # the same bound in interactions
        self.assertEqual(len(self.counts(1, stop_frame=2)), 3)


def _exitsNonZero():
    """A worker that dies, for testing that the parent notices."""

    import sys
    sys.exit(1)


def _exitsCleanly():
    """A worker that does nothing at all, successfully."""

    return


class TestJoinProcesses(unittest.TestCase):
    """A Process reports failure only through its exitcode: the exception is raised in
    the child, whose traceback goes to stderr, and a worker that writes into a shared
    list leaves its slot untouched.  Without checking the exit code, a crashed worker
    is indistinguishable from a genuinely empty result."""

    def testRaisesWhenAWorkerDies(self):

        import multiprocessing as mp
        from prody.utilities import joinProcesses

        processes = [mp.Process(target=_exitsNonZero) for _ in range(2)]
        for process in processes:
            process.start()

        self.assertRaises(RuntimeError, joinProcesses, processes, 'frame')

    def testSilentWhenEveryWorkerSucceeds(self):

        import multiprocessing as mp
        from prody.utilities import joinProcesses

        processes = [mp.Process(target=_exitsCleanly) for _ in range(2)]
        for process in processes:
            process.start()

        joinProcesses(processes, 'frame')       # must not raise
        for process in processes:
            self.assertEqual(process.exitcode, 0)
