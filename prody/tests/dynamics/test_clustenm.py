"""This module contains unit tests for :mod:`~prody.dynamics.clustenm`."""

import os
import shutil
import sys
import tempfile
import types
import warnings
import numpy as np
from numpy.testing import *
from prody.utilities import importDec
dec = importDec()

import prody
from prody import *
from prody import LOGGER
from prody.dynamics.clustenm import ClustENM, ClustRTB, ClustImANM
from prody.tests import unittest
from prody.tests.datafiles import *

# Import mock tools
try:
    from unittest.mock import MagicMock, patch
except ImportError:
    from mock import MagicMock, patch

# Prevent threading hangs on remote servers
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

LOGGER.verbosity = 'none'

# The title a ClustENM gets from setAtoms is the title of the atoms plus
# '_clustenm', and writePDBFixed strips those 8 characters back off to build
# its default filename. Setting the title directly gives the same result
# without needing PDBFixer to fix a real structure.
TITLE = 'my prot'
CLUSTENM_TITLE = TITLE + '_clustenm'
SPACED_NAME = TITLE + '_fixed.pdb'
UNDERSCORED_NAME = TITLE.replace(' ', '_') + '_fixed.pdb'

# Fixing and minimising a structure needs both of these. They are optional
# dependencies, so the tests that use them are skipped rather than failed when
# they are missing, with a warning so that the gap is not silent. The rest of
# the tests either stub out OpenMM or do not get as far as needing it, and so
# run either way.
OPENMM_SKIP_MSG = 'PDBFixer and OpenMM are needed to fix a structure'
try:
    import openmm
    import pdbfixer
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False
    warnings.warn(OPENMM_SKIP_MSG + ', so the ClustENM tests that need them '
                  'are being skipped')

# Which platforms OpenMM can run a simulation on depends on how it was built
# and on the hardware it finds, so the tests of each are skipped rather than
# failed when it cannot use that one, with a warning so that not having tried
# them is not silent.
def usablePlatform(name):
    """Returns whether OpenMM can make a context on the platform called *name*.

    Being one of the platforms it lists is not enough. The OpenCL plugin loads
    on machines that have no device it can use, where the platform is offered
    but making a context on it raises, so each one is tried rather than looked
    up."""

    if not HAVE_OPENMM:
        return False

    from openmm import Context, Platform, System, VerletIntegrator

    listed = [Platform.getPlatform(i).getName()
              for i in range(Platform.getNumPlatforms())]
    if name not in listed:
        return False

    # a single particle is enough to make a context and find out
    system = System()
    system.addParticle(1.)
    try:
        Context(system, VerletIntegrator(0.001),
                Platform.getPlatformByName(name))
    except Exception:
        return False

    return True


PLATFORMS = {name: usablePlatform(name) for name in ('CPU', 'CUDA', 'OpenCL')}
_unusable = sorted(name for name, usable in PLATFORMS.items() if not usable)
if HAVE_OPENMM and _unusable:
    warnings.warn('OpenMM cannot run on the %s platform here, so the ClustENM '
                  'tests that run on it are being skipped'
                  % ' or '.join(_unusable))


def platformSkipMsg(platform):
    """Returns the reason for skipping a test of *platform*."""

    return 'OpenMM cannot run on the %s platform here' % platform


def mockOpenMMApp():
    """Returns a ``patch.dict`` of :mod:`sys.modules` and a mock standing in
    for :class:`openmm.app.PDBFile`, so that tests of code which only writes
    files can run without OpenMM installed."""

    pdbfile = MagicMock()
    app = types.ModuleType('openmm.app')
    app.PDBFile = pdbfile
    openmm = types.ModuleType('openmm')
    openmm.app = app
    patcher = patch.dict(sys.modules, {'openmm': openmm, 'openmm.app': app})
    return patcher, pdbfile


class TestClustENM(unittest.TestCase):
    """Tests of the arguments of :meth:`.ClustENM.run`, which are checked
    before any structure or force field is touched.

    Every one of these but the first raises ValueError, so each is matched
    against its message to make sure it is not passing for another reason."""

    def testParallelWrongType(self):
        """Test response to wrong type *parallel* argument."""
        if prody.PY3K:
            with self.assertRaises(
                    TypeError,
                    msg='ClustENM run with a wrong type parallel failed to '
                        'raise a TypeError'):
                ClustENM().run(parallel='nogood')

    def testNoMaxclustOrThreshold(self):
        """Test response to neither *maxclust* nor *threshold* being set.

        Only generations that sample conformers need to be clustered, so this
        is not required when there are none, as TestClustENMMinimiseOnly
        covers."""
        if prody.PY3K:
            with self.assertRaisesRegex(
                    ValueError, 'maxclust or threshold',
                    msg='ClustENM run without maxclust or threshold failed to '
                        'say that one of them is needed'):
                ClustENM().run()

    def testMaxclustSizeMismatch(self):
        """Test response to *maxclust* not covering every generation."""
        if prody.PY3K:
            with self.assertRaisesRegex(
                    ValueError, 'maxclusts were given',
                    msg='ClustENM run with too few maxclusts failed to say so'):
                ClustENM().run(n_gens=5, maxclust=(10, 30))

    def testThresholdSizeMismatch(self):
        """Test response to *threshold* not covering every generation."""
        if prody.PY3K:
            with self.assertRaisesRegex(
                    ValueError, 'thresholds were given',
                    msg='ClustENM run with too few thresholds failed to say so'):
                ClustENM().run(n_gens=5, threshold=(1.5, 2.0))

    def testAtomsUnset(self):
        """Test response to running without having set atoms."""
        if prody.PY3K:
            with self.assertRaisesRegex(
                    ValueError, 'atoms are not set',
                    msg='ClustENM run without atoms failed to say so'):
                ClustENM().run(maxclust=2)

    def testRunWhenBuilt(self):
        """Test response to running an ensemble that is already built."""
        if prody.PY3K:
            clustenm = ClustENM()
            # standing in for a run, which _isBuilt tests for by its conformers
            clustenm._confs = np.zeros((1, 3, 3))
            with self.assertRaisesRegex(
                    ValueError, 'has been built',
                    msg='ClustENM run on a built ensemble failed to say so'):
                clustenm.run(maxclust=2)

    def testEmptyIndex(self):
        """Test response to indexing with an empty tuple."""
        if prody.PY3K:
            with self.assertRaisesRegex(
                    ValueError, 'cannot be empty',
                    msg='ClustENM indexed with an empty tuple failed to say so'):
                ClustENM()[()]


class TestClustENMUnbuilt(unittest.TestCase):
    """Tests of a ClustENM that has not been given atoms or run yet."""

    def setUp(self):
        if prody.PY3K:
            self.clustenm = ClustENM()

    def testNotBuilt(self):
        """Test that a new ClustENM is not built"""
        if prody.PY3K:
            self.assertFalse(self.clustenm._isBuilt(),
                             'a new ClustENM claimed to be built')

    def testNoAtoms(self):
        """Test that a new ClustENM has no atoms"""
        if prody.PY3K:
            self.assertIsNone(self.clustenm.getAtoms(),
                              'a new ClustENM had atoms')

    def testNoConfs(self):
        """Test that a new ClustENM has no conformers"""
        if prody.PY3K:
            assert_equal(self.clustenm.numConfs(), 0,
                         'a new ClustENM had conformers')

    def testNoKeys(self):
        """Test that a new ClustENM has no generation keys"""
        if prody.PY3K:
            self.assertIsNone(self.clustenm.getKeys(),
                              'a new ClustENM had generation keys')

    def testDefaultGenerations(self):
        """Test the default number of generations"""
        if prody.PY3K:
            assert_equal(self.clustenm.numGenerations(), 5,
                         'ClustENM failed to default to 5 generations')


class TestClustENMTitle(unittest.TestCase):

    def testTitleWrongType(self):
        """Test response to wrong type *title* argument."""
        if prody.PY3K:
            with self.assertRaises(
                    TypeError,
                    msg='ClustENM given a wrong type title failed to raise a '
                        'TypeError'):
                ClustENM().setTitle(1)

    def testTitleUnset(self):
        """Test title of a ClustENM without atoms or an explicit title"""
        if prody.PY3K:
            assert_equal(ClustENM().getTitle(), 'Unknown',
                         'ClustENM without atoms failed to give a dummy title')

    def testTitleRoundTrip(self):
        """Test that a title survives being set"""
        if prody.PY3K:
            clustenm = ClustENM()
            clustenm.setTitle(CLUSTENM_TITLE)
            assert_equal(clustenm.getTitle(), CLUSTENM_TITLE,
                         'ClustENM failed to give back the title it was set')


class TestClustENMBlocks(unittest.TestCase):
    """The block based subclasses cannot build their ANM without blocks, and
    check for them before anything else. That check must not stand in for the
    ones their parent makes, so each of those is tested again here with blocks
    set, where it would otherwise be masked."""

    BLOCK_CLASSES = (ClustRTB, ClustImANM)

    def withBlocks(self, cls):
        """Returns an instance of *cls* whose blocks are set, so that its own
        check passes and those of :meth:`.ClustENM.run` are reached. The blocks
        themselves are never used, because no run here gets as far as building
        an ANM."""

        clustenm = cls()
        clustenm.setBlocks(np.zeros(10, dtype=int))
        return clustenm

    def testBlocksUnset(self):
        """Test response to running without blocks."""
        if prody.PY3K:
            for cls in self.BLOCK_CLASSES:
                with self.assertRaisesRegex(
                        ValueError, 'blocks are not set',
                        msg='%s run without blocks failed to say so'
                            % cls.__name__):
                    cls().run(maxclust=2)

    def testBlocksSetReachesAtomsCheck(self):
        """Test that blocks being set does not mask unset atoms."""
        if prody.PY3K:
            for cls in self.BLOCK_CLASSES:
                with self.assertRaisesRegex(
                        ValueError, 'atoms are not set',
                        msg='%s run with blocks but no atoms failed to say '
                            'that the atoms are missing' % cls.__name__):
                    self.withBlocks(cls).run(maxclust=2)

    def testBlocksSetReachesParallelCheck(self):
        """Test that blocks being set does not mask a wrong type *parallel*."""
        if prody.PY3K:
            for cls in self.BLOCK_CLASSES:
                with self.assertRaises(
                        TypeError,
                        msg='%s run with blocks and a wrong type parallel '
                            'failed to raise a TypeError' % cls.__name__):
                    self.withBlocks(cls).run(maxclust=2, parallel='nogood')

    def testBlocksSetReachesClusteringCheck(self):
        """Test that blocks being set does not mask unset clustering."""
        if prody.PY3K:
            for cls in self.BLOCK_CLASSES:
                with self.assertRaisesRegex(
                        ValueError, 'maxclust or threshold',
                        msg='%s run with blocks but no maxclust or threshold '
                            'failed to say that one of them is needed'
                            % cls.__name__):
                    self.withBlocks(cls).run()


class TestClustENMBlocksMatchAtoms(unittest.TestCase):
    """The blocks of the block based subclasses have to describe the coarse
    grained nodes of the structure, whose number is only known once atoms are
    set, so these need a really fixed structure to check against. What counts
    as valid blocks is checked without one in TestCheckBlocks, and against RTB
    and imANM directly in TestRTBBlocks."""

    BLOCK_CLASSES = (ClustRTB, ClustImANM)

    # how many conformers each generation samples, and the most it can keep
    N_CONFS = 2

    def setUp(self):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        self.cwd = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

        self.ATOMS = parseDatafile('1ubi')

    def tearDown(self):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        os.chdir(self.cwd)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def runWithBlocks(self, cls, blocks):
        """Runs *cls* over the fixed structure with *blocks*, which is given the
        number of coarse grained nodes to size itself against."""

        clustenm = cls()
        clustenm.setAtoms(self.ATOMS)
        clustenm.setBlocks(blocks(clustenm._n_cg))
        clustenm.run(n_confs=self.N_CONFS, n_gens=1, n_modes=2,
                     maxclust=self.N_CONFS, sim=False, outlier=False)
        return clustenm

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testBlocksWrongLength(self):
        """Test response to blocks that do not cover every node.

        The nodes the blocks are counted against are the coarse grained ones,
        not every atom of the fixed structure, so the message says so."""
        if prody.PY3K:
            for cls in self.BLOCK_CLASSES:
                with self.assertRaisesRegex(
                        ValueError, 'must match number of coarse grained nodes',
                        msg='%s run with too few blocks failed to say so'
                            % cls.__name__):
                    self.runWithBlocks(
                        cls, lambda n: np.zeros(n - 5, dtype=int))

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testBlocksAllOneValue(self):
        """Test response to blocks that are all the same value.

        One block makes the structure a single rigid body, which has no
        internal modes to sample along, so this has to be rejected rather than
        left to fail in the eigendecomposition downstream."""
        if prody.PY3K:
            for cls in self.BLOCK_CLASSES:
                with self.assertRaisesRegex(
                        ValueError, 'degrees of freedom of a rigid body',
                        msg='%s run with a single block failed to say so'
                            % cls.__name__):
                    self.runWithBlocks(cls, lambda n: np.zeros(n, dtype=int))

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testBlocksTwoValues(self):
        """Test that two blocks are enough to sample along.

        How many of the conformers sampled from them survive being filtered
        and clustered depends on the conformers themselves, and so on the
        modes the eigensolver of the machine gives, so this counts them rather
        than expecting a number of its own."""
        if prody.PY3K:
            for cls in self.BLOCK_CLASSES:
                clustenm = self.runWithBlocks(
                    cls, lambda n: np.arange(n) // (n // 2))

                sampled = clustenm.numConfs(1)
                assert_equal(clustenm.numConfs(0), 1,
                             '%s run with two blocks failed to give the '
                             'minimised starting structure' % cls.__name__)
                self.assertGreaterEqual(sampled, 1,
                                        '%s run with two blocks failed to '
                                        'sample along them' % cls.__name__)
                self.assertLessEqual(sampled, self.N_CONFS,
                                     '%s run with two blocks gave more '
                                     'conformers than it sampled'
                                     % cls.__name__)
                assert_equal(clustenm.numConfs(), 1 + sampled,
                             '%s run with two blocks failed to give the '
                             'conformers of both generations' % cls.__name__)


class TestClustENMWritePDBFixed(unittest.TestCase):
    """Tests for the filename handling of :meth:`.ClustENM.writePDBFixed`,
    which :func:`.runANMD` depends on to find the fixed structure it wrote.
    A title containing spaces used to make the two disagree, because ANMD
    replaced the spaces with underscores and writePDBFixed did not."""

    def setUp(self):
        if not prody.PY3K:
            return

        # writePDBFixed writes into the current directory
        self.cwd = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

        self.patcher, self.PDBFile = mockOpenMMApp()
        self.patcher.start()

        self.clustenm = ClustENM()
        self.clustenm.setTitle(CLUSTENM_TITLE)

    def tearDown(self):
        if not prody.PY3K:
            return

        self.patcher.stop()
        os.chdir(self.cwd)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def testDefaultFilename(self):
        """Test that the default filename is derived from the title"""
        if prody.PY3K:
            filename = self.clustenm.writePDBFixed()
            assert_equal(filename, SPACED_NAME,
                         'writePDBFixed failed to derive the default filename '
                         'from the title')
            self.assertTrue(os.path.exists(SPACED_NAME),
                            'writePDBFixed failed to write the file it named')

    def testReturnsFilenameWritten(self):
        """Test that the returned filename is the one written to"""
        if prody.PY3K:
            filename = self.clustenm.writePDBFixed(UNDERSCORED_NAME)
            self.assertIsNotNone(filename,
                                 'writePDBFixed failed to return a filename')
            self.assertTrue(os.path.exists(filename),
                            'writePDBFixed returned a filename it did not write')

    def testExplicitFilename(self):
        """Test that an explicit filename is used verbatim"""
        if prody.PY3K:
            filename = self.clustenm.writePDBFixed(UNDERSCORED_NAME)
            assert_equal(filename, UNDERSCORED_NAME,
                         'writePDBFixed failed to use the filename it was given')
            self.assertTrue(os.path.exists(UNDERSCORED_NAME),
                            'writePDBFixed failed to write to the filename it '
                            'was given')

    def testExplicitFilenameOverridesTitle(self):
        """Test that an explicit filename suppresses the title-derived one.

        This is the case that made runANMD raise FileNotFoundError: it asked
        for the underscored name and the spaced one was written instead."""
        if prody.PY3K:
            self.clustenm.writePDBFixed(UNDERSCORED_NAME)
            self.assertFalse(os.path.exists(SPACED_NAME),
                             'writePDBFixed wrote the title-derived filename '
                             'instead of the one it was given')

    def testReplaceSpaces(self):
        """Test that *replace_spaces* underscores the default filename"""
        if prody.PY3K:
            filename = self.clustenm.writePDBFixed(replace_spaces=True)
            assert_equal(filename, UNDERSCORED_NAME,
                         'writePDBFixed with replace_spaces failed to replace '
                         'the spaces of the default filename')
            self.assertTrue(os.path.exists(UNDERSCORED_NAME),
                            'writePDBFixed with replace_spaces failed to write '
                            'the file it named')
            self.assertFalse(os.path.exists(SPACED_NAME),
                             'writePDBFixed with replace_spaces wrote a '
                             'filename containing spaces')

    def testReplaceSpacesWithFilename(self):
        """Test that *replace_spaces* also applies to an explicit filename"""
        if prody.PY3K:
            filename = self.clustenm.writePDBFixed(SPACED_NAME,
                                                   replace_spaces=True)
            assert_equal(filename, UNDERSCORED_NAME,
                         'writePDBFixed with replace_spaces failed to replace '
                         'the spaces of the given filename')

    def testWriteFileArguments(self):
        """Test that the fixed topology and positions are the ones written"""
        if prody.PY3K:
            self.clustenm.writePDBFixed(UNDERSCORED_NAME)
            assert_equal(self.PDBFile.writeFile.call_count, 1,
                         'writePDBFixed failed to write the file once')

            args, kwargs = self.PDBFile.writeFile.call_args
            assert_equal(args[0], self.clustenm._topology,
                         'writePDBFixed failed to write the fixed topology')
            assert_equal(args[1], self.clustenm._positions,
                         'writePDBFixed failed to write the fixed positions')
            assert_equal(kwargs.get('keepIds'), True,
                         'writePDBFixed failed to keep the original ids')

    def testStreamClosed(self):
        """Test that the file written to is closed afterwards"""
        if prody.PY3K:
            self.clustenm.writePDBFixed(UNDERSCORED_NAME)
            stream = self.PDBFile.writeFile.call_args[0][2]
            self.assertTrue(stream.closed,
                            'writePDBFixed failed to close the file it wrote')


class TestClustENMFixed(unittest.TestCase):
    """Tests against a really fixed structure, rather than the stubbed OpenMM
    the tests above use, so that fixing itself and the file it is written to
    are both checked."""

    def setUp(self):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        # writePDBFixed writes into the current directory
        self.cwd = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

        self.ATOMS = parseDatafile('1ubi')
        self.ATOMS.setTitle(TITLE)

    def tearDown(self):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        os.chdir(self.cwd)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testSetAtomsFixes(self):
        """Test that setAtoms fixes the structure it is given"""
        if prody.PY3K:
            clustenm = ClustENM()
            clustenm.setAtoms(self.ATOMS)

            assert_equal(clustenm.getTitle(), CLUSTENM_TITLE,
                         'setAtoms failed to take the title of the atoms')
            # fixing adds the hydrogens that 1ubi does not have
            self.assertGreater(clustenm.getAtoms().numAtoms(),
                               self.ATOMS.numAtoms(),
                               'setAtoms failed to add the missing atoms')
            self.assertIsNotNone(clustenm.getAtoms().hydrogen,
                                 'setAtoms failed to add the missing hydrogens')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testWritePDBFixedSpacedTitle(self):
        """Test writing a really fixed structure whose title has a space"""
        if prody.PY3K:
            clustenm = ClustENM()
            clustenm.setAtoms(self.ATOMS)

            filename = clustenm.writePDBFixed(UNDERSCORED_NAME)
            assert_equal(filename, UNDERSCORED_NAME,
                         'writePDBFixed failed to use the filename it was given')
            self.assertTrue(os.path.exists(filename),
                            'writePDBFixed returned a filename it did not write')

            fixed = parsePDB(filename, compressed=False)
            self.assertIsNotNone(fixed,
                                 'writePDBFixed failed to write a PDB file')
            assert_equal(fixed.numAtoms(), clustenm.getAtoms().numAtoms(),
                         'writePDBFixed failed to write every fixed atom')


class TestClustENMRun(unittest.TestCase):
    """Tests of a small ClustENM run. Sampling and clustering are cheap, so
    *sim* is turned off to skip the molecular dynamics, leaving the
    minimisation of each conformer as the only expensive part."""

    N_CONFS = 2
    N_GENS = 1

    @classmethod
    def setUpClass(cls):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        cls.cwd = os.getcwd()
        cls.tmpdir = tempfile.mkdtemp()
        os.chdir(cls.tmpdir)

        cls.ATOMS = parseDatafile('1ubi')
        cls.CLUSTENM = ClustENM()
        cls.CLUSTENM.setAtoms(cls.ATOMS)
        cls.CLUSTENM.run(n_confs=cls.N_CONFS, n_gens=cls.N_GENS, n_modes=2,
                         maxclust=2, sim=False, outlier=False)

    @classmethod
    def tearDownClass(cls):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        os.chdir(cls.cwd)
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testBuilt(self):
        """Test that a run builds the ensemble"""
        if prody.PY3K:
            self.assertTrue(self.CLUSTENM._isBuilt(),
                            'ClustENM run failed to build the ensemble')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testGenerations(self):
        """Test that a run gives the number of generations asked for"""
        if prody.PY3K:
            assert_equal(self.CLUSTENM.numGenerations(), self.N_GENS,
                         'ClustENM run failed to give the number of '
                         'generations it was asked for')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testNumConfs(self):
        """Test the conformers of each generation of a run.

        How many of the conformers sampled survive being filtered and
        clustered depends on the conformers themselves, and so on the modes
        the eigensolver of the machine gives, so generation 1 is counted
        rather than expected to keep every one of them."""
        if prody.PY3K:
            # generation 0 is the minimised starting structure
            assert_equal(self.CLUSTENM.numConfs(0), 1,
                         'ClustENM run failed to give 1 conformer for '
                         'generation 0')
            sampled = self.CLUSTENM.numConfs(1)
            self.assertGreaterEqual(sampled, 1,
                                    'ClustENM run failed to give any '
                                    'conformers for generation 1')
            self.assertLessEqual(sampled, self.N_CONFS,
                                 'ClustENM run gave more conformers for '
                                 'generation 1 than it sampled')
            assert_equal(self.CLUSTENM.numConfs(), 1 + sampled,
                         'ClustENM run failed to give the conformers of every '
                         'generation')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testKeysAndLabels(self):
        """Test that every conformer is keyed and labelled by generation"""
        if prody.PY3K:
            # the starting structure, then however many generation 1 kept
            expect = [[0, 0]] + [[1, i]
                                 for i in range(self.CLUSTENM.numConfs(1))]

            keys = [list(key) for key in self.CLUSTENM.getKeys()]
            assert_equal(keys, expect,
                         'ClustENM run failed to key the conformers by '
                         'generation')
            assert_equal(self.CLUSTENM.getLabels(),
                         ['%d_%d' % tuple(key) for key in expect],
                         'ClustENM run failed to label the conformers by '
                         'generation')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testCoordsets(self):
        """Test that a run gives coordinates for every fixed atom"""
        if prody.PY3K:
            coordsets = self.CLUSTENM.getCoordsets()
            assert_equal(coordsets.shape,
                         (self.CLUSTENM.numConfs(),
                          self.CLUSTENM.getAtoms().numAtoms(), 3),
                         'ClustENM run failed to give coordinates for every '
                         'conformer and fixed atom')
            self.assertTrue(np.isfinite(coordsets).all(),
                            'ClustENM run gave coordinates that are not finite')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testPotentials(self):
        """Test that a run gives a potential energy for every conformer"""
        if prody.PY3K:
            potentials = self.CLUSTENM.getPotentials()
            assert_equal(len(potentials), self.CLUSTENM.numConfs(),
                         'ClustENM run failed to give a potential energy for '
                         'every conformer')
            self.assertTrue(np.isfinite(potentials).all(),
                            'ClustENM run gave potential energies that are '
                            'not finite')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testRerun(self):
        """Test response to running an ensemble that has been run already"""
        if prody.PY3K:
            with self.assertRaisesRegex(
                    ValueError, 'has been built',
                    msg='ClustENM run again after a run failed to say that it '
                        'has been built already'):
                self.CLUSTENM.run(maxclust=2)


class TestClustENMMinimiseOnly(unittest.TestCase):
    """Zero generations is a supported way to minimise the starting structure
    without generating any conformers along combinations of normal modes. No
    generation is clustered, so *maxclust* and *threshold* are not needed
    either, which is why this run gives neither."""

    @classmethod
    def setUpClass(cls):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        cls.cwd = os.getcwd()
        cls.tmpdir = tempfile.mkdtemp()
        os.chdir(cls.tmpdir)

        cls.ATOMS = parseDatafile('1ubi')
        cls.CLUSTENM = ClustENM()
        cls.CLUSTENM.setAtoms(cls.ATOMS)
        cls.CLUSTENM.run(n_gens=0, n_modes=2, sim=False, outlier=False)

    @classmethod
    def tearDownClass(cls):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        os.chdir(cls.cwd)
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testBuilt(self):
        """Test that a run without generations still builds the ensemble"""
        if prody.PY3K:
            self.assertTrue(self.CLUSTENM._isBuilt(),
                            'ClustENM run with 0 generations failed to build '
                            'the ensemble')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testNoGenerations(self):
        """Test that no generations are sampled"""
        if prody.PY3K:
            assert_equal(self.CLUSTENM.numGenerations(), 0,
                         'ClustENM run with 0 generations failed to give 0 '
                         'generations')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testOnlyStartingStructure(self):
        """Test that only the minimised starting structure is given"""
        if prody.PY3K:
            assert_equal(self.CLUSTENM.numConfs(), 1,
                         'ClustENM run with 0 generations failed to give only '
                         'the minimised starting structure')

            keys = [list(key) for key in self.CLUSTENM.getKeys()]
            assert_equal(keys, [[0, 0]],
                         'ClustENM run with 0 generations failed to key the '
                         'starting structure as generation 0')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testMinimised(self):
        """Test that the starting structure was minimised"""
        if prody.PY3K:
            coordsets = self.CLUSTENM.getCoordsets()
            assert_equal(coordsets.shape,
                         (1, self.CLUSTENM.getAtoms().numAtoms(), 3),
                         'ClustENM run with 0 generations failed to give '
                         'coordinates for every fixed atom')
            self.assertTrue(np.isfinite(coordsets).all(),
                            'ClustENM run with 0 generations gave coordinates '
                            'that are not finite')

            potentials = self.CLUSTENM.getPotentials()
            assert_equal(len(potentials), 1,
                         'ClustENM run with 0 generations failed to give a '
                         'potential energy for the starting structure')
            self.assertTrue(np.isfinite(potentials).all(),
                            'ClustENM run with 0 generations gave a potential '
                            'energy that is not finite')


class TestClustENMPlatform(unittest.TestCase):
    """Properties can only be given to OpenMM for a platform that was named,
    so each platform is run to make sure ClustENM asks for a combination of the
    two that OpenMM accepts. Without a platform it picks the fastest it has,
    and asking for properties as well used to be silently ignored but now
    raises, so the run without one is the case that matters most here.

    These are minimisation only runs, which is enough to build a simulation."""

    def setUp(self):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        self.cwd = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

        self.ATOMS = parseDatafile('1ubi')

    def tearDown(self):
        if not prody.PY3K or not HAVE_OPENMM:
            return

        os.chdir(self.cwd)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def runOnPlatform(self, platform):
        """Runs ClustENM over the fixed structure on *platform*, or on whichever
        one OpenMM picks if it is **None**, and checks the minimised starting
        structure that comes back."""

        clustenm = ClustENM()
        clustenm.setAtoms(self.ATOMS)

        kwargs = {} if platform is None else {'platform': platform}
        clustenm.run(n_gens=0, n_modes=2, sim=False, outlier=False, **kwargs)

        named = 'no platform' if platform is None else platform
        assert_equal(clustenm.numConfs(), 1,
                     'ClustENM run with %s failed to give the minimised '
                     'starting structure' % named)
        self.assertTrue(np.isfinite(clustenm.getCoordsets()).all(),
                        'ClustENM run with %s gave coordinates that are not '
                        'finite' % named)

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    def testNoPlatform(self):
        """Test running without naming a platform"""
        if prody.PY3K:
            self.runOnPlatform(None)

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    @unittest.skipUnless(PLATFORMS['CPU'], platformSkipMsg('CPU'))
    def testCPUPlatform(self):
        """Test running on the CPU platform, which is given a thread count"""
        if prody.PY3K:
            self.runOnPlatform('CPU')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    @unittest.skipUnless(PLATFORMS['CUDA'], platformSkipMsg('CUDA'))
    def testCUDAPlatform(self):
        """Test running on the CUDA platform, which is given a precision"""
        if prody.PY3K:
            self.runOnPlatform('CUDA')

    @dec.slow
    @unittest.skipUnless(HAVE_OPENMM, OPENMM_SKIP_MSG)
    @unittest.skipUnless(PLATFORMS['OpenCL'], platformSkipMsg('OpenCL'))
    def testOpenCLPlatform(self):
        """Test running on the OpenCL platform, which is given a precision"""
        if prody.PY3K:
            self.runOnPlatform('OpenCL')


if __name__ == '__main__':
    unittest.main()
