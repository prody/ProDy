from os import remove
import shlex
from os.path import isfile, join, split, splitext
from prody.tests import TestCase, skipIf, skipUnless

from numpy.testing import *

from prody.utilities import importDec
dec = importDec()

from prody import parsePDB, DCDFile, parseDCD

from prody.tests.datafiles import TEMPDIR, pathDatafile

from prody.apps import prody_parser

from prody.tests import MATPLOTLIB, NOPRODYCMD, WINDOWS

class TestCatdcdCommand(TestCase):

    @classmethod
    def setUpClass(cls):
        # Optimization: Perform the heavy PDB parsing only ONCE.
        cls.dcdpath = pathDatafile('dcd')
        cls.pdbpath = pathDatafile('multi_model_truncated')
        cls.ag = parsePDB(cls.pdbpath, model=1)

    def setUp(self):
        # Safety: Open a fresh DCD handle for every test. 
        self.dcd = DCDFile(self.dcdpath)
        
        self.output = join(TEMPDIR, 'test_prody_catdcd.dcd')
        self.command = 'catdcd -o ' + self.output

        if isfile(self.output):
            remove(self.output)

    def tearDown(self):
        # CRITICAL FIX: Close the main DCD handle before teardown
        if hasattr(self, 'dcd') and self.dcd is not None:
            self.dcd.close()

        # Remove output only after ensuring input is closed
        if isfile(self.output): 
            try:
                remove(self.output)
            except OSError:
                pass

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSimpleConcat(self):

        command = self.command + ' {0:s} {0:s} {0:s}'.format(self.dcdpath)

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        coords = self.dcd[:]._getCoordsets()
        
        # parseDCD returns an Ensemble, which auto-closes the file.
        # We do not need to call .close() on 'concat'.
        concat = parseDCD(self.output)._getCoordsets()
        
        assert_equal(coords, concat[:3])
        assert_equal(coords, concat[3:6])
        assert_equal(coords, concat[6:])

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSelectConcat(self):

        command = self.command + ' -s ca --pdb {1:s} {0:s} {0:s}'.format(
                                self.dcdpath, self.pdbpath)

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        select = self.ag.ca
        assert_equal(select.numAtoms(), 10)

        coords = self.dcd[:]
        coords.setAtoms(select)
        coords = coords._getCoordsets()

        concat = parseDCD(self.output)
        assert_equal(concat.numAtoms(), select.numAtoms())
        concat = concat._getCoordsets()

        assert_equal(select.numAtoms(), coords.shape[1])
        assert_equal(select.numAtoms(), concat.shape[1])
        assert_equal(coords, concat[:3])
        assert_equal(coords, concat[3:])

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testAlignConcat(self):

        command = self.command + ' --align ca --pdb {1:s} {0:s} {0:s}'.format(
                                self.dcdpath, self.pdbpath)

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        select = self.ag.ca

        coords = self.dcd[:]
        
        concat = parseDCD(self.output)
        assert_equal(concat.numAtoms(), coords.numAtoms())

        coords.setCoords(self.ag.getCoords())
        coords.setAtoms(select)
        coords.superpose()
        coords.setAtoms(None)
        coords = coords._getCoordsets()

        concat = concat._getCoordsets()

        assert_equal(coords, concat[:3])
        assert_equal(coords, concat[3:])

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSelectException(self):

        command = self.command + ' -s ca {0:s} {0:s}'.format(
                                self.dcdpath)
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(ValueError, namespace.func, namespace)

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testAlignException(self):

        command = self.command + ' --align ca {0:s} {0:s}'.format(
                                self.dcdpath)
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(ValueError, namespace.func, namespace)

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testIOException(self):

        command = self.command + ' {0:s} {0:s}'.format('deneme.dcd')
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(IOError, namespace.func, namespace)

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSelectException2(self):

        command = self.command + ' -s None {0:s} {0:s}'.format(self.dcdpath)
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(ValueError, namespace.func, namespace)
