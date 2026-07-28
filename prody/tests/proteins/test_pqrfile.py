"""This module contains unit tests for :mod:`~prody.proteins`."""

import os

import numpy as np
from numpy.testing import *

from prody.utilities import importDec
dec = importDec()

from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

class TestParsePRQ(unittest.TestCase):

    def setUp(self):
        """Set PQR file data."""

        self.pqr1 = DATA_FILES['pqrUnknown']
        self.pqr2 = DATA_FILES['pqrTranscomp']
        self.pqr3 = DATA_FILES['pqrFpocket']
        self.pqr4 = DATA_FILES['pqrPymol']

    def testExamplePQR(self):
        """Test the outcome of a simple parsing scenario with example 1 (unnamed)."""

        path = pathDatafile(self.pqr1['file'])
        ag = parsePQR(path)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePQR failed to return an AtomGroup instance')

        self.assertEqual(ag.numAtoms(), self.pqr1['atoms'],
            'parsePQR failed to parse correct number of atoms')

        self.assertEqual(ag.numCoordsets(), self.pqr1['models'],
            'parsePQR failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(ag.getTitle(),
             os.path.splitext(self.pqr1['file'])[0],
            'failed to set AtomGroup title based on filename')

    def testTranscompPQR(self):
        """Test the outcome of a simple parsing scenario with transcomp example."""

        path = pathDatafile(self.pqr2['file'])
        ag = parsePQR(path)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePQR failed to return an AtomGroup instance')

        self.assertEqual(ag.numAtoms(), self.pqr2['atoms'],
            'parsePQR failed to parse correct number of atoms')

        self.assertEqual(ag.numCoordsets(), self.pqr2['models'],
            'parsePQR failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(ag.getTitle(),
             os.path.splitext(self.pqr2['file'])[0],
            'failed to set AtomGroup title based on filename')

    def testFpocketPQR(self):
        """Test the outcome of a simple parsing scenario with example 1 (unnamed)."""

        path = pathDatafile(self.pqr3['file'])
        ag = parsePQR(path)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePQR failed to return an AtomGroup instance')

        self.assertEqual(ag.numAtoms(), self.pqr3['atoms'],
            'parsePQR failed to parse correct number of atoms')

        self.assertEqual(ag.numCoordsets(), self.pqr3['models'],
            'parsePQR failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(ag.getTitle(),
             os.path.splitext(self.pqr3['file'])[0],
            'failed to set AtomGroup title based on filename')

    def testPymolPQR(self):
        """Test the outcome of a simple parsing scenario with example 1 (unnamed)."""

        path = pathDatafile(self.pqr4['file'])
        ag = parsePQR(path)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePQR failed to return an AtomGroup instance')

        self.assertEqual(ag.numAtoms(), self.pqr4['atoms'],
            'parsePQR failed to parse correct number of atoms')

        self.assertEqual(ag.numCoordsets(), self.pqr4['models'],
            'parsePQR failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(ag.getTitle(),
             os.path.splitext(self.pqr4['file'])[0],
            'failed to set AtomGroup title based on filename')

    def testTitleArgument(self):
        """Test outcome of *title* argument."""

        path = pathDatafile(self.pqr1['file'])
        title = 'small protein'
        self.assertEqual(parsePQR(path, title=title).getTitle(),
             title, 'parsePQR failed to set user given title')

    def testSubsetArgument(self):
        """Test outcome of valid and invalid *subset* arguments."""

        path = pathDatafile(self.pqr2['file'])
        self.assertRaises(TypeError, parsePQR, path, subset=['A'])
        self.assertEqual(parsePQR(path, subset='ca').numAtoms(), 1,
                        'failed to parse correct number of "ca" atoms')
        self.assertEqual(parsePQR(path, subset='bb').numAtoms(), 2,
                        'failed to parse correct number of "bb" atoms')

    def testAgArgument(self):
        """Test outcome of 2 invalid and 2 valid *ag* arguments."""

        path = pathDatafile(self.pqr1['file'])
        self.assertRaises(TypeError, parsePQR, path, ag='AtomGroup')

        ag = prody.AtomGroup('One atom')
        ag.setCoords(np.array([[0., 0., 0.]]))
        self.assertRaises(ValueError, parsePQR, path, ag=ag)

        ag = prody.AtomGroup('Test')
        self.assertEqual(parsePQR(path, ag=ag).numAtoms(),
            self.pqr1['atoms'],
            'parsePQR failed to parse correct number of atoms')

        ag = prody.AtomGroup('Test')
        ag.setCoords(np.array([[0., 0., 0.]]*5))
        self.assertEqual(parsePQR(path, ag=ag).numAtoms(),
            self.pqr1['atoms'],
            'parsePQR failed to parse correct number of atoms')