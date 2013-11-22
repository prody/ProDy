"""This module contains unit tests for fragmenting function and methods."""

from prody.tests import TestCase

from numpy import arange
from numpy.random import shuffle

from prody import *
from prody.tests.datafiles import pathDatafile

AG = prody.parsePDB(pathDatafile('pdb3mht.pdb'))
SHUFFLED = arange(len(AG))
shuffle(SHUFFLED)
SHUFFLED = AtomMap(AG, SHUFFLED)

RTER = prody.parsePDB(pathDatafile('pdbRTER.pdb'))


class TestShuffled(TestCase):

    def testCA(self):

        self.assertEqual(AG.ca.numAtoms(), SHUFFLED.ca.numAtoms())

    def testProtein(self):

        self.assertEqual(AG.protein.numAtoms(), SHUFFLED.protein.numAtoms())


class TestTerRecord(TestCase):

    def testNumResidues(self):

        self.assertEqual(RTER.getHierView(ter=True).numResidues(), 9)

    def testResidueIndexing(self):

        self.assertEqual(len(RTER.getHierView(ter=True)['A', 864]), 2)

    def testSelectionResidueIndexing(self):

        residues = RTER[:32].getHierView(ter=True)['A', 864]
        self.assertEqual(len(residues), 2)
        self.assertEqual((residues[0] + residues[1]).numAtoms(), 5)

    def testSelectionResidueIndexing2(self):

        self.assertEqual(len(RTER[20:].getHierView(ter=True)['A', 866]), 3)