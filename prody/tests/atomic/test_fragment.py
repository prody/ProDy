"""This module contains unit tests for fragmenting function and methods."""

from prody.tests import TestCase

from prody import *
from prody.tests.datafiles import pathDatafile

WHOLE = fetchPDBLigand(pathDatafile('sti'))['ideal']
WHOLE_ALL = WHOLE.all
WHOLE_COPY = WHOLE.copy()
WHOLE_NOH = WHOLE.noh
WHOLE_NOH_COPY = WHOLE_NOH.copy()
SPLIT = WHOLE.select('not name N13 N21 C46 C22')
SPLIT_COPY = SPLIT.copy()
SPLIT_NOH = WHOLE.select('noh and not name N13 N21 C46 C22')
SPLIT_NOH_COPY = SPLIT_NOH.copy()

class TestFragment(TestCase):


    def testWhole(self):

        self.assertEqual(WHOLE.numFragments(), 1)

    def testWhole2(self):

        self.assertEqual(len(list(iterFragments(WHOLE))), 1)

    def testWholeAll(self):

        self.assertEqual(len(list(iterFragments(WHOLE_ALL))), 1)

    def testWholeCopy(self):

        self.assertEqual(WHOLE_COPY.numFragments(), 1)

    def testWholeNoh(self):

        self.assertEqual(len(list(iterFragments(WHOLE_NOH))), 1)

    def testWholeNohCopy(self):

        self.assertEqual(WHOLE_NOH_COPY.numFragments(), 1)

    def testSplit(self):

        self.assertEqual(len(list(iterFragments(SPLIT))), 9)

    def testSplitCopy(self):

        self.assertEqual(SPLIT_COPY.numFragments(), 9)

    def testSplitNoh(self):

        self.assertEqual(len(list(iterFragments(SPLIT_NOH))), 5)

    def testSplitNohCopy(self):

        self.assertEqual(SPLIT_NOH_COPY.numFragments(), 5)
