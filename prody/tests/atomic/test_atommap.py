"""This module contains unit tests for :mod:`.atommap` module."""

from prody import *
from prody import LOGGER, SELECT
from prody.atomic import atommap
from prody.tests import unittest
from prody.tests.datafiles import *

prody.atomic.select.DEBUG = False
LOGGER.verbosity = 'none'

AG = parseDatafile('1ubi')

AM = AtomMap(AG, list(range(9, -1, -1)) + list(range(11, 20)),
             mapping=list(range(10)) + list(range(11, 20)), dummies=[10])


class TestInstantiation(unittest.TestCase):


    def testInstantiation(self):

        am = AtomMap(AG, list(range(9, -1, -1)) + [atommap.DUMMY] +
                     list(range(11, 20)),
                     dummies=True)

        self.assertEqual(AM, am)


class TestSelection(unittest.TestCase):

    def testAllSelection(self):

        self.assertEqual(AM, AM.all)

    def testChainSelection(self):

        self.assertEqual(AM, AM.select('chain A _'))

    def testNameSelection(self):

        self.assertEqual((AM.getNames() == 'CA').sum(),
                          len(AM.select('name CA')))

    def testGetBoolArray(self):

        self.assertEqual(len(AM), len(SELECT.getBoolArray(AM, 'name CA')))
