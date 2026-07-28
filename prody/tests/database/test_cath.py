"""This module contains unit tests for :mod:`.cath` module."""

from prody import LOGGER
from prody.database.cath import CATHDB, CATHCollection, CATHElement
from prody.tests.datafiles import *
from prody.tests import unittest

LOGGER.verbosity = 'none'

ATOL = 1e-5
RTOL = 0

try:
    cath = CATHDB()
except:
    cath = CATHDB(pathDatafile('cath'))

class TestCATH(unittest.TestCase):

    def setUp(self):

        self.cath = cath

    def testGetChildren(self):
        """Testing return type of :meth:`~.GNMBase.getCutoff`."""

        root = cath.root
        self.assertIsInstance(root, CATHElement,
                              'cath.root failed to return a CATHElement')

        node = root.getchildren()

        self.assertIsInstance(node, CATHCollection,
                              'cath.root.getchildren failed to return a CATHCollection')
        
        self.assertEqual(len(node), 5,
                         'cath.root.getchildren failed to return 5 elements')
