"""This module contains unit tests for :mod:`~prody.proteins`."""

from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'


class TestParsePDBHeaderOnly(unittest.TestCase):

    def setUp(self):
        self.header = parsePDB(pathDatafile('pdb2k39_truncated.pdb'),
                               header=True, model=0)

    def testHeaderType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')

    def testHeaderContent(self):
        self.assertEqual(self.header.get('classification'),
            'SIGNALING PROTEIN',
            'failed to get expected value for classification from header')
        self.assertEqual(self.header.get('experiment'), 'SOLUTION NMR',
            'failed to get expected value for experiment from header')
        self.assertEqual(self.header.get('deposition_date'), '25-APR-08',
            'failed to get expected value for deposition_date from header')
        self.assertEqual(self.header.get('identifier'), '2K39',
            'failed to get expected value for identifier from header')
        self.assertEqual(self.header.get('title'),
            'RECOGNITION DYNAMICS UP TO MICROSECONDS REVEALED FROM RDC '
            'DERIVED UBIQUITIN ENSEMBLE IN SOLUTION',
            'failed to get expected value for title from header dictionary')

    def tearDown(self):

        self.header = None





if __name__ == '__main__':
    unittest.main()
