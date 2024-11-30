"""This module contains unit tests for :mod:`~prody.proteins`."""
import numpy as np
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

class TestParsePDBHeaderCRYST1andSCALE(unittest.TestCase):

    def setUp(self):
        self.header = parsePDB(pathDatafile('1pwc.pdb'),
                               header=True, model=0)

    def testHeaderType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')

    def testHeaderContent(self):
        # CRYST1
        cryst1 = self.header.get('CRYST1', None)
        self.assertIsInstance(cryst1, dict,
            'CRYST1 type is incorrect in the header')
        self.assertEqual(cryst1.get('Z value'), 4,
            'failed to get expected value for Z value')
        self.assertEqual(cryst1.get('spaceGroup'), 'P 21 21 21',
            'failed to get expected value for spaceGroup from header')
        self.assertEqual(cryst1.get('cellAngles', None), (90.0, 90.0, 90.0),
            'failed to get expected value for cellAngles from header')
        self.assertEqual(cryst1.get('cellLength'), (50.9, 66.7, 99.6),
            'failed to get expected value for cellLength from header')

        # SCALE
        scale = self.header.get('SCALE', None)
        self.assertIsInstance(scale, dict,
            'SCALE type is incorrect in the header')
        self.assertIsInstance(scale.get('ctof', None), np.ndarray)
        self.assertEqual(scale.get('ctof').shape, (4,4), 
            'failed to get expected value for scale.ctof value')
        _ctof = np.array([[0.019646, 0.      , 0.      , 0.      ],
                          [0.      , 0.014993, 0.      , 0.      ],
                          [0.      , 0.      , 0.01004 , 0.      ],
                          [0.      , 0.      , 0.      , 1.      ]])
              
        assert_array_equal(scale.get('ctof'), _ctof, 
            'failed to get expected value for sacle.ctof value')

        self.assertIsInstance(scale.get('ftoc', None), np.ndarray)
        self.assertEqual(scale.get('ftoc').shape, (4,4), 
            'failed to get expected value for scale.ftoc value')
        _ftoc = np.array([[50.90094676,  0.        ,  0.        ,  0.        ],
                          [ 0.        , 66.6977923 ,  0.        ,  0.        ],
                          [ 0.        ,  0.        , 99.60159363,  0.        ],
                          [ 0.        ,  0.        ,  0.        ,  1.        ]])
        assert_array_almost_equal(scale.get('ftoc'), _ftoc, decimal=7,
            err_msg='failed to get expected value for scale.ftoc value')

        # REMARK 365
        self.assertIsInstance(self.header.get('missing_residues', None), list,
            'missing_residues type is incorrect in the header')
        mis_res = self.header.get('missing_residues')
        
        self.assertEqual(len(mis_res), 4, 
                         'failed to get expected value for missing_residues value')        
        self.assertEqual( mis_res[0], ('A', 'ALA', 1, '', ''),
            'failed to get expected value for missing_residues[0]')
        self.assertEqual( mis_res[1], ('A', 'ASP', 2, '', ''),
            'failed to get expected value for missing_residues[1]')
        self.assertEqual( mis_res[2], ('A', 'THR', 348, '', ''),
            'failed to get expected value for missing_residues[2]')
        self.assertEqual( mis_res[3], ('A', 'THR', 349, '', ''),
            'failed to get expected value for missing_residues[3]')
        
    def tearDown(self):

        self.header = None

class TestParsePDBHeaderMissingAtoms(unittest.TestCase):

    def setUp(self):
        self.header = parsePDB(pathDatafile('pdb3hsy.pdb'),
                               header=True, model=0)

    def testHeaderType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')

    def testHeaderContent(self):
        # CRYST1
        self.assertIsInstance(self.header.get('missing_atoms', None), list,
            'missing_atoms type is incorrect in the header')
        mis_atoms = self.header.get('missing_atoms')
        
        self.assertEqual(len(mis_atoms), 5, 
                         'failed to get expected value for missing_atoms value')
        self.assertEqual( mis_atoms[0],
            ('A', 'ASN', 305, '', ['CG', 'OD1', 'ND2']), 
            'failed to get expected value for missing_atoms[0]')
        self.assertEqual( mis_atoms[1],
            ('B', 'GLN', 296, '', ['CG', 'CD', 'OE1', 'NE2']),
            'failed to get expected value for missing_atoms[1]')
        self.assertEqual( mis_atoms[2],
            ('B', 'ARG', 297, '', ['CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2']),
            'failed to get expected value for missing_atoms[2]')
        self.assertEqual( mis_atoms[3],
            ('B', 'ILE', 298, '', ['CG1', 'CG2', 'CD1']),
            'failed to get expected value for missing_atoms[3]')
        self.assertEqual( mis_atoms[4],
            ('B', 'GLU', 299, '', ['CG', 'CD', 'OE1', 'OE2']),
            'failed to get expected value for missing_atoms[4]')
        
    def tearDown(self):

        self.header = None


if __name__ == '__main__':
    unittest.main()
