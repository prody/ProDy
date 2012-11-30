#!/usr/bin/python
# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module contains unit tests for :mod:`~prody.proteins`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os

import unittest
import numpy as np
from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR
from prody.tests.test_datafiles import *

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
