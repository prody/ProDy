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

"""This module contains unit tests for :mod:`prody.measure.measure` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import unittest
import numpy as np
from numpy.testing import assert_almost_equal

from prody.tests.test_datafiles import parseDatafile

from prody.measure import calcDistance, buildDistMatrix
from prody.measure import calcAngle, calcPsi, calcPhi 
ATOL = 1e-5
RTOL = 0


UBI = parseDatafile('1ubi')
UBI_NTER = UBI['A', 1]
UBI_GLY10 = UBI['A', 10]
UBI_CTER = UBI['A', 76]

class TestDihedrals(unittest.TestCase):
    
    """Test functions that calculate unittests."""
    
    
    def testCalcPsi(self):
        
        assert_almost_equal(14.386, calcPsi(UBI_GLY10), 2)
        
    def testCalcPhi(self):
        
        assert_almost_equal(87.723, calcPhi(UBI_GLY10), 2)

    def testNTerPsi(self):
        
        assert_almost_equal(153.553, calcPsi(UBI_NTER), 2)
        
    def testCTerPhi(self):
        
        assert_almost_equal(174.160, calcPhi(UBI_CTER), 2)
    
    def testNTerPhi(self):
        
        self.assertRaises(ValueError, calcPhi, (UBI_NTER))
        
    def testCTerPsi(self):
        
        self.assertRaises(ValueError, calcPsi, (UBI_CTER))
