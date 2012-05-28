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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import unittest
from numpy import array, concatenate, unique
from numpy.testing import assert_array_equal

from prody.tests.test_datafiles import parseDatafile, getDatafilePath

from prody.measure import Contacts


UBI = parseDatafile('1ubi')
XYZ = UBI._getCoords()
MIN = XYZ.min(0)
MAX = XYZ.max(0)
UC = MAX - MIN
N0P = [-1, 0, 1]
EDGES = MIN + UC * array([(i, j, k) for i in N0P for j in N0P for k in N0P])
RADIUS = 20 

CONTACTS = Contacts(XYZ)
CONTACTS_PBC = Contacts(XYZ, UC)


class TestContacts(unittest.TestCase):
    
    """Test functions that calculate dihedral angles."""
    
    def testPBCvsNONE(self):
        
        wout_pbc = unique(concatenate([CONTACTS(RADIUS, EDGES)]))
        with_pbc = CONTACTS_PBC(RADIUS, MIN)
        assert_array_equal(wout_pbc, with_pbc)
        
    def testMINvsMAX(self):
        
        assert_array_equal(CONTACTS_PBC(RADIUS, MIN), 
                           CONTACTS_PBC(RADIUS, MAX))
        
        
