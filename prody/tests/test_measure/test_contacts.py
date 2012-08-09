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
from numpy.testing import assert_array_equal, assert_equal

from prody.tests.test_datafiles import parseDatafile, pathDatafile

from prody.measure import Contacts, findNeighbors, iterNeighbors
from prody.measure import buildDistMatrix, calcDistance


UBI = parseDatafile('1ubi')
UBI_XYZ = UBI._getCoords()
UBI_MIN = UBI_XYZ.min(0)
UBI_MAX = UBI_XYZ.max(0)
UBI_UC = UBI_MAX - UBI_MIN
N0P = [-1, 0, 1]
UBI_EDGES = UBI_MIN + UBI_UC * array([(i, j, k) for i in N0P 
                                                for j in N0P for k in N0P])
UBI_RADIUS = 20 

UBI_CONTACTS = Contacts(UBI_XYZ)
UBI_CONTACTS_PBC = Contacts(UBI_XYZ, UBI_UC)


class TestContacts(unittest.TestCase):
    
    """Test functions that calculate dihedral angles."""
    
    def testPBCvsNONE(self):
        
        wout_pbc = unique(concatenate([UBI_CONTACTS(UBI_RADIUS, UBI_EDGES)]))
        with_pbc = UBI_CONTACTS_PBC(UBI_RADIUS, UBI_MIN)
        assert_array_equal(wout_pbc, with_pbc)
        
    def testMINvsMAX(self):
        
        assert_array_equal(UBI_CONTACTS_PBC(UBI_RADIUS, UBI_MIN), 
                           UBI_CONTACTS_PBC(UBI_RADIUS, UBI_MAX))
        

UCA = UBI.ca.copy()[::2].copy()
UCA_XYZ = UCA._getCoords()
UCA_MIN = UCA_XYZ.min(0)
UCA_MAX = UCA_XYZ.max(0)
UCA_UC = UCA_MAX - UCA_MIN
N0P = [-1, 0, 1]
UCA_EDGES = UCA_MIN + UCA_UC * array([(i, j, k) for i in N0P 
                                                for j in N0P for k in N0P])
UCA_RADIUS = 20 

UCA_CONTACTS = Contacts(UCA_XYZ)
UCA_CONTACTS_PBC = Contacts(UCA_XYZ, UCA_UC)
        
class TestNeighbors(unittest.TestCase):
    
    def testPBC(self):
        
        neighbors = findNeighbors(UCA_XYZ, UCA_RADIUS, unitcell=UCA_UC)
        n_neighbors = (buildDistMatrix(UCA_XYZ, unitcell=UCA_UC, 
                                       format='arr') <= UCA_RADIUS).sum()
        assert_equal(len(neighbors), n_neighbors)
        
    def testNoPBC(self):
        
        neighbors = findNeighbors(UCA_XYZ, UCA_RADIUS)
        n_neighbors = (buildDistMatrix(UCA_XYZ, 
                                       format='arr') <= UCA_RADIUS).sum()
        assert_equal(len(neighbors), n_neighbors)

    def testCoordArgumentSwitching(self):
        
        dist = 12.
        neighbors1 = findNeighbors(UCA_XYZ, dist, UCA_XYZ[1])
        neighbors2 = findNeighbors(UCA_XYZ[1], dist, UCA_XYZ)
        n_neighbors = (calcDistance(UCA_XYZ, UCA_XYZ[1]) <= dist).sum()
        assert_equal(len(neighbors1), n_neighbors)
        neighbors1.sort()
        neighbors2.sort()
        assert_array_equal(array(neighbors1)[:,-1], 
                           array(neighbors2)[:,-1])
        
    def testAtomicArgumentSwitching(self):
        
        dist = 12.
        neighbors1 = [(a.getIndex(), d) 
                      for a, b, d in iterNeighbors(UCA, dist, UCA[1])]
        neighbors2 = [(b.getIndex(), d)
                      for a, b, d in iterNeighbors(UCA[1], dist, UCA)]
        n_neighbors = (calcDistance(UCA, UCA[1]) <= dist).sum()
        assert_equal(len(neighbors1), n_neighbors)
        neighbors1.sort()
        neighbors2.sort()
        self.assertEqual(neighbors1, neighbors2)
        
    def testPBCCoordArgumentSwitching(self):
        
        dist = 12.
        neighbors1 = findNeighbors(UCA_XYZ, dist, UCA_XYZ[1], unitcell=UCA_UC)
        neighbors2 = findNeighbors(UCA_XYZ[1], dist, UCA_XYZ, unitcell=UCA_UC)
        n_neighbors = (calcDistance(UCA_XYZ, UCA_XYZ[1], unitcell=UCA_UC) <= 
                       dist).sum()
        assert_equal(len(neighbors1), n_neighbors)
        neighbors1.sort()
        neighbors2.sort()
        assert_array_equal(array(neighbors1)[:,-1], 
                           array(neighbors2)[:,-1])
        
    def testPBCAtomicArgumentSwitching(self):
        
        dist = 12.
        neighbors1 = [(a.getIndex(), d) for a, b, d in 
                            iterNeighbors(UCA, dist, UCA[1], unitcell=UCA_UC)]
        neighbors2 = [(b.getIndex(), d) for a, b, d in 
                            iterNeighbors(UCA[1], dist, UCA, unitcell=UCA_UC)]
        n_neighbors = (calcDistance(UCA, UCA[1], unitcell=UCA_UC) <= 
                       dist).sum()
        assert_equal(len(neighbors1), n_neighbors)
        neighbors1.sort()
        neighbors2.sort()
        self.assertEqual(neighbors1, neighbors2)
