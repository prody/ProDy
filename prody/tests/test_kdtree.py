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

"""This module contains unit tests for :mod:`~prody.KDTree` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import unittest
import numpy as np
from numpy.testing import *
from prody.KDTree import KDTree
ATOL = 1e-5
RTOL = 0

class TestKDTree(unittest.TestCase):
    
    
    def setUp(self):
        
        self.ndim = 3
        self.kdtree = KDTree(self.ndim)
        self.coords = np.tile(np.arange(10), (3,1)).T.astype(float)
        self.kdtree.set_coords(self.coords)
        
        
    def testSearch(self):
        
        kdtree = self.kdtree 
        coords = self.coords
        kdtree.search(np.array([0,0,0]), 10000)
        indices = kdtree.get_indices()
        dist = (coords**2).sum(1)**0.5
        radii = kdtree.get_radii()
        assert_allclose(radii, dist[indices], rtol=RTOL, atol=ATOL,
                        err_msg='KDTree search failed')

    def testAllSearch(self):
        kdtree = self.kdtree 
        coords = self.coords
        kdtree.all_search(1.75)
        radii = kdtree.all_get_radii()
        indices = kdtree.all_get_indices()
        self.assertEqual(len(indices), 9, 'KDTree all search failed')
        for pair, radius in zip(indices, radii):
            x, y = coords[pair] 
            assert_allclose(((x - y)**2).sum()**0.5, radius, 
                            rtol=RTOL, atol=ATOL,
                            err_msg='KDTree all search failed')

        
