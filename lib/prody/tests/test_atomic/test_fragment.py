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

"""This module contains unit tests for fragmenting function and methods."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from unittest import TestCase

from prody import *
from prody.tests.test_datafiles import pathDatafile

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
