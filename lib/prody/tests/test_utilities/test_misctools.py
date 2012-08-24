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

from unittest import TestCase

from prody.utilities import rangeString


class TestRangeString(TestCase):
    
    def testContinuous(self):
        
        self.assertEqual(rangeString(range(10)), '0 to 9')
        
    def testNegative(self):
        
        self.assertEqual(rangeString(range(-5, 10), pos=False), '-5 to 9')

    def testGapped(self):
        
        self.assertEqual(rangeString(range(-5, 10) + range(15, 20) + 
                                     range(25, 30), pos=False), 
                                     '-5 to 9 15 to 19 25 to 29')

    def testRepeated(self):
        
        self.assertEqual(rangeString(range(10, 20) + range(15, 20) + 
                                     range(30)), '0 to 29')
