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

"""This module contains unit tests for :mod:`.atommap` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from code import interact 

import unittest

from prody import *
from prody.atomic import atommap 
from prody.tests.test_datafiles import *


prody.atomic.select.DEBUG = False
prody.setVerbosity('none')

AG = parseDatafile('1ubi')

AM = AtomMap(AG, range(9, -1, -1) + range(11, 20), 
             mapping=range(10) + range(11, 20), dummies=[10])


class TestInstantiation(unittest.TestCase):
    

    def testInstantiation(self):
        
        am = AtomMap(AG, range(9, -1, -1) + [atommap.DUMMY] + range(11, 20), 
                     dummies=True)

        self.assertEqual(AM, am)


class TestSelection(unittest.TestCase):
    
    def testAllSelection(self):
        
        self.assertEqual(AM, AM.all)
    
    def testChainSelection(self):
        
        self.assertEqual(AM, AM.select('chain A _'))
