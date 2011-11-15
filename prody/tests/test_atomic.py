#!/usr/bin/python
# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

"""This module contains unit tests for :mod:`~prody.atomic`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import os.path

import unittest
import numpy as np
from numpy.testing import *

from prody import *
from prody.tools import *
from test_datafiles import *

prody.changeVerbosity('none')

ATOL = 1e-5
RTOL = 0

ATOMS = parseDatafile('multi_model_truncated', subset='ca')


class TestCopying(unittest.TestCase):
    
    def TestCopyAtomGroup(self):
        
        atoms = ATOMS.copy()
        
        assert_equal(atoms.getCoordsets(), ATOMS.getCoordsets())
        for label in ATOMS.getDataLabels():
            assert_equal(atoms.getData(label), ATOMS.getData(label))
            
    def TestCopyAtomChain(self):
        
        CHAIN = ATOMS['A']
        chain = CHAIN.copy()
        
        assert_equal(chain.getCoordsets(), CHAIN.getCoordsets())
        for label in ATOMS.getDataLabels():
            print label
            assert_equal(chain.getData(label), CHAIN.getData(label),
                         'failed to copy ' + label)
