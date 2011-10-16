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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import sys
import os.path
import unittest
import test_proteins
import test_select
import test_ensemble

def test(verbosity=2, descriptions=True, stream=sys.stderr):
    testrunner = unittest.TextTestRunner(stream, descriptions, verbosity)
    for module in [test_select, test_ensemble, test_proteins]:
        testrunner.run(unittest.defaultTestLoader.loadTestsFromModule(module))
    
if __name__ == '__main__':
    test()
