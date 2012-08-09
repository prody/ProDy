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

"""This module contains unit tests for :mod:`~prody.ensemble`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from os.path import join
from unittest import TestCase

from numpy.testing import assert_equal, assert_allclose

from prody import DCDFile, writeDCD, parseDCD

from prody.tests import TEMPDIR
from prody.tests.test_ensemble import ALLATOMS, ENSEMBLE, RTOL, ATOL, DCD

class TestDCDFile(TestCase):
    
    def setUp(self):
        
        self.dcd = join(TEMPDIR, 'temp.dcd')
    
    def testWriteDCD(self):
        dcd = writeDCD(self.dcd, ALLATOMS)
        self.assertEqual(dcd, self.dcd, 'failed to write DCD file')
        
    def testParseDCD(self):
        e = parseDCD(writeDCD(self.dcd, ALLATOMS))
        assert_equal(e._getCoordsets(), DCD._getCoordsets(),
                     err_msg='failed to parse DCD file correctly')

    def testWrite(self):
        dcd = DCDFile(self.dcd, 'w')
        dcd.write(ENSEMBLE.getCoordsets())
        dcd.close()
        e = parseDCD(self.dcd)
        assert_allclose(e._getCoordsets(), ENSEMBLE._getCoordsets(),
                        rtol=RTOL, atol=ATOL,
                        err_msg='failed to parse DCD file correctly')

    def testWriteModeAppend(self):
        dcd = DCDFile(writeDCD(self.dcd, ENSEMBLE), 'a')
        dcd.write(ENSEMBLE.getCoordsets())
        dcd.close()
        e = parseDCD(self.dcd)
        n_csets = len(ENSEMBLE)
        coordsets = e._getCoordsets()
        assert_equal(coordsets, coordsets, 
                     'failed to parse DCD file correctly')
        assert_allclose(coordsets[:n_csets], ENSEMBLE._getCoordsets(),
                        rtol=RTOL, atol=ATOL,
                        err_msg='failed to parse DCD file correctly')
