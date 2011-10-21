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

"""This module contains unit tests for :mod:`~prody.ensemble` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import unittest
import numpy as np
from numpy.testing import *

from prody import *
from test_datafiles import *

prody.changeVerbosity('none')

ATOL = 1e-5
RTOL = 0


ATOMS = parseDatafile('1ubi')

ANM_HESSIAN = parseDatafile('anm1ubi_hessian', symmetric=True)
ANM_EVALUES = parseDatafile('anm1ubi_evalues')[:,1].flatten()
ANM_EVECTORS = parseDatafile('anm1ubi_vectors')[:,1:]

GNM_KIRCHHOFF = parseDatafile('gnm1ubi_kirchhoff', symmetric=True, skiprows=1)
GNM_EVALUES = parseDatafile('gnm1ubi_evalues')[:,1].flatten()
GNM_EVECTORS = parseDatafile('gnm1ubi_vectors', usecols=range(3,23))

 
anm = ANM() 
anm.buildHessian(ATOMS)
anm.calcModes(n_modes=30, zeros=True)

gnm = GNM() 
gnm.buildKirchhoff(ATOMS)
gnm.calcModes(n_modes=None, zeros=True)

class TestANM(unittest.TestCase):

    @unittest.skipIf(NONOSE, NONOSE_MSG)    
    def testEigenvalues(self):
        
        assert_allclose(anm.getEigenvalues(), ANM_EVALUES, 
                        rtol=RTOL, atol=ATOL*10,
                        err_msg='failed when comparing eigenvalues')


    @unittest.skipIf(NONOSE, NONOSE_MSG)
    def testEigenvectors(self):
        _temp = np.abs(np.dot(anm[6:26].getEigenvectors().T, ANM_EVECTORS))
        assert_allclose(_temp, np.eye(20), rtol=RTOL, atol=ATOL,
                        err_msg='comparing eigenvectors')

    @unittest.skipIf(NONOSE, NONOSE_MSG)
    def testHessian(self):
        assert_allclose(anm.getHessian(), ANM_HESSIAN, rtol=0, atol=ATOL,
                        err_msg='comparing Hessian')

    @unittest.skipIf(NONOSE, NONOSE_MSG)
    def testVariances(self):
        assert_allclose(anm[6:].getVariances(), 1/ANM_EVALUES[6:], 
                        rtol=0, atol=ATOL*100,
                        err_msg='comparing variances')

class TestGNM(unittest.TestCase):
    
    @unittest.skipIf(NONOSE, NONOSE_MSG)
    def testEigenvalues(self):
        assert_allclose(gnm[:21].getEigenvalues(), GNM_EVALUES[:21], 
                        rtol=RTOL, atol=ATOL*100,
                        err_msg='failed when comparing gnm slow eigenvalues')
                        
        assert_allclose(gnm[-21:].getEigenvalues(), GNM_EVALUES[21:], 
                        rtol=RTOL, atol=ATOL*100,
                        err_msg='failed when comparing gnm fast eigenvalues')

    @unittest.skipIf(NONOSE, NONOSE_MSG)
    def testEigenvectors(self):
        _temp = np.abs(np.dot(gnm[1:21].getEigenvectors().T, GNM_EVECTORS))
        assert_allclose(_temp, np.eye(20), rtol=RTOL, atol=ATOL*10,
                       err_msg='comparing gnm eigenvectors')

    @unittest.skipIf(NONOSE, NONOSE_MSG)
    def testKirchhoff(self):
        assert_allclose(gnm.getKirchhoff(), GNM_KIRCHHOFF, 
                        rtol=0, atol=ATOL,
                        err_msg='comparing Kirchhoff')

if __name__ == '__main__':
    unittest.main()
