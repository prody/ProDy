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

"""This module contains unit tests for :mod:`prody.measure.measure` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import array, ones, arange
from numpy.testing import assert_approx_equal, assert_equal
from numpy.testing import assert_array_almost_equal

from prody.tests import unittest
from prody.tests.test_datafiles import parseDatafile, pathDatafile

from prody.trajectory import DCDFile
from prody.measure import calcDistance, buildDistMatrix
from prody.measure import calcAngle, calcPsi, calcPhi
from prody.measure import calcCenter
from prody.measure import calcMSF
from prody import LOGGER
LOGGER.verbosity = None

ATOL = 1e-5
RTOL = 0


UBI = parseDatafile('1ubi')
UBI_NTER = UBI['A', 1]
UBI_GLY10 = UBI['A', 10]
UBI_CTER = UBI['A', 76]

class TestDihedrals(unittest.TestCase):

    """Test functions that calculate dihedral angles."""

    def testCalcPsi(self):

        assert_approx_equal(14.386, calcPsi(UBI_GLY10), 2)

    def testCalcPhi(self):

        assert_approx_equal(87.723, calcPhi(UBI_GLY10), 2)

    def testNTerPsi(self):

        assert_approx_equal(153.553, calcPsi(UBI_NTER), 2)

    def testCTerPhi(self):

        assert_approx_equal(174.160, calcPhi(UBI_CTER), 2)

    def testNTerPhi(self):

        self.assertRaises(ValueError, calcPhi, (UBI_NTER))

    def testCTerPsi(self):

        self.assertRaises(ValueError, calcPsi, (UBI_CTER))


SYM_ONE = ones((5,3)) * arange(5).reshape((5,1))
SYM_TWO = ones((5,3)) * arange(5).reshape((5,1)) * 3
SYM_UC = ones(3) * 5

PBC_ONE = array([7., 0., 0.])
PBC_TWO = array([1., 0., 0.])
PBC_UC = array([5., 5., 5.])
PBC_DIST = 1.

class TestDistances(unittest.TestCase):


    def testPBCSymmetry(self):

        dist1 = calcDistance(SYM_ONE, SYM_TWO, unitcell=SYM_UC)
        dist2 = calcDistance(SYM_TWO, SYM_ONE, unitcell=SYM_UC)
        assert_equal(dist1, dist2)

    def testPBC(self):

        assert_equal(PBC_DIST, calcDistance(PBC_ONE, PBC_TWO, unitcell=PBC_UC))
        assert_equal(PBC_DIST, calcDistance(PBC_TWO, PBC_ONE, unitcell=PBC_UC))

ATOMS = parseDatafile('multi_model_truncated')
CENTERS = ATOMS.getCoordsets().mean(-2)

class TestCenter(unittest.TestCase):

    def testCenter(self):

        assert_equal(calcCenter(ATOMS), CENTERS[0])

    def testCenterWithWeights(self):

        assert_equal(calcCenter(ATOMS, weights=ones(len(ATOMS))),
                     CENTERS[0])

    def testMultiCoordsets(self):

        assert_equal(calcCenter(ATOMS.getCoordsets()), CENTERS)

    def testMultiCoordsetsWithWeights(self):

        assert_equal(calcCenter(ATOMS.getCoordsets(),
                                weights=ones(len(ATOMS))), CENTERS)

class TestMSF(unittest.TestCase):

    def testMSF(self):

        dcd = DCDFile(pathDatafile('dcd'))
        ens = parseDatafile('dcd')
        ens.superpose()
        assert_array_almost_equal(calcMSF(dcd), calcMSF(ens), 4)

    def testMSFfloat(self):

        dcd = DCDFile(pathDatafile('dcd'), astype=float)
        ens = parseDatafile('dcd', astype=float)
        ens.superpose()
        assert_array_almost_equal(calcMSF(dcd), calcMSF(ens), 10)
