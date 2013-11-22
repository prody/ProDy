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

from numpy.testing import assert_array_equal

from numpy import concatenate

from prody.dynamics import calcGNM, calcANM
from prody.dynamics import extendModel, extendMode, extendVector

from prody.tests import unittest
from prody.tests.datafiles import parseDatafile

from prody import _PY3K

if not _PY3K:
    range = xrange

ATOMS = parseDatafile('1ubi').protein.copy()
ANM, NODES = calcANM(ATOMS)
GNM = calcGNM(ATOMS)[0]

EXT1D = concatenate([[i] * len(res)
                     for i, res in enumerate(ATOMS.iterResidues())])

EXT3D = concatenate([list(range(i*3, (i+1)*3)) * len(res)
                     for i, res in enumerate(ATOMS.iterResidues())])


class TestExtending(unittest.TestCase):

    def testModel3d(self):

        ext = extendModel(ANM, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), ANM._getArray()[EXT3D, :])

    def testModel1d(self):

        ext = extendModel(GNM, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), GNM._getArray()[EXT1D, :])

    def testMode3d(self):

        mode = ANM[0]
        ext = extendMode(mode, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), mode._getArray()[EXT3D] *
                                            mode.getVariance()**0.5)

    def testMode1d(self):

        mode = GNM[0]
        ext = extendMode(mode, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), mode._getArray()[EXT1D] *
                                            mode.getVariance()**0.5)

    def testVector3d(self):

        vector = ANM[0] * 1
        ext = extendVector(vector, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), vector._getArray()[EXT3D])

    def testVector1d(self):

        vector = GNM[0] * 1
        ext = extendVector(vector, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), vector._getArray()[EXT1D])

