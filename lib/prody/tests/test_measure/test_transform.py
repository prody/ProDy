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

"""This module contains unit tests for :mod:`prody.measure.transform` module.
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import zeros, ones, eye, all
from numpy.testing import assert_equal

from prody.tests import unittest
from prody.tests.test_datafiles import parseDatafile

from prody.measure import moveAtoms, wrapAtoms

UBI = parseDatafile('1ubi')


class TestMoveAtoms(unittest.TestCase):

    def testToArgument(self):

        atoms = UBI.ca

        coords = UBI.getCoords()
        center = atoms._getCoords().mean(0)
        moveAtoms(atoms, to=zeros(3), ag=True)
        assert_equal(UBI._getCoords(), coords - center)

    def testByArgument(self):

        atoms = UBI
        offset = ones(3) * 10.
        coords = atoms.getCoords()
        moveAtoms(atoms, by=offset)
        assert_equal(atoms._getCoords(), coords + offset)

    def testTransformation(self):

        atoms = UBI.ca
        coords = UBI.getCoords()
        matrix = eye(4)
        moveAtoms(atoms, by=matrix, ag=True)
        assert_equal(UBI._getCoords(), coords)


class TestWrapAtoms(unittest.TestCase):

    def testWrap(self):

        xyz = UBI.getCoords()
        xyzmin = xyz.min(0) - 0.5
        unitcell = xyz.max(0) - xyzmin + 0.5
        center = xyzmin - unitcell / 2
        wrapAtoms(UBI, unitcell, center)
        diff = xyz - UBI.getCoords()
        self.assertTrue(all(diff == unitcell))

