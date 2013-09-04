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

from prody.tests import TestCase

from numpy.random import random
from numpy.testing import assert_equal

from prody.utilities import checkCoords, checkTypes

COORDS = random((10, 3))*10
COORDSET = random((2, 10, 3))*10

class TestCheckCoords(TestCase):

    def testInvalidCoords(self):

        self.assertRaises(TypeError, checkCoords, [None])

    def testCoords(self):

        self.assertTrue(checkCoords(COORDS))

    def testCoordset(self):

        self.assertTrue(checkCoords(COORDSET, csets=True))


    def testCoordsetNatoms(self):

        self.assertRaises(ValueError, checkCoords, COORDSET, csets=True,
                          natoms=20)


class testCheckTypes(TestCase):

    def testCorrectMonotypeOneArg(self):

        self.assertTrue(checkTypes({'i': 1}, i=int))

    def testCorrectMonotypeTwoArgs(self):

        self.assertTrue(checkTypes({'i': 1, 'n': 10}, i=int, n=int))

    def testCorrectMultitypeOneArg(self):

        self.assertTrue(checkTypes({'i': 1.}, i=(int, float)))

    def testCorrectMonotypeTwoArgs(self):

        self.assertTrue(checkTypes({'i': 1, 'n': 10.}, i=(int, float),
                                                       n=(int, float)))

    def testWrongMonotypeOneArg(self):

        self.assertRaises(TypeError, checkTypes, {'i': 1.0}, i=int)

    def testWrongMonotypeTwoArgs(self):

        self.assertRaises(TypeError, checkTypes, {'i': 1, 'n': 10},
                                                  i=int, n=int)

    def testWrongMultitypeOneArg(self):

        self.assertRaises(TypeError, checkTypes, {'i': '1.'},
                                      i=(int, float))

    def testWrongMonotypeTwoArgs(self):

        self.assertRaises(TypeError, checkTypes, {'i': 1, 'n': '10.'},
                            i=(int, float), n=(int, float))
