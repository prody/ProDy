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

from prody.tests import TestCase

from numpy.testing import assert_equal

from prody import calcOccupancies, PDBEnsemble
from . import PDBENSEMBLE, WEIGHTS, ENSEMBLE


class TestCalcOccupancies(TestCase):

    def testResults(self):
        assert_equal(calcOccupancies(PDBENSEMBLE), WEIGHTS.sum(0).flatten(),
                     'calcOccupancies failed')

    def testInvalidType(self):

        self.assertRaises(TypeError, calcOccupancies, ENSEMBLE)

    def testWeightsNone(self):

        self.assertRaises(ValueError, calcOccupancies, PDBEnsemble())

