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

"""This module contains unit tests for :mod:`.frame` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from code import interact
from prody.tests import TestCase

from numpy import array
from numpy.testing import assert_allclose

from prody.trajectory import Trajectory
from prody.tests.test_datafiles import parseDatafile, pathDatafile

DCD = Trajectory(pathDatafile('dcd'))
PDB = parseDatafile('multi_model_truncated', model=1)
RMSD_ALL = array([0.0, 1.380, 1.745])
RMSD_CARBON = array([0.0, 0.964, 1.148])


class TestSuperpose(TestCase):

    def setUp(self):

        DCD.setCoords(PDB.getCoords())
        DCD.reset()

    def testAll(self):

        rmsd = []
        for frame in DCD:
            frame.superpose()
            rmsd.append(frame.getRMSD())
        assert_allclose(rmsd, RMSD_ALL, atol=0.001)

    def testCarbon(self):

        DCD.setAtoms(PDB.carbon)
        rmsd = []
        for frame in DCD:
            frame.superpose()
            rmsd.append(frame.getRMSD())
        assert_allclose(rmsd, RMSD_CARBON, atol=0.001)
