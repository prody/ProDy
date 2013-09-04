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

"""This module contains unit tests for fragmenting function and methods."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from prody.tests import TestCase

from numpy import arange
from numpy.random import shuffle

from prody import *
from prody.tests.test_datafiles import pathDatafile

AG = prody.parsePDB(pathDatafile('pdb3mht.pdb'))
SHUFFLED = arange(len(AG))
SHUFFLED = AtomMap(AG, SHUFFLED)



class TestShuffled(TestCase):


    def testCA(self):

        self.assertEqual(AG.ca.numAtoms(), SHUFFLED.ca.numAtoms())

    def testProtein(self):

        self.assertEqual(AG.protein.numAtoms(), SHUFFLED.protein.numAtoms())
