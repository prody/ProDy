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

"""This module contains unit tests for :mod:`.atommap` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from prody import *
from prody import LOGGER
from prody.atomic import atommap
from prody.tests import unittest
from prody.tests.test_datafiles import *

try:
    range = xrange
except NameError:
    pass

prody.atomic.select.DEBUG = False
LOGGER.verbosity = 'none'

AG = parseDatafile('1ubi')

AM = AtomMap(AG, list(range(9, -1, -1)) + list(range(11, 20)),
             mapping=list(range(10)) + list(range(11, 20)), dummies=[10])


class TestInstantiation(unittest.TestCase):


    def testInstantiation(self):

        am = AtomMap(AG, list(range(9, -1, -1)) + [atommap.DUMMY] +
                     list(range(11, 20)),
                     dummies=True)

        self.assertEqual(AM, am)


class TestSelection(unittest.TestCase):

    def testAllSelection(self):

        self.assertEqual(AM, AM.all)

    def testChainSelection(self):

        self.assertEqual(AM, AM.select('chain A _'))
