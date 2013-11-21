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

"""ProDy test suite.  Usage::

  from prody import *
  prody.test()

Testing will use :mod:`nose` if it is available, otherwise it will use
:mod:`unittest`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'


from glob import glob
from os.path import abspath, split, join, relpath, splitext
from os.path import sep as dirsep
import inspect
import tempfile

try:
    import unittest2 as unittest
    from unittest2 import TestCase, skipIf, skipUnless
except ImportError:
    import unittest
    from unittest import TestCase, skipIf, skipUnless

from prody.utilities import PLATFORM
from prody import LOGGER

from prody.utilities import which
NOPRODYCMD = which('prody') is None

WINDOWS = PLATFORM == 'Windows'

try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    MATPLOTLIB = False
else:
    try:
        from matplotlib import pyplot
    except ImportError:
        MATPLOTLIB = False
    else:
        MATPLOTLIB = True

TESTDIR = abspath(split(inspect.getfile(inspect.currentframe()))[0])
TEMPDIR = tempfile.gettempdir()