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

"""This module contains unit tests for :mod:`~prody.proteins`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os

import numpy as np
from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

class TestFetchPDB(unittest.TestCase):

    """Test :func:`~.fetchPDB` function."""

    @dec.slow
    def setUp(self):
        """Instantiate a list for storing downloaded file names."""

        self.filenames = []

    @dec.slow
    def testFetchingSingleFile(self):
        """Test the outcome of fetching a single PDB file."""

        fn = fetchPDB('1p38', folder=TEMPDIR, copy=True)
        self.assertTrue(os.path.isfile(fn),
            'fetching a single PDB file failed')
        self.filenames.append(fn)

    @dec.slow
    def testFetchingMultipleFiles(self):
        """Test the outcome of fetching multiple PDB files."""

        fns = fetchPDB(['1p38', '1r39'], folder=TEMPDIR, copy=True)
        self.assertIsInstance(fns, list,
            'failed to return a list of filenames')
        self.assertTrue(all([os.path.isfile(fn) for fn in fns]),
            'fetching multiple PDB files failed')
        self.filenames.extend(fns)

    @dec.slow
    def testCompressedArgument(self):
        """Test decompressing fetched PDB files."""

        fns = fetchPDB(['1p38', '1r39'], folder=TEMPDIR, compressed=False)
        self.assertTrue(all([os.path.isfile(fn) for fn in fns]),
            'fetching decompressed PDB files failed')
        self.assertTrue(all([os.path.splitext(fn)[1] != '.gz' for fn in fns]),
            'decompressing PDB files failed')
        self.filenames.extend(fns)


    def testInvalidPDBIdentifier(self):
        """Test outcome of passing invalid PDB identifiers."""

        self.assertIsNone(fetchPDB('XXXXX', folder=TEMPDIR),
            'failed to return None for an invalid PDB identifier')

        self.assertFalse(all(fetchPDB(
                    ['XXXXX', '654654', '-/-*/+', ''], folder=TEMPDIR)),
            'failed to return None for invalid PDB identifiers')
    @dec.slow
    def tearDown(self):
        """Remove downloaded files from disk."""

        for fn in self.filenames:
            try:
                os.remove(fn)
            except:
                pass
