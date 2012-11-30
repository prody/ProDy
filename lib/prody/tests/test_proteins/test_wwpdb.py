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

import unittest
import numpy as np
from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR
from prody.tests.test_datafiles import *

LOGGER.verbosity = 'none'

        
class TestFTP(unittest.TestCase):
    
    def setUp(self):
        
        self.pdb = ['1ubi', '1aar', 'arg', 1234]
        self.fns = []
        self.len = [683, 1218, None, None]
        self.fetch = fetchPDBviaFTP
        self.protocol = 'FTP'

        
    @dec.slow
    def testCompressed(self):
        
        self.fns = self.fetch(*self.pdb, folder=TEMPDIR)

        for fn, pdb, n_atoms in zip(self.fns, self.pdb, self.len):
            if fn is None:
                continue
            self.assertTrue(os.path.isfile(fn), 
                'failed to fetch PDB file via ' + self.protocol)
            atoms = parsePDB(fn)
            self.assertEqual(n_atoms, len(atoms))


    @dec.slow
    def testDecompressed(self):
        
        self.fns = self.fetch(*self.pdb, folder=TEMPDIR, compressed=False)

        for fn, pdb, n_atoms in zip(self.fns, self.pdb, self.len):
            if fn is None:
                continue
            self.assertTrue(os.path.isfile(fn), 
                'failed to fetch PDB file via ' + self.protocol)
            atoms = parsePDB(fn)
            self.assertEqual(n_atoms, len(atoms))


    def tearDown(self):

        
        for fn in self.fns:
            if fn is None:
                continue
            try:
                os.remove(fn)
            except:
                pass

class TestHTTP(TestFTP):
    
    def setUp(self):
        
        TestFTP.setUp(self)
        self.fetch = fetchPDBviaHTTP
        self.protocol = 'HTTP'

