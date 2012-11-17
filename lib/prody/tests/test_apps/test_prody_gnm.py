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

from os import remove
import shlex
from os.path import isfile, join, split, splitext
from unittest import TestCase, skipIf

from numpy.testing import *

from prody.tests.test_datafiles import TEMPDIR, pathDatafile

from prody.apps import prody_parser 

from . import NOPRODYCMD

class TestGNMCommand(TestCase): 
    
    def setUp(self):

        self.command = ('gnm -e -r -o {outdir} -v -z -t all '
                        '-f %8g -d , -x .dat '
                        '-R -Q '
                        '-F png -D 120 -W 5 -H 4 ').format(outdir=TEMPDIR)
        self.suffixes = [
            '_gnm_cc.png',
            '_gnm.gnm.npz',
            '_gnm_covariance.dat',
            '_gnm_cross-correlations.dat',
            '_gnm_evalues.dat',
            '_gnm_evectors.dat',
            '_gnm_sf.png',
            '_gnm_extended_all.nmd',
            '_gnm.nmd',                     
        ]
        
        self.tearDown()

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    def testGNMCommand(self):
            
        pdb = pathDatafile('multi_model_truncated')
        command = self.command + pdb
        prefix = splitext(split(pdb)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')
            
            
    def tearDown(self):
        
        for pdb in [pathDatafile('multi_model_truncated')]:
            prefix = splitext(split(pdb)[1])[0]
            for suffix in self.suffixes:
                fn = join(TEMPDIR, prefix + suffix)
                if isfile(fn): remove(fn)
