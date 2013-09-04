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
from prody.tests import TestCase, skipIf, skipUnless

from numpy.testing import *

from prody.tests.test_datafiles import TEMPDIR, pathDatafile

from prody.apps import prody_parser

from prody.tests import MATPLOTLIB, NOPRODYCMD, WINDOWS

class TestPCACommand(TestCase):

    def setUp(self):

        self.command = ('pca --pdb {pdb} '
                        '-e -r -o {outdir} -v -z -t all -j '
                        '-f %8g -d , -x .dat '
                        '-R -Q -J 1,2 '
                        '-F png -D 120 -W 5 -H 4 ').format(outdir=TEMPDIR,
                        pdb=pathDatafile('multi_model_truncated'))
        self.suffixes = [
            '_pca_cc.png',
            '_pca.pca.npz',
            '_pca_covariance.dat',
            '_pca_cross-correlations.dat',
            '_pca_proj_1_2.png',
            '_pca_evalues.dat',
            '_pca_proj.dat',
            '_pca_evectors.dat',
            '_pca_sf.png',
            '_pca_extended_all.nmd',
            '_pca.nmd',
        ]

        self.tearDown()

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipUnless(MATPLOTLIB, 'matplotlib not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testPCACommandDCD(self):

        dcd = pathDatafile('dcd')
        command = self.command + dcd
        prefix = splitext(split(dcd)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipUnless(MATPLOTLIB, 'matplotlib not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testPCACommandPDB(self):

        dcd = pathDatafile('multi_model_truncated')
        command = self.command + dcd
        prefix = splitext(split(dcd)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')


    def tearDown(self):

        for dcd in [pathDatafile('multi_model_truncated'),
                    pathDatafile('dcd')]:
            prefix = splitext(split(dcd)[1])[0]
            for suffix in self.suffixes:
                fn = join(TEMPDIR, prefix + suffix)
                if isfile(fn): remove(fn)
