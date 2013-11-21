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

from prody import parsePDB, DCDFile, parseDCD

from prody.tests.datafiles import TEMPDIR, pathDatafile

from prody.apps import prody_parser

from prody.tests import MATPLOTLIB, NOPRODYCMD, WINDOWS

class TestCatdcdCommand(TestCase):

    def setUp(self):

        self.output = join(TEMPDIR, 'test_prody_catdcd.dcd')

        self.dcdpath = pathDatafile('dcd')
        self.pdbpath = pathDatafile('multi_model_truncated')

        self.dcd = DCDFile(self.dcdpath)
        self.ag = parsePDB(self.pdbpath, model=1)

        self.command = 'catdcd -o ' + self.output

        self.tearDown()

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSimpleConcat(self):

        command = self.command + ' {0:s} {0:s} {0:s}'.format(self.dcdpath)

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        coords = self.dcd[:]._getCoordsets()
        concat = parseDCD(self.output)._getCoordsets()
        assert_equal(coords, concat[:3])
        assert_equal(coords, concat[3:6])
        assert_equal(coords, concat[6:])

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSelectConcat(self):

        command = self.command + ' -s ca --pdb {1:s} {0:s} {0:s}'.format(
                                self.dcdpath, self.pdbpath)

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        select = self.ag.ca
        assert_equal(select.numAtoms(), 10)

        coords = self.dcd[:]
        coords.setAtoms(select)
        coords = coords._getCoordsets()

        concat = parseDCD(self.output)
        assert_equal(concat.numAtoms(), select.numAtoms())
        concat = concat._getCoordsets()

        assert_equal(select.numAtoms(), coords.shape[1])
        assert_equal(select.numAtoms(), concat.shape[1])
        assert_equal(coords, concat[:3])
        assert_equal(coords, concat[3:])

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testAlignConcat(self):

        command = self.command + ' --align ca --pdb {1:s} {0:s} {0:s}'.format(
                                self.dcdpath, self.pdbpath)

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        select = self.ag.ca

        coords = self.dcd[:]
        concat = parseDCD(self.output)

        assert_equal(concat.numAtoms(), coords.numAtoms())

        coords.setCoords(self.ag.getCoords())
        coords.setAtoms(select)
        coords.superpose()
        coords.setAtoms(None)
        coords = coords._getCoordsets()

        concat = concat._getCoordsets()

        assert_equal(coords, concat[:3])
        assert_equal(coords, concat[3:])

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSelectException(self):

        command = self.command + ' -s ca {0:s} {0:s}'.format(
                                self.dcdpath)
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(ValueError, namespace.func, namespace)

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testAlignException(self):

        command = self.command + ' --align ca {0:s} {0:s}'.format(
                                self.dcdpath)
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(ValueError, namespace.func, namespace)

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testIOException(self):

        command = self.command + ' {0:s} {0:s}'.format('deneme.dcd')
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(IOError, namespace.func, namespace)

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command is not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testSelectException2(self):

        command = self.command + ' -s None {0:s} {0:s}'.format(self.dcdpath)
        namespace = prody_parser.parse_args(shlex.split(command))
        self.assertRaises(ValueError, namespace.func, namespace)


    def tearDown(self):

        if isfile(self.output): remove(self.output)
