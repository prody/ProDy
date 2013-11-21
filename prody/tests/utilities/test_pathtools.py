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

from prody.tests import TestCase

from glob import glob
from os import remove
from os.path import join

from prody.tests.datafiles import TEMPDIR
from prody.utilities import gunzip, openFile


class TestGunzip(TestCase):

    def setUp(self):

        self.pref = join(TEMPDIR, 'compressed.txt')
        self.gzfn = self.pref + '.gz'
        self.text = ''.join(['some random text '] * 100)
        try:
            self.bytes = bytes(self.text, encoding='utf-8')
        except TypeError:
            self.bytes = self.text
        out = openFile(self.gzfn, 'wt')
        out.write(self.text)
        out.close()

    def testFile(self):

        fn = gunzip(self.gzfn)
        text = open(fn).read()
        self.assertEqual(text, self.text)

    def testBuffer(self):

        buff = open(self.gzfn, 'rb').read()
        text = gunzip(buff)
        self.assertEqual(text, self.bytes)

    def tearDown(self):

        for fn in glob(self.pref + '*'):
            remove(fn)
