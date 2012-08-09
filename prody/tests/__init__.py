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
  
or::

  import prody.tests
  prody.tests.test()


Testing will use :mod:`nose` if it is available, otherwise it will use 
:mod:`unittest`.

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from sys import stderr
from glob import glob
from os.path import abspath, split, join, relpath, splitext
from os.path import sep as dirsep
import inspect
import tempfile

from prody import LOGGER

TESTDIR = abspath(split(inspect.getfile(inspect.currentframe()))[0])
TEMPDIR = tempfile.gettempdir()

try:
    import nose
    
except ImportError:

    LOGGER.warning('Failed to import nose, using unittest for testing.')
    LOGGER.info('nose is available at http://readthedocs.org/docs/nose/')
    import unittest
    def test(verbosity=2, descriptions=True, stream=stderr):
        testrunner = unittest.TextTestRunner(stream, descriptions, 
                                             verbosity)
        for pyfile in glob(join(TESTDIR, '*', '*.py')):
            pyfile = splitext(relpath(pyfile, TESTDIR))[0]
            module = '.'.join(pyfile.split(dirsep))
            print module
            testrunner.run(unittest.defaultTestLoader.
                           loadTestsFromName('prody.tests.' + module))
else:
    from numpy.testing import Tester
    test = Tester().test

if __name__ == '__main__':
    test()
