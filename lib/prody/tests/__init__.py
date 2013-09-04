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
except ImportError:
    MATPLOTLIB = False
else:
    MATPLOTLIB = True
    matplotlib.use('Agg')

TESTDIR = abspath(split(inspect.getfile(inspect.currentframe()))[0])
TEMPDIR = tempfile.gettempdir()

MODULES = dict()
PREFIX = 'prody.tests.'

for pyfile in glob(join(TESTDIR, '*', '*.py')):
    pyfile = splitext(relpath(pyfile, TESTDIR))[0]

    items = pyfile.split(dirsep)
    if items[-1] == '__init__':
        items = items[:-1]

    if items[-1].startswith('test_'):
        MODULES['.'.join([i[5:] for i in items])] = PREFIX + '.'.join(items)


def runTests(*mods, **kwargs):

    if mods:
        modules = []
        for mod in mods:
            try:
                modules.append(MODULES[mod])
            except KeyError:
                raise ValueError(mod + ' is not a valid test module name')
    else:
        modules = MODULES.values() # PY3K: OK

    try:
        import nose

    except ImportError:

        LOGGER.warning('Failed to import nose, using unittest for testing.')
        LOGGER.info('nose is available at http://readthedocs.org/docs/nose/')
        from sys import stderr

        verbosity = kwargs.get('verbose', 2)
        descriptions = kwargs.get('descriptions', True)
        stream = kwargs.get('stream', stderr)

        testrunner = unittest.TextTestRunner(stream, descriptions, verbosity)

        for module in modules:
            testrunner.run(unittest.defaultTestLoader.
                           loadTestsFromName(module))
    else:
        from numpy.testing import Tester
        verbose = kwargs.get('verbose', 1)
        label = kwargs.get('label', 'fast')

        if mods:
            for module in modules:
                Tester(module).test(label=label, verbose=verbose)
        else:
            Tester('prody.tests').test(label=label, verbose=verbose)

if __name__ == '__main__':
    runTests()
