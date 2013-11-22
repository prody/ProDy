# -*- coding: utf-8 -*-
"""ProDy test suite.  Usage::

  from prody import *
  prody.test()

Testing will use :mod:`nose` if it is available, otherwise it will use
:mod:`unittest`."""

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