# -*- coding: utf-8 -*-
"""This module provides utility functions and classes for handling files,
logging, type checking, etc.  Contents of this module are not included in
ProDy namespace, as it is not safe to import them all due to name conflicts.
Required or classes should be imported explicitly, e.g.
``from prody.utilities import PackageLogger, openFile``.


Package utilities
===============================================================================

  * :class:`.PackageLogger`
  * :class:`.PackageSettings`
  * :func:`.getPackagePath`
  * :func:`.setPackagePath`

Type/Value checkers
===============================================================================

  * :func:`.checkCoords`
  * :func:`.checkWeights`
  * :func:`.checkTypes`

Path/file handling
===============================================================================

  * :func:`.gunzip`
  * :func:`.openFile`
  * :func:`.openDB`
  * :func:`.openSQLite`
  * :func:`.openURL`
  * :func:`.copyFile`
  * :func:`.isExecutable`
  * :func:`.isReadable`
  * :func:`.isWritable`
  * :func:`.makePath`
  * :func:`.relpath`
  * :func:`.which`
  * :func:`.pickle`
  * :func:`.unpickle`
  * :func:`.glob`


Documentation tools
===============================================================================

  * :func:`.joinRepr`
  * :func:`.joinRepr`
  * :func:`.joinTerms`
  * :func:`.tabulate`
  * :func:`.wrapText`


Miscellaneous tools
===============================================================================

  * :func:`.rangeString`
  * :func:`.alnum`
  * :func:`.importLA`
  * :func:`.dictElement`

"""

__all__ = []

from .checkers import *
from .logger import *
from .settings import *
from .pathtools import *
from .misctools import *
from .doctools import *
