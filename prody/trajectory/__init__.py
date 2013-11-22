# -*- coding: utf-8 -*-
"""This module defines classes for handling trajectory files in DCD format.


Parse/write DCD files
===============================================================================

  * :class:`.DCDFile`
  * :func:`.parseDCD`
  * :func:`.writeDCD`

Parse structure files
===============================================================================

  * :func:`.parsePSF`

Handle multiple files
===============================================================================

  * :class:`.Trajectory`

Handle frame data
===============================================================================

  * :class:`.Frame`

Examples
===============================================================================

Following examples show how to use trajectory classes and functions:

  * :ref:`trajectory`
  * :ref:`trajectory2`
  * :ref:`eda`"""

from os.path import splitext

__all__ = []

def openTrajFile(filename, *args, **kwargs):

    ext = splitext(filename)[1][1:].lower()
    try:
        return TRAJFILE[ext](filename, *args, **kwargs)
    except KeyError:
        raise ValueError('trajectory file type {0} is not recognized'
                         .format(repr(ext)))

from . import trajbase
from .trajbase import *
__all__.extend(trajbase.__all__)

from . import trajectory
from .trajectory import *
__all__.extend(trajectory.__all__)

from . import dcdfile
from .dcdfile import *
__all__.extend(dcdfile.__all__)

from . import frame
from .frame import *
__all__.extend(frame.__all__)

from . import psffile
from .psffile import *
__all__.extend(psffile.__all__)

TRAJFILE = {'dcd': DCDFile}

