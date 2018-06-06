# -*- coding: utf-8 -*-
"""This module defines classes and functions to parse and write Hi-C data files, 
visualize, and perform elastic network analysis on Hi-C data.

Parse/write Hi-C files
=====================

Following ProDy functions are for parsing and writing Hi-C files:

  * :func:`.parseHiC` - parse Hi-C data file
  * :func:`.parseHiCStream` - parse Hi-C data stream
  * :func:`.writeMap` - write Hi-C data to text file

Visualize Hi-C data
=====================

Following ProDy functions are for visualizing Hi-C data:

  * :func:`.showMap` - show Hi-C contact map
  * :func:`.showDomains` - show Hi-C structural domains

Save/load HiC class
-------------------

    * :func:`.saveHiC`
    * :func:`.loadHiC`

"""

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = []

from . import hic
from .hic import *
__all__.extend(hic.__all__)

from . import norm
from .norm import *
__all__.extend(norm.__all__)

from . import cluster
from .cluster import *
__all__.extend(cluster.__all__)

from . import functions
from .functions import *
__all__.extend(functions.__all__)
