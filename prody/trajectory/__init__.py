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

"""This module defines classes for handling trajectory files in DCD format.


Parse/write DCD files
===============================================================================

  * :class:`~.DCDFile`
  * :func:`~.parseDCD`
  * :func:`~.writeDCD`

Parse structure files
===============================================================================

  * :func:`~.parsePSF`  

Handle multiple files
===============================================================================
  
  * :class:`~.Trajectory`

Handle frame data
===============================================================================
  
  * :class:`~.Frame`

Examples
===============================================================================

Following examples show how to use trajectory classes and functions:
    
  * :ref:`trajectory`
  * :ref:`trajectory2`
  * :ref:`eda`"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = []

def openTrajFile(filename):
    
    ext = os.path.splitext(filename)[1]
    if ext.lower() == '.dcd':
        traj = DCDFile(filename)
    else: 
        raise ValueError('Trajectory file type {0:s} is not recognized.'
                         .format(repr(ext)))

    return traj 
    
    
import trajbase
from trajbase import *
__all__.extend(trajbase.__all__)

import trajectory
from trajectory import *
__all__.extend(trajectory.__all__)

import dcdfile
from dcdfile import *
__all__.extend(dcdfile.__all__)

import frame
from frame import *
__all__.extend(frame.__all__)
   
import psffile
from psffile import *
__all__.extend(psffile.__all__)
   
