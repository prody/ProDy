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

"""This module contains features for accessing databases containing protein 
related data. 

Pfam_
=====

Following functions can be used to search and retrieve Pfam data:

  * :func:`.fetchPfamMSA` - download MSA files
  * :func:`.searchPfam` - search families of a protein
  
    
.. _Pfam: http://pfam.sanger.ac.uk/"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

__all__ = []

from . import pfam
from .pfam import *
__all__.extend(pfam.__all__)
