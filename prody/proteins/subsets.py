# ProDy: A Python Package for Protein Structural Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
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
""":mod:`subsets` module defines classes for manipulating subsets of atomic 
data.

Classes:
    
    * :class:`AtomSubset`
    * :class:`Chain`
    * :class:`Residue`
    * :class:`Selection`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np

import prody
from prody import ProDyLogger as LOGGER
from . import Atom
from .select import ProDyAtomSelect as SELECT

__all__ = ['AtomSubset', 'Selection', 'Chain', 'Residue']

