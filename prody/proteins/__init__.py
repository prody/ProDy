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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np

DTYPES = {'coordinates': np.float64, 

          'atomnames': '|S6', 
          'altlocs': '|S1',
          'anisou': np.float64, 
          'chainids': '|S1', 
          'elements': '|S2', 
          'hetero': np.bool, 
          'occupancies': np.float64, 
          'resnames': '|S6',
          'resnums': np.int64, 
          'secondary': '|S1',
          'segnames': '|S6', 
          'siguij': np.float64,
          'bfactors': np.float64,
          'icodes': '|S1',
          
          'atomtypes': '|S6',
          'charges': np.float64, 
          'masses': np.float64, 
          'radii': np.float64, 
}
 
__all__ = []

import atom
from atom import *
__all__ += atom.__all__

import prodb
from prodb import *  
__all__ += prodb.__all__

from . import keywords
from keywords import *
__all__ += keywords.__all__

from . import measure
from measure import *
__all__ += measure.__all__

import select

from . import subsets 
from subsets import *
__all__ += subsets.__all__

from . import hierview
from hierview import *
__all__ += hierview.__all__

from . import atomgroup 
from atomgroup import *
__all__ += atomgroup.__all__

from . import atommap 
from atommap import *
__all__ += atommap.__all__

from . import hierview 
from hierview import *
__all__ += hierview.__all__

from . import compare
from compare import *
__all__ += compare.__all__
