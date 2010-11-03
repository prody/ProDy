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

from prody import ProDyLogger as LOGGER

__all__ = []

from . import nma
from .nma import *
__all__ += nma.__all__

from . import modes
from .modes import *
__all__ += modes.__all__

from . import ensemble
from .ensemble import *
__all__ += ensemble.__all__

from . import functions
from .functions import *
__all__ += functions.__all__

try:
    from . import plotting
    from .plotting import *
    __all__ += plotting.__all__
except ImportError:
    LOGGER.warning('Matplotlib package could not be imported. Plotting functions will be unavailable.')    
