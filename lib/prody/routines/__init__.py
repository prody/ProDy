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

"""This module contains routines that are used as command line programs.


Dynamics analysis routines
===============================================================================

  * :func:`.prody_anm`
  * :func:`.prody_gnm`
  * :func:`.prody_pca`

Structure analysis routines
===============================================================================

  * :func:`.prody_align`
  * :func:`.prody_biomol`
  * :func:`.prody_blast`
  * :func:`.prody_catdcd`
  * :func:`.prody_contacts`
  * :func:`.prody_fetch`
  * :func:`.prody_select`
  
Sequence analysis routines
===============================================================================

  * :func:`.evol_fetch`
  * :func:`.evol_search`
  * :func:`.evol_refine`
  * :func:`.evol_entropy`
  * :func:`.evol_mutinfo`
  * :func:`.evol_dynamics`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from prody_routines import prody_main, prody_parser, PRODY_COMMANDS
from evol_routines import evol_main, evol_parser, EVOL_COMMANDS
