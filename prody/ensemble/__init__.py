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

"""This module defines classes for handling conformational ensembles.

Conformational ensembles
------------------------

The following two classes are implemented for handling arbitrary but uniform
conformational ensembles, e.g. NMR models, MD snapshots:

  * :class:`.Ensemble`
  * :class:`.Conformation`

See usage examples in :ref:`pca-nmr` and :ref:`eda`.

PDB ensembles
-------------

PDB ensembles, such as multiple structures of the same protein, are in general
heterogeneous.  This just means that different residues in different structures
are missing.  The following classes extend above to support this heterogeneity:

  * :class:`.PDBEnsemble`
  * :class:`.PDBConformation`

Following functions are for editing PDB ensembles, e.g. finding and removing
residues that are missing in too many structures:

  * :func:`.alignPDBEnsemble`
  * :func:`.calcOccupancies`
  * :func:`.showOccupancies`
  * :func:`.trimPDBEnsemble`

See usage examples in :ref:`pca-xray`, :ref:`pca-dimer`, :ref:`pca-blast`.

Save/load ensembles
-------------------

    * :func:`.saveEnsemble`
    * :func:`.loadEnsemble`"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

__all__ = []

from . import functions
from .functions import *
__all__.extend(functions.__all__)

from . import ensemble
from .ensemble import *
__all__.extend(ensemble.__all__)

from . import pdbensemble
from .pdbensemble import *
__all__.extend(pdbensemble.__all__)

from . import conformation
from .conformation import *
__all__.extend(conformation.__all__)

from .functions import checkWeights
ensemble.checkWeights = checkWeights
pdbensemble.checkWeights = checkWeights
