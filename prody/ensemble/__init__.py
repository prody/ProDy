# -*- coding: utf-8 -*-
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

Following functions are for creating or editing PDB ensembles, e.g. finding and 
removing residues that are missing in too many structures:
  
  * :func:`.buildPDBEnsemble`
  * :func:`.alignPDBEnsemble`
  * :func:`.calcOccupancies`
  * :func:`.showOccupancies`
  * :func:`.trimPDBEnsemble`

See usage examples in :ref:`pca-xray`, :ref:`pca-dimer`, :ref:`pca-blast`.

Save/load ensembles
-------------------

    * :func:`.saveEnsemble`
    * :func:`.loadEnsemble`"""

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

from . import dali
from .dali import *
__all__.extend(dali.__all__)
