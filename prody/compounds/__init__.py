# -*- coding: utf-8 -*-
"""This module defines classes and functions to fetch, parse, and write
structural data files associated with compounds and to access and search 
structural databases, e.g. `ProteinDataBank <http://wwpdb.org>`_.

Ligand data
===========

The following function can be used to fetch metadata on PDB ligands:

  * :func:`.fetchPDBLigand` - retrieve ligand from Ligand-Expo


BIRD data
===========

The following functions can be used to fetch metadata from the 
Biologically Interesting Reference Dictionary (BIRD) for peptide-like 
ligands:

  * :func:`.fetchBIRDviaFTP` - retrieve BIRD data from PDB
  * :func:`.parseBIRD` - parse BIRD data for a particular compound or family


CCD data
===========

The following function can be used to fetch metadata from the 
Chemical Component Dictionary (CCD) for residues within ligands:

  * :func:`.parseCCD` - parse CCD data for a particular component


Functions
===========

The following functions can be used to calculate 2D similarity using Morgan Fingerprints:

  * :func:`.calc2DSimilarity` - calculate 2D similarity for a pair of SMILES
  * :func:`.calc2DSimilarityMatrix` - calculate 2D similarity matrix for a collection of SMILES

"""

__all__ = []

from . import pdbligands
from .pdbligands import *
__all__.extend(pdbligands.__all__)

from . import bird
from .bird import *
__all__.extend(bird.__all__)

from . import ccd
from .ccd import *
__all__.extend(ccd.__all__)

from . import functions
from .functions import *
__all__.extend(functions.__all__)


