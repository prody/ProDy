# -*- coding: utf-8 -*-
"""This module defines classes and functions to fetch, parse, and write
structural data files, execute structural analysis programs, and to access
and search structural databases, e.g. `ProteinDataBank <http://wwpdb.org>`_.

PDB resources
=============

  * :func:`.fetchPDB` - retrieve PDB files
  * :func:`.fetchPDBviaFTP` - download PDB/PDBML/mmCIF files
  * :func:`.fetchPDBviaHTTP` - download PDB files

You can use following functions to manage PDB file resources:

  * :func:`.pathPDBFolder` - local folder for storing PDB files
  * :func:`.pathPDBMirror` - local PDB mirror path
  * :func:`.wwPDBServer` - set wwPDB FTP/HTTP server for downloads

The following functions can be used to handle local PDB files:

  * :func:`.findPDBFiles` - return a dictionary containing files in a path
  * :func:`.iterPDBFilenames` - yield file names in a path or local PDB mirror


Blast search PDB
================

The following are for blast searching PDB content.

  * :func:`.blastPDB` - blast search NCBI PDB database
  * :class:`.PDBBlastRecord` - store/evaluate NCBI PDB blast search results

PDB clusters biopolymer chains using blast weekly.  These clusters can be
retrieved using the following functions.  Using cluster data is as good
as blast searching PDB most of the time and incredibly faster always.

  * :func:`.listPDBCluster` - get list of identifiers in a PDB sequence cluster
  * :func:`.loadPDBClusters` - load PDB clusters into memory
  * :func:`.fetchPDBClusters` - retrieve PDB sequence cluster data from wwPDB

Parse/write PDB files
=====================

Following ProDy functions are for parsing and writing :file:`.pdb` files:

  * :func:`.parsePDB` - parse :file:`.pdb` formated file
  * :func:`.parsePDBStream` - parse :file:`.pdb` formated stream
  * :func:`.writePDB` - write :file:`.pdb` formatted file
  * :func:`.writePDBStream`  write :file:`.pdb` formated stream

Since :file:`.pqr` format is similar to :file:`.pdb` format, following
functions come as bonus features:

  * :func:`.writePQR` - write atomic data to a file in :file:`.pqr` format
  * :func:`.parsePQR` - parse atomic data from files in :file:`.pqr` format


.. seealso::

   Atom data (coordinates, atom names, residue names, etc.) parsed from
   PDB/PSF/PQR/mmCIF files are stored in :class:`~.AtomGroup` instances.
   See :mod:`~prody.atomic` module documentation for more details.

Parse mmCIF files
=====================

Following ProDy functions are for parsing :file:`.cif` files:

  * :func:`.parseMMCIF` - parse :file:`.cif` formated file
  * :func:`.parseMMCIFStream` - parse :file:`.cif` formated stream

.. seealso::

   Atom data (coordinates, atom names, residue names, etc.) parsed from
   PDB/PSF/PQR/mmCIF files are stored in :class:`~.AtomGroup` instances.
   See :mod:`~prody.atomic` module documentation for more details.

Quick visualization
===================

:func:`.showProtein` function can be used to take a quick look at protein
structures.

Edit structures
===============

The following functions allow editing structures using structural data from PDB
header records:

  * :func:`.assignSecstr` - add secondary structure data from header to atoms
  * :func:`.buildBiomolecules` - build biomolecule from header records


PDB header data
===============

Use the following to parse and access header data in PDB files:

  * :func:`.parsePDBHeader` - parse header data from :file:`.pdb` files
  * :class:`.Chemical` - store PDB chemical (heterogen) component data
  * :class:`.Polymer` - store PDB polymer (macromolecule) component data
  * :class:`.DBRef` - store polymer sequence database reference records

Analyze interactions and stability with InSty and find water bridges with WatFinder
====================

Use the following to analyze interactions within protein structure or
between protein and ligand structure in single PDB file or in trajectory:

  * :func:`.addHydrogens` - add missing hydrogens to :file:`.pdb` files
  * :func:`.calcHydrogenBonds` - compute hydrogen bonds in proteins
  * :func:`.calcSaltBridges` - compute salt bridges in proteins
  * :func:`.calcRepulsiveIonicBonding` - compute repulsive ionic bonding in proteins
  * :func:`.calcPiStacking` - compute Pi-stacking interactions in proteins
  * :func:`.calcPiCation` - compute Pi-cation interactions in proteins
  * :func:`.calcHydrophobic` - compute hydrophobic interactions in proteins
  * :func:`.calcProteinInteractions` - compute all above interaction types at once
  * :func:`.showProteinInteractions_VMD` - return TCL file for visualization in VMD 

  * :func:`.calcHydrogenBondsDCD` - compute hydrogen bonds in a trajectory for proteins
  * :func:`.calcSaltBridgesDCD` - ompute salt bridges in a trajectory for proteins
  * :func:`.calcRepulsiveIonicBondingDCD` - compute repulsive ionic bonding in a trajectory for proteins    
  * :func:`.calcPiStackingDCD` - compute Pi-stacking interactions in a trajectory for proteins
  * :func:`.calcPiCationDCD` - compute Pi-cation interactions in a trajectory for proteins
  * :func:`.calcHydrophobicDCD` - compute hydrophobic interactions in a trajectory for proteins
  * :func:`.calcStatisticsInteractions` - return statistical information for each interaction type 

  * :func:`.calcLigandInteractions` - compute all type of interactions between protein and ligand
  * :func:`.listLigandInteractions` - return list of interactions between protein and ligand
  * :func:`.showLigandInteraction_VMD` - return TCL file for visualization of interactions for VMD

  * :class:`.Interactions` - store inteactions for a single PDB structure
  * :class:`.InteractionsDCD` - store interactions for a trajectory

Compare/align chains
====================

The following functions can be used to match, align, and map polypeptide chains:

  * :func:`.matchChains` - finds matching chains in two protein structures
  * :func:`.matchAlign` - finds best matching chains and aligns structures
  * :func:`.mapOntoChain` - maps chains in a structure onto a reference chain

The following functions can be used to adjust alignment parameters:

  * :func:`.getAlignmentMethod`, :func:`.setAlignmentMethod`
  * :func:`.getMatchScore`, :func:`.setMatchScore`
  * :func:`.getMismatchScore`, :func:`.setMismatchScore`
  * :func:`.getGapPenalty`, :func:`.setGapPenalty`
  * :func:`.getGapExtPenalty`, :func:`.setGapExtPenalty`

Execute DSSP
============

The following functions can be used to execute DSSP structural analysis program
and/or parse results:

  * :func:`.execDSSP` - execute :program:`dssp`
  * :func:`.performDSSP` - execute :program:`dssp` and parse results
  * :func:`.parseDSSP` - parse structural data from :program:`dssp` output

Execute STRIDE
==============

The following functions can be used to execute STRIDE structural analysis program
and/or parse results:

  * :func:`.execSTRIDE` - execute :program:`stride`
  * :func:`.parseSTRIDE` - parse structural data from :program:`stride` output
  * :func:`.performSTRIDE` - execute :program:`stride` and parse results

Handle EMD Map Files and Build Pseudoatoms into them
===========

Use the following to parse and access header data in EMD files:

  * :func:`.parseEMD` - parse structural data from :file:`.emd` files
  * :class:`.EMDMAP` - access structural data from :file:`.emd` files
  * :class:`.TRNET` - fit pseudoatoms to EM density maps using the TRN algorithm

Add missing atoms including hydrogens
===========

Use the following to add missing atoms

  * :func:`.addMissingAtoms` - add missing atoms with separately installed OpenBabel or PDBFixer

Missing residues can also be added if a PDB or mmCIF file with SEQRES entries is provided.
"""

__all__ = []

from . import compare
from .compare import *
__all__.extend(compare.__all__)

from . import localpdb
from .localpdb import *
__all__.extend(localpdb.__all__)

from . import wwpdb
from .wwpdb import *
__all__.extend(wwpdb.__all__)

from . import pdbclusters
from .pdbclusters import *
__all__.extend(pdbclusters.__all__)

from . import blastpdb
from .blastpdb import *
__all__.extend(blastpdb.__all__)

from . import functions
from .functions import *
__all__.extend(functions.__all__)

from . import header
from .header import *
__all__.extend(header.__all__)

from . import cifheader
from .cifheader import *
__all__.extend(cifheader.__all__)

from . import dssp
from .dssp import *
__all__.extend(dssp.__all__)

from . import stride
from .stride import *
__all__.extend(stride.__all__)

from . import pdbfile
from .pdbfile import *
__all__.extend(pdbfile.__all__)

from . import mmtffile
from .mmtffile import *
__all__.extend(mmtffile.__all__)

from . import emdfile
from .emdfile import *
__all__.extend(emdfile.__all__)

from . import ciffile
from .ciffile import *
__all__.extend(ciffile.__all__)

from . import starfile
from .starfile import *
__all__.extend(starfile.__all__)

from . import interactions
from .interactions import *
__all__.extend(interactions.__all__)

try:
    from . import waterbridges
    from .waterbridges import *
except SyntaxError:
    import logging
    logger = logging.getLogger()
    logger.warn("Cannot import waterbridges")
else:
    __all__.extend(waterbridges.__all__)

from . import fixer
from .fixer import *
__all__.extend(fixer.__all__)

from .pdbfile import PDBParseError
from .ciffile import MMCIFParseError

from . import opm
from .opm import *
__all__.extend(opm.__all__)
