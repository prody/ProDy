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

Following functions can be used to handle local PDB files:

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
   PDB/PSF/PQR files are stored in :class:`~.AtomGroup` instances.
   See :mod:`~prody.atomic` module documentation for more details.

Quick visualization
===================

:func:`.showProtein` function can be used to take a quick look at protein
structures.

Edit structures
===============

Following functions allow editing structures using structural data from PDB
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

Ligand data
===========

Following function can be used to fetch meta data on PDB ligands:

  * :func:`.fetchPDBLigand` - retrieve ligand from Ligand-Expo

Compare/align chains
====================

Following functions can be used to match, align, and map polypeptide chains:

  * :func:`.matchChains` - finds matching chains in two protein structures
  * :func:`.matchAlign` - finds best matching chains and aligns structures
  * :func:`.mapOntoChain` - maps chains in a structure onto a reference chain

Following functions can be used to adjust alignment parameters:

  * :func:`.getAlignmentMethod`, :func:`.setAlignmentMethod`
  * :func:`.getMatchScore`, :func:`.setMatchScore`
  * :func:`.getMismatchScore`, :func:`.setMismatchScore`
  * :func:`.getGapPenalty`, :func:`.setGapPenalty`
  * :func:`.getGapExtPenalty`, :func:`.setGapExtPenalty`



Execute DSSP
============

Following functions can be used to execute DSSP structural analysis program
and/or parse results:

  * :func:`.execDSSP` - execute :program:`dssp`
  * :func:`.performDSSP` - execute :program:`dssp` and parse results
  * :func:`.parseDSSP` - parse structural data from :program:`dssp` output

Execute STRIDE
==============

Following functions can be used to execute STRIDE structural analysis program
and/or parse results:

  * :func:`.execSTRIDE` - execute :program:`stride`
  * :func:`.performSTRIDE` - execute :program:`stride` and parse results
  * :func:`.parseSTRIDE` - parse structural data from :program:`stride` output
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

from . import pdbligands
from .pdbligands import *
__all__.extend(pdbligands.__all__)

from . import functions
from .functions import *
__all__.extend(functions.__all__)

from . import header
from .header import *
__all__.extend(header.__all__)

from . import dssp
from .dssp import *
__all__.extend(dssp.__all__)

from . import stride
from .stride import *
__all__.extend(stride.__all__)

from . import pdbfile
from .pdbfile import *
__all__.extend(pdbfile.__all__)

from .pdbfile import PDBParseError
