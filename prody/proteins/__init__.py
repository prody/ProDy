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

"""This module defines classes and functions to fetch, parse, and write 
structural data files, execute structural analysis programs, and to access 
and search structural databases, e.g. `ProteinDataBank <http://wwpdb.org>`_.

Retrieve PDB files
==================

Following functions provide access to wwPDB FTP servers:

  * :func:`~.fetchPDB` - retrieve PDB/PDBML/mmCIF files from wwPDB
  * :func:`~.setWWPDBFTPServer` - set a wwPDB FTP server for downloads 
  * :func:`~.getWWPDBFTPServer` - get preset wwPDB FTP server


Local PDB folder/mirror 
=======================

You can use following functions to manage local PDB file resources:

  * :func:`~.setPDBLocalFolder` - set a local folder for storing PDB files
  * :func:`~.setPDBMirrorPath` - set a local PDB mirror path
  * :func:`~.getPDBLocalFolder` - get preset local PDB folder
  * :func:`~.getPDBMirrorPath` - get preset local PDB mirror path

Blast search PDB 
================

The following are for blast searching PDB content.

  * :func:`~.blastPDB` - blast search NCBI PDB database
  * :class:`~.PDBBlastRecord` - store/evaluate NCBI PDB blast search results 

PDB clusters biopolymer chains using blast weekly.  These clusters can be
retrieved using the following functions.  Using cluster data is as good
as blast searching PDB most of the time and incredibly faster always. 

  * :func:`~.getPDBCluster` - access PDB sequence cluster
  * :func:`~.loadPDBClusters` - load PDB clusters into memory
  * :func:`~.fetchPDBClusters` - retrieve PDB sequence cluster data from wwPDB

Parse/write PDB files
=====================

Following ProDy functions are for parsing and writing :file:`.pdb` files:

  * :func:`~.parsePDB` - parse :file:`.pdb` formated file
  * :func:`~.parsePDBStream` - parse :file:`.pdb` formated stream
  * :func:`~.writePDB` - write :file:`.pdb` formatted file
  * :func:`~.writePDBStream`  write :file:`.pdb` formated stream

Since :file:`.pqr` format is similar to :file:`.pdb` format, following 
functions come as bonus features:
 
  * :func:`~.writePQR`        write atomic data to a file in :file:`.pqr` format
  * :func:`~.parsePQR`        parse atomic data from files in :file:`.pqr` format
  

.. seealso::
   
   Atom data (coordinates, atom names, residue names, etc.) parsed from 
   PDB/PSF/PQR files are stored in :class:`~.AtomGroup` instances.  
   See :mod:`~prody.atomic` module documentation for more details. 

Quick visualization
===================

:func:`~.showProtein` function can be used to take a quick look at protein 
structures. 

Edit structures
===============

Following functions allow editing structures using structural data from PDB 
header records:

:func:`~.assignSecstr`       add secondary structure data from header to atoms
:func:`~.buildBiomolecules`  build biomolecule data based on header records


PDB header data
===============

Use the following to parse and access header data in PDB files:
  
  * :func:`~.parsePDBHeader` - parse header data from :file:`.pdb` files 
  * :class:`~.Chemical` - store PDB chemical (heterogen) component data
  * :class:`~.Polymer` - store PDB polymer (macromolecule) component data
  * :class:`~.DBRef` - store polymer sequence database reference records

  * :func:`~.parsePSF`        parse atomic data from files in :file:`.psf` format

Ligand data
===========

Following function can be used to fetch meta data on PDB ligands:

  * :func:`~.fetchLigandData` - retrieve ligand from Ligand-Expo 

Compare/align chains
====================

Following functions can be used to match, align, and map polypeptide chains:

  * :func:`~.matchChains` - finds matching chains in two protein structures
  * :func:`~.matchAlign` - finds best matching chains and aligns structures
  * :func:`~.mapOntoChain` - maps chains in a structure onto a reference chain
        

The following functions can be used to adjust alignment parameters: 
        
  * :func:`~.getAlignmentMethod`, :func:`~.setAlignmentMethod`
  * :func:`~.getMatchScore`, :func:`~.setMatchScore`
  * :func:`~.getMismatchScore`, :func:`~.setMismatchScore`
  * :func:`~.getGapPenalty`, :func:`~.setGapPenalty`
  * :func:`~.getGapExtPenalty`, :func:`~.setGapExtPenalty`
    


Execute DSSP
============

Following functions can be used to execute DSSP structural analysis program
and/or parse results:

  * :func:`~.execDSSP` - execute :program:`dssp`
  * :func:`~.performDSSP` - execute :program:`dssp` and parse results
  * :func:`~.parseDSSP` - parse structural data from :program:`dssp` output

Execute STRIDE
==============

Following functions can be used to execute STRIDE structural analysis program
and/or parse results:

  * :func:`~.execSTRIDE` - execute :program:`stride`
  * :func:`~.performSTRIDE` - execute :program:`stride` and parse results
  * :func:`~.parseSTRIDE` - parse structural data from :program:`stride` output

.. doctest::
    :hide:
        
    >>> from prody import *
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS


__all__ = []

import compare
from compare import *
__all__.extend(compare.__all__)

import databank
from databank import *
__all__.extend(databank.__all__)

import functions
from functions import *
__all__.extend(functions.__all__)

import header
from header import *
__all__.extend(header.__all__)

import programs
from programs import *
__all__.extend(programs.__all__)

import pdbfile
from pdbfile import *
__all__.extend(pdbfile.__all__)

from pdbfile import PDBParseError
