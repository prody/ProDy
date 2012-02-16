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

Protein structural data
===============================================================================

Access databases 
-------------------------------------------------------------------------------

Following functions are provided for access to protein structural data:

===========================  ==================================================
Function                     Description
===========================  ==================================================
:func:`~.blastPDB`           blast search NCBI PDB database
:func:`~.fetchPDB`           retrieve PDB/PDBML/mmCIF files from wwPDB
:func:`~.fetchLigandData`    retrieve ligand from Ligand-Expo 
:func:`~.fetchPDBClusters`   retrieve PDB sequence cluster data from wwPDB
:func:`~.getPDBCluster`      access PDB sequence clusters
:func:`~.setPDBLocalFolder`  set a local folder for storing PDB files
:func:`~.setPDBMirrorPath`   set a local PDB mirror path
:func:`~.setWWPDBFTPServer`  set a wwPDB FTP server for downloads 
:func:`~.getPDBLocalFolder`  get preset local PDB folder
:func:`~.getPDBMirrorPath`   get preset local PDB mirror path
:func:`~.getWWPDBFTPServer`  get preset wwPDB FTP server
===========================  ==================================================


Parse/write files
-------------------------------------------------------------------------------

Following ProDy functions are for parsing and writing atomic data:

========================  =====================================================
Function                  Description
========================  =====================================================
:func:`~.parsePDB`        parse atomic data from files in :file:`.pdb` format
:func:`~.parsePSF`        parse atomic data from files in :file:`.psf` format
:func:`~.parsePQR`        parse atomic data from files in :file:`.pqr` format
:func:`~.parsePDBHeader`  parse header data from :file:`.pdb` files 
:func:`~.parsePDBStream`  parse atomic data from a :file:`.pdb` formated stream
:func:`~.parseDSSP`       parse structural data from :program:`dssp` output
:func:`~.parseSTRIDE`     parse structural data from :program:`stride` output
:func:`~.writePDB`        write atomic data to a file in :file:`.pdb` format
:func:`~.writePQR`        write atomic data to a file in :file:`.pqr` format
:func:`~.writePDBStream`  write atomic data in :file:`.pdb` format to a stream
========================  =====================================================

.. seealso::
   
   Atom data (coordinates, atom names, residue names, etc.) parsed from 
   PDB/PSF/PQR files are stored in :class:`~prody.atomic.AtomGroup` instances.  
   See :mod:`~prody.atomic` module documentation for more details. 


Store data
-------------------------------------------------------------------------------

Following classes are for storing meta, structural, and/or search data: 

=========================  ====================================================
Function                   Description
=========================  ====================================================
:class:`~.Chemical`        store PDB chemical (heterogen) component data
:class:`~.Polymer`         store PDB polymer (macromolecule) component data
:class:`~.DBRef`           store polymer sequence database reference records
:class:`~.PDBBlastRecord`  store and evaluate NCBI PDB blast search results 
=========================  ====================================================

Execute programs
-------------------------------------------------------------------------------

Following functions can be used to execute structural analysis programs from 
within Python:

==========================  ===================================================
Function                    Description
==========================  ===================================================
:func:`~.execDSSP`          execute :program:`dssp`
:func:`~.execSTRIDE`        execute :program:`stride`
:func:`~.performDSSP`       execute :program:`dssp` and parse results
:func:`~.performSTRIDE`     execute :program:`stride` and parse results
==========================  ===================================================


Edit structures
-------------------------------------------------------------------------------

Following functions allow editing structures using structural data from PDB 
header records:

=========================  ====================================================
Function                   Description
=========================  ====================================================
:func:`assignSecstr`       add secondary structure data from header to atoms
:func:`buildBiomolecules`  build biomolecule data based on header records
=========================  ====================================================


Visualization
-------------------------------------------------------------------------------

:func:`showProtein` function can be used to take a quick look at protein 
structures. 

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
