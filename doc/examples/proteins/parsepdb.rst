.. currentmodule:: prody.proteins

.. _parsepdb:

*******************************************************************************
Parse PDB files
*******************************************************************************

Synopsis
===============================================================================

ProDy offers a fast and flexible PDB parser, :func:`parsePDB`. 
Parser can be used to read well defined subsets of atoms, specific chains or 
models (in NMR structures) to boost the performance. This example shows how to 
use the flexible parsing options. 

|more| For performance benchmarks of ProDy PDB parser see 
:ref:`pdbparser-performance`.

Input
-------------------------------------------------------------------------------

Three types of input are accepted from user:

  * PDB file path, e.g. ``"../1MKP.pdb"``
  * compressed (gzipped) PDB file path, e.g. ``"1p38.pdb.gz"`` 
  * PDB identifier, e.g. ``"2k39"``
 
Output
-------------------------------------------------------------------------------
 
Output is an :class:`~prody.atomic.AtomGroup` instance that stores atomic data
and can be used as input to functions and classes for dynamics analysis.  
 
ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Parse a file
-------------------------------------------------------------------------------

You can parse PDB files by passing a filename (gzipped files are handled).
We do so after downloading a PDB file (see :ref:`fetchpdb` for more 
information): 

>>> fetchPDB('1p38') # doctest: +SKIP
'1p38.pdb.gz'
>>> atoms = parsePDB('1p38.pdb.gz')
>>> atoms
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>

Parser returns an :class:`~prody.atomic.AtomGroup` instance.

Also note that the time it took to parse the file is printed on
the screen. This includes the time that it takes to evaluate 
coordinate lines and build an :class:`~prody.atomic.AtomGroup` instance and 
excludes the time spent on reading the file from disk.

Use an identifier
-------------------------------------------------------------------------------

PDB files can be parsed by passing simply an identifier. arser will look for a 
PDB file that matches the given identifier in the current working directory. 
If a matching file is not found, ProDy will downloaded it from PDB FTP server 
automatically and saved it in the current working directory.

>>> atoms = parsePDB('1mkp')
>>> atoms
<AtomGroup: 1mkp (1183 atoms; 1 coordinate sets, active set index: 0)>


Subsets of atoms
-------------------------------------------------------------------------------

Parser can be used to parse backbone or Cα atoms:

>>> backbone = parsePDB('1mkp', subset='bb')
>>> backbone
<AtomGroup: 1mkp (616 atoms; 1 coordinate sets, active set index: 0)>
>>> calpha = parsePDB('1mkp', subset='ca')
>>> calpha
<AtomGroup: 1mkp (144 atoms; 1 coordinate sets, active set index: 0)>


Specific chains
-------------------------------------------------------------------------------

Parser can be used to parse a specific chain from a PDB file:

>>> chA = parsePDB('3mkb', chain='A')
>>> chA
<AtomGroup: 3mkb (1198 atoms; 1 coordinate sets, active set index: 0)>
>>> chC = parsePDB('3mkb', chain='C')
>>> chC
<AtomGroup: 3mkb (1189 atoms; 1 coordinate sets, active set index: 0)>

Multiple chains can also be parsed in the same way:

>>> chAC = parsePDB('3mkb', chain='AC')
>>> chAC
<AtomGroup: 3mkb (2387 atoms; 1 coordinate sets, active set index: 0)>

Specific models
-------------------------------------------------------------------------------

Parser can be used to parse a specific model from a file:

>>> model1 = parsePDB('2k39', model=10)
>>> model1
<AtomGroup: 2k39 (1231 atoms; 1 coordinate sets, active set index: 0)>

Alternate locations
-------------------------------------------------------------------------------

When a PDB file contains alternate locations for some of the atoms, by default
alternate locations with indicator ``A`` are parsed. 

>>> altlocA = parsePDB('1ejg')
>>> altlocA
<AtomGroup: 1ejg (637 atoms; 1 coordinate sets, active set index: 0)>

Specific alternate locations can be parsed as follows:

>>> altlocB = parsePDB('1ejg', altloc='B')
>>> altlocB
<AtomGroup: 1ejg (634 atoms; 1 coordinate sets, active set index: 0)>

Note that in this case number of atoms are different between the two atom 
groups. This is because the residue types of atoms with alternate locations
are different.

Also, all alternate locations can be parsed as follows:

>>> all_altlocs = parsePDB('1ejg', altloc=True)
>>> all_altlocs
<AtomGroup: 1ejg (637 atoms; 3 coordinate sets, active set index: 0)>

Note that this time parser returned three coordinate sets. One for each 
alternate location indicator found in this file (A, B, C). When parsing
multiple alternate locations, parser will expect for the same residue type
for each atom with an alternate location. If residue names differ, a warning
message will be printed.

Composite arguments
-------------------------------------------------------------------------------

Parser can be used to parse coordinates from a specific model for a subset of 
atoms of a specific chain:

>>> composite = parsePDB('2k39', model=10, chain='A', subset='ca')
>>> composite
<AtomGroup: 2k39 (76 atoms; 1 coordinate sets, active set index: 0)>

Header data
-------------------------------------------------------------------------------

PDB parser can be used to extract header data from PDB files as follows:

>>> atoms, header = parsePDB('1mkp', header=True)

Header data is returned in a :class:`dict`. Printing its keys will show what
was parsed.

>>> header['experiment']
'X-RAY DIFFRACTION'
>>> header['resolution']
'2.35 ANGSTROMS'
>>> print header.keys() # doctest: +SKIP
['biomolecular_transformations', 'reference', 'classification', 'compounds', 
'resolution', 'title', 'source', 'experiment', 'helix', 'authors', 'sheet', 
'identifier', 'deposition_date']

It is also possible to parse only header data by passing `model=0` as an 
argument:

>>> header = parsePDB('1mkp', header=True, model=0)
>>> print header.keys() # doctest: +SKIP
['biomolecular_transformations', 'reference', 'classification', 'compounds', 
'resolution', 'title', 'source', 'experiment', 'helix', 'authors', 'sheet', 
'identifier', 'deposition_date']


|questions|

|suggestions|
