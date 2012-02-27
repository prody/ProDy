.. _writepdb:

*******************************************************************************
Write PDB file
*******************************************************************************

Synopsis
===============================================================================

PDB files can be written using :func:`~.writePDB` function. This
example shows how to write PDB files for :class:`~.AtomGroup` 
instances and subsets of atoms. 

Parse a PDB file
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *
 
You can parse PDB files by passing a PDB identifier:

>>> atoms = parsePDB('1p38')
>>> atoms
<AtomGroup: 1p38 (2962 atoms)>

:func:`~.parsePDB` function returns atomic data in an 
:class:`~.AtomGroup`.


Write all atoms
===============================================================================

All atoms in an :class:`~.AtomGroup` can be written in PDB format
as follows:

>>> writePDB('1p38.pdb', atoms)
'1p38.pdb'

Upon successful writing of PDB file, filename is returned.

Write a subset
===============================================================================

It is also possible to write subsets of atoms in PDB format:

>>> alpha_carbons = atoms.select('calpha')
>>> writePDB('1p38_ca.pdb', alpha_carbons)
'1p38_ca.pdb'

>>> backbone = atoms.select('backbone')
>>> writePDB('1p38_bb.pdb', backbone)
'1p38_bb.pdb'

|questions|

|suggestions|
