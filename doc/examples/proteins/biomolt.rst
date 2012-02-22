.. _biomolt:

*******************************************************************************
Building biomolecules
*******************************************************************************

Synopsis
===============================================================================

Some PDB files contain coordinates for a monomer from a functional/biological 
multimer (biomolecule).  ProDy offers functions to build structures of 
biomolecules using the header data from the PDB file.

Input
-------------------------------------------------------------------------------

A PDB file that contains the coordinates for a monomer of a biological 
multimeric protein and the transformations in the header section to
generate the multimer coordinates.

Output
-------------------------------------------------------------------------------

An :class:`~.AtomGroup` instance that contains the multimer 
coordinates.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Parse a PDB file
-------------------------------------------------------------------------------

We parse a PDB file that contains the coordinates for a monomer of a dimeric
protein:

>>> monomer, header = parsePDB('3enl', header=True)
>>> monomer
<AtomGroup: 3enl (3647 atoms)>

Note that we passed ``header=True`` argument to parse header data in addition
to coordinates.

Build multimer
-------------------------------------------------------------------------------

Let's get the dimer coordinates using :func:`~.buildBiomolecules` function:

>>> dimer = buildBiomolecules(header, monomer)
>>> dimer
<AtomGroup: 3enl biomolecule 1 (7294 atoms)>

This function takes biomolecular tarnsformations from the *header* dictionary
(item with key ``'biomolecular_transformations'``) and applies them to the 
*monomer*.  

Iterate monomers
-------------------------------------------------------------------------------

The *dimer* object now has two chains. Let's see by iterating over the chains 
in the dimer:

>>> for chain in dimer.getHierView(): print chain
Chain A
Chain B

|questions|

|suggestions|
