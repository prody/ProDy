.. currentmodule:: prody.proteins

.. _biomolt:

*******************************************************************************
Biomolecular transformations
*******************************************************************************

Synopsis
===============================================================================

Some PDB files contain coordinates for monomer from a functional/biological 
multimer. ProDy offers functions to build the structure of the multimer
using the header data from the PDB file.

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
<AtomGroup: 3enl (3647 atoms; 1 coordinate sets, active set index: 0)>

Note that we passed ``header=True`` argument to parse header data in addition
to coordinates.

Build multimer
-------------------------------------------------------------------------------

Let's get the dimer coordinates:

>>> dimer = applyBiomolecularTransformations(header, monomer)
>>> dimer
<AtomGroup: 3enl biomolecule 1 (7294 atoms; 1 coordinate sets, active set index: 0)>

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
