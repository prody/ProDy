.. currentmodule:: prody.proteins

.. _biomolt:

*******************************************************************************
Biomolecular transformations
*******************************************************************************

Let's say the PDB file what you have contains coordinates for monomer of 
a functional/biological dimer, and you would like to have coordinates for the dimer.
There is a function for that: :func:`applyBiomolecularTransformations`

>>> from prody import *
>>> monomer, header = parsePDB('3enl', header=True)
@> 3enl downloaded (./3enl.pdb.gz)
@> 3647 atoms and 1 coordinate sets were parsed in 0.05s.
>>> print monomer
3enl (3647 atoms; 1 coordinate sets, active set index: 0)

Let's get the dimer coordinates:

>>> dimer = applyBiomolecularTransformations(header, monomer)
>>> print dimer
Biomolecule 3enl (7294 atoms; 1 coordinate sets, active set index: 0)

Let's iterate over chains in the dimer:

>>> for chain in dimer.getHierView(): print chain
Chain A from Biomolecule 3enl (3647 atoms; 1 coordinate sets, active set index: 0)
Chain B from Biomolecule 3enl (3647 atoms; 1 coordinate sets, active set index: 0)

You can use :class:`prody.atomic.AtomGroup` instance *dimer* in anywhere you can use
*monomer*. 
