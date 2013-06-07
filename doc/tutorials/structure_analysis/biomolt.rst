.. _biomolt:

Building Biomolecules
===============================================================================

Some PDB files contain coordinates for a monomer of a functional/biological
multimer (biomolecule).  ProDy offers functions to build structures of
biomolecules using the header data from the PDB file.  We will use PDB file
that contains the coordinates for a monomer of a biological
multimeric protein and the transformations in the header section to
generate the multimer coordinates.  Output will be an :class:`~.AtomGroup`
instance that contains the multimer coordinates.

Parse a PDB file
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

We parse a PDB file that contains the coordinates for a monomer of a dimeric
protein:

.. ipython:: python

   monomer, header = parsePDB('3enl', header=True)
   repr(monomer)

Note that we passed ``header=True`` argument to parse header data in addition
to coordinates.

.. ipython:: python

   showProtein(monomer);
   @savefig structure_analysis_biomolt_monomer.png width=4in
   legend();

Build multimer
-------------------------------------------------------------------------------

Let's get the dimer coordinates using :func:`~.buildBiomolecules` function:

.. ipython:: python

   dimer = buildBiomolecules(header, monomer)
   repr(dimer)


This function takes biomolecular tarnsformations from the *header* dictionary
(item with key ``'biomoltrans'``) and applies them to the
*monomer*.

.. ipython:: python

   showProtein(dimer);
   @savefig structure_analysis_biomolt_dimer.png width=4in
   legend();

Iterate monomers
-------------------------------------------------------------------------------

The *dimer* object now has two chains. Let's see by iterating over the chains
in the dimer:

.. ipython:: python

   for chain in dimer.iterChains():
       repr(chain)
