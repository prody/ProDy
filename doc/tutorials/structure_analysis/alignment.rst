.. _aligncoordsets:


Align coordinate sets
===============================================================================

Synopsis
-------------------------------------------------------------------------------

:class:`.AtomGroup` instances can store multiple coordinate sets,
i.e. multiple models from an NMR structure. This example shows how to align
such coordinate sets using :func:`.alignCoordsets` function.

Resulting :class:`.AtomGroup` will have its coordinate sets superposed onto
the active coordinate set selected by the user.

Parse an NMR structure
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

We use :pdb:`1joy` that contains 21 models homodimeric domain of EnvZ protein
from E. coli.

.. ipython:: python

   pdb = parsePDB('1joy')
   pdb.numCoordsets()

Calculate RMSD
-------------------------------------------------------------------------------

.. ipython:: python

   rmsds = calcRMSD(pdb)
   rmsds.mean()

This function calculates RMSDs with respect to the active coordinate set,
which is the first model in this case.


.. ipython:: python

   showProtein(pdb);
   pdb.setACSIndex(1) # model 2 in PDB is now the active coordinate set
   showProtein(pdb);
   @savefig structure_analysis_alignment_unaligned.png width=4in
   legend();


Align coordinate sets
-------------------------------------------------------------------------------

We will superpose all models onto the first model in the file using
based on Cα atom positions:

.. ipython:: python

   pdb.setACSIndex(0)
   alignCoordsets(pdb.calpha);

To use all backbone atoms, ``pdb.backbone`` can be passed as argument. See
:ref:`selections` for more information on making selections.

Coordinate sets are superposed onto the first model (the active coordinate
set).

.. ipython:: python

   rmsds = calcRMSD(pdb)
   rmsds.mean()

.. ipython:: python

   showProtein(pdb);
   pdb.setACSIndex(1) # model 2 in PDB is now the active coordinate set
   showProtein(pdb);
   @savefig structure_analysis_alignment_aligned.png width=4in
   legend();

Write aligned coordinates
-------------------------------------------------------------------------------

Using :func:`.writePDB` function, we can write the aligned
coordinate sets in PDB format:

.. ipython:: python

   writePDB('1joy_aligned.pdb', pdb)

