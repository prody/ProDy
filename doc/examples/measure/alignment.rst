.. _aligncoordsets:

*******************************************************************************
Align coordinate sets
*******************************************************************************

Synopsis
===============================================================================

:class:`~.AtomGroup` instances can store multiple coordinate sets,
i.e. multiple models from an NMR structure. This example shows how to align
such coordinate sets using :func:`~.alignCoordsets` function. 

Resulting :class:`~.AtomGroup` will have its coordinate sets superposed onto 
the active coordinate set selected by the user.

Parse an NMR structure
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

We use 1joy that contains 21 models homodimeric domain of EnvZ protein 
from E. coli.

>>> pdb = parsePDB('1joy')
>>> print( pdb.numCoordsets() )
21

Calculate RMSD
===============================================================================
   
>>> rmsds = calcRMSD( pdb )
>>> print( rmsds.mean().round(2) )
37.51

This function calculates RMSDs with respect to the active coordinate set,
which is the first model in this case.

Align coordinate sets
===============================================================================

We will superpose all models onto the first model in the file using
based on Cα atom positions:
   
>>> alignCoordsets(pdb.calpha)
<Selection: 'calpha' from 1joy (134 atoms; active #0 of 21 coordsets)>

To use all backbone atoms, ``pdb.backbone`` can be passed as argument. See 
:ref:`selections` for more information on making selections.

Coordinate sets are superposed onto the first model (the active coordinate 
set).
   
>>> rmsds = calcRMSD(pdb)
>>> print(rmsds.mean().round(2))
3.28

Write aligned coordinates
===============================================================================

Using :func:`~.writePDB` function, we can write the aligned
coordinate sets in PDB format: 

>>> writePDB('1joy_aligned.pdb', pdb)
'1joy_aligned.pdb'

|questions|

|suggestions|
