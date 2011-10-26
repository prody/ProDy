.. currentmodule:: prody.measure

.. _aligncoordsets:

*******************************************************************************
Align coordinate sets
*******************************************************************************

Synopsis
===============================================================================

:class:`~prody.atomic.AtomGroup` instances can store multiple coordinate sets,
i.e. multiple models from an NMR structure. This example shows how to align
such coordinate sets using :func:`alignCoordsets` function. 

Input
-------------------------------------------------------------------------------

Multiple structures/models of the same protein in a single PDB file.

Output
-------------------------------------------------------------------------------

Output is a :class:`~prody.atomic.AtomGroup` instance whose coordinate
sets are superposed on the the active coordinate set selected by the user.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Parse an NMR structure
-------------------------------------------------------------------------------
   
We use 1joy that contains 21 models homodimeric domain of EnvZ protein 
from E. coli.

>>> pdb = parsePDB('1joy')
>>> print( pdb.numCoordsets() )
21

Calculate RMSD
-------------------------------------------------------------------------------
   
>>> rmsds = calcRMSD( pdb )
>>> print( rmsds.mean().round(2) )
37.51

This function calculates RMSDs with respect to the active coordinate set,
which is the first model in this case.

Align coordsets
-------------------------------------------------------------------------------

We will superpose all models onto the first model in the file using
based on Cα atom positions:
   
>>> alignCoordsets( pdb, selstr='calpha')

``selstr='calpha'`` argument selects Cα atoms. To use all backbone
atoms, ``selstr='backbone'`` can be passed as argument. See :ref:`selections`
for more information on making selections.

Coordinate sets are superposed onto the first model (the active coordinate 
set).
   
>>> rmsds = calcRMSD( pdb )
>>> print( rmsds.mean().round(2) )
3.28

Write aligned coordinates
-------------------------------------------------------------------------------

Using :func:`~prody.proteins.writePDB` function, we can write the aligned
coordinate sets in PDB format: 

>>> writePDB('1joy_aligned.pdb', pdb)
'1joy_aligned.pdb'

|questions|

|suggestions|
