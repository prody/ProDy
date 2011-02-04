.. currentmodule:: prody.proteins

.. _aligncoordsets:

*******************************************************************************
Align coordinate sets
*******************************************************************************

Synopsis
===============================================================================

:class:`~prody.atomic.AtomGroup` instances can store multiple coordinate sets,
i.e. multiple models from an NMR structure. This example shows how to align
such coordinate sets. 

User input
-------------------------------------------------------------------------------

A PDB file that contains multiple NMR models.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Parse an NMR structure
-------------------------------------------------------------------------------
   
We use 1joy that contains 21 models homodimeric domain of EnvZ protein 
from E. coli.

>>> pdb = parsePDB('1joy')
>>> print pdb.getNumOfCoordsets()
21

Calculate RMSD
-------------------------------------------------------------------------------
   
>>> rmsds = calcRMSD( pdb )
>>> rmsds.mean() # doctest: +SKIP
37.5069116784

This function calculates RMSDs with respect to the active coordinate set,
which is the first model in this case.



Align coordsets
-------------------------------------------------------------------------------

We will superpose all models onto the first model in the file using
based on alpha carbon atom positions:
   
>>> alignCoordsets( pdb, selstr='calpha')

``selstr='calpha'`` argument selects alpha carbon atoms. To use all backbone
atoms, ``selstr='backbone'`` can be passed as argument. See :ref:`selections`
for more information on making selections.

Coordinate sets are superposed onto the first model (the active coordinate 
set).
   
>>> rmsds = calcRMSD( pdb )
>>> rmsds.mean() # doctest: +SKIP
3.27689121518

|questions|

|suggestions|
