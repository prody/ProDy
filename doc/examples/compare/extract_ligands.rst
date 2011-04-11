.. currentmodule:: prody.compare

.. _extract-ligands:

*******************************************************************************
Extract ligands
*******************************************************************************

Synopsis
=============================================================================

:func:`matchAlign` function can be used for aligning protein structures.
This example shows to to use it to extract ligands from multiple PDB structures 
after superposing the structures onto a reference.

.. versionadded:: 0.7.1

Input
-------------------------------------------------------------------------------

A reference PDB file.

Output
-------------------------------------------------------------------------------

Output is PDB files that contain ligands superposed onto the reference
structure.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Parse reference and blast PDB
-------------------------------------------------------------------------------

First, we parse the reference structure and blast search PDB for similar 
structure::

  p = parsePDB('1p38')
  seq = p.getHierView()['A'].getSequence()
  b = blastPDB(seq)

Align structures and extract ligands
-------------------------------------------------------------------------------

Then, we parse the hits one-by-one, superpose them onto the reference 
structure, and extract ligands::

  for pdb in b.getHits().keys():
      # blast search may return PDB identifiers of deprecated structures,
      # so we parse structures within a try statement
      try:
          m = parsePDB(pdb)
          m = matchAlign(m,p)[0] 
      except:
          continue
      s = m.select('not protein and not water')
      if s is not None:
          writePDB(pdb+'_ligand.pdb', s)

This code will write PDB files that contain non-protein and non-water atoms.
The output can be loaded into a molecular visualization tool for analysis.

|questions|

|suggestions|
