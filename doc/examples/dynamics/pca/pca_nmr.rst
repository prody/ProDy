.. currentmodule:: prody.dynamics

.. _pca-nmr:

*******************************************************************************
PCA of NMR Structures
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform PCA of an NMR structure with multiple models. 
The protein of interest is Ubiquitin, and example will repeat the calculations 
for ubiquitin that was published in [AB09]_.

User input
-------------------------------------------------------------------------------

A PDB identifier for an NMR structure.
  
ProDy Code
===============================================================================
  
We start by importing everything from the ProDy package:

>>> from prody import *

Prepare protein
-------------------------------------------------------------------------------

We parse only CA atoms using :func:`~prody.proteins.parsePDB` 
(note that it is possible to repeat this calculation for all atoms):
 
>>> ubi = parsePDB('2k39', subset='calpha')

We use residues 1 to 70, as residues 71 to 76 are very mobile and including
them skews the results.

>>> ubi = ubi.copy('resnum < 71')


Align models
-------------------------------------------------------------------------------

We use :func:`~prody.measure.alignCoordsets` function to superimpose all
models onto the first one.

>>> alignCoordsets(ubi)

PCA calculations
-------------------------------------------------------------------------------

Performing :class:`PCA` is only three lines of code:

>>> pca = PCA('Ubiquitin')
>>> pca.buildCovariance(ubi)
>>> pca.calcModes()
>>> pca
<PCA: Ubiquitin (20 modes, 70 atoms)>

.. note::
   Note than in this PCA example we did not use :class:`~prody.ensemble.Ensemble`
   class. This is because the models in an NMR structure file contains 
   coordinates of all atoms, which makes :class:`~prody.ensemble.Ensemble`
   not so useful. In PCA calculations for heterogeneous structural datasets, 
   on the other hand, we deal with structures with missing residues. In that
   case :class:`~prody.ensemble.Ensemble` class becomes handy to assemble
   the coordinate data and perform structural alignment despite missing atoms.
   

Write NMD file
-------------------------------------------------------------------------------

Write principal modes into an NMD file for NMWiz using :func:`writeNMD` 
function:

>>> writeNMD('ubi_pca.nmd', pca[:3], ubi)
'ubi_pca.nmd'

Print data
-------------------------------------------------------------------------------
Let's print fraction of variance for top raking 4 PCs (listed in the Table S3):

>>> for mode in pca[:4]:
...     print mode.getFractOfVariance() # doctest: +SKIP
0.299016803492
0.0959780950608
0.0647918823066
0.058247703612


Compare with ANM results
-------------------------------------------------------------------------------

We set the active coordinate set to 79, which is the one that is closest 
to the mean structure (note that indices start from 0 in Python).
Then, we perform ANM calculations using :func:`calcANM` for the active coordset:

>>> ubi.setActiveCoordsetIndex(78)
>>> anm = calcANM(ubi)
>>> anm.setName('Ubiquitin')

We calculate overlaps between ANM and PCA modes (presented in Table 1).
:func:`printOverlapTable` function is handy to print a formatted overlap table:

>>> printOverlapTable(pca[:4], anm[:4]) # doctest: +SKIP
Overlap Table
                         ANM Ubiquitin
                     #1     #2     #3     #4
PCA Ubiquitin #1   +0.02  -0.16  -0.12  -0.10
PCA Ubiquitin #2   +0.19  +0.35  -0.20  +0.63
PCA Ubiquitin #3   -0.20  +0.65  +0.24  -0.26
PCA Ubiquitin #4   -0.27  -0.16  +0.14  +0.21

