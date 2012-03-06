.. _pca-nmr:

*******************************************************************************
PCA of NMR Structures
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform PCA of an NMR structure with multiple models. 
The protein of interest is Ubiquitin, and example will repeat the calculations 
for ubiquitin that was published in [AB09]_.  A :class:`~.PCA` instance that 
stores covariance matrix and principal modes that describes the dominant 
changes in the dataset will be obtained. :class:`~.PCA` instance
and principal modes (:class:`~.Mode`) can be used as input to functions in 
:mod:`~.dynamics` module for further analysis.


Notes
-------------------------------------------------------------------------------

Note that this example is slightly different from that in the :ref:`tutorial`.
This example uses the :class:`~.Ensemble` class which has a method for 
performing iterative superposition.

Also, note that this example applies to any PDB file that contains multiple 
models. 
  
Prepare ensemble
===============================================================================
  
We start by importing everything from the ProDy package:

>>> from prody import *

We parse only Cα atoms using :func:`~.parsePDB` (note that it is possible to 
repeat this calculation for all atoms):
 
>>> ubi = parsePDB('2k39', subset='calpha')

We use residues 1 to 70, as residues 71 to 76 are very mobile and including
them skews the results.

>>> ubi = ubi.copy('resnum < 71')

>>> ensemble = Ensemble('Ubiquitin NMR ensemble')
>>> ensemble.setCoords( ubi.getCoords() )
	
Then, we add all of the coordinate sets to the ensemble, and perform an
iterative superposition: 
	
>>> ensemble.addCoordset( ubi.getCoordsets() ) 
>>> ensemble.iterpose()


PCA calculations
===============================================================================

Performing :class:`~.PCA` is only three lines of code:

>>> pca = PCA('Ubiquitin')
>>> pca.buildCovariance(ensemble)
>>> pca.calcModes()
>>> pca
<PCA: Ubiquitin (20 modes; 70 atoms)>


**Faster method**

Principal modes can be calculated faster using singular value decomposition:

>>> svd = PCA('Ubiquitin')
>>> svd.performSVD(ensemble)

For heterogeneous NMR datasets, both methods yields identical results:

>>> '%.3f' % abs(svd.getEigenvalues()[:20] - pca.getEigenvalues()).max()
'0.000'
>>> '%.3f' % abs(calcOverlap(pca, svd).diagonal()[:20]).min()
'1.000'

Write NMD file
===============================================================================

Write principal modes into an :term:`NMD` file for NMWiz using 
:func:`~.writeNMD` function:

>>> writeNMD('ubi_pca.nmd', pca[:3], ubi)
'ubi_pca.nmd'

Print data
===============================================================================
Let's print fraction of variance for top raking 4 PCs (listed in the Table S3):

>>> for mode in pca[:4]:
...     print calcFractVariance(mode).round(3) # doctest: +ELLIPSIS
0.134
0.094
0.083
0.065

Compare with ANM results
===============================================================================

We set the active coordinate set to 79, which is the one that is closest 
to the mean structure (note that indices start from 0 in Python).
Then, we perform ANM calculations using :func:`~.calcANM` for the active 
coordset:

>>> ubi.setACSIndex(78)
>>> anm, temp = calcANM(ubi)
>>> anm.setTitle('Ubiquitin')

We calculate overlaps between ANM and PCA modes (presented in Table 1).
:func:`~.printOverlapTable` function is handy to print a formatted overlap 
table:

>>> printOverlapTable(pca[:4], anm[:4])
Overlap Table
                         ANM Ubiquitin
                     #1     #2     #3     #4
PCA Ubiquitin #1   -0.19  -0.30  +0.22  -0.62
PCA Ubiquitin #2   +0.09  -0.72  -0.16  +0.16
PCA Ubiquitin #3   +0.31  -0.06  -0.23   0.00
PCA Ubiquitin #4   +0.11  +0.02  +0.16  -0.31
<BLANKLINE>

See Also
===============================================================================
   
User is referred to other examples in :ref:`pca-xray` for illustration of 
comparative analysis of theoretical and computational data.

|questions|

|suggestions|
