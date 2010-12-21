.. currentmodule:: prody.dynamics

.. _pca-nmr:

*******************************************************************************
PCA of Ubiquitin NMR Structure
*******************************************************************************

This example repeats the calculations for ubiquitin that was 
published in [AB09]_.

Following classes and functions are used in the example:

Classes:
  * :class:`PCA`
  * :class:`~prody.ensemble.Ensemble`
Functions:
  * :func:`~prody.proteins.parsePDB`
  * :func:`calcANM`
  * :func:`printOverlapTable`
  * :func:`writeNMD`
  
>>> from prody import *

We start with parsing only CA atoms:
 
>>> ubi = parsePDB('2k39', subset='calpha')

We use residues 1 to 70. 71 and above are disordered.

>>> ubi = ubi.copy('resnum < 71')

We instantiate an ensemble and set the reference coordinates 
(coordinates from model 1).

>>> ensemble = Ensemble('Ubiquitin NMR ensemble')
>>> ensemble.setCoordinates( ubi.getCoordinates() )

Then, we add all of the coordinate sets to the ensemble, and perform an
iterative superposition: 

>>> ensemble.addCoordset( ubi.getCoordsets() ) 
>>> ensemble.iterpose()

PCA calculations:

>>> pca = PCA('ubi PCA')
>>> pca.buildCovariance(ensemble)
>>> pca.calcModes()
>>> # Write principal modes into an NMD file for NMWiz
>>> writeNMD('ubi_pca.nmd', pca[:3], ubi)
'ubi_pca.nmd'

Let's print fraction of variance for top raking 4 PCs (listed in the Table S3):

>>> for mode in pca[:4]:
...     print mode.getFractOfVariance() # doctest: +SKIP
0.133691677009
0.0942276000043
0.0833640062736
0.0654647139302

We set the active coordinate set to 79, which is the one that is closest 
to the mean structure (note that indices start from 0 in Python).
Then, we perform ANM calculations for the active coordset:

>>> ubi.setActiveCoordsetIndex(78)
>>> anm = calcANM(ubi)
>>> anm.setName('ubi ANM')

We calculate overlaps between ANM and PCA modes (presented in Table 1):

>>> printOverlapTable(pca[:4], anm[:4]) # doctest: +SKIP
Overlap Table
                        ANM ubi ANM
                   #1     #2     #3     #4
PCA ubi PCA #1   -0.19  -0.30  +0.22  -0.62
PCA ubi PCA #2   +0.09  -0.72  -0.16  +0.16
PCA ubi PCA #3   +0.31  -0.06  -0.23  -0.00
PCA ubi PCA #4   +0.11  +0.02  +0.16  -0.31

Finally, we do some cleaing:

>>> import os
>>> os.remove('ubi_pca.nmd')
>>> os.remove('2k39.pdb.gz')
