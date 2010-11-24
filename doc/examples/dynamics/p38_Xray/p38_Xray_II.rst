.. currentmodule:: prody.dynamics

.. _p38-Xray-analysis:

*******************************************************************************
p38 X-ray Ensemble - Part II: Analysis
*******************************************************************************

|more| Continued from :ref:`p38-PCA-ANM-I`.

The aim of this part is to print/save the numerical data that was presented
in our paper [AB09]_.

If you have left your session in the previous part, you will need to
load saved data::

  from prody import *

  pca = loadModel('p38_xray.pca.npz')
  anm = loadModel('1p38.anm.npz')

Variance along PCs
===============================================================================

Of interest is the fraction of variance that is explained by principal 
components, which are the dominant modes of variability in the dataset.
We can print this information to screen for top 6 PCs as follows::

  for mode in pca[:3]:    # Print % variance explained by top PCs
      print '{0:s}  % variance = {1:.2f}'.format( mode, mode.getFractOfVariance()*100 )
  
This prints::

  Mode 1 from PCA p38 xray  % variance = 29.19
  Mode 2 from PCA p38 xray  % variance = 16.51
  Mode 3 from PCA p38 xray  % variance = 10.63

  
This data was included in Table 1 in [AB09]_.

Collectivity of modes 
===============================================================================

Collectivity of a normal mode ([BR95]_) can be obtained using 
:meth:`Mode.getCollectivity`::

  for mode in pca[:3]:    # Print PCA mode collectivity
      print '{0:s}  collectivity = {1:.2f}'.format( mode, mode.getCollectivity() )

This prints::

  Mode 1 from PCA p38 xray  collectivity = 0.50
  Mode 2 from PCA p38 xray  collectivity = 0.50
  Mode 3 from PCA p38 xray  collectivity = 0.32

Similarly, we can get collectivity of ANM modes::

  for mode in anm[:3]:    # Print ANM mode collectivity
      print '{0:s}  collectivity = {1:.2f}'.format( mode, mode.getCollectivity() )

This prints::

  Mode 1 from ANM 1p38  collectivity = 0.65
  Mode 2 from ANM 1p38  collectivity = 0.55
  Mode 3 from ANM 1p38  collectivity = 0.68

This show that top PCA modes and slow ANM modes are highly collective.

PCA - ANM overlap  
===============================================================================

We also calculated overlap between PCA and ANM modes to see whether 
structural changes observed upon inhibitor binding correlated with 
the intrinsic fluctuations of the p38 MAP kinase (Table 1 in [AB09]_).

There are a number of functions to calculate or show overlaps between modes 
(see list of them in :ref:`dynamics`). In an interactive session, most useful is 
:func:`printOverlapTable`::

  printOverlapTable(pca[:3], anm[:3])   # Compare top 3 PCs and slowest 3 ANM modes

This prints::

  Overlap Table
                          ANM 1p38
                      #1     #2     #3
  PCA p38 xray #1   -0.39  +0.04  -0.71
  PCA p38 xray #2   -0.78  -0.20  +0.22
  PCA p38 xray #3   +0.05  -0.57  +0.06

This formatted table can also be written into a file using 
:func:`writeOverlapTable` function. 

Save numeric data
===============================================================================

:class:`ANM` and :class:`PCA` instances store calculated numeric data. Their 
class documentation lists methods that return eigenvalue, eigenvector, 
covariance matrix etc. data to the user. Such data can easily be written into
text files for analysis using external software; The function is to use is 
:func:`writeArray`::

  writeArray( 'p38_PCA_eigvectors.txt', pca.getEigenvectors() ) # PCA eigenvectors
  
  writeModes( 'p38_ANM_modes.txt', anm ) # ANM eigenvectors, same as using above function
  
  # It is also possible to write arbitrary arrays
  overlap = getOverlap(pca[:3], anm[:3])
  writeArray( 'p38_PCA_ANM_overlap.txt', abs(overlap), format='.2f')
  
|more| See all functions and class methods in :ref:`dynamics`.

|more| This example continues in :ref:`p38-PCA-ANM-III`
