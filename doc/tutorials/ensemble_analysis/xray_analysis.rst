.. _pca-xray-analysis:


Analysis
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example is continued from :ref:`pca-xray-calculations`.  The aim of this
part is to perform a quantitative comparison of experimental and theoretical
data and to print/save the numerical data that were presented in [AB09]_.


We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Then we load data saved in the previous part:

.. ipython:: python

   pca = loadModel('p38_xray.pca.npz')
   anm = loadModel('1p38.anm.npz')

Variance along PCs
-------------------------------------------------------------------------------

Of interest is the fraction of variance that is explained by principal
components, which are the dominant modes of variability in the dataset.
We can print this information to screen for top 6 PCs as follows:

.. ipython:: python

   for mode in pca[:3]:    # Print % variance explained by top PCs
       var = calcFractVariance(mode)*100
       print('{0:s}  % variance = {1:.2f}'.format(mode, var))


These data were included in Table 1 in [AB09]_.

Collectivity of modes
-------------------------------------------------------------------------------

Collectivity of a normal mode ([BR95]_) can be obtained using
:func:`.calcCollectivity`:

.. ipython:: python

   for mode in pca[:3]:    # Print PCA mode collectivity
       coll = calcCollectivity(mode)
       print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))

Similarly, we can get collectivity of ANM modes:

.. ipython:: python

   for mode in anm[:3]:    # Print ANM mode collectivity
       coll = calcCollectivity(mode)
       print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))

This shows that top PCA modes and slow ANM modes are highly collective.

PCA - ANM overlap
-------------------------------------------------------------------------------

We calculate the overlap between PCA and ANM modes in order to see whether
structural changes observed upon inhibitor binding correlate with
the intrinsic fluctuations of the p38 MAP kinase (Table 1 in [AB09]_).

There are a number of functions to calculate or show overlaps between modes
(see list of them in :ref:`dynamics`). In an interactive session, most useful
is :func:`.printOverlapTable`:

.. ipython:: python

   printOverlapTable(pca[:3], anm[:3]) # Top 3 PCs vs slowest 3 ANM modes


This formatted table can also be written into a file using
:func:`.writeOverlapTable` function.

Save numeric data
-------------------------------------------------------------------------------

:class:`.ANM` and :class:`.PCA` instances store calculated numeric data.
Their class documentation lists methods that return eigenvalue, eigenvector,
covariance matrix etc. data to the user. Such data can easily be written into
text files for analysis using external software. The function is to use is
:func:`.writeArray`:

.. ipython:: python

   writeArray('p38_PCA_eigvecs.txt', pca.getEigvecs() ) # PCA eigenvectors
   writeModes('p38_ANM_modes.txt', anm) # ANM modes, same as using above func


It is also possible to write arbitrary arrays:

.. ipython:: python

   overlap = calcOverlap(pca[:3], anm[:3])
   writeArray('p38_PCA_ANM_overlap.txt', abs(overlap), format='%.2f')


