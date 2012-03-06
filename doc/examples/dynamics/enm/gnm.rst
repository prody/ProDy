.. _gnm:

*******************************************************************************
GNM analysis
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform GNM calculations using an X-ray structure 
of ubiquitin.  A :class:`~.GNM` instance that stores Kirchhoff matrix and 
normal mode data describing intrinsic dynamics of the protein structure will
be obtained.  :class:`~.GNM` instances and individual normal modes
(:class:`~.Mode`) can be used as input to functions in :mod:`~prody.dynamics` 
module.

Parse structure
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

First we parse a PDB file by passing its identifier to 
:func:`~.parsePDB` function. Note that if file is not found in 
the current working directory, it will be downloaded.


>>> ubi = parsePDB('1aar')
>>> ubi 
<AtomGroup: 1aar (1218 atoms)>

This file contains 2 chains, and a flexible C-terminal (residues 71-76).
We only want to use Cα atoms of first 70 residues from chain A, 
so we select them:

>>> calphas = ubi.select('calpha and chain A and resnum < 71')
>>> calphas 
<Selection: 'calpha and chai...and resnum < 71' from 1aar (70 atoms)>

|more| See definition of "calpha", "chain", and other selection 
keywords in :ref:`selections`.

Note that, flexible design of classes allows users to select atoms other than 
alpha carbons to be used in GNM calculations.

Build Kirchoff matrix
===============================================================================

    
Instantiate a :class:`~.GNM` instance:

>>> gnm = GNM('Ubiquitin')

We can build Kirchhoff matrix using selected atoms and 
:meth:`.GNM.buildKirchhoff` method:

>>> gnm.buildKirchhoff(calphas)


We can get a copy of the Kirchhoff matrix using :meth:`.GNM.getKirchhoff` 
method:

>>> print( gnm.getKirchhoff() ) # doctest: +ELLIPSIS
[[ 11.  -1.  -1. ...,   0.   0.   0.]
 [ -1.  15.  -1. ...,   0.   0.   0.]
 [ -1.  -1.  20. ...,   0.   0.   0.]
 ...
 [  0.   0.   0. ...,  20.  -1.  -1.]
 [  0.   0.   0. ...,  -1.  21.  -1.]
 [  0.   0.   0. ...,  -1.  -1.  12.]]

Parameters
===============================================================================

We didn't pass any parameters, but :meth:`.GNM.buildKirchhoff` method accepts 
two of them, which by default are ``cutoff=10.0`` and ``gamma=1.0``, i.e.
``buildKirchhoff(calphas, cutoff=10., gamma=1.)`` 


>>> print( gnm.getCutoff() )
10.0
>>> print( gnm.getGamma() )
1.0

Note that it is also possible to use an externally calculated Kirchhoff 
matrix. Just pass it to the GNM instance using :meth:`.GNM.setKirchhoff` method.

Calculate normal modes
===============================================================================

# calculate modes (by default slowest 20 will be calculated)
   
>>> gnm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate different number
of non-zero modes or to keep zero modes, try ``gnm.calcModes(50, zeros=True)``

Normal mode data
===============================================================================

Get eigenvalues and eigenvectors:

>>> print( gnm.getEigenvalues().round(3) )
[  2.502   2.812   4.366   5.05    7.184   7.65    7.877   9.08    9.713
  10.132  10.502  10.644  10.888  11.157  11.285  11.632  11.78   11.936
  12.006  12.218]
>>> print( gnm.getEigenvectors().round(3) ) # doctest: +ELLIPSIS
[[-0.064 -0.131 -0.245 ..., -0.256  0.538 -0.   ]
 [-0.073 -0.085 -0.19  ...,  0.006 -0.069  0.032]
 [-0.076 -0.043 -0.135 ...,  0.017 -0.047  0.018]
 ...
 [-0.092  0.064  0.105 ...,  0.032 -0.042  0.006]
 [-0.07   0.099  0.054 ...,  0.031  0.024 -0.014]
 [-0.081  0.135  0.124 ...,  0.013 -0.04  -0.018]]

Get covariance matrix:

>>> print( gnm.getCovariance().round(2) ) # doctest: +ELLIPSIS
[[ 0.08  0.02  0.01 ..., -0.01 -0.01 -0.01]
 [ 0.02  0.02  0.01 ..., -0.   -0.   -0.01]
 [ 0.01  0.01  0.01 ...,  0.   -0.   -0.  ]
 ...
 [-0.01 -0.    0.   ...,  0.01  0.01  0.01]
 [-0.01 -0.   -0.   ...,  0.01  0.01  0.02]
 [-0.01 -0.01 -0.   ...,  0.01  0.02  0.05]]
              
Note that covariance matrices are calculated using available modes 
in the model, which is slowest 20 modes in this case. 
If user calculated M slowest modes, only they will be used in the 
calculation of covariance.

Individual modes
===============================================================================

Normal mode indices start from 0, so slowest mode has index 0. 

>>> slowest_mode = gnm[0]
>>> print( slowest_mode.getEigenvalue().round(3) )
2.502
>>> print( slowest_mode.getEigenvector().round(3) ) # doctest: +ELLIPSIS
[-0.064 -0.073 -0.076 -0.112 -0.092 -0.143 -0.164 -0.205 -0.24  -0.313
 -0.192 -0.152 -0.066 -0.07  -0.025 -0.031  0.001 -0.006 -0.015  0.027
  0.042  0.055  0.063  0.09   0.09   0.069  0.132  0.175  0.145  0.121
  0.195  0.218  0.158  0.217  0.245  0.214  0.225  0.171  0.2    0.151
  0.102  0.043 -0.029 -0.064 -0.072 -0.086 -0.09  -0.078 -0.057 -0.011
  0.016  0.061  0.058  0.043  0.029  0.013  0.004  0.011 -0.013 -0.037
 -0.05  -0.059 -0.07  -0.094 -0.094 -0.099 -0.097 -0.092 -0.07  -0.081]

By default, modes with 0 eigenvalue are excluded. If they were retained, 
slowest non-trivial mode would have index 6.

Plot results
===============================================================================

.. plot::
   :context:
   :nofigs:

   from prody import *
   ubi = parsePDB('1aar')
   calphas = ubi.select('calpha and chain A and resnum < 71')
   gnm = GNM('Ubiquitin')
   gnm.buildKirchhoff(calphas)
   gnm.calcModes()

ProDy plotting functions are prefixed with ``show``. Let's use some of them
to plot data:

**Contact Map**

.. plot::
   :context:
   :include-source:
   
   # We import plotting library
   import matplotlib.pyplot as plt
   
   plt.figure(figsize=(5,4))
   showContactMap(gnm)
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  

**Cross-correlations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showCrossCorr(gnm)
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  

**Slow mode shape**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showMode(gnm[0])
   plt.grid()
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  

**Square fluctuations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showSqFlucts(gnm[0])
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  

|questions|

|suggestions|
