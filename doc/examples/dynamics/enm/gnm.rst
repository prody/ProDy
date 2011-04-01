.. currentmodule:: prody.dynamics

.. _gnm:

*******************************************************************************
GNM analysis
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform GNM calculations using an X-ray structure 
of ubiquitin.  

Input
-------------------------------------------------------------------------------

Protein structure data in PDB file format or specified by a PDB identifier.

Output
-------------------------------------------------------------------------------

A :class:`GNM` instance that stores Kirchhoff matrix and normal mode data 
describing intrinsic dynamics of the protein structure.
:class:`GNM` instances and individual normal modes (:class:`Mode`) can be
used as input to functions in :mod:`~prody.dynamics` module.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Prepare protein
-------------------------------------------------------------------------------

    
First we parse a PDB file by passing its identifier to 
:func:`~prody.proteins.parsePDB` function. Note that if file is not found in 
the current working directory, it will be downloaded.


>>> ubi = parsePDB('1aar')
>>> ubi 
<AtomGroup: 1aar (1218 atoms; 1 coordinate sets, active set index: 0)>

This file contains 2 chains, and a flexible C-terminal (residues 71-76).
We only want to use Cα atoms of first 70 residues from chain A, 
so we select them:

>>> calphas = ubi.select('calpha and chain A and resnum < 71')
>>> calphas 
<Selection: "calpha and chai...and resnum < 71" from 1aar (70 atoms; 1 coordinate sets, active set index: 0)>

|more| See definition of "calpha", "chain", and other selection 
keywords in :ref:`selections`.

Note that, flexible design of classes allows users to select atoms other than 
alpha carbons to be used in GNM calculations.

Build Kirchoff matrix
------------------------------------------------------------------------------- 

    
Instantiate a :class:`GNM` instance:

>>> gnm = GNM('Ubiquitin')

We can build Kirchhoff matrix using selected atoms and :meth:`GNM.buildKirchhoff`
method:

>>> gnm.buildKirchhoff(calphas)


We can get a copy of the Kirchhoff matrix using :meth:`GNM.getKirchhoff` method:

>>> gnm.getKirchhoff() # doctest: +SKIP
array([[ 11.,  -1.,  -1., ...,   0.,   0.,   0.],
       [ -1.,  15.,  -1., ...,   0.,   0.,   0.],
       [ -1.,  -1.,  20., ...,   0.,   0.,   0.],
       ..., 
       [  0.,   0.,   0., ...,  20.,  -1.,  -1.],
       [  0.,   0.,   0., ...,  -1.,  21.,  -1.],
       [  0.,   0.,   0., ...,  -1.,  -1.,  12.]])

Change cutoff distance and force constant
-------------------------------------------------------------------------------


We didn't pass any parameters, but :meth:`~GNM.buildKirchhoff` method accepts 
two of them, which by default are ``cutoff=10.0`` and ``gamma=1.0``, i.e.
``buildKirchhoff(calphas, cutoff=10., gamma=1.)`` 


>>> gnm.getCutoff()
10.0
>>> gnm.getGamma()
1.0

Note that it is also possible to use an externally calculated Kirchhoff 
matrix. Just pass it to the GNM instance using :meth:`~GNM.setKirchhoff` method.

Calculate normal modes
-------------------------------------------------------------------------------

# calculate modes (by default slowest 20 will be calculated)
   
>>> gnm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate different number
of non-zero modes or to keep zero modes, try ``gnm.calcModes(50, zeros=True)``

Access calculated data
-------------------------------------------------------------------------------

Get eigenvalues and eigenvectors:

>>> gnm.getEigenvalues() # doctest: +SKIP
array([  2.50159996,   2.81198884,   4.36615755,   5.04996398,
         7.18368407,   7.64988589,   7.87722993,   9.08034538,
         9.71281297,  10.13238988,  10.50197833,  10.64403971,
        10.88838721,  11.15731079,  11.2850227 ,  11.63219978,
        11.7801197 ,  11.93585376,  12.00584891,  12.2183852 ])
>>> gnm.getEigenvectors() # doctest: +SKIP
array([[ -6.38051740e-02,  -1.30544492e-01,  -2.45334248e-01, ...,
         -2.56196713e-01,   5.37771418e-01,  -1.08369816e-04],
       [ -7.25770215e-02,  -8.45109261e-02,  -1.90052460e-01, ...,
          6.24201652e-03,  -6.93922518e-02,   3.22234776e-02],
       [ -7.58827175e-02,  -4.30040725e-02,  -1.34783179e-01, ...,
          1.66397155e-02,  -4.70341843e-02,   1.77800008e-02],
       ..., 
       [ -9.21438212e-02,   6.36544251e-02,   1.04633107e-01, ...,
          3.21057282e-02,  -4.19619569e-02,   6.20810516e-03],
       [ -7.00439044e-02,   9.90390676e-02,   5.42931923e-02, ...,
          3.08618840e-02,   2.37876216e-02,  -1.36026733e-02],
       [ -8.08918894e-02,   1.34634053e-01,   1.24148042e-01, ...,
          1.25198059e-02,  -3.99479375e-02,  -1.84106396e-02]])


Get covariance matrix:

>>> gnm.getCovariance() # doctest: +SKIP
array([[ 0.08038538,  0.0214721 ,  0.01154276, ..., -0.00810262,
        -0.0068922 , -0.00519627],
       [ 0.0214721 ,  0.01841361,  0.01294682, ..., -0.00164253,
        -0.0030955 , -0.00590236],
       [ 0.01154276,  0.01294682,  0.01015144, ...,  0.00079905,
        -0.00034443, -0.00350823],
       ..., 
       [-0.00810262, -0.00164253,  0.00079905, ...,  0.01120008,
         0.00810447,  0.01423867],
       [-0.0068922 , -0.0030955 , -0.00034443, ...,  0.00810447,
         0.00920142,  0.01602405],
       [-0.00519627, -0.00590236, -0.00350823, ...,  0.01423867,
         0.01602405,  0.05269755]])


   
Note that covariance matrices are calculated using available modes 
in the model, which is slowest 20 modes in this case. 
If user calculated M slowest modes, only they will be used in the 
calculation of covariance.

Individual modes
-------------------------------------------------------------------------------

Normal mode indices start from 0, so slowest mode has index 0. 

>>> slowest_mode = gnm[0]
>>> print slowest_mode.getEigenvalue() # doctest: +SKIP
2.50159995996
>>> print slowest_mode.getEigenvector() # doctest: +SKIP
[-0.06380517 -0.07257702 -0.07588272 -0.11243047 -0.09221726 -0.14296277
 -0.16396628 -0.2053729  -0.24020989 -0.3134546  -0.19166723 -0.15241052
 -0.0657869  -0.07000933 -0.0253338  -0.03129034  0.0007266  -0.0058239
 -0.01543903  0.02655975  0.04238382  0.05465168  0.06268858  0.08997896
  0.08987359  0.06877047  0.13222327  0.1748625   0.14497343  0.12132971
  0.19509205  0.21771503  0.15807374  0.21746306  0.24530097  0.21361762
  0.22544686  0.1707908   0.1998145   0.15089032  0.10228462  0.04310205
 -0.028911   -0.06447926 -0.0721635  -0.08602533 -0.09012947 -0.07761776
 -0.05653068 -0.01097932  0.0162001   0.06091741  0.05773185  0.04274781
  0.02903195  0.01313377  0.00376044  0.01133737 -0.01299074 -0.03725209
 -0.05023782 -0.05862394 -0.07037726 -0.0936061  -0.09401438 -0.09896954
 -0.09684669 -0.09214382 -0.0700439  -0.08089189]


By default,
modes with 0 eigenvalue are excluded. If they were retained, slowest 
non-trivial mode would have index 6.

Plot results
-------------------------------------------------------------------------------

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
   showCrossCorrelations(gnm)
   
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

Reduce the model
-------------------------------------------------------------------------------


Slice the model
-------------------------------------------------------------------------------

If you want to use analysis and plotting functions such as :func:`showSqFlucts`
for a subset of atoms, you can take a slice of the GNM model and pass
it to the functions:

>>> gnm_res1_40, calphas_res1_40 = sliceModel(gnm, calphas, 'resnum 1 to 40')
>>> gnm_res1_40
<GNM: Ubiquitin slice "resnum 1 to 40" (20 modes, 40 nodes)>


|questions|

|suggestions|
