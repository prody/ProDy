.. currentmodule:: prody.dynamics

.. _anm:

*******************************************************************************
ANM analysis
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform ANM calculations, and retrieve
normal mode data.  

Input
-------------------------------------------------------------------------------

Protein structure data in PDB file format or specified by a PDB identifier.

Output
-------------------------------------------------------------------------------

An :class:`ANM` instance that stores Hessian and Kirchhoff matrices and 
normal mode data describing intrinsic dynamics of the protein structure. 
:class:`ANM` instances and individual normal modes 
(:class:`Mode`) can be used as input to functions in :mod:`~prody.dynamics` 
module.


ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Prepare protein
-------------------------------------------------------------------------------

We start with parsing a PDB file by passing an identifier.
Note that if a file is not found in the current working directory, it will be 
downloaded.

>>> p38 = parsePDB('1p38')
>>> p38
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>

We want to use only Cα atoms, so we select them:

>>> calphas = p38.select('protein and name CA')
>>> calphas
<Selection: "protein and name CA" from 1p38 (351 atoms; 1 coordinate sets, active set index: 0)>

We can also make the same selection like this:

>>> calphas2 = p38.select('calpha')
>>> calphas2
<Selection: "calpha" from 1p38 (351 atoms; 1 coordinate sets, active set index: 0)>

To check whether the selections are the same, we can try:

>>> calphas == calphas2
True

Note that, ProDy atom selector gives the flexibility to select any set of atoms 
to be used in ANM  calculations.

Build Hessian matrix
-------------------------------------------------------------------------------

We instantiate an :class:`ANM` instance:

>>> anm = ANM('p38 ANM analysis')

Then, build the Hessian matrix by passing selected atoms (351 Cα's)
to :meth:`ANM.buildHessian` method:

>>> anm.buildHessian(calphas)

We can get a copy of the Hessian matrix using :meth:`ANM.getHessian` method:

>>> anm.getHessian() # doctest: +SKIP
array([[ 9.959, -3.788,  0.624, ...,  0.   ,  0.   ,  0.   ],
       [-3.788,  7.581,  1.051, ...,  0.   ,  0.   ,  0.   ],
       [ 0.624,  1.051,  5.46 , ...,  0.   ,  0.   ,  0.   ],
       ..., 
       [ 0.   ,  0.   ,  0.   , ...,  1.002, -0.282,  0.607],
       [ 0.   ,  0.   ,  0.   , ..., -0.282,  3.785, -2.504],
       [ 0.   ,  0.   ,  0.   , ...,  0.607, -2.504,  4.214]])

Change cutoff distance and force constant
-------------------------------------------------------------------------------

We didn't pass any parameters to :meth:`ANM.buildHessian` method, but it 
accepts *cutoff* and *gamma* parameters, for which  default values are
``cutoff=15.0`` and ``gamma=1.0``.
 
>>> anm.getCutoff()
15.0
>>> anm.getGamma()
1.0

Note that it is also possible to use an externally calculated Hessian 
matrix. Just pass it to the ANM instance using :meth:`ANM.setHessian` method.

Calculate normal modes
-------------------------------------------------------------------------------

Calculate modes using :meth:`ANM.calcModes` method: 

>>> anm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate different number
of non-zero modes or to keep zero modes, try: ``anm.calcModes(50, zeros=True)``

Access calculated data
-------------------------------------------------------------------------------

>>> anm.getEigenvalues() # doctest: +SKIP
array([ 0.179,  0.334,  0.346,  0.791,  0.942,  1.012,  1.188,  1.304,
        1.469,  1.546,  1.608,  1.811,  1.925,  1.983,  2.14 ,  2.298,
        2.33 ,  2.364,  2.69 ,  2.794])
>>> anm.getEigenvectors() # doctest: +SKIP
array([[ 0.039, -0.045,  0.007, ...,  0.105,  0.032, -0.038],
       [ 0.009, -0.096, -0.044, ...,  0.091,  0.036, -0.037],
       [ 0.058, -0.009,  0.08 , ..., -0.188, -0.08 , -0.063],
       ..., 
       [ 0.046, -0.093, -0.131, ...,  0.018, -0.008,  0.006],
       [ 0.042, -0.018, -0.023, ...,  0.014, -0.043,  0.037],
       [ 0.08 , -0.002, -0.023, ...,  0.024, -0.023, -0.009]])


You can get the covariance matrix as follows:

>>> anm.getCovariance() # doctest: +SKIP
    array([[  2.99535077e-02,   3.21530229e-02,  -4.58288186e-03, ...,
              3.51513585e-03,   1.48465694e-03,   1.29723420e-02],
           [  3.21530229e-02,   5.92090638e-02,  -2.83705957e-02, ...,
              1.28889834e-02,  -1.19675523e-03,   6.94071253e-03],
           [ -4.58288186e-03,  -2.83705957e-02,   8.51327574e-02, ...,
             -8.65973670e-03,  -1.31024621e-03,   9.64403685e-03],
           ..., 
           [  3.51513585e-03,   1.28889834e-02,  -8.65973670e-03, ...,
              1.21472334e+00,   4.11272480e-04,  -1.74884224e-01],
           [  1.48465694e-03,  -1.19675523e-03,  -1.31024621e-03, ...,
              4.11272480e-04,   4.05594282e-01,   3.76288771e-01],
           [  1.29723420e-02,   6.94071253e-03,   9.64403685e-03, ...,
             -1.74884224e-01,   3.76288771e-01,   3.99781786e-01]])

Covariance matrices are calculated using available modes (slowest 20 modes in
this case). If user calculates M slowest modes, only they will be used in the 
calculation of covariance.

Individual modes
-------------------------------------------------------------------------------

Normal mode indices in Python start from 0, so slowest mode has index 0. 
By default, modes with zero eigenvalues are excluded. If they were retained, 
slowest non-trivial mode would have index 6.

Get the slowest mode by indexing :class:`ANM` instance as follows:

>>> slowest_mode = anm[0]
>>> slowest_mode.getEigenvalue() # doctest: +SKIP
0.17886928858119394
>>> slowest_mode.getEigenvector() # doctest: +SKIP
array([ 0.039,  0.009,  0.058, ...,  0.046,  0.042,  0.08 ])

Write NMD file
-------------------------------------------------------------------------------

ANM results in NMD format can be visualized using NMWiz VMD plugin |nmwiz|.


>>> # write slowest 3 ANM modes into an NMD file
>>> writeNMD('p38_anm_modes.nmd', anm[:3], calphas)
'p38_anm_modes.nmd'

Note that slicing an ANM (or GNM, EDA) instances returns a list of modes.
In this case, slowest 3 ANM modes were written into NMD file.

View modes in VMD
-------------------------------------------------------------------------------

First make sure that the VMD path is correct

>>> print getVMDpath()
/usr/local/bin/vmd

::

   # if this is incorrect use setVMDpath to correct it
   viewNMDinVMD('p38_anm_modes.nmd')

This will show the slowest 3 modes in VMD using NMWiz. This concludes ANM
example. Many of these apply to other NMA models, such as GNM and EDA instances.

Reduce the model
-------------------------------------------------------------------------------

Slice the model
-------------------------------------------------------------------------------

If you want to use analysis and plotting functions such as :func:`showSqFlucts`
for a subset of atoms, you can take a slice of the GNM model and pass
it to the functions:

>>> anm_res1_300, calphas_res1_300 = sliceModel(anm, calphas, 'resnum 1 to 300')
>>> anm_res1_300
<ANM: p38 ANM analysis slice "resnum 1 to 300" (20 modes, 297 nodes)>

|questions|

|suggestions|
