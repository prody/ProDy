.. currentmodule:: prody.dynamics

.. _anm:

*******************************************************************************
Anisotropic Network Model analysis
*******************************************************************************

This example shows how to perform ANM calculations.  

Prepare protein
===============================================================================

>>> from prody import *

>>> # parse pdb file by identifier
>>> # if file is not found in current working directory, it will be downloaded
>>> p38 = parsePDB('1p38')
>>> p38
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>

>>> # we only want to use alpha carbons, so we select them
>>> calphas = p38.select('protein and name CA')
>>> # we can also make the same selection like this
>>> calphas = p38.select('calpha')
>>> calphas 
<Selection: "calpha" from 1p38 (351 atoms; 1 coordinate sets, active set index: 0)>

Note that, atoms other than alpha carbons can be selected and used in ANM 
calculations

Perform calculations
===============================================================================

**Build Hessian**

>>> # instantiate ANM instance
>>> anm = ANM('p38 ANM analysis')

>>> # build Hessian matrix using selected atoms (351 alpha carbons)
>>> anm.buildHessian(calphas)

>>> # get a copy of the hessian matrix
>>> anm.getHessian() # doctest: +SKIP
array([[ 9.959, -3.788,  0.624, ...,  0.   ,  0.   ,  0.   ],
       [-3.788,  7.581,  1.051, ...,  0.   ,  0.   ,  0.   ],
       [ 0.624,  1.051,  5.46 , ...,  0.   ,  0.   ,  0.   ],
       ..., 
       [ 0.   ,  0.   ,  0.   , ...,  1.002, -0.282,  0.607],
       [ 0.   ,  0.   ,  0.   , ..., -0.282,  3.785, -2.504],
       [ 0.   ,  0.   ,  0.   , ...,  0.607, -2.504,  4.214]])

**How to change cutoff and gamma parameters?**

We didn't pass any parameters, but ``buildHessian`` method accepts two of them
by default ``cutoff=15.0`` and ``gamma=1.0`` is passed
that is, ``buildHessian(calphas, cutoff=15., gamma=1.)``.
 
>>> anm.getCutoff()
15.0
>>> anm.getGamma()
1.0

Note that it is also possible to use an externally calculated Hessian 
matrix. Just pass it to the ANM instance using :meth:`~ANM.setHessian` method.

**Calculate normal modes**

>>> # calculate modes (by default slowest 20 will be calculated)
>>> anm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate different number
of non-zero modes or to keep zero modes, try: ``anm.calcModes(50, zeros=True)``

Access calculated data
===============================================================================

You can get the covariance matrix as follows:
>>> anm.getCovariance()

Note that covariance calculated using 20 modes when you call this method.

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

Covariance matrices are calculated using available modes. If user calculates
M slowest modes, only they will be used in the calculation of covariance.

Individual modes
===============================================================================

Normal mode indices in Python start from 0, so slowest mode has index 0. By default,
modes with zero eigenvalues are excluded. If they were retained, slowest 
non-trivial mode would have index 6.

>>> # investigate individual modes
>>> slowest_mode = anm[0]

>>> slowest_mode.getEigenvalue() # doctest: +SKIP
0.17886928858119394
>>> slowest_mode.getEigenvector() # doctest: +SKIP
array([ 0.039,  0.009,  0.058, ...,  0.046,  0.042,  0.08 ])

Write NMD file
===============================================================================

ANM results in NMD format can be visualized using |nmwiz|. NMWiz is
a |vmd| plugin.


>>> # write slowest 3 ANM modes into an NMD file
>>> writeNMD('p38_anm_modes.nmd', anm[:3], calphas)
'p38_anm_modes.nmd'

Note that slicing an ANM (or GNM, EDA) instances returns a list of modes.
In this case, slowest 3 ANM modes were written into NMD file.

View modes in VMD
===============================================================================

>>> # first make sure that the VMD path is correct
>>> print getVMDpath()
/usr/local/bin/vmd

::

   # if this is incorrect use setVMDpath to correct it
   viewNMDinVMD('p38_anm_modes.nmd')

This will show the slowest 3 modes in VMD using NMWiz. This concludes ANM
example. Many of these apply to other NMA models, such as GNM and EDA instances.



