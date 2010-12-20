.. currentmodule:: prody.dynamics

.. _anm:

*******************************************************************************
Anisotropic Network Model analysis
*******************************************************************************

This example shows how to perform ANM calculations.  

Prepare protein
===============================================================================

.. plot::
   :context:
   :nofigs:
   :include-source:

   from prody import *

   # parse pdb file by identifier
   # if file is not found in current working directory, it will be downloaded
   p38 = parsePDB('1p38')
   print p38 # print p38 to see number of atoms and number of coordinate sets

   # we only want to use alpha carbons, so we select them
   calphas = p38.select('protein and name CA')
   # we can also make the same selection like this
   calphas = p38.select('calpha')
   print calphas # this will show the number of atoms in the selection 

Note that, atoms other than alpha carbons can be selected and used in ANM 
calculations

Perform calculations
===============================================================================

**Build Hessian**

.. plot::
   :context:
   :nofigs:
   :include-source: 

   # instantiate ANM instance
   anm = ANM('p38 ANM analysis')

   # build Hessian matrix using selected atoms (351 alpha carbons)
   anm.buildHessian(calphas)

   # get a copy of the hessian matrix
   anm.getHessian()

**How to change cutoff and gamma parameters?**

We didn't pass any parameters, but ``buildHessian`` method accepts two of them
by default ``cutoff=15.0`` and ``gamma=1.0`` is passed
that is, ``buildHessian(calphas, cutoff=15., gamma=1.)``.
 
.. plot::
   :context:
   :nofigs:
   :include-source:
   
   anm.getCutoff()
   # this prints 15.0
   anm.getGamma()
   # this prints 1.0

Note that it is also possible to use an externally calculated Hessian 
matrix. Just pass it to the ANM instance using :meth:`~ANM.setHessian` method.

**Calculate normal modes**

.. plot::
   :context:
   :nofigs:
   :include-source:

   # calculate modes (by default slowest 20 will be calculated)
   anm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate different number
of non-zero modes or to keep zero modes, try: ``anm.calcModes(50, zeros=True)``

Access calculated data
===============================================================================

::

   # get covariance matrix
   # note that covariance calculated using 20 modes when you call this method
   anm.getCovariance()

   anm.getEigenvalues()
   anm.getEigenvectors()

Covariance matrices are calculated using available modes. If user calculates
M slowest modes, only they will be used in the calculation of covariance.

Individual modes
===============================================================================

Normal mode indices in Python start from 0, so slowest mode has index 0. By default,
modes with zero eigenvalues are excluded. If they were retained, slowest 
non-trivial mode would have index 6.

::

   # investigate individual modes
   slowest_mode = anm[0]

   slowest_mode.getEigenvalue()
   slowest_mode.getEigenvector()

Write NMD file
===============================================================================

ANM results in NMD format can be visualized using |nmwiz|. NMWiz is
a |vmd| plugin.

::

   # write slowest 3 ANM modes into an NMD file
   writeNMD('p38_anm_modes.nmd', anm[:3], calphas)
   # Note that slicing an ANM (or GNM, EDA) instances returns a list of modes.
   # In this case, slowest 3 ANM modes were written into NMD file.

View modes in VMD
===============================================================================

::

   # first make sure that the VMD path is correct
   print getVMDpath()

   # if this is incorrect use setVMDpath to correct it
   viewNMDinVMD('p38_anm_modes.nmd')

This will show the slowest 3 modes in VMD using NMWiz. This concludes ANM
example. Many of these apply to other NMA models, such as GNM and EDA instances.

.. plot::
   :context:
   :nofigs:

   import os
   os.remove('1p38.pdb.gz')

