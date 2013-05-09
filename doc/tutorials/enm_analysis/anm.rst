.. _anm:

ANM analysis
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example shows how to perform ANM calculations, and retrieve
normal mode data.  An :class:`~.ANM` instance that stores Hessian and Kirchhoff
matrices and normal mode data describing intrinsic dynamics of the protein
structure will be obtained.  :class:`~.ANM` instances and individual normal
modes (:class:`~.Mode`) can be used as input to functions in
:mod:`~prody.dynamics` module.


Parse structure
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *

.. ipython:: python
   :suppress:

   from prody import *; from pylab import *; ion()

We start with parsing a PDB file by passing an identifier.
Note that if a file is not found in the current working directory, it will be
downloaded.

.. ipython:: python

   p38 = parsePDB('1p38')
   repr(p38)

We want to use only Cα atoms, so we select them:

.. ipython:: python

   calphas = p38.select('protein and name CA')
   repr(calphas)

We can also make the same selection like this:

.. ipython:: python

   calphas2 = p38.select('calpha')
   repr(calphas2)


To check whether the selections are the same, we can try:

.. ipython:: python

   calphas == calphas2


Note that, ProDy atom selector gives the flexibility to select any set of atoms
to be used in ANM  calculations.

Build Hessian
-------------------------------------------------------------------------------

We instantiate an :class:`~.ANM` instance:

.. ipython:: python

   anm = ANM('p38 ANM analysis')

Then, build the Hessian matrix by passing selected atoms (351 Cα's)
to :meth:`.ANM.buildHessian` method:

.. ipython:: python

   anm.buildHessian(calphas)

We can get a copy of the Hessian matrix using :meth:`.ANM.getHessian` method:

.. ipython:: python

   anm.getHessian().round(3)


Parameters
-------------------------------------------------------------------------------

We didn't pass any parameters to :meth:`.ANM.buildHessian` method, but it
accepts *cutoff* and *gamma* parameters, for which  default values are
``cutoff=15.0`` and ``gamma=1.0``.

.. ipython:: python

   anm.getCutoff()
   anm.getGamma()


Note that it is also possible to use an externally calculated Hessian
matrix. Just pass it to the ANM instance using :meth:`.ANM.setHessian` method.

Calculate normal modes
-------------------------------------------------------------------------------

Calculate modes using :meth:`.ANM.calcModes` method:

.. ipython:: python

   anm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate different number
of non-zero modes or to keep zero modes, try: ``anm.calcModes(50, zeros=True)``

Normal mode data
-------------------------------------------------------------------------------

.. ipython:: python

   anm.getEigvals().round(3)
   anm.getEigvecs().round(3)


You can get the covariance matrix as follows:

.. ipython:: python

   anm.getCovariance().round(2)

Covariance matrices are calculated using available modes (slowest 20 modes in
this case). If user calculates M slowest modes, only they will be used in the
calculation of covariance.

Individual modes
-------------------------------------------------------------------------------

Normal mode indices in Python start from 0, so slowest mode has index 0.
By default, modes with zero eigenvalues are excluded. If they were retained,
slowest non-trivial mode would have index 6.

Get the slowest mode by indexing :class:`~.ANM` instance as follows:

.. ipython:: python

   slowest_mode = anm[0]
   slowest_mode.getEigval().round(3)
   slowest_mode.getEigvec().round(3)


Write NMD file
-------------------------------------------------------------------------------

ANM results in NMD format can be visualized using NMWiz VMD plugin |nmwiz|.
Following statement writes slowest 3 ANM modes into an NMD file:

.. ipython:: python

   writeNMD('p38_anm_modes.nmd', anm[:3], calphas)


Note that slicing an ANM (or GNM, EDA) instances returns a list of modes.
In this case, slowest 3 ANM modes were written into NMD file.

View modes in VMD
-------------------------------------------------------------------------------

First make sure that the VMD path is correct

.. ipython:: python

   getVMDpath()


::

   # if this is incorrect use setVMDpath to correct it
   viewNMDinVMD('p38_anm_modes.nmd')

This will show the slowest 3 modes in VMD using NMWiz. This concludes ANM
example. Many of these apply to other NMA models, such as GNM and EDA instances.

