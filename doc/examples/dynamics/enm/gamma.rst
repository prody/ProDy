.. _gamma:

*******************************************************************************
Developing Gamma Functions
*******************************************************************************

Synopsis
===============================================================================

This example shows how to develop custom force constant functions for
:class:`~.ANM` (or :class:`~.GNM`) calculations. 

We will use the relation shown in the figure below. For Cα atoms that are
10 to 15 Å apart from each other, we use a unit force constant. For those
that are 4 to 10 Å apart, we use a 2 times stronger force constant. 
For those that are within 4 Å of each other (i.e. those from connected 
residue pairs), we use a 10 times stronger force constant.  

.. plot::
   :context:
   
   import matplotlib.pyplot as plt
   plt.figure(figsize=(5,4))
   plt.plot([0,4,4,10,10,15,15,20], [10,10,2,2,1,1,0,0], lw=4)
   plt.axis([0, 20, 0, 12])
   plt.xlabel('Distance')
   plt.ylabel('Force constant')

.. plot::
   :context:
   :nofigs:

   plt.close('all')  


Input
-------------------------------------------------------------------------------

Protein structure data in PDB file format or specified by a PDB identifier.

Output
-------------------------------------------------------------------------------

An :class:`~.ANM` instance that stores Hessian and Kirchhoff matrices and 
normal mode data describing intrinsic dynamics of the protein structure. 
:class:`~.ANM` instances and individual normal modes 
(:class:`~.Mode`) can be used as input to functions in :mod:`~prody.dynamics` 
module.


ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Prepare protein
-------------------------------------------------------------------------------

We start with parsing a PDB file by passing an identifier.

>>> p38 = parsePDB('1p38')
>>> p38
<AtomGroup: 1p38 (2962 atoms)>

We want to use only Cα atoms, so we select them:

>>> calphas = p38.select('protein and name CA')
>>> calphas
<Selection: "protein and name CA" from 1p38 (351 atoms)>

Define Force Constant Function
-------------------------------------------------------------------------------

We define the aformentioned function as follows:

>>> def gammaDistanceDependent(dist2, *args):
...     """Return a force constant based on the given square distance."""
...
...     if dist2 <= 16:
...         return 10 
...     elif dist2 <= 100:
...         return 2
...     elif dist2 <= 225:
...         return 1
...     else:
...         return 0

Note that the input to this function from :class:`~.ANM` or :class:`~.GNM` 
is the square of the distance. In addition, node (atom or residue) indices
are passed to this function, that's why we used ``*args`` in the function
definition.

Let's test how it works:

>>> gammaDistanceDependent(3.8**2)
10
>>> gammaDistanceDependent(10**2)
2
>>> gammaDistanceDependent(10.1**2)
1
>>> gammaDistanceDependent(25**2)
0

ANM calculations
-------------------------------------------------------------------------------

We use selected atoms (351 Cα's) and ``gammaDistanceDependent`` function
for ANM calculations as follows:

>>> anm = ANM('1p38')
>>> anm.buildHessian(calphas, cutoff=15, gamma=gammaDistanceDependent)
>>> anm.calcModes()


For more detailed examples see :ref:`anm` or :ref:`gnm`.

|questions|

|suggestions|
