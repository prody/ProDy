.. _gamma:

Custom Gamma Functions
===============================================================================

This example shows how to develop custom force constant functions for
:class:`.ANM` (or :class:`.GNM`) calculations.

We will use the relation shown in the figure below. For Cα atoms that are
10 to 15 Å apart from each other, we use a unit force constant. For those
that are 4 to 10 Å apart, we use a 2 times stronger force constant.
For those that are within 4 Å of each other (i.e. those from connected
residue pairs), we use a 10 times stronger force constant.

We will obtain an :class:`.ANM` instance that stores Hessian and Kirchhoff
matrices and normal mode data describing the intrinsic dynamics of the protein
structure. :class:`.ANM` instances and individual normal modes
(:class:`.Mode`) can be used as input to functions in :mod:`prody.dynamics`
module.


Parse structure
-------------------------------------------------------------------------------

We start by importing everything from ProDy, Numpy, and Matplotlib packages:

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion() # turn interactive mode on


We start with parsing a PDB file by passing an identifier.

.. ipython:: python

   p38 = parsePDB('1p38')
   p38


We want to use only Cα atoms, so we select them:

.. ipython:: python

   calphas = p38.select('protein and name CA')
   calphas

Force Constant Function
-------------------------------------------------------------------------------

We define the aformentioned function as follows:

.. ipython:: python

   def gammaDistanceDependent(dist2, *args):
       """Return a force constant based on the given square distance."""
       if dist2 <= 16:
           return 10
       elif dist2 <= 100:
           return 2
       elif dist2 <= 225:
           return 1
       else:
           return 0

Note that the input to this function from :class:`.ANM` or :class:`.GNM`
is the square of the distance. In addition, node (atom or residue) indices
are passed to this function, that's why we used ``*args`` in the function
definition.

Let's test how it works:


.. ipython:: python

   dist = arange(0, 20, 0.1)
   gamma = map(gammaDistanceDependent, dist ** 2)
   plot(dist, gamma, lw=4);
   axis([0, 20, 0, 12]);
   xlabel('Distance (A)');
   ylabel('Force constant');
   @savefig enm_analysis_gamma.png width=4in
   grid();


ANM calculations
-------------------------------------------------------------------------------

We use selected atoms (351 Cα's) and ``gammaDistanceDependent`` function
for ANM calculations as follows:

.. ipython:: python

   anm = ANM('1p38')
   anm.buildHessian(calphas, cutoff=15, gamma=gammaDistanceDependent)
   anm.calcModes()


For more detailed examples see :ref:`anm` or :ref:`gnm`.
