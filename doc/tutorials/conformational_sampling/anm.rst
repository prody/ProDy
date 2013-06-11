ANM Calculations
===============================================================================

Required imports:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Atom Selection
-------------------------------------------------------------------------------

First, we parse :download:`p38 structure
<conformational_sampling_files/p38.pdb>`:

.. ipython:: python

   p38 = parsePDB('conformational_sampling_files/p38.pdb')
   p38

Let's take a look at the structure:

.. ipython:: python

   showProtein(p38);
   @savefig conformational_sampling_p38.png width=4in
   legend();

Note that this structure has hydrogen atoms which are added using PSFGEN that
comes with NAMD:

.. ipython:: python

   p38.numAtoms('hydrogen')


We will perform ANM calculations for 351 Cα atoms of the structure:


.. ipython:: python

   p38_ca = p38.ca
   p38_ca


ANM Calculation
-------------------------------------------------------------------------------

First, let's instantiate an :class:`.ANM` object:

.. ipython:: python

   p38_anm = ANM('p38 ca')
   p38_anm

Now, we can build Hessian matrix, simply by calling :meth:`.ANM.buildHessian`
method:


.. ipython:: python

   p38_anm.buildHessian(p38_ca)
   p38_anm

We see that :class:`.ANM` object contains 351 nodes, which correspond to the
Cα atoms.

We will calculate only top ranking three ANM modes, since we are going to
use only that many in sampling:

.. ipython:: python

   p38_anm.calcModes(n_modes=3)
   p38_anm


Analysis & Plotting
-------------------------------------------------------------------------------

Let's plot mobility of residues along ANM modes:

.. ipython:: python

   @savefig conformational_sampling_sqflucts.png width=4in
   showSqFlucts(p38_anm);

We can also calculate collectivity of these modes as follows:

.. ipython:: python

   for mode in p38_anm:
       print('{}\tcollectivity: {}'.format(str(mode), calcCollectivity(mode)))


Visualization
-------------------------------------------------------------------------------

You can visualize ANM modes using :ref:`nmwiz`. You need to write an
:file:`.nmd` file using :func:`writeNMD` and open it using VMD:

.. ipython:: python

   writeNMD('p38_anm.nmd', p38_anm, p38_ca)

For visualization, you can use :func:`viewNMDinVMD`, i.e.
``viewNMDinVMD('p38_anm.nmd')``

Extend Model
-------------------------------------------------------------------------------

We want to use ANM model to sample all atoms conformations of p38 MAPK, but
we have a coarse-grained model. We will use :func:`.extendModel` function
for this purpose:


.. ipython:: python

   p38_anm_ext, p38_all = extendModel(p38_anm, p38_ca, p38, norm=True)
   p38_anm_ext
   p38_all


Note ``p38_anm_ext`` is an :class:`.NMA` model, which has similar features as
an :class:`.ANM` object. Extended model has 3 modes, but 5668 atoms as opposed
to 351 nodes in the original :class:`.ANM` model.

Let's plot mobility of residues again to help understand what extending a
model does:

.. ipython:: python

   @savefig conformational_sampling_sqflucts_ext.png width=4in
   showSqFlucts(p38_anm_ext);

As you see, shape of the mobility plot is identical.  In the extended model,
each in the same direction as the Cα atoms of the residues that they belong to.
The mobility profile is scaled down, however, due to renormalization of
the mode vectors.

Save Results
-------------------------------------------------------------------------------

Now let's save the original and extended model, and atoms:

.. ipython:: python

   saveAtoms(p38)
   saveModel(p38_anm)
   saveModel(p38_anm_ext, 'p38_ext')

More Examples
-------------------------------------------------------------------------------

We have performed a quick ANM calculation and extended the resulting model
to all atoms of of the structure. You can see more examples on this
in :ref:`enm-analysis` tutorial.