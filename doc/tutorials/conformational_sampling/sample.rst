Sample Conformations
===============================================================================

In this part, we will sample conformations along ANM modes.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Load results
-------------------------------------------------------------------------------

First, we load results produced in the previous part. If you are in the
same Python session, you don't need to do this.

.. ipython:: python

   p38 = loadAtoms('p38.ag.npz')
   p38_anm = loadModel('p38_ca.anm.npz')
   p38_anm_ext = loadModel('p38_ext.nma.npz')


Sampling
-------------------------------------------------------------------------------

We will use :func:`.sampleModes` function:

.. ipython:: python

   ens = sampleModes(p38_anm_ext, atoms=p38.protein, n_confs=40, rmsd=1.0)
   ens

We passed extended model, p38 structure, and two other parameters.
This will produce 40 (``n_confs``) conformations.  The conformations
will have an average 1.0 Å RMSD from the input structure.

We can write this ensemble in :file:`.dcd` for visualization in VMD:

.. ipython:: python

   writeDCD('p38all.dcd', ens)


Analysis
-------------------------------------------------------------------------------

Let's analyze the :class:`.Ensemble` by plotting RMSD of conformations
to the input structure:

.. ipython:: python

   rmsd = ens.getRMSDs()
   hist(rmsd, normed=False);
   @savefig conformational_sampling_ens_rmsd.png width=4in
   xlabel('RMSD');

This histogram might look like a flat distribution  due to the small size
of the ensemble. For larger numbers of conformations it will get closer to
a normal distribution. Let's calculate average and extremum RMSD values:

.. ipython:: python

   rmsd.mean()
   rmsd.max()
   rmsd.min()

Let's see the projection of these conformations in the ANM slow mode space:


.. ipython:: python

   @savefig conformational_ensemble_sampling_projection.png width=4in
   showProjection(ens, p38_anm_ext[:3], rmsd=True);
   proj = calcProjection(ens, p38_anm_ext[:3])


Write conformations
-------------------------------------------------------------------------------

We will write them in :file:`p38_ensemble` folder:

.. ipython::

   In [1]: mkdir -p p38_ensemble

Let's add the conformations to the :class:`.AtomGroup` object and set
:term:`beta` values of Cα atoms to 1 and of other atoms to 0:

.. ipython:: python

   p38.addCoordset(ens.getCoordsets())
   p38
   p38.all.setBetas(0)
   p38.ca.setBetas(1)

In the next step, we will place a harmonic constraint on atoms with beta
values 1. The optimization is aims for refining covalent geometry of atoms.
We do not want the new Cα to change much to keep the refined ensemble
diverse. We can easily verify that only Cα atoms have beta values set to 1:

.. ipython:: python

   p38.ca == p38.beta_1


Now we write these conformations out:

.. ipython:: python

   import os
   for i in range(1, p38.numCoordsets()):  # skipping 0th coordinate set
       fn = os.path.join('p38_ensemble', 'p38_' + str(i) + '.pdb')
       writePDB(fn, p38, csets=i)


Visualization
-------------------------------------------------------------------------------

You can visualize all of these conformations using VMD as follows::

  $ vmd -m p38_ensemble/*pdb