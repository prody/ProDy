Analyze Conformations
===============================================================================

First, necessary imports:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()
   import os, glob

Parse conformations
-------------------------------------------------------------------------------

Now, let's read initial and refined conformations:

.. ipython:: python

   initial = AtomGroup('p38 initial')
   refined = AtomGroup('p38 refined')
   for pdb in glob.glob('p38_ensemble/*pdb'):
       fn = os.path.splitext(os.path.split(pdb)[1])[0]
       opt = os.path.join('p38_optimize', fn + '.coor')
       parsePDB(pdb, ag=initial)
       parsePDB(opt, ag=refined)
   initial
   refined


Calculate RMSD change
-------------------------------------------------------------------------------

We can plot RMSD change after refinement as follows:

.. ipython:: python

   rmsd_ca = []
   rmsd_all = []

   initial_ca = initial.ca
   refined_ca = refined.ca
   for i in range(initial.numCoordsets()):
       initial.setACSIndex(i)
       refined.setACSIndex(i)
       initial_ca.setACSIndex(i)
       refined_ca.setACSIndex(i)
       rmsd_ca.append(calcRMSD(initial_ca, refined_ca))
       rmsd_all.append(calcRMSD(initial, refined))

   plot(rmsd_all, label='all');
   plot(rmsd_ca, label='ca');
   xlabel('Conformation index');
   ylabel('RMSD');
   @savefig conformational_sampling_comparison.png width=4in
   legend();


Select a diverse set
-------------------------------------------------------------------------------

To select a diverse set of refined conformations, let's calculate average RMSD
for each conformation to others:

.. ipython:: python

   rmsd_mean = []
   for i in range(refined.numCoordsets()):
       refined.setACSIndex(i)
       alignCoordsets(refined)
       rmsd = calcRMSD(refined)
       rmsd_mean.append(rmsd.sum() / (len(rmsd) - 1))

   bar(arange(1, len(rmsd_mean) + 1), rmsd_mean);
   xlabel('Conformation index');
   @savefig conformational_sampling_mean_rmsd.png width=4in
   ylabel('Mean RMSD');

Let's select conformations that are 1.2 Å away from other on average:

.. ipython:: python

   rmsd_mean = array(rmsd_mean)
   selected = (rmsd_mean >= 1.2).nonzero()[0] + 1
   selected
   len(selected)


Visualization
-------------------------------------------------------------------------------

When you visualize the refined ensemble, you should see something similar to
this:

.. image:: _static/p38_sampling.png
   :width: 3in