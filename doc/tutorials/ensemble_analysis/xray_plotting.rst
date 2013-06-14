.. _pca-xray-plotting:


Plotting
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example is continued from :ref:`pca-xray-analysis`. The aim of this part
is to produce graphical comparison of experimental and theoretical data.
We will reproduce the plots that was presented in our paper [AB09]_.

Load data
-------------------------------------------------------------------------------


First, we load data saved in :ref:`pca-xray-calculations`:


.. ipython:: python

   from prody import *
   from pylab import *
   ion()

.. ipython:: python

   pca = loadModel('p38_xray.pca.npz')
   anm = loadModel('1p38.anm.npz')
   ensemble = loadEnsemble('p38_X-ray.ens.npz')
   ref_chain = parsePDB('p38_ref_chain.pdb')


PCA - ANM overlap
-------------------------------------------------------------------------------

In previous page, we compared PCA and ANM modes to get some numbers. In this
case, we will use plotting functions to make similar comparisons:

.. ipython:: python

   showOverlapTable(pca[:6], anm[:6]);
   # Let's change the title of the figure
   @savefig ensemble_analysis_xray_overlap_table.png width=4in
   title('PCA - ANM Overlap Table');


It is also possible to plot overlap of a single mode from one model with
multiple modes from another:

.. ipython:: python

   @savefig ensemble_analysis_xray_overlap.png width=4in
   showOverlap(pca[0], anm);
   @savefig ensemble_analysis_xray_cumulative_overlap.png width=4in
   showCumulOverlap(pca[0], anm);

Let's also plot the cumulative overlap in the same figure:


Square fluctuations
-------------------------------------------------------------------------------

.. ipython:: python

   @savefig ensemble_analysis_xray_pca_sqflucts.png width=4in
   showSqFlucts(pca[:3]);
   @savefig ensemble_analysis_xray_anm_sqflucts.png width=4in
   showSqFlucts(anm[:3]);


Now let's plot square fluctuations along PCA and ANM modes in the same plot:

.. ipython:: python

   showScaledSqFlucts(pca[0], anm[2]);
   @savefig ensemble_analysis_xray_pca_anm_sqflucts_1.png width=4in
   legend();


.. ipython:: python

   showScaledSqFlucts(pca[1], anm[0]);
   @savefig ensemble_analysis_xray_pca_anm_sqflucts_2.png width=4in
   legend();


In above example, ANM modes are scaled to have the same mean as PCA modes.
Alternatively, we could plot normalized square fluctuations:

.. ipython:: python

   showNormedSqFlucts(pca[0], anm[1]);
   @savefig ensemble_analysis_xray_pca_anm_sqflucts_3.png width=4in
   legend();



Projections
-------------------------------------------------------------------------------

Now we will project the ensemble onto PC 1 and 2 using
:func:`.showProjection`:

.. ipython:: python

   showProjection(ensemble, pca[:2]);
   @savefig ensemble_analysis_xray_pca_projection.png width=4in
   axis([-0.8, 0.8, -0.8, 0.8]);


Now we will do a little more work, and get a colorful picture:

======  =====================
red     unbound
blue    inhibitor bound
yellow  glucoside bound
purple  peptide/protein bound
======  =====================


.. ipython:: python

   color_list = ['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue',
                 'blue', 'purple', 'purple', 'blue', 'blue', 'blue',
                 'blue', 'blue', 'red', 'red', 'red', 'blue', 'blue',
                 'blue', 'blue', 'blue','blue', 'blue', 'blue', 'blue',
                 'blue', 'red', 'blue', 'blue','blue', 'blue', 'blue',
                 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'yellow',
                 'yellow', 'yellow', 'yellow', 'blue', 'blue','blue',
                 'blue', 'blue', 'blue', 'yellow', 'purple', 'purple',
                 'blue', 'yellow', 'yellow', 'yellow', 'blue', 'yellow',
                 'yellow', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue',
                 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue',
                 'blue', 'purple']
   color2label = {'red': 'Unbound', 'blue': 'Inhibitor bound',
                  'yellow': 'Glucoside bound',
                  'purple': 'Peptide/protein bound'}
   label_list = [color2label[color] for color in color_list]
   showProjection(ensemble, pca[:2], color=color_list,
                  label=label_list);
   axis([-0.8, 0.8, -0.8, 0.8]);
   @savefig ensemble_analysis_xray_pca_projection_2.png width=4in
   legend();


Now let's project conformations onto 3d principal space and label conformations
using ``text`` keyword argument and :meth:`.PDBEnsemble.getLabels` method:

.. ipython:: python

   @savefig ensemble_analysis_xray_pca_projection_3.png width=4in
   showProjection(ensemble, pca[:3], color=color_list, label=label_list,
                  text=ensemble.getLabels(), fontsize=10);

The figure with all conformation labels is crowded, but in an interactive
session you can zoom in and out to make text readable.


Cross-projections
-------------------------------------------------------------------------------

Finally, we will make a cross-projection plot using
:func:`.showCrossProjection`. We will pass ``scale='y'`` argument, which will
scale the width of the projection along ANM mode:


.. ipython:: python

   showCrossProjection(ensemble, pca[0], anm[2], scale="y",
                        color=color_list, label=label_list);
   plot([-0.8, 0.8], [-0.8, 0.8], 'k');
   axis([-0.8, 0.8, -0.8, 0.8]);
   @savefig ensemble_analysis_xray_cross_projection_1.png width=4in
   legend(loc='upper left');


.. ipython:: python

   showCrossProjection(ensemble, pca[1], anm[0], scale="y",
                       color=color_list, label=label_list);
   plot([-0.8, 0.8], [-0.8, 0.8], 'k');
   @savefig ensemble_analysis_xray_cross_projection_2.png width=4in
   axis([-0.8, 0.8, -0.8, 0.8]);

It is also possible to find the correlation between these projections:

.. ipython:: python

   pca_coords, anm_coords = calcCrossProjection(ensemble, pca[0], anm[2])
   print(np.corrcoef(pca_coords, anm_coords))


This is going to print 0.95 for PC 1 and ANM mode 2 pair.


Finally, it is also possible to label conformations in cross projection plots
too:

.. ipython:: python

   showCrossProjection(ensemble, pca[1], anm[0], scale="y",
       color=color_list, label=label_list, text=ensemble.getLabels(),
       fontsize=10);
   plot([-0.8, 0.8], [-0.8, 0.8], 'k');
   @savefig ensemble_analysis_xray_cross_projection_2.png width=4in
   axis([-0.8, 0.8, -0.8, 0.8]);
