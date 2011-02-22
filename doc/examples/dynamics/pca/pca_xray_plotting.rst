.. currentmodule:: prody.dynamics

.. _pca-xray-plotting:

*******************************************************************************
PCA of X-ray structures: Plotting
*******************************************************************************

Synopsis
===============================================================================

This example is continued from :ref:`pca-xray-analysis`. The aim of this part 
is to make the plots that was presented in our paper [AB09]_.

Input
-------------------------------------------------------------------------------

PCA and ANM data calculated in :ref:`pca-xray-calculations`.

Output
-------------------------------------------------------------------------------

Graphical quantitative output comparing experimental and theoretical data.

ProDy Code
===============================================================================


We start by importing everything from the ProDy package and the matplotlib 
library:

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   from prody import *
   import matplotlib.pyplot as plt


Then we load data saved in :ref:`pca-xray-calculations`:

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   pca = loadModel('p38_xray.pca.npz')
   anm = loadModel('1p38.anm.npz')
   ensemble = loadEnsemble('p38_X-ray.ens.npz')
   ref_chain = parsePDB('p38_ref_chain.pdb')
   
 
 
PCA - ANM overlap  
-------------------------------------------------------------------------------

In previous page, we compared PCA and ANM modes to get some numbers. In this
case, we will use plotting functions to make similar comparisons:

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4)) # This opens a new empty figure
   showOverlapTable(pca[:6], anm[:6])

   # Let's change the title of the figure
   plt.title('PCA - ANM Overlap Table')

.. plot::
   :context:
   :nofigs:

   plt.close('all')
   

It is also possible to plot overlap of a single mode from one model with
multiple modes from another:

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showOverlap(pca[0], anm)

Let's also plot the cumulative overlap in the same figure:

.. plot::
   :context:
   :include-source:
   
   # plt.figure(figsize=(5,4)) # Note that we don't want to call this function in this case
   showCumulativeOverlap(pca[0], anm)

.. plot::
   :context:
   :nofigs:

   plt.close('all')  

Square fluctuations  
-------------------------------------------------------------------------------

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showSqFlucts(pca[:3])


   plt.figure(figsize=(5,4))
   showSqFlucts(anm[:3])

.. plot::
   :context:
   :nofigs:

   plt.close('all')
   
Now let's plot square fluctuations along PCA and ANM modes in the same plot:

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showScaledSqFlucts(pca[0], anm[2])
   plt.legend()

   plt.figure(figsize=(5,4))
   showScaledSqFlucts(pca[1], anm[0])
   plt.legend()

.. plot::
   :context:
   :nofigs:

   plt.close('all')

In above example, ANM modes are scaled to have the same mean as PCA modes. 
Alternatively, we could plot normalized square fluctuations:

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showNormedSqFlucts(pca[0], anm[1])
   plt.legend()

.. plot::
   :context:
   :nofigs:

   plt.close('all')


Projections  
-------------------------------------------------------------------------------

Now we will project the ensemble onto PC 1 and 2 using :func:`showProjection`:

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showProjection(ensemble, pca[:2])

.. plot::
   :context:
   :nofigs:

   plt.close('all')

Now we will do a little more work, and get a colorful picture:

.. plot::
   :context:
   :include-source:

   # list of colors, 
   #   red for unbound
   #   blue for inhibitor bound
   #   yellow for glucoside bound
   #   purple for peptide/protein bound
   # the order of 
   color_list = ['purple', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 
                 'purple', 'purple', 'blue', 'blue', 'blue', 'blue', 'blue', 
                 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 
                 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'blue', 'blue', 
                 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 
                 'blue', 'yellow', 'yellow', 'yellow', 'yellow', 'blue', 'blue', 
                 'blue', 'blue', 'blue', 'blue', 'yellow', 'purple', 'purple', 
                 'blue', 'yellow', 'yellow', 'yellow', 'blue', 'yellow', 'yellow', 
                 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 
                 'blue', 'blue', 'blue', 'blue', 'blue', 'blue'] 
   import numpy as np
   color_array = np.array(color_list) # Having an array will be handier  
   color_assignments = [('Unbound', 'red'), ('Inhibitor bound', 'blue'), ('Glucoside bound', 'yellow'), ('Peptide/protein bound', 'purple')]
   
   plt.figure(figsize=(5,4))
   for lbl, clr in color_assignments:
       showProjection(ensemble[ color_array == clr], pca[:2], color=clr, label=lbl)
   
It is possible to show the legend for this plot, but the figure gets crowded:
   
.. plot::
   :context:
   :include-source:

   plt.legend()

.. plot::
   :context:
   :nofigs:

   plt.close('all')

Cross-projections
-------------------------------------------------------------------------------

Finally, we will make a cross-projection plot using :func:`showCrossProjection`.
We will pass ``scale='y'`` argument, which will scale the width of the projection
along ANM mode:


.. plot::
   :context:
   :include-source:

   plt.figure(figsize=(5,4))
   for lbl, clr in color_assignments:
       showCrossProjection(ensemble[color_array == clr], pca[0], anm[2], scale="y", scalar=-0.78, color=clr, label=lbl)
   plt.plot([-25, 25], [-25, 25], 'k')
   plt.axis([-25, 25, -25, 25])

   plt.figure(figsize=(5,4))
   for lbl, clr in color_assignments:
       showCrossProjection(ensemble[color_array == clr], pca[1], anm[0], scale="y", scalar=-0.94, color=clr, label=lbl)
   plt.plot([-15, 15], [-15, 15], 'k')
   plt.axis([-15, 15, -15, 15])

It is also possible to find the correlation between these projections:

.. plot::
   :context:
   :include-source:
   :nofigs:
   
   pca_coords = calcProjection(ensemble, pca[0])
   anm_coords = calcProjection(ensemble, anm[2])
   
   print np.corrcoef(pca_coords, anm_coords)
   
This is going to print 0.95 for PC 1 and ANM mode 2 pair.

.. plot::
   :context:
   :nofigs:

   plt.close('all')

See Also
===============================================================================

This example is continued in :ref:`pca-xray-visualization`.

|more| See a list of all plotting functions in :ref:`dynamics`.

|questions|

|suggestions|
