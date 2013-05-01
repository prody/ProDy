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

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   >>> from prody import *
   >>> from matplotlib import pyplot as plt
   >>> pca = loadModel('p38_xray.pca.npz')
   >>> anm = loadModel('1p38.anm.npz')
   >>> ensemble = loadEnsemble('p38_X-ray.ens.npz')
   >>> ref_chain = parsePDB('p38_ref_chain.pdb')
   
 
PCA - ANM overlap  
-------------------------------------------------------------------------------

In previous page, we compared PCA and ANM modes to get some numbers. In this
case, we will use plotting functions to make similar comparisons:

.. plot::
   :context:
   :include-source:
   
   >>> showOverlapTable(pca[:6], anm[:6]) # doctest: +SKIP
   >>> # Let's change the title of the figure
   >>> plt.title('PCA - ANM Overlap Table') # doctest: +SKIP

   

It is also possible to plot overlap of a single mode from one model with
multiple modes from another:

.. plot::
   :context:
   :include-source:
   
   >>> showOverlap(pca[0], anm) # doctest: +SKIP
   >>> showCumulOverlap(pca[0], anm) # doctest: +SKIP

Let's also plot the cumulative overlap in the same figure:


Square fluctuations  
-------------------------------------------------------------------------------

.. plot::
   :context:
   :include-source:
   
   >>> showSqFlucts(pca[:3]) # doctest: +SKIP

.. plot::
   :context:
   :include-source:

   >>> showSqFlucts(anm[:3]) # doctest: +SKIP

   
Now let's plot square fluctuations along PCA and ANM modes in the same plot:

.. plot::
   :context:
   :include-source:
   
   >>> showScaledSqFlucts(pca[0], anm[2]) # doctest: +SKIP
   >>> plt.legend()


.. plot::
   :context:
   :include-source:

   >>> showScaledSqFlucts(pca[1], anm[0]) # doctest: +SKIP
   >>> plt.legend()


In above example, ANM modes are scaled to have the same mean as PCA modes. 
Alternatively, we could plot normalized square fluctuations:

.. plot::
   :context:
   :include-source:
   
   >>> showNormedSqFlucts(pca[0], anm[1]) # doctest: +SKIP
   >>> plt.legend()



Projections  
-------------------------------------------------------------------------------

Now we will project the ensemble onto PC 1 and 2 using 
:func:`.showProjection`:

.. plot::
   :context:
   :include-source:
   
   >>> showProjection(ensemble, pca[:2]) # doctest: +SKIP
   >>> plt.axis([-0.8, 0.8, -0.8, 0.8]) # doctest: +SKIP


Now we will do a little more work, and get a colorful picture:

======  =====================
red     unbound
blue    inhibitor bound
yellow  glucoside bound
purple  peptide/protein bound
======  =====================


.. plot::
   :context:
   :include-source:

   >>> color_list = ['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 
   ...               'blue', 'purple', 'purple', 'blue', 'blue', 'blue', 
   ...               'blue', 'blue', 'red', 'red', 'red', 'blue', 'blue',  
   ...               'blue', 'blue', 'blue','blue', 'blue', 'blue', 'blue', 
   ...               'blue', 'red', 'blue', 'blue','blue', 'blue', 'blue',  
   ...               'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'yellow', 
   ...               'yellow', 'yellow', 'yellow', 'blue', 'blue','blue', 
   ...               'blue', 'blue', 'blue', 'yellow', 'purple', 'purple', 
   ...               'blue', 'yellow', 'yellow', 'yellow', 'blue', 'yellow', 
   ...               'yellow', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue',
   ...               'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 
   ...               'blue', 'purple'] 
   >>> color2label = {'red': 'Unbound', 'blue': 'Inhibitor bound', 
   ...               'yellow': 'Glucoside bound', 
   ...               'purple': 'Peptide/protein bound'}
   >>> label_list = [color2label[color] for color in color_list]
   >>> showProjection(ensemble, pca[:2], color=color_list, 
   ...                label=label_list) # doctest: +SKIP
   >>> plt.axis([-0.8, 0.8, -0.8, 0.8]) # doctest: +SKIP
   >>> plt.legend()

   
Now let's project conformations onto 3d principal space and label conformations 
using ``text`` keyword argument and :meth:`.PDBEnsemble.getLabels` method:
 
.. plot::
   :context:
   :include-source:

   >>> showProjection(ensemble, pca[:3], color=color_list, label=label_list,  
   ...                text=ensemble.getLabels(), fontsize=10) # doctest: +SKIP

The figure with all conformation labels is crowded, but in an interactive 
session you can zoom in and out to make text readable.

   
Cross-projections
-------------------------------------------------------------------------------

Finally, we will make a cross-projection plot using 
:func:`.showCrossProjection`. We will pass ``scale='y'`` argument, which will 
scale the width of the projection along ANM mode:


.. plot::
   :context:
   :include-source:

   >>> showCrossProjection(ensemble, pca[0], anm[2], scale="y", 
   ...                     color=color_list, label=label_list) # doctest: +SKIP
   >>> plt.plot([-0.8, 0.8], [-0.8, 0.8], 'k') # doctest: +SKIP
   >>> plt.axis([-0.8, 0.8, -0.8, 0.8]) # doctest: +SKIP
   >>> plt.legend(loc='upper left') # doctest: +SKIP
   

.. plot::
   :context:
   :include-source:
   
   >>> showCrossProjection(ensemble, pca[1], anm[0], scale="y", 
   ...                     color=color_list, label=label_list) # doctest: +SKIP
   >>> plt.plot([-0.8, 0.8], [-0.8, 0.8], 'k') # doctest: +SKIP
   >>> plt.axis([-0.8, 0.8, -0.8, 0.8]) # doctest: +SKIP

It is also possible to find the correlation between these projections:

.. plot::
   :context:
   :include-source:
   
   >>> import numpy as np 
   >>> pca_coords, anm_coords = calcCrossProjection(ensemble, pca[0], anm[2])
   >>> print(np.corrcoef(pca_coords, anm_coords))
   [[ 1.         -0.94621454]
    [-0.94621454  1.        ]]
    
    
This is going to print 0.95 for PC 1 and ANM mode 2 pair.


Finally, it is also possible to label conformations in cross projection plots 
too:

.. plot::
   :context:
   :include-source:

   >>> showCrossProjection(ensemble, pca[1], anm[0], scale="y", 
   ... color=color_list, label=label_list, text=ensemble.getLabels(), 
   ... fontsize=10) # doctest: +SKIP
   >>> plt.plot([-0.8, 0.8], [-0.8, 0.8], 'k') # doctest: +SKIP
   >>> plt.axis([-0.8, 0.8, -0.8, 0.8]) # doctest: +SKIP
