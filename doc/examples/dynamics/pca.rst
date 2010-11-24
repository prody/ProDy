.. currentmodule:: prody.dynamics

.. _howtopca:

*******************************************************************************
Principal Component Analysis
*******************************************************************************

ProDy can be used to perform principal component analysis of:

  * X-ray structures
  * NMR models
  * MD snapshots 
  
X-ray structures
===============================================================================

X-ray structures of the same protein are typically not uniform. That is termini 
or loop residues may be missing in some structures. This situation hardens 
construction of a covariance matrix for PC analysis. 

ProDy implements flexible classes to enable a forgiving analysis of an X-ray 
ensemble. That is, despite missing residues in some structures most of a protein 
can be analyzed. An example is provided for p38 MAP kinase [AB09]_.

.. toctree::
   :maxdepth: 2

   p38_Xray/p38_Xray_I
   p38_Xray/p38_Xray_II
   p38_Xray/p38_Xray_III
   p38_Xray/p38_Xray_IV


NMR models
===============================================================================

Analysis of NMR models is fairly easier, as it requires parsing a single 
structure. An example is provided for Ubiquitin [AB09]_.

.. toctree::
   :maxdepth: 2

   pca_nmr_anm_ubi


MD snapshots
===============================================================================

Snapshots from an MD trajectory can also be analyzed in the
same way. Replacing the PDB file in the above example with a file that 
contains snaphots from a simulation will do the trick. However, if a large
number of frames is to be analyzed, the following example should be followed.

.. toctree::
   :maxdepth: 2

   eda_md


Functions and Plotting
===============================================================================

ProDy comes with several functions to compare PCA/ENM models or to plot data
from such models. These functions are listed in :ref:`dynamics`.

