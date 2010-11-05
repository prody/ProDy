.. module:: prody.dynamics

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
can be analyzed. An example script is provided for p38 MAP kinase [AB09]_.

.. seealso:: :ref:`egp38`


NMR models
===============================================================================

Analysis of NMR models is fairly easier, as it requires parsing a single 
structure. An example is provided for Ubiquitin [AB09]_.

.. seealso:: :ref:`egubi`



MD snapshots
===============================================================================

Snapshots from an MD trajectory can also be analyzed in the
same way. Replacing the PDB file in the above example with a file that 
contains snaphots from a simulation will do the trick. However, if a large
number of frames is to be analyzed, the following example should be followed.

.. seealso:: :ref:`egeda`

Functions and Plotting
===============================================================================

ProDy comes with several functions to compare PCA/ENM models or to plot data
from such models. Interested user should see the following pages:

  * :ref:`dyfunctions`
  * :ref:`dyplotting`
