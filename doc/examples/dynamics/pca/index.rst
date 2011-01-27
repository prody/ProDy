.. currentmodule:: prody.dynamics

.. _pca:

*******************************************************************************
Principal Component Analysis
*******************************************************************************

ProDy can be used to perform principal component analysis of:

  * X-ray structures
  * NMR models
  * Mixed datasets from Blast 
  * MD snapshots 
  
X-ray structures
===============================================================================

X-ray structures of the same protein usually do not contain uniform data. 
That is termini or loop residues are missing non-uniformly among X-ray structures. 
This situation hardens construction of a covariance matrix for PC analysis. 

ProDy implements flexible classes that enable a forgiving analysis of an X-ray 
structure dataset. That is, despite missing residues in some structures most of 
a protein can be analyzed. An example is provided for p38 MAP kinase.
p38 MAP kinase structures were studied in [AB09]_ and [AB11]_. The following
example will reproduce the results presented in [AB09]_.

.. toctree::
   :maxdepth: 2

   pca_xray/index

NMR models
===============================================================================

Analysis of NMR models is fairly easier, as it requires parsing a single 
structure. An example is provided for Ubiquitin [AB09]_.

.. toctree::
   :maxdepth: 2

   pca_nmr

Mixed datasets
===============================================================================

Mixed structural datasets can be analyzed as well. Below example shows
how to retrieve and analyze cytochrome c structures by a Blast searching. 

.. toctree::
   :maxdepth: 2
   
   pca_blast


Protein Dimer
===============================================================================

Above examples deal with protein monomers. It is also possible to analyze
dimer structures.   

.. toctree::
   :maxdepth: 2
   
   pca_dimer


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

