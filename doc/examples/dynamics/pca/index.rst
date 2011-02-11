.. currentmodule:: prody.dynamics

.. _pca:

*******************************************************************************
PCA examples
*******************************************************************************

ProDy can be used to perform principal component analysis for:

  * X-ray structures
  * NMR models
  * Mixed datasets
  * Multimeric proteins
  * MD snapshots 
  
X-ray structures
===============================================================================

X-ray structures of the same protein usually do not contain uniform data. 
For example, termini or loop residues are can be unevenly represented in the 
X-ray structures, which may make it difficult to construct the covariance 
matrix for PC analysis. 

To address the above, we implemented flexible classes in ProDy that facilitate 
the analysis of heterogeneous structural datasets. An example is provided for 
p38 MAP kinase. The p38 MAP kinase structures were studied in [AB09]_ and 
[AB11]_. The following example reproduces the results presented in [AB09]_.

.. toctree::
   :maxdepth: 2

   pca_xray
   
NMR models
===============================================================================

The analysis of NMR models is fairly easy with ProDy and requires parsing a 
single structure file. An example is provided for Ubiquitin [AB09]_.

.. toctree::
   :maxdepth: 1

   pca_nmr

Mixed datasets
===============================================================================

Mixed structural datasets can also be analyzed. The example below shows
how to retrieve and analyze cytochrome c structures by Blast searching. 

.. toctree::
   :maxdepth: 1
   
   pca_blast


Multimeric proteins
===============================================================================

The previous examples have dealt with monomeric proteins. It is also possible 
to analyze multimeric structures.   

.. toctree::
   :maxdepth: 1
   
   pca_dimer


MD snapshots
===============================================================================

Snapshots from an MD trajectory can also be analyzed by replacing the PDB file 
in the above examples with a file that contains snaphots from a simulation. 
However, when analyzing a large number of frames is to be analyzed, the 
following example should be used:

.. toctree::
   :maxdepth: 1

   eda_md


Functions and Plotting
===============================================================================

ProDy comes with several functions to compare PCA/ENM models or to plot data
from such models. These functions can be found in :ref:`dynamics`.

