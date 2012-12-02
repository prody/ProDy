.. _msa-analysis:

*******************************************************************************
Conservation and Co-evolution Analysis
*******************************************************************************

Synopsis
===============================================================================

This example follows from :ref:`msafiles`. The aim of this part
is to show how to:

  * entropy calculations for conservation analsis
  * mutual information for co-evolution analsis
  * refine mutual info for comparison with dynamics
  * calculate entropy and mutual information
  * write and plot the calculated data

Get MSA data
===============================================================================

First, we parse an MSA file. See also :ref:`msafiles`:

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   >>> from prody import *
   >>> from matplotlib import pyplot as plt
   >>> searchPfam('1K2A').keys()
   ['PF00074']
   >>> msa = parseMSA(fetchPfamMSA('PF00074'))
   
 
Refine MSA
===============================================================================

Here, we refine MSA such that columns in the MSA that have gaps for a given 
sequence can be eliminated. We want to refine the MSA such that the sequence 
corresponding to the PDB has no gaps. Then we refine the msa so that sequences 
in the refined msa that have more than 20% gaps can be eliminated. We use the 
:func:`.refineMSA` to do this. It returns an :class:`.MSA` object.  

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   >>> msa_refine = msa[:,'RNAS2_HUMAN']
   >>> msa_refine
   <MSA: PF00074_full' (525 sequences, 128 residues)>
   >>> msa_refine = refineMSA(msa_refine, row_occ=0.8)
   >>> msa_refine
   <MSA: PF00074_full' refined (row_occ>=0.8) (466 sequences, 128 residues)>


Plotting occupancy
===============================================================================   

We can plot the ocuupancy for each column to see if there are any positions in 
the MSA that have a lot of gaps. We use the function :func:`.showMSAOccupancy` 
that uses :func:`.calcMSAOccupancy` to calculate occupancy for MSA. 

.. plot::
   :context:
   :include-source:
   
   >>> showMSAOccupancy(msa_refine, occ='res') # doctest: +SKIP
   
We can also specify indices based on the PDB.

.. plot::
   :context:
   :include-source:
	
	>>> indices = list(range(4,132))
	>>> showMSAOccupancy(msa_refine, occ='res', indices=indices) # doctest: +SKIP

We can further refine the MSA to remove positions that have low occupancy, but 
that will change the start-end positions of the labels in the MSA that is not 
corrected automatically on refinement. We can also plot occupancy based on rows 
for the seqeunces in the MSA.


Calculating and Plotting Entropy
===============================================================================

Here, we show how to calculate Shannon Entropy and plot entropy. Entropy for 
each position in the MSA is calculated using :func:`.calcShannonEntropy`. It 
takes :class:`.MSA` object or a numpy 2D array containg MSA as input. Returns
a 1D numpy arrauy. Plotting is done using :func:`.showShannonEntropy`. 

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   >>> entropy = calcShannonEntropy(msa_refine)

*entropy* is a 1D numpy array. 

.. plot::
   :context:
   :include-source:

   >>> showShannonEntropy(entropy, indices) # doctest: +SKIP
  
 
Calculating and Plotting Mutual Information
===============================================================================

Here, we show how to calculate mutual information between the positions of the 
MSA using :func:`.buildMutinfoMatrix` which also takes  :class:`.MSA` object 
or a numpy 2D array containg MSA as input. We can also apply normalization 
using :func:`.applyMutinfoNorm` and correction using :func:`.applyMutinfoCorr` 
to the mutual information matrix based on references [MLC05]_ and [DSD08]_ 
respectively. Returns a numpy 2D array.

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   >>> mutinfo = buildMutinfoMatrix(msa_refine)
   >>> mutinfo_norm = applyMutinfoNorm(mutinfo, entropy, norm='minent')
   >>> mutinfo_corr = applyMutinfoCorr(mutinfo, corr='apc')

Note that by default ``norm="sument"`` normalization is applied in 
``applyMutinfoNorm`` and ``corr="prod"`` is applied in ``applyMutinfoCorr``. 

Now we plot the mutual information matrices that we obtained above and see
the effects of different corrections and normalizations. 

.. plot::
   :context:
   :include-source:

   >>> showMutinfoMatrix(mutinfo) # doctest: +SKIP

.. plot::
   :context:
   :include-source:
   
   >>> showMutinfoMatrix(mutinfo_corr, clim=[-1,1], xlabel='1KA2: 4-131') # doctest: +SKIP
   
Note ylabel does not need to be set, since xlabel = ylabel
   
   
Writing Mutual Information and Entropy
===============================================================================

Here we show how to write the mutual information and entropy array. We use the
:func:`.writeArray` to write numpy array data. 

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   >>> writeArray('1KA2_MI.txt', mutinfo)
   '1KA2_MI.txt'

This can be later loaded using :func:`numpy.loadtxt`. Further analysis can also
be done by rank ordering the matrix and analyzing the pairs with highest mutual
information or the most co-evolving residues. This is done using  
:func:`.calcRankorder`. A zscore normalization can also be applied to select 
coevolving pairs based on a zscore cutoff.

.. plot::
   :context:
   :nofigs:
   :include-source:
   
   >>> import numpy
   >>> rank_row, rank_col, zscore_sort = calcRankorder(mutinfo, zscore=True)
   >>> print(numpy.asarray(indices)[rank_row[:5]])
   [128 129 130 130 130]
   >>> print(numpy.asarray(indices)[rank_col[:5]])
   [127 127 127 129 128]
   >>> print(zscore_sort[:5])
   [ 4.73041929  4.32016678  4.1165174   3.62089428  3.10104779]
   
   
See Also
===============================================================================

See :mod:`.sequence` module for all sequence analysis functions.

|questions|

|suggestions|

.. sectionauthor:: Anindita Dutta
