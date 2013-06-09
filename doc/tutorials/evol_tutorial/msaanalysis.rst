.. _msaanalysis:

Evolution Analysis
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example follows from :ref:`msafiles`. The aim of this part
is to show how to:

  * entropy calculations for conservation analsis
  * mutual information for co-evolution analsis
  * refine mutual info for comparison with dynamics
  * calculate entropy and mutual information
  * write and plot the calculated data

Get MSA data
-------------------------------------------------------------------------------

First, we import everything from the ProDy package.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()  # turn interactive mode on

Then, we parse an MSA file for protein family :pfam:`PF00074`.
See also :ref:`msafiles`:

.. ipython:: python

   searchPfam('1K2A').keys()
   msa = parseMSA(fetchPfamMSA('PF00074'))


Refine MSA
-------------------------------------------------------------------------------

Here, we refine MSA such that columns in the MSA that have gaps for a given
sequence can be eliminated. We want to refine the MSA such that the sequence
corresponding to the PDB has no gaps. Then we refine the msa so that sequences
in the refined msa that have more than 20% gaps can be eliminated. We use the
:func:`.refineMSA` to do this. It returns an :class:`.MSA` object.

.. ipython:: python

   msa_refine = refineMSA(msa, label='RNAS2_HUMAN', rowocc=0.8, seqid=0.98)
   msa_refine

MSA is refined based on the sequence of :uniprot:`RNAS2_HUMAN`.

Plotting occupancy
-------------------------------------------------------------------------------

Evol plotting functions are prefixed with ``show``. We can plot the ocuupancy
for each column to see if there are any positions in the MSA that have a lot of
gaps. We use the function :func:`.showMSAOccupancy` that uses
:func:`.calcMSAOccupancy` to calculate occupancy for MSA.

.. ipython:: python

   @savefig msa_analysis_occ_res.png width=4in
   showMSAOccupancy(msa_refine, occ='res');

We can also specify indices based on the PDB.

.. ipython:: python

   indices = list(range(4,132))
   @savefig msa_analysis_occ_res_indices.png width=4in
   showMSAOccupancy(msa_refine, occ='res', indices=indices);

We can further refine the MSA to remove positions that have low occupancy, but
that will change the start-end positions of the labels in the MSA that is not
corrected automatically on refinement. We can also plot occupancy based on rows
for the seqeunces in the MSA.

Entropy Calculation
-------------------------------------------------------------------------------

Here, we show how to calculate Shannon Entropy and plot entropy. Entropy for
each position in the MSA is calculated using :func:`.calcShannonEntropy`. It
takes :class:`.MSA` object or a numpy 2D array containg MSA as input. Returns
a 1D numpy arrauy. Plotting is done using :func:`.showShannonEntropy`.

.. ipython:: python

   entropy = calcShannonEntropy(msa_refine)

*entropy* is a 1D numpy array.

.. ipython:: python

   @savefig msa_analysis_entropy.png width=6in
   showShannonEntropy(entropy, indices);


Mutual Information
-------------------------------------------------------------------------------

Here, we show how to calculate mutual information between the positions of the
MSA using :func:`.buildMutinfoMatrix` which also takes  :class:`.MSA` object
or a numpy 2D array containg MSA as input. We can also apply normalization
using :func:`.applyMutinfoNorm` and correction using :func:`.applyMutinfoCorr`
to the mutual information matrix based on references [MLC05]_ and [DSD08]_
respectively. Returns a numpy 2D array.

.. ipython:: python

   mutinfo = buildMutinfoMatrix(msa_refine)
   mutinfo_norm = applyMutinfoNorm(mutinfo, entropy, norm='minent')
   mutinfo_corr = applyMutinfoCorr(mutinfo, corr='apc')

Note that by default ``norm="sument"`` normalization is applied in
``applyMutinfoNorm`` and ``corr="prod"`` is applied in ``applyMutinfoCorr``.

Now we plot the mutual information matrices that we obtained above and see
the effects of different corrections and normalizations.

.. ipython:: python

   @savefig msa_analysis_mutinfo.png width=4in
   showMutinfoMatrix(mutinfo);

.. ipython:: python

   @savefig msa_analysis_mutinfo_corr.png width=4in
   showMutinfoMatrix(mutinfo_corr, clim=[0, mutinfo_corr.max()],
      xlabel='1K2A: 4-131');

Note ylabel does not need to be set, since xlabel = ylabel


Output Results
-------------------------------------------------------------------------------

Here we show how to write the mutual information and entropy array. We use the
:func:`.writeArray` to write numpy array data.

.. ipython:: python

   writeArray('1K2A_MI.txt', mutinfo)


This can be later loaded using :func:`numpy.loadtxt`. Further analysis can also
be done by rank ordering the matrix and analyzing the pairs with highest mutual
information or the most co-evolving residues. This is done using
:func:`.calcRankorder`. A zscore normalization can also be applied to select
coevolving pairs based on a zscore cutoff.

.. ipython:: python

   rank_row, rank_col, zscore_sort = calcRankorder(mutinfo, zscore=True)
   asarray(indices)[rank_row[:5]]
   asarray(indices)[rank_col[:5]]
   zscore_sort[:5]