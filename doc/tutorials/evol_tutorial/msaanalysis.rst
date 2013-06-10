.. _msaanalysis:

Evolution Analysis
===============================================================================

This part follows from :ref:`msafiles`. The aim of this part is to show how to
calculate residue conservation and coevolution properties based on multiple
sequence alignments (MSAs). MSA


First, we import everything from the ProDy package.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()  # turn interactive mode on

Get MSA data
-------------------------------------------------------------------------------

Let's download full MSA file for protein family :pfam:`RnaseA`.
We can do this by specifying the PDB ID of a protein in this family.

.. ipython::
   :verbatim:

   In [1]: searchPfam('1K2A').keys()
   Out[1]: ['PF00074']
   In [2]: fetchPfamMSA('PF00074')
   Out[2]: 'PF00074_full.sth'

Let's parse the downloaded file:

.. ipython:: python

   msa = parseMSA('PF00074_full.sth')


Refine MSA
-------------------------------------------------------------------------------

Here, we refine the MSA to decrease the number of gaps.  We will remove any
columns in the alignment for which there is a gap in the specified PDB file,
and then remove any rows that have more than 20% gaps.  :func:`.refineMSA`
does all of this and returns an :class:`.MSA` object.

.. ipython:: python

   msa_refine = refineMSA(msa, label='RNAS2_HUMAN', rowocc=0.8, seqid=0.98)
   msa_refine

MSA is refined based on the sequence of :uniprot:`RNAS2_HUMAN`, corresponding
to :pdb:`1K2A`.

Occupancy calculation
-------------------------------------------------------------------------------

Evol plotting functions are prefixed with ``show``. We can plot the occupancy
for each column to see if there are any positions in the MSA that have a lot of
gaps. We use the function :func:`.showMSAOccupancy` that uses
:func:`.calcMSAOccupancy` to calculate occupancy for MSA.

.. ipython:: python

   @savefig msa_analysis_occ_res.png width=4in
   showMSAOccupancy(msa_refine, occ='res');

Let's find the minimum:

.. ipython:: python

   calcMSAOccupancy(msa_refine, occ='res').min();

We can also specify indices based on the PDB.

.. ipython:: python

   indices = list(range(4,132))
   @savefig msa_analysis_occ_res_indices.png width=4in
   showMSAOccupancy(msa_refine, occ='res', indices=indices);

Further refining the MSA to remove positions that have low occupancy will
change the start and end positions of the labels in the MSA. This is not
corrected automatically on refinement. We can also plot occupancy based on
rows for the sequences in the MSA.

Entropy Calculation
-------------------------------------------------------------------------------

Here, we show how to calculate and plot Shannon Entropy. Entropy for
each position in the MSA is calculated using :func:`.calcShannonEntropy`. It
takes :class:`.MSA` object or a numpy 2D array containg MSA as input and returns
a 1D numpy array.

.. ipython:: python

   entropy = calcShannonEntropy(msa_refine)
   entropy


*entropy* is a 1D Numpy array. Plotting is done using
:func:`.showShannonEntropy`.

.. ipython:: python

   @savefig msa_analysis_entropy.png width=6in
   showShannonEntropy(entropy, indices);


Mutual Information
-------------------------------------------------------------------------------

We can calculate mutual information between the positions of the MSA using
:func:`.buildMutinfoMatrix` which also takes an :class:`.MSA` object
or a numpy 2D array containing MSA as input.

.. ipython:: python

   mutinfo = buildMutinfoMatrix(msa_refine)
   mutinfo

Result is a 2D Numpy array.

We can also apply normalization using :func:`.applyMutinfoNorm` and
correction using :func:`.applyMutinfoCorr` to the mutual information matrix
based on references [MLC05]_ and [DSD08]_, respectively.

.. ipython:: python

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


Output Results
-------------------------------------------------------------------------------

Here we show how to write the mutual information and entropy arrays to file.
We use the :func:`.writeArray` to write Numpy array data.

.. ipython:: python

   writeArray('1K2A_MI.txt', mutinfo)


This can be later loaded using :func:`.parseArray`.

Rank-ordering
-------------------------------------------------------------------------------

Further analysis can also be done by rank ordering the matrix and analyzing
the pairs with highest mutual information or the most co-evolving residues.
This is done using :func:`.calcRankorder`. A z-score normalization can also
be applied to select coevolving pairs based on a z score cutoff.

.. ipython:: python

   rank_row, rank_col, zscore_sort = calcRankorder(mutinfo, zscore=True)
   asarray(indices)[rank_row[:5]]
   asarray(indices)[rank_col[:5]]
   zscore_sort[:5]
