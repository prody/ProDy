.. _comparison:

Sequence-Structure Comparison
===============================================================================

The part shows how to compare sequence conservation properties with
structural mobility obtained from Gaussian network model (GNM) calculations.

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()  # turn interactive mode on


Entropy Calculation
-------------------------------------------------------------------------------

First, we retrieve MSA for protein for protein family :pfam:`PF00074`:

.. ipython::
   :verbatim:


   In [1]: fetchPfamMSA('PF00074')
   Out[1]: 'PF00074_full.sth'

We parse the MSA file:

.. ipython:: python

   msa = parseMSA('PF00074_full.sth')

Then, we refine it using :func:`.refineMSA` based on the sequence of
:uniprot:`RNAS1_BOVIN`:

.. ipython:: python

   msa_refine = refineMSA(msa, label='RNAS1_BOVIN', rowocc=0.8, seqid=0.98)


We calculate the entropy for the refined MSA:

.. ipython:: python

   entropy = calcShannonEntropy(msa_refine)


Mobility Calculation
-------------------------------------------------------------------------------

Next, we obtain residue fluctuations or mobility for protein member of the
above family. We will use chain B of :pdb:`2W5I`.


.. ipython:: python

   pdb = parsePDB('2W5I', chain='B')
   chB_ca = pdb.select('protein and name CA and resid 1 to 121')

We perform GNM as follows:

.. ipython:: python

   gnm = GNM('2W5I')
   gnm.buildKirchhoff(chB_ca)
   gnm.calcModes(n_modes=None)  # calculate all modes

Now, let's obtain residue mobility using slowest mode, slowest 8 modes,
and all modes:


.. ipython:: python

   mobility_1 = calcSqFlucts(gnm[0])
   mobility_1to8 = calcSqFlucts(gnm[:8])
   mobility_all = calcSqFlucts(gnm[:])


See :ref:`gnm` for details.

Comparison of mobility and conservation
-------------------------------------------------------------------------------

We use the above data to compare structural mobility and degree of
conservation. We can calculate a correlation coefficient between the two
quantities:

.. ipython:: python

   result = corrcoef(mobility_all, entropy)
   result.round(3)[0,1]

We can plot the two curves simultaneously to visualize the correlation.
We have to scale the values of mobility to display them in the same plot.

Plotting
^^^^^^^^

.. ipython:: python

   indices = range(1,122)
   bar(indices, entropy, width=1.2, color='grey', hold='True');
   xlim(min(indices)-1, max(indices)+1);
   @savefig entropy_mobility.png width=4in
   plot(indices, mobility_all*(max(entropy)/mean(mobility_all)), color='b',
   linewidth=2);


Writing PDB files
-------------------------------------------------------------------------------

We can also write PDB with b-factor column replaced by entropy and mobility
values respectively. We can then load the PDB structure in VMD or PyMol to
see the distribution of entropy and mobility on the structure.

.. ipython:: python

   selprot = pdb.select('protein and resid 1 to 121')
   resindex = selprot.getResindices()
   index = unique(resindex)
   count = 0
   entropy_prot = []
   mobility_prot = []
   for ind in index:
       while(count < len(resindex)):
           if(ind == resindex[count]):
               entropy_prot.append(entropy[ind])
               mobility_prot.append(mobility_all[ind]*100)
           count = count + 1
   selprot.setBetas(entropy_prot)
   writePDB('2W5I_entropy.pdb', selprot)
   selprot.setBetas(mobility_prot)
   writePDB('2W5I_mobility.pdb', selprot)