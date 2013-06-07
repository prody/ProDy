.. _comparison:

Sequence-Structure Comparison
===============================================================================

Synopsis
-------------------------------------------------------------------------------

The part shows how to:

  * Compare sequence conservation properties with mobility obtained from GNM.

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()  # turn interactive mode on

   
Entropy Calculation
-------------------------------------------------------------------------------

Fisrt, we obtain conservation for protein for :pfam:`PF00074` family.
See also :ref:`msaanalysis`: 

.. ipython:: python

   msa_refine = refineMSA(parseMSA(fetchPfamMSA('PF00074', timeout=40)),
   label='RNAS1_BOVIN', rowocc=0.8, seqid=0.98)
   entropy = calcShannonEntropy(msa_refine)
   
   
GNM Mobility Calculation
-------------------------------------------------------------------------------

Next, we obtain residue fluctuations or mobility for a PDB member of the above 
family like :pdb:`2W5I` or specifically ``2W5IB``.
See also :ref:`gnm`: 

.. ipython:: python

   pdb = parsePDB('2W5I', chain='B')
   selca = pdb.select('protein and name CA and resid 1 to 121')
   gnm = GNM('2W5I')
   gnm.buildKirchhoff(selca)
   gnm.calcModes(n_modes=None)  # calculate all modes
   mobility_1 = calcSqFlucts(gnm[0])
   mobility_1to8 = calcSqFlucts(gnm[:8])
   mobility_all = calcSqFlucts(gnm[:])
   indices = range(1,122)
  

Comparison of mobility and conservation
-------------------------------------------------------------------------------

Use the above data to compare structural mobility and degree of conservation.
We can calculate a correlation coefficient between the two quantities. 

.. ipython:: python

   result = corrcoef(mobility_all, entropy)
   print(result.round(3)[0,1])
   
We can plot the two curves simulataneously to visualize the correlation. We have
to scale the values of mobility to display them in the same plot. 

Plotting
^^^^^^^^^
.. ipython:: python

   bar(indices,entropy,width=1.2,color='grey', hold='True'); 
   xlim(min(indices)-1, max(indices)+1);
   @savefig entropy_mobility.png width=4in
   plot(indices, mobility_all*(max(entropy)/mean(mobility_all)), color='b', 
   linewidth=2);
   

Writing PDB files
-------------------------------------------------------------------------------

We can also write PDB with bfactor column replaced by entropy and mobility 
values respectively. We can then load the PDB structure in vmd or pymol to
see the distribution of entropy and mobility on the structure. 

.. ipython:: python

   selprot = pdb.select('protein and resid 1 to 121')
   resindex = selprot.getResindices()
   index = unique(resindex)
   count = 0; entropy_prot = []; mobility_prot = []
   for ind in index:
       while(count < len(resindex)):
           if(ind == resindex[count]):
               entropy_prot.append(entropy[ind])
               mobility_prot.append(mobility_all[ind]*100)
           count = count + 1
   selprot.setBetas(entropy_prot)
   writePDB('2W5I_Entropy.pdb', selprot)
   selprot.setBetas(mobility_prot)
   writePDB('2W5I_Mobility.pdb', selprot)