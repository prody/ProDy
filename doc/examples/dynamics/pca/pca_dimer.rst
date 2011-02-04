.. currentmodule:: prody.dynamics

.. _pca-dimer:

*******************************************************************************
PCA of a dimeric protein 
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform PCA of structural dataset for a dimeric
protein. The protein of interest is HIV reverse transcriptase (RT). 
Dataset will be obtained by blast searching PDB.

.. versionadded:: 0.6

User input
-------------------------------------------------------------------------------
 
* Amino acid sequence of the protein
* A reference PDB structure
* List of PDB files to be excluded from the analysis, if any 

Parameters
-------------------------------------------------------------------------------

* Percent sequence identity used for selecting blast hits (PDB structures)
* Selection of the RT chains and residues to be considered in analysis

How to Use
-------------------------------------------------------------------------------

This example can be user in the following ways:

  * in an interactive Python session (start a Python interpreter and insert
    code lines one by one).
  * in an Python script to perform calculations all at once (copy the code
    into a Python file and run it). 

Notes
-------------------------------------------------------------------------------

* This example needs internet connectivity for blast searching PDB and 
  retrieving files from PDB FTP server.

* Also note that this example will attempt to download well over 100 structure 
  files, which make take several minutes depending on connection speed.  

* For plotting results, |matplotlib| library is required.

 
ProDy Code
===============================================================================

Imports
-------------------------------------------------------------------------------

Import ProDy into the current namespace.

>>> from prody import *

Definitions
-------------------------------------------------------------------------------

Set the name of the protein/dataset (a name without a white space is preferred): 

>>> name = 'HIV-RT'

>>> # Amino acid sequence of the protein
>>> sequence = '''PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFREL
... NKRTQDFWEVQLGIPHPAGLKKNKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIF
... QSSMTKILEPFKKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWT
... VQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLSKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYD
... PSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWET
... WWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNKGRQKVVPLTNTTNQKTELQ
... AIYLALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKIL'''

Set the reference PDB file:    

>>> ref_pdb = '1dlo'


Parameters
-------------------------------------------------------------------------------

Set the minimum sequence identity of hits:

>>> sequence_identity = 94

Set reference chain identifiers. HIV-RT 

>>> ref_chids = 'AB'

Set the sequence coverage:

>>> sequence_coverage = 85

Step 1: Blast and download
-------------------------------------------------------------------------------

>>> blast_record = blastPDB(sequence)
>>> pdb_hits = blast_record.getHits(sequence_identity).keys()
>>> pdb_files = fetchPDB(pdb_hits, folder='pdbfiles')

Step 2: Set reference
-------------------------------------------------------------------------------

>>> # Parse reference structure
>>> reference_structure = parsePDB('pdbfiles/'+ref_pdb+'.pdb.gz', subset='calpha')
>>> # Get the reference chain from this structure
>>> reference_hierview = reference_structure.getHierView() 
>>> reference_chains = [reference_hierview[chid] for chid in ref_chids]
>>> reference_chains
[<Chain: A from 1dlo (556 atoms; 1 coordinate sets, active set index: 0)>, <Chain: B from 1dlo (415 atoms; 1 coordinate sets, active set index: 0)>]
 
Chain A is the p53 domain, and chain B is the p38 domain of HIV-RT.
 
Step 3: Prepare ensemble
-------------------------------------------------------------------------------
 
Instantiate an :class:`~prody.ensemble.Ensemble`

>>> ensemble = Ensemble(name)

We now combine the reference chains and set the reference coordinates 
of the ensemble.

>>> reference_chain = reference_chains[0] + reference_chains[1]
>>> ensemble.setCoordinates(reference_chain.getCoordinates())

We also start a log file using :func:`prody.startLogfile`. 
Screen output will be save in this file, and can be
used to check if structures are added to the ensemble as expected.

>>> startLogfile(name)

We defined a list to keep track of PDB files that are not added to the ensemble:

>>> unmapped = []

Now, we parse the PDB files one by one and add them to the ensemble: 

>>> for pdb_hit, pdb_file in zip(pdb_hits, pdb_files):
...     # Parse the PDB file   
...     structure = parsePDB(pdb_file, subset='calpha', model=1)
...     atommaps = []
...     for reference_chain in reference_chains:
...         # Map current PDB file to the reference chain
...         mappings = mapOntoChain(structure, reference_chain, seqid=sequence_identity, coverage=sequence_coverage)
...         if len(mappings) == 0:
...             print 'Failed to map', pdb_hit
...             break
...         atommaps.append(mappings[0][0])
...         # Make sure all chains are mapped
...     if len(atommaps) != len(reference_chains):
...         unmapped.append(pdb_hit)
...         continue
...     atommap = atommaps[0] + atommaps[1]
...     ensemble.addCoordset(atommap, weights=atommap.getMappedFlags()) 
>>> ensemble.iterpose()
>>> saveEnsemble(ensemble)
'HIV-RT.ens.npz'

We can now close the logfile using :func:`~prody.closeLogfile`:

>>> closeLogfile(name)

Let's check which structures, if any, are not mapped (added to the ensemble):

>>> print unmapped
[]

We can write the aligned conformations into a PDB file as follows:

>>> reference_structure.addCoordset(ensemble.getCoordsets())
>>> writePDB(name+'.pdb', reference_structure)
'HIV-RT.pdb'

This file can be used to visualize the aligned conformations in modeling 
software.

This is a heterogeneous dataset, i.e. many structures had missing residues.
We want to make sure that we include residues in PCA analysis if they
are resolved in more than 94% of the time.

We can show this by using :func:`~prody.ensemble.showSumOfWeights` function:

.. plot::
   :context:
   :nofigs:
   
   from prody import *
   ensemble = loadEnsemble('HIV-RT.ens.npz')

.. plot::
   :context:
   :include-source:
   
   from matplotlib import pyplot as plt
   plt.figure(figsize=(5,4))
   showSumOfWeights(ensemble)

.. plot::
   :context:
   :nofigs:
   
   plt.close('all')

This shows that some residues were resolved only in around 40 structures,
which about 30% of the dataset.
We trim the ensemble to contain residues resolved in more than 94% of the 
ensemble:

>>> ensemble = trimEnsemble(ensemble, occupancy=0.94)

After trimmin, another round of iterative superposition may be useful:

>>> ensemble.iterpose() 

Step 4: Perform PCA
-------------------------------------------------------------------------------

Once the ensemble is ready, performing :class:`PCA` is 3 easy steps:

>>> pca = PCA(name)
>>> pca.buildCovariance(ensemble)
>>> pca.calcModes()
   
The calculated data can be saved as a compressed file using :func:`saveModel`

>>> saveModel(pca) 
'HIV-RT.pca.npz'

Step 5: Plot data and results
-------------------------------------------------------------------------------

Let's plot RMSD to the average structure:

.. plot::
   :context:
   :nofigs:
   
   plt.close('all')
   pca = loadModel('HIV-RT.pca.npz')

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   plt.plot(calcRMSD(ensemble))
   plt.xlabel('Conformation index')
   plt.ylabel('RMSD (A)')

Let's show a projection of the ensemble onto PC1 and PC2:

.. plot::
   :context:
   :nofigs:
   
   plt.close('all')

.. plot::
   :context:
   :include-source:

   plt.figure(figsize=(5,4))
   showProjection(ensemble, pca[:2])

Only some of the ProDy plotting functions are shown here. A complete list
can be found in :ref:`dynamics` module. 

See Also
===============================================================================
   
User is referred to other examples in :ref:`pca` for illustration of 
comparative analysis of theoretical and computational data.

|questions|

|suggestions|
