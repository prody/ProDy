.. currentmodule:: prody.dynamics

.. _pca-dimer:

*******************************************************************************
PCA of a protein dimer (HIV-RT) dataset 
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform PCA of structural dataset for a dimeric
protein. The protein of interest is HIV reverse transcriptase (RT). 
Dataset will be obtained by blast searching PDB.

User Inputs
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

>>> # Name of the protein (a name without a white space is preferred) 
>>> name = 'HIV-RT'

>>> # Amino acid sequence of the protein
>>> sequence = '''PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFREL
... NKRTQDFWEVQLGIPHPAGLKKNKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIF
... QSSMTKILEPFKKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWT
... VQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLSKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYD
... PSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWET
... WWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNKGRQKVVPLTNTTNQKTELQ
... AIYLALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKIL'''

>>> # Reference PDB file   
>>> ref_pdb = '1hrc'


Parameters
-------------------------------------------------------------------------------

>>> # Minimum sequence identity of hits
>>> sequence_identity = 44

>>> # Reference chain identifier
>>> ref_chid = 'A'

>>> # Selection string ("all" can be used if all of the chain is to be analyzed) 
>>> selstr = 'resnum 1 to 103'
    

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
>>> reference_chain = reference_hierview[ref_chid]
 
Step 3: Prepare ensemble
-------------------------------------------------------------------------------
 
>>> # Instantiate an ensemble
>>> ensemble = Ensemble(name)
>>> # Set reference coordinates
>>> ensemble.setCoordinates(reference_chain.getCoordinates())
   
>>> # Parse hits 
>>> for pdb_hit, pdb_file in zip(pdb_hits, pdb_files):
...     # Skip the PDB file if its in the exclude list
...     if pdb_hit in exclude:
...         continue
...     
...     # Parse the current PDB file   
...     current_structure = parsePDB(pdb_file, subset='calpha', model=1)
...     # Map current PDB file to the reference chain
...     current_mapping = mapOntoChain(current_structure, reference_chain, seqid=sequence_identity)
...     current_atommap = current_mapping[0][0]
...     ensemble.addCoordset(current_atommap, weights=current_atommap.getMappedFlags()) 
>>> ensemble.iterpose()
>>> saveEnsemble(ensemble)
'CytC.ens.npz'

Write aligned conformations into a PDB file as follows:

>>> reference_structure.addCoordset(ensemble.getCoordsets())
>>> writePDB(name+'.pdb', reference_structure)
'CytC.pdb'

This file can be used to visualize the aligned conformations in modeling 
software.


>>> plt.close('all')
>>> import os
>>> os.remove(name+'.pdb')

Step 4: Perform PCA
-------------------------------------------------------------------------------

Once the ensemble is ready, performing PCA is 3 easy steps:

>>> # Instantiate a PCA
>>> pca = PCA(name)
>>> # Build covariance matrix
>>> pca.buildCovariance(ensemble)
>>> # Calculate modes
>>> pca.calcModes()
   
The calculated data can be saved as a compressed file using saveModel 

Step 5: Plot data and results
-------------------------------------------------------------------------------

Let's plot RMSD to the average structure:

>>> plt.figure(figsize=(5,4))
>>> plt.plot(calcRMSD(ensemble))
>>> plt.xlabel('Conformation index')
>>> plt.ylabel('RMSD (A)')

Let's show a projection of the ensemble onto PC1 and PC2:

>>> plt.figure(figsize=(5,4))
>>> showProjection(ensemble, pca[:2])


See Also
===============================================================================
   
User is referred to other examples in :ref:`pca` for illustration of 
comparative analysis of theoretical and computational data.
