.. currentmodule:: prody.dynamics

.. _pca-blast:

*******************************************************************************
PCA from a Blast search
*******************************************************************************

In this example we will perform PCA of Cytochrome C (cyt C) structures. We will
start with a sequence and a reference structure PDB identifier. 
Blast search will be performed to identify cyt C structures. 
We will set the sequence identity to 44%.
The dataset will contain homologous/orthologous cyt C structures.

**Imports and definitions**

In this part we import ProDy and Matplotlib libraries, and define
the sequence and the reference chain for analysis. This example can
be applied to other proteins by changing the definitions in this part. 

.. plot::
   :context:
   :include-source:
   :nofigs:

   # Import libraries
   from prody import *
   from matplotlib import pyplot as plt

   # Name of the protein (a name without a white space is preferred) 
   name = 'CytC'
   # Amino acid sequence of the protein
   sequence = '''GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETL
   MEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE'''
   # Minimum sequence identity of hits
   sequence_identity = 44
   # Reference PDB file   
   ref_pdb = '1hrc'
   # Reference chain identifier
   ref_chid = 'A'
   # Selection string ("all" can be used if all of the chain is to be analyzed) 
   selstr = 'resnum 1 to 103'
   # Optionally, a list of PDB files to be excluded from analysis can be provided
   # In this case dimeric Cyt C structures are excluded from the analysis
   # If all PDB hits will be used, provide an empty list
   exclude = ['3nbt', '3nbs']
    

**Step 1: Blast search and downloads**

.. plot::
   :context:
   :include-source:
   :nofigs:
   
   blast_record = blastPDB(sequence)
   pdb_hits = blast_record.getHits(sequence_identity).keys()
   pdb_files = fetchPDB(pdb_hits, folder='pdbfiles')

**Step 2: Prepare the reference structure and the ensemble**

.. plot::
   :context:
   :include-source:
   :nofigs:

   # Parse reference structure
   reference_structure = parsePDB('pdbfiles/'+ref_pdb+'.pdb.gz', subset='calpha')
   # Get the reference chain from this structure
   reference_hierview = reference_structure.getHierView() 
   reference_chain = reference_hierview[ref_chid]
 
   # Instantiate an ensemble
   ensemble = Ensemble(name)
   # Set reference coordinates
   ensemble.setCoordinates(reference_chain.getCoordinates())
   
**Step 3: Parse PDB files and populate the ensemble**
 
.. plot::
   :context:
   :include-source:

   # Parse hits 
   for pdb_hit, pdb_file in zip(pdb_hits, pdb_files):
       # Skip the PDB file if its in the exclude list
       skip = False
       for pdb in exclude:
           if pdb in pdb_file:
               skip = True
       if skip: continue
       
       # Parse the current PDB file   
       current_structure = parsePDB(pdb_file, subset='calpha')
       # Map current PDB file to the reference chain
       current_mapping = mapAtomsToChain(current_structure, reference_chain, seqid=sequence_identity)
       # 
       current_atommap = current_mapping[0][0]
       ensemble.addCoordset(current_atommap, weights=current_atommap.getMappedFlags()) 
   ensemble.iterpose()

   # Let's plot RMSD to the average structure
   plt.figure(figsize=(5,4))
   plt.plot(calcRMSD(ensemble))
   plt.xlabel('Conformation index')
   plt.ylabel('RMSD (A)')
   
   # Write aligned conformations into a PDB file
   reference_structure.addCoordset(ensemble.getCoordsets())
   writePDB(name+'.pdb', reference_structure)
   # This file can be used to visualize the aligned conformations in modeling software

.. plot::
   :context:
   :nofigs:
   
   plt.close('all')
   import os
   os.remove(name+'.pdb')

**Step 4: Perform PCA**

.. plot::
   :context:
   :include-source:

   # Instantiate a PCA
   pca = PCA(name)
   # Build covariance matrix
   pca.buildCovariance(ensemble)
   # Calculate modes
   pca.calcModes()
   
   # Let's show a projection of the ensemble onto PC1 and PC2
   plt.figure(figsize=(5,4))
   showProjection(ensemble, pca[:2])

**Further reading**
   
User is referred to other examples in :ref:`pca` for illustration of comparative 
analysis of theoretical and computational data. 
