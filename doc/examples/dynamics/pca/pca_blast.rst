.. currentmodule:: prody.dynamics

.. _pca-blast:

*******************************************************************************
PCA from a Blast search
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform PCA of a structural dataset obtained by blast
searching PDB. The protein of interest is cytochrome C (*cyt C*). 
Dataset will contain structures sharing 44% or more 
sequence identity with human *cyt C*, i.e. its homologs and/or orthologs.

Input
-------------------------------------------------------------------------------
 
* Amino acid sequence of the protein
* A reference PDB structure
* Optionally, a list of PDB files to be excluded from the analysis 
* Percent sequence identity used for selecting blast hits (PDB structures)
* Selection of the protein chains and residues to be considered in analysis

Output
-------------------------------------------------------------------------------

A :class:`PCA` instance that stores covariance matrix and principal modes
that describes the dominant changes in the dataset. :class:`PCA` instance
and principal modes (:class:`Mode`) can be used as input to functions in 
:mod:`~prody.dynamics` module for further analysis.
 
ProDy Code
===============================================================================

Imports
-------------------------------------------------------------------------------

Import ProDy and matplotlib into the current namespace.

>>> from prody import *
>>> from matplotlib import pyplot as plt

Definitions
-------------------------------------------------------------------------------

>>> # Name of the protein (a name without a white space is preferred) 
>>> name = 'CytC'

>>> # Amino acid sequence of the protein
>>> sequence = '''GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETL
... MEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE'''

>>> # Reference PDB file   
>>> ref_pdb = '1hrc'

>>> # Optionally, a list of PDB files to be excluded from analysis can be provided
>>> # In this case dimeric Cyt C structures are excluded from the analysis
>>> # If all PDB hits will be used, provide an empty list
>>> exclude = ['3nbt', '3nbs']

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
>>> pdb_hits = blast_record.getHits(sequence_identity)
>>> pdb_files = fetchPDB(pdb_hits.keys(), folder='pdbfiles')

Let's check number of downloaded files:

>>> len(pdb_files)
82

Step 2: Set reference
-------------------------------------------------------------------------------

We first parse the reference structure. Note that we parse only Cα atoms from
chain A. The analysis will be performed for a single chain (monomeric) protein.
For analysis of a dimeric protein see :ref:`pca-dimer`

>>> reference_structure = parsePDB('pdbfiles/'+ref_pdb+'.pdb.gz', 
...                                subset='calpha', chain=ref_chid)
>>> # Get the reference chain from this structure
>>> reference_hierview = reference_structure.getHierView() 
>>> reference_chain = reference_hierview[ref_chid]
 
Step 3: Prepare ensemble
-------------------------------------------------------------------------------
 
>>> # Start a log file
>>> startLogfile('pca_blast') 
>>> # Instantiate a PDB ensemble
>>> ensemble = PDBEnsemble(name)
>>> # Set reference coordinates
>>> ensemble.setCoordinates(reference_chain.getCoordinates())
   
>>> # Parse hits 
>>> for pdb_hit, pdb_file in zip(pdb_hits.keys(), pdb_files):
...     # Skip the PDB file if its in the exclude list
...     if pdb_hit in exclude:
...         continue
...     
...     # Parse the current PDB file   
...     structure = parsePDB(pdb_file, subset='calpha', 
...                          chain=pdb_hits[pdb_hit]['chain_id'])
...     if structure is None:
...         plog('Failed to parse ' + pdb_file)
...         continue
...     # Map current PDB file to the reference chain
...     mappings = mapOntoChain(structure, reference_chain, seqid=sequence_identity)
...     if len(mappings) == 0:
...         plog('Failed to map', pdb_hit)
...         continue  
...     atommap = mappings[0][0]
...     ensemble.addCoordset(atommap, weights=atommap.getMappedFlags())
>>> ensemble.iterpose()
>>> saveEnsemble(ensemble)
'CytC.ens.npz'

Let's check how many conformations are extracted from PDB files:

>>> len(ensemble)
349

Note that number of conformations are more than the number of PDB structures
we evaluated. This is because some of the PDB files contained NMR structures
with multiple models. Each model in NMR structures are added to the ensemble
as individual conformations.

Write aligned conformations into a PDB file as follows:

>>> reference_structure.addCoordset(ensemble.getCoordsets())
>>> writePDB(name+'.pdb', reference_structure)
'CytC.pdb'

This file can be used to visualize the aligned conformations in a modeling 
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
   
The calculated data can be saved as a compressed file using :func:`saveModel`
function:

>>> saveModel(pca)
'CytC.pca.npz'

Step 5: Plot data and results
-------------------------------------------------------------------------------

.. plot::
   :context:
   :nofigs:
   
   from prody import *
   from matplotlib import pyplot as plt
   ensemble = loadEnsemble('CytC.ens.npz')
   pca = loadModel('CytC.pca.npz')

Let's plot RMSD to the average structure:


.. plot::
   :context:
   :include-source:

   rmsd = calcRMSD(ensemble)

   plt.figure(figsize=(5,4))
   plt.plot( rmsd )
   plt.xlabel('Conformation index')
   plt.ylabel('RMSD (A)')

Let's show a projection of the ensemble onto PC1 and PC2:

.. plot::
   :context:
   :include-source:

   plt.figure(figsize=(5,4))
   showProjection(ensemble, pca[:2])

See Also
===============================================================================
   
User is referred to other examples in :ref:`pca-xray` for illustration of 
comparative analysis of theoretical and computational data.

|questions|

|suggestions|
