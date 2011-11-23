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

The results are displayed for following list of structures:

>>> pdb_hits = [('1ccr', 'A'), ('2frc', 'A'), ('1fhb', 'A'), ('1chh', 'A'), 
...  ('1chj', 'A'), ('1lc2', 'A'), ('1cie', 'A'), ('1kyo', 'W'), ('1m60', 'A'),
...  ('1fi7', 'A'), ('2bgv', 'X'), ('2hv4', 'A'), ('2aiu', 'A'), ('2bcn', 'B'),
...  ('1csu', 'A'), ('1csw', 'A'), ('1csv', 'A'), ('1cyc', 'A'), ('2giw', 'A'),
...  ('1fi9', 'A'), ('1i55', 'A'), ('1i54', 'A'), ('3nbs', 'A'), ('1nmi', 'A'),
...  ('1cry', 'A'), ('1ocd', 'A'), ('1ycc', 'A'), ('1cih', 'A'), ('1lms', 'A'),
...  ('2ycc', 'A'), ('1hrc', 'A'), ('2gb8', 'B'), ('2b12', 'B'), ('2b11', 'B'),
...  ('3cx5', 'W'), ('1io3', 'A'), ('1wej', 'F'), ('2jqr', 'A'), ('2b0z', 'B'),
...  ('1csx', 'A'), ('1raq', 'A'), ('1u75', 'B'), ('3cyt', 'O'), ('1crc', 'A'),
...  ('1cig', 'A'), ('3nwv', 'A'), ('1crg', 'A'), ('2b4z', 'A'), ('1yfc', 'A'),
...  ('1crj', 'A'), ('1cif', 'A'), ('1crh', 'A'), ('1cri', 'A'), ('2w9k', 'A'),
...  ('1cty', 'A'), ('1ctz', 'A'), ('2orl', 'A'), ('2pcc', 'B'), ('2pcb', 'B'),
...  ('5cyt', 'R'), ('1akk', 'A'), ('1ytc', 'A'), ('1co6', 'A'), ('1rap', 'A'),
...  ('3nbt', 'A'), ('1i5t', 'A'), ('1yic', 'A'), ('1irv', 'A'), ('1irw', 'A'),
...  ('2yk3', 'A'), ('1giw', 'A'), ('2jti', 'B'), ('3cxh', 'W'), ('1lc1', 'A'),
...  ('1s6v', 'B'), ('1yea', 'A'), ('1yeb', 'A'), ('1lfm', 'A'), ('2b10', 'B'),
...  ('1u74', 'B'), ('1j3s', 'A'), ('1chi', 'A'),]

List of PDB structures can be updated using :func:`~prody.proteins.blastPDB` 
as follows::

  blast_record = blastPDB(sequence)
  pdb_hits = []
  for key, item blast_record.getHits(sequence_identity).iteritems():
      pdb_hits.append((key, item['chain_id']))

>>> pdb_files = fetchPDB([pdb for pdb, ch in pdb_hits], folder='pdbfiles', compressed=False)

Let's check number of downloaded files:

>>> len(pdb_files)
82

Step 2: Set reference
-------------------------------------------------------------------------------

We first parse the reference structure. Note that we parse only Cα atoms from
chain A. The analysis will be performed for a single chain (monomeric) protein.
For analysis of a dimeric protein see :ref:`pca-dimer`

>>> reference_structure = parsePDB('pdbfiles/'+ref_pdb+'.pdb', 
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
>>> ensemble.setCoords(reference_chain.getCoords())
   
>>> # Parse hits 
>>> for pdb_hit, pdb_file in zip(pdb_hits, pdb_files):
...     pdb_id, chain_id = pdb_hit
...     # Skip the PDB file if its in the exclude list
...     if pdb_id in exclude:
...         continue
...     
...     # Parse the current PDB file   
...     structure = parsePDB(pdb_file, subset='calpha', chain=chain_id)
...     if structure is None:
...         plog('Failed to parse ' + pdb_file)
...         continue
...     # Map current PDB file to the reference chain
...     mappings = mapOntoChain(structure, reference_chain, seqid=sequence_identity)
...     if len(mappings) == 0:
...         plog('Failed to map', pdb_id)
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
   plt.close('all')

Let's plot RMSD to the average structure:


.. plot::
   :context:
   :include-source:

   rmsd = calcRMSD(ensemble)

   plt.figure(figsize=(5,4))
   plt.plot( rmsd )
   plt.xlabel('Conformation index')
   plt.ylabel('RMSD (A)')


.. plot::
   :context:
   :nofigs:

   plt.close('all')
   
   
Let's show a projection of the ensemble onto PC1 and PC2:

.. plot::
   :context:
   :include-source:

   plt.figure(figsize=(5,4))
   showProjection(ensemble, pca[:2])


.. plot::
   :context:
   :nofigs:

   plt.close('all')
   

See Also
===============================================================================
   
User is referred to other examples in :ref:`pca-xray` for illustration of 
comparative analysis of theoretical and computational data.

|questions|

|suggestions|
