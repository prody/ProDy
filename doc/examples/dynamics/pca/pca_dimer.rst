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

.. versionadded:: 0.5.3

User input
-------------------------------------------------------------------------------
 
* List of PDB identifiers
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

Import ProDy into the current namespace.

>>> from prody import *

Definitions
-------------------------------------------------------------------------------

Set the name of the protein/dataset (a name without a white space is preferred): 

>>> name = 'HIV-RT'

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

Structures
-------------------------------------------------------------------------------

>>> pdb_ids = ['3kk3', '3kk2', '3kk1', '1suq', '1dtt', '3i0r', '3i0s', '3m8p',
...            '3m8q', '1jlq', '3nbp', '1klm', '2ops', '2opr', '1s9g', '2jle', 
...            '1s9e', '1jla', '1jlc', '1jlb', '1jle', '1jlg', '1jlf', '3drs', 
...            '3e01', '3drp', '1hpz', '3ith', '1s1v', '1s1u', '1s1t', '1ep4', 
...            '3klf', '2wom', '2won', '1s1x', '2zd1', '3kle', '1hqe', '1n5y', 
...            '1fko', '1hnv', '1hni', '1hqu', '1iky', '1ikx', '1t03', '1ikw', 
...            '1ikv', '1t05', '3qip', '3jsm', '1c0t', '1c0u', '2ze2', '1hys', 
...            '1rev', '3dle', '1uwb', '3dlg', '3qo9', '1tv6', '2i5j', '3meg', 
...            '3mee', '3med', '3mec', '3dya', '2be2', '2opp', '3di6', '1tl3', 
...            '1jkh', '1sv5', '1tl1', '1n6q', '2rki', '1tvr', '3klh', '3kli', 
...            '1dtq', '1bqn', '3klg', '1bqm', '3ig1', '2b5j', '1r0a', '3dol', 
...            '1fk9', '2ykm', '1rtd', '1hmv', '3dok', '1rti', '1rth', '1rtj', 
...            '1dlo', '1fkp', '3bgr', '1c1c', '1c1b', '3lan', '3lal', '3lam', 
...            '3lak', '3drr', '2rf2', '1rt1', '1j5o', '1rt3', '1rt2', '1rt5', 
...            '1rt4', '1rt7', '1rt6', '3lp1', '3lp0', '2iaj', '3lp2', '1qe1', 
...            '3dlk', '1s1w', '3isn', '3kjv', '3jyt', '2ban', '3dmj', '2vg5', 
...            '1vru', '1vrt', '1lw2', '1lw0', '2ic3', '3c6t', '3c6u', '3is9', 
...            '2ykn', '1hvu', '3irx', '2b6a', '3hvt', '1tkz', '1eet', '1tkx', 
...            '2vg7', '2hmi', '1lwf', '1tkt', '2vg6', '1s6p', '1s6q', '3dm2', 
...            '1lwc', '3ffi', '1lwe']

Step 1: Retrieve structures
-------------------------------------------------------------------------------

>>> pdb_files = fetchPDB(pdb_ids, folder='pdbfiles')

Step 2: Set reference
-------------------------------------------------------------------------------

>>> # Parse reference structure
>>> reference_structure = parsePDB('pdbfiles/'+ref_pdb+'.pdb.gz', subset='calpha')
>>> # Get the reference chain from this structure
>>> reference_hierview = reference_structure.getHierView() 
>>> reference_chains = [reference_hierview[chid] for chid in ref_chids]
>>> reference_chains
[<Chain: A from 1dlo_ca (556 atoms; 1 coordinate sets, active set index: 0)>, <Chain: B from 1dlo_ca (415 atoms; 1 coordinate sets, active set index: 0)>]
 
Chain A is the p66 domain, and chain B is the p51 domain of HIV-RT.
 
Step 3: Prepare ensemble
-------------------------------------------------------------------------------
 
Instantiate an :class:`~prody.ensemble.PDBEnsemble`

>>> ensemble = PDBEnsemble(name)

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

>>> for pdb_hit, pdb_file in zip(pdb_ids, pdb_files):
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
>>> ensemble
<PDB Ensemble: HIV-RT (155 conformations, 971 atoms, 971 selected)>
>>> ensemble.iterpose()
>>> saveEnsemble(ensemble)
'HIV-RT.ens.npz'

We can now close the logfile using :func:`~prody.closeLogfile`:

>>> closeLogfile(name)

Let's check which structures, if any, are not mapped (added to the ensemble):

>>> print( unmapped )
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

We can find out this using :func:`~prody.ensemble.calcSumOfWeights` function:

>>> print( (calcSumOfWeights(ensemble).min() / len(ensemble)).round(3) )
0.252


This shows that some residues were resolved in only 24% of the dataset.
We trim the ensemble to contain residues resolved in more than 94% of the 
ensemble:

>>> ensemble = trimEnsemble(ensemble, occupancy=0.94)

After trimmin, another round of iterative superposition may be useful:

>>> ensemble.iterpose()
>>> saveEnsemble(ensemble)
'HIV-RT.ens.npz'

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

.. plot::
   :context:
   :nofigs:
   
   from prody import *
   ensemble = loadEnsemble('HIV-RT.ens.npz')
   pca = loadModel('HIV-RT.pca.npz')

Let's plot RMSD to the average structure:

.. plot::
   :context:
   :include-source:
   
   from matplotlib import pyplot as plt
   
   plt.figure(figsize=(5,4))
   plt.plot(calcRMSD(ensemble))
   plt.xlabel('Conformation')
   plt.ylabel('RMSD (A)')
   plt.title(ensemble)

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
   plt.title(ensemble)
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')
   
Only some of the ProDy plotting functions are shown here. A complete list
can be found in :ref:`dynamics` module. 

   
|more| See also other examples in :ref:`pca-xray` for illustration of 
comparative analysis of theoretical and computational data.

|questions|

|suggestions|
