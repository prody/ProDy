.. _pca-dimer:


Multimeric Structures
===============================================================================


Synopsis
-------------------------------------------------------------------------------

In this part, we perform PCA of HIV :wiki:`Reverse Transcriptase` (RT), which
is a dimeric protein.


Input and Parameters
-------------------------------------------------------------------------------

First, we make necessary imports:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Reference structure
^^^^^^^^^^^^^^^^^^^

We set the name of the protein/dataset (a name without a white space is
preferred) and also reference structure id and chain identifiers:

.. ipython:: python


   name = 'HIV-RT'  # dataset name
   ref_pdb = '1dlo'  # reference PDB file
   ref_chids = 'AB'  # reference chain identifiers


Parameters
^^^^^^^^^^

Following parameters are for comparing two structures to determine matching
chain.

.. ipython:: python

   sequence_identity = 94
   sequence_coverage = 85

Chains from two different structures will be paired if they share
94% sequence identity and the aligned part of the sequences cover
85% of the longer sequence.

Structures
^^^^^^^^^^

We are going to use the following list of structures:

.. ipython:: python

   pdb_ids = ['3kk3', '3kk2', '3kk1', '1suq', '1dtt', '3i0r', '3i0s', '3m8p',
              '3m8q', '1jlq', '3nbp', '1klm', '2ops', '2opr', '1s9g', '2jle',
              '1s9e', '1jla', '1jlc', '1jlb', '1jle', '1jlg', '1jlf', '3drs',
              '3e01', '3drp', '1hpz', '3ith', '1s1v', '1s1u', '1s1t', '1ep4',
              '3klf', '2wom', '2won', '1s1x', '2zd1', '3kle', '1hqe', '1n5y',
              '1fko', '1hnv', '1hni', '1hqu', '1iky', '1ikx', '1t03', '1ikw',
              '1ikv', '1t05', '3qip', '3jsm', '1c0t', '1c0u', '2ze2', '1hys',
              '1rev', '3dle', '1uwb', '3dlg', '3qo9', '1tv6', '2i5j', '3meg',
              '3mee', '3med', '3mec', '3dya', '2be2', '2opp', '3di6', '1tl3',
              '1jkh', '1sv5', '1tl1', '1n6q', '2rki', '1tvr', '3klh', '3kli',
              '1dtq', '1bqn', '3klg', '1bqm', '3ig1', '2b5j', '1r0a', '3dol',
              '1fk9', '2ykm', '1rtd', '1hmv', '3dok', '1rti', '1rth', '1rtj',
              '1dlo', '1fkp', '3bgr', '1c1c', '1c1b', '3lan', '3lal', '3lam',
              '3lak', '3drr', '2rf2', '1rt1', '1j5o', '1rt3', '1rt2', '1rt5',
              '1rt4', '1rt7', '1rt6', '3lp1', '3lp0', '2iaj', '3lp2', '1qe1',
              '3dlk', '1s1w', '3isn', '3kjv', '3jyt', '2ban', '3dmj', '2vg5',
              '1vru', '1vrt', '1lw2', '1lw0', '2ic3', '3c6t', '3c6u', '3is9',
              '2ykn', '1hvu', '3irx', '2b6a', '3hvt', '1tkz', '1eet', '1tkx',
              '2vg7', '2hmi', '1lwf', '1tkt', '2vg6', '1s6p', '1s6q', '3dm2',
              '1lwc', '3ffi', '1lwe']

A predefined set of structures will be used, but an up-to-date list can be
obtained by blast searching PDB. See :ref:`pca-blast` and :ref:`blastpdb`
examples.

Set reference
^^^^^^^^^^^^^

Now we set the reference chains that will be used for compared to the
structures in the ensemble and will be basis of the structural alignment.

.. ipython:: python

   # Parse reference structure
   reference_structure = parsePDB(ref_pdb + '.pdb', subset='calpha')
   # Get the reference chain from this structure
   reference_hierview = reference_structure.getHierView()
   reference_chains = [reference_hierview[chid] for chid in ref_chids]
   reference_chains

Chain A is the p66 domain, and chain B is the p51 domain of HIV-RT.
Let's take a quick look at that:

.. ipython:: python

   showProtein(reference_structure);
   @savefig ensemble_analysis_dimer_protein.png width=4in
   legend();

Prepare Ensemble
-------------------------------------------------------------------------------

We handle an ensembles of heterogeneous conformations using
:class:`.PDBEnsemble` objects, so let's instantiate one:

.. ipython:: python

   ensemble = PDBEnsemble(name)

We now combine the reference chains and set the reference coordinates
of the ensemble.

.. ipython:: python

   reference_chain = reference_chains[0] + reference_chains[1]
   ensemble.setAtoms(reference_chain)
   ensemble.setCoords(reference_chain.getCoords())

Coordinates of the reference structure are set as the coordinates of the
ensemble onto which other conformations will be superposed.


We can also start a log file using :func:`.startLogfile`.
Screen output will be save in this file, and can be
used to check if structures are added to the ensemble as expected.

.. ipython:: python

   startLogfile(name)

Let's also start a list to keep track of PDB files that are not added to the
ensemble:

.. ipython:: python

   unmapped = []

Now, we parse the PDB files one by one and add them to the ensemble:

.. ipython:: python

   for pdb in pdb_ids:
       # Parse the PDB file
       structure = parsePDB(pdb, subset='calpha', model=1)
       atommaps = []
       for reference_chain in reference_chains:
           # Map current PDB file to the reference chain
           mappings = mapOntoChain(structure, reference_chain,
                                   seqid=sequence_identity,
                                   coverage=sequence_coverage)
           if len(mappings) == 0:
               print 'Failed to map', pdb
               break
           atommaps.append(mappings[0][0])
           # Make sure all chains are mapped
       if len(atommaps) != len(reference_chains):
           unmapped.append(pdb)
           continue
       atommap = atommaps[0] + atommaps[1]
       ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
   ensemble
   ensemble.iterpose()
   saveEnsemble(ensemble)

We can now close the logfile using :func:`.closeLogfile`:

.. ipython:: python

   closeLogfile(name)

Let's check which structures, if any, are not mapped (added to the ensemble):

.. ipython:: python

   unmapped

We can write the aligned conformations into a PDB file as follows:

.. ipython:: python

   writePDB(name + '.pdb', ensemble)

This file can be used to visualize the aligned conformations in modeling
software.

This is a heterogeneous dataset, i.e. many structures had missing residues.
We want to make sure that we include residues in PCA analysis if they
are resolved in more than 94% of the time.

We can find out this using :func:`.calcOccupancies` function:

.. ipython:: python

   calcOccupancies(ensemble, normed=True).min()


This shows that some residues were resolved in only 24% of the dataset.
We trim the ensemble to contain residues resolved in more than 94% of the
ensemble:

.. ipython:: python

   ensemble = trimPDBEnsemble(ensemble, occupancy=0.94)

After trimmin, another round of iterative superposition may be useful:

.. ipython:: python

   ensemble.iterpose()
   saveEnsemble(ensemble)


Perform PCA
-------------------------------------------------------------------------------

Once the ensemble is ready, performing :class:`.PCA` is 3 easy steps:

.. ipython:: python

   pca = PCA(name)
   pca.buildCovariance(ensemble)
   pca.calcModes()

The calculated data can be saved as a compressed file using :func:`.saveModel`

.. ipython:: python

   saveModel(pca)

Plot results
-------------------------------------------------------------------------------


Let's plot RMSD to the average structure:

.. ipython:: python

   plot(calcRMSD(ensemble));
   xlabel('Conformation');
   ylabel('RMSD (A)');
   @savefig ensemble_analysis_dimer_rmsd.png width=4in
   title(ensemble);

Let's show a projection of the ensemble onto PC1 and PC2:

.. ipython:: python

   showProjection(ensemble, pca[:2]);
   @savefig ensemble_analysis_dimer_proj.png width=4in
   title(ensemble);


Only some of the ProDy plotting functions are shown here. A complete list
can be found in :ref:`dynamics` module.
