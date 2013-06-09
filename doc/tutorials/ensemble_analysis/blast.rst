.. _pca-blast:

Homologous Proteins
===============================================================================

This example shows how to perform PCA of a structural dataset obtained by blast
searching PDB. The protein of interest is :wiki:`cytochrome c` (cyt *c*).
Dataset will contain structures sharing 44% or more sequence identity with
human *cyt c*, i.e. its homologs and/or orthologs.

A :class:`.PCA` instance that stores covariance matrix and principal modes that
describe the dominant changes in the dataset will be obtained. :class:`.PCA`
instance and principal modes (:class:`.Mode`) can be used as input to functions
in :mod:`.dynamics` module for further analysis.

Input is amino acid sequence of the protein, a reference PDB identifier,
and some parameters.

Setup
-------------------------------------------------------------------------------

Import ProDy and matplotlib into the current namespace.


.. ipython:: python

   from prody import *
   from pylab import *
   ion()



Name of the protein (a name without a white space is preferred)

.. ipython:: python

   name = 'cyt_c'
   sequence = '''GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEE
   TLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE'''
   ref_pdb = '1hrc'

Optionally, a list of PDB files to be excluded from analysis can be provided.
In this case dimeric Cyt c structures are excluded from the analysis. To use
all PDB hits, provide an empty list.

.. ipython:: python

   exclude = ['3nbt', '3nbs']

Parameters
-------------------------------------------------------------------------------



.. ipython:: python

   # Minimum sequence identity of hits
   seqid = 44
   # Reference chain identifier
   ref_chid = 'A'

.. ipython:: python

   # Selection string ("all" can be used if all of the chain is to be analyzed)
   selstr = 'resnum 1 to 103'

Blast and download
-------------------------------------------------------------------------------

The results are displayed for following list of structures:



List of PDB structures can be updated using :func:`.blastPDB`
as follows:

.. ipython:: python

   blast_record = blastPDB(sequence)
   pdb_hits = []
   for key, item in blast_record.getHits(seqid).iteritems():
       pdb_hits.append((key, item['chain_id']))

Let's fetch PDB files and see how many there are:

.. ipython:: python

   pdb_files = fetchPDB(*[pdb for pdb, ch in pdb_hits], compressed=False)
   len(pdb_files)


Set reference
-------------------------------------------------------------------------------

We first parse the reference structure. Note that we parse only Cα atoms from
chain A. The analysis will be performed for a single chain (monomeric) protein.
For analysis of a dimeric protein see :ref:`pca-dimer`

.. ipython:: python

   reference_structure = parsePDB(ref_pdb, subset='ca', chain=ref_chid)
   # Get the reference chain from this structure
   reference_hierview = reference_structure.getHierView()
   reference_chain = reference_hierview[ref_chid]

Prepare ensemble
-------------------------------------------------------------------------------

.. ipython:: python

   # Start a log file
   startLogfile('pca_blast')
   # Instantiate a PDB ensemble
   ensemble = PDBEnsemble(name)
   # Set ensemble atoms
   ensemble.setAtoms(reference_chain)
   # Set reference coordinates
   ensemble.setCoords(reference_chain.getCoords())

.. ipython:: python

   for (pdb_id, chain_id), pdb_file in zip(pdb_hits, pdb_files):
       if pdb_id in exclude:
           continue
       structure = parsePDB(pdb_file, subset='calpha', chain=chain_id)
       if structure is None:
           plog('Failed to parse ' + pdb_file)
           continue
       mappings = mapOntoChain(structure, reference_chain, seqid=seqid)
       if len(mappings) == 0:
           plog('Failed to map', pdb_id)
           continue
       atommap = mappings[0][0]
       ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
   ensemble.iterpose()
   saveEnsemble(ensemble)


Let's check how many conformations are extracted from PDB files:

.. ipython:: python

   len(ensemble)

Note that number of conformations is larger than the number of PDB structures
we retrieved. This is because some of the PDB files contained NMR structures
with multiple models. Each model in NMR structures are added to the ensemble
as individual conformations.

Write aligned conformations into a PDB file as follows:

.. ipython:: python

   writePDB(name+'.pdb', ensemble)


This file can be used to visualize the aligned conformations in a modeling
software.



Align PDB files
-------------------------------------------------------------------------------

:func:`.alignPDBEnsemble` function can be used to align all PDB structures used
in the analysis, e.g. ``alignPDBEnsemble(ensemble)``.  Outputted files will
contain intact structures and can be used for visualization purposes in other
software.  In this case, we will align only select PDB files:

.. ipython:: python

   conf1_alinged = alignPDBEnsemble(ensemble[0])
   conf2_alinged = alignPDBEnsemble(ensemble[1])


Let's take a quick look at the aligned structures:

.. ipython:: python


   showProtein(parsePDB(conf1_alinged), parsePDB(conf2_alinged));
   @savefig ensemble_analysis_blast_aligned.png width=4in
   legend();


Perform PCA
-------------------------------------------------------------------------------

Once the ensemble is ready, performing PCA is 3 easy steps:

.. ipython:: python

   # Instantiate a PCA
   pca = PCA(name)
   # Build covariance matrix
   pca.buildCovariance(ensemble)
   # Calculate modes
   pca.calcModes()

The calculated data can be saved as a compressed file using :func:`.saveModel`
function:

.. ipython:: python

   saveModel(pca)


Plot results
-------------------------------------------------------------------------------


Let's plot RMSDs of all conformations from the average conformation:


.. ipython:: python

   rmsd = calcRMSD(ensemble)
   plot(rmsd);
   xlabel('Conformation index');
   @savefig ensemble_analysis_blast_rmsd.png width=4in
   ylabel('RMSD (A)');


Let's show a projection of the ensemble onto PC1 and PC2:

.. ipython:: python

   @savefig ensemble_analysis_blast_projection.png width=4in
   showProjection(ensemble, pca[:2]);
