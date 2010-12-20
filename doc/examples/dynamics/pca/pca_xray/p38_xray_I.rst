.. currentmodule:: prody.dynamics

.. _p38-xray-calculations:

*******************************************************************************
p38 X-ray Ensemble - Part I: Calculations
*******************************************************************************

This is the first part of a lengthy ProDy example. The aim is to repeat the 
calculations for p38 MAP kinase (MAPK) that was published in [AB09]_. 
In this part, we perform the calculations, using the same p38 MAPK dataset.


We start with importing ProDy into the current namespace, which may be
an interactive Python session::

  from prody import *

Gather dataset
===============================================================================

We need a list of PDB identifiers that we want to include in our analysis.
We start with a list of identifiers which is slightly different from that we 
have listed in supporting material of our paper [AB09]_. 
Since the paper was published, ProteinDataBank had refined some
structures. Below is an updated list that contains the same structures::
  
  pdbids = ['1A9U', '1BL6', '1BL7', '1BMK', '1DI9', '1IAN', '1KV1', '1KV2', '1LEW', '1LEZ', 
            '1M7Q', '1OUK', '1OUY', '1OVE', '1OZ1', '1P38', '1R39', '1R3C', '1W7H', '1W82', 
            '1W83', '1W84', '1WBN', '1WBO', '1WBS', '1WBT', '1WBV', '1WBW', '1WFC', '1YQJ', 
            '1YW2', '1YWR', '1ZYJ', '1ZZ2', '1ZZL', '2BAJ', '2BAK', '2BAL', '2BAQ', '2EWA', 
            '2FSL', '2FSM', '2FSO', '2FST', '2GFS', '2GHL', '2GHM', '2GTM', '2GTN', '2I0H', 
            '2NPQ', '2OKR', '2OZA', '3HVC', '3MH0', '3MH3', '3MH2', '2PUU', '3MGY', '3MH1', 
            '2QD9', '2RG5', '2RG6', '2ZAZ', '2ZB0', '2ZB1', '3BV2', '3BV3', '3BX5', '3C5U', 
            '3L8X', '3CTQ', '3D7Z', '3D83', '2ONL']

  pdbfiles = fetchPDB(pdbids, folder='p38')
  
After we set the list of PDB identifiers, we fetched them using 
:func:`~prody.proteins.fetchPDB` function. ``pdbfile`` variable now 
contains list of fetched PDB filenames. 

Set reference chain
===============================================================================

The next step is setting one of these p38 structures as the reference
structure. We will use 1p38 chain A as the reference. Note that we won't use
all resolved residues in this structure. We are going to select residues
that are resolved at least in 90% of the dataset. 

:: 

  ref_structure = parsePDB('p38/1p38.pdb.gz')
  ref_structure = ref_structure.copy('resnum 5 to 31 36 to 114 122 to 169 185 to 351 and calpha')

  # Rename the reference structure
  ref_structure.setName('p38 reference')

  # Select chain A of the reference structure
  ref_chain = ref_structure.getHierView().getChain('A')

We used :func:`~prody.proteins.parsePDB` function to parse a PDB file.
This returned a :class:`~prody.atomic.AtomGroup` instance. We made a copy
of alpha carbons of select residues for analysis.   

|more| See :ref:`selections` for reference on making selections.

Prepare ensemble
===============================================================================

X-ray structural ensembles are not uniform, i.e. different structures
have different sets of unresolved residues. Hence, it is not straightforward
to analyzed them as it would be for NMR models. 

ProDy has special functions and classes for facilitating efficient analysis
of PDB X-ray data. In this example we use :func:`~prody.compare.mapAtomsToChain` 
which returns a instances :class:`~prody.atomic.AtomMap`.

|more| See :ref:`atommaps` for more details on how they work.   

::

  # Start a logfile to save screen output 
  ProDyStartLogfile('p38_pca_anm_calculations') 

  # Start an ensemble instance
  ensemble = Ensemble('p38 X-ray')
  
  # Set the reference coordinates
  ensemble.setCoordinates(ref_chain) 
      
  # For each PDB file we find matching chain and add it to the ensemble
  for pdbfile in pdbfiles:
      # Parse next PDB file. (only alpha carbons, since it's faster)
      current_structure = parsePDB(pdbfile, subset='calpha')
      
      # Get mapping to the reference chain
      current_mapping = mapAtomsToChain(current_structure, ref_chain)
      
      # Get the atom mapping
      current_atommap = current_mapping[0][0]
      
      # Rename the atom mapping
      #current_atommap.setName(current_structure.getName())

      # Add the atommap (mapped coordinates) to the ensemble
      # Note that some structures do not completely map (missing residues)
      # so we pass weights (1 for mapped atoms, 0 for unmapped atoms)
      ensemble.addCoordset(current_atommap, weights=current_atommap.getMappedFlags())    

  # Perform an iterative superimposition
  ensemble.iterpose()

  # Close the logfile (contents of the file shows how chains were paired/mapped)
  ProDyCloseLogfile('p38_pca_anm_calculations')

Save coordinates
===============================================================================

We used :class:`~prody.ensemble.Ensemble` to store coordinates of the X-ray 
structures. Its instances do not store any other atomic data. So, if we want
to write aligned coordinates into a file, we need to pass coordinates
to an :class:`~prody.atomic.AtomGroup` instance::

  xray_coords = ref_structure.copy()
  xray_coords.delCoordset(0) # Delete existing coordinate set
  xray_coords.addCoordset( ensemble.getCoordsets() )
  writePDB('p38_xray_coors.pdb', xray_coords)

We used :func:`~prody.proteins.writePDB` function to save coordinates.

Perform PCA
===============================================================================

Now we have the coordinates ready, it is straightforward to perform PCA
calculations:: 

  pca = PCA('p38 xray')           # Instantiate a PCA instance
  pca.buildCovariance(ensemble)   # Build covariance for the ensemble
  pca.calcModes()                 # Calculate modes (20 of the by default)

.. _p38-xray-calculations-anm:

Perform ANM
===============================================================================

So it is to perform ANM calculations:: 

  anm = ANM('1p38')             # Instantiate a ANM instance
  anm.buildHessian(ref_chain)   # Build Hessian for the reference chain  
  anm.calcModes()               # Calculate slowest non-trivial 20 modes 

Save your work
===============================================================================

Now, we will save calculated data and use it in the next part of the example.

If you are in an interactive Python session, and wish to continue without
leaving your session, you do not need to save the data. Saving data is useful
if you want to use it in another session or at a later time.

::

  saveModel(pca)
  saveModel(anm)
  saveEnsemble(ensemble)
  writePDB('p38_ref_chain.pdb', ref_chain)

We used :func:`saveModel` to save calculated data. In :ref:`p38-II`, we will
used :func:`loadModel` to load them.

|more| This example continues in :ref:`p38-xray-analysis` 
