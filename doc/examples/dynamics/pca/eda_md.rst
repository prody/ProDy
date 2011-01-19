.. currentmodule:: prody.dynamics

.. _eda-md:

*******************************************************************************
EDA of MD Trajectories
*******************************************************************************

This example shows how to perform essential dynamics analysis of a molecular
dynamics trajector.

Following ProDy classes and functions are used in the example:

Classes:
  * :class:`EDA`
  * :class:`~prody.ensemble.Ensemble`
Functions:
  * :func:`~prody.proteins.parsePDB`

Also required:

  * PSF, PDB, and DCD files from a simulation
  * MDAnalysis Python package is required for parsing coordinates from DCD
    trajectory files.
   
.. note::
   |mdanalysis| is a very useful Python package for analyzing MD trajectories.
   It can be used to parse trajectories in several different file formats. 
   For more information interested user should consult the wiki pages of the 
   code.

::
  
  from prody import *

  # Import MDAnalysis
  import MDAnalysis

  # Instantiate a Universe for the simulated system
  universe = MDAnalysis.Universe('protein.psf', 'protein.dcd')

  # Select atoms of interest
  universe_ca = universe.selectAtoms('name CA')

  # Get coordinates of CA atoms
  ca_coords = universe.dcd.timeseries(universe_ca, format='fac')

  # Instantiate an emseble
  ensemble = Ensemble('MD snapshots')
  # Add all coordinate sets to ensemble
  ensemble.addCoordset(ca_coords)
  # Set reference coordinates 
  ensemble.setCoordinates(ca_coords[0])
  # Perform iterative sueprimposition
  ensemble.iterpose()

  # Instantiate EDA and perform calculations
  eda = EDA('EDA')
  eda.buildCovariance(ensemble)
  eda.calcModes()



  # Write essential modes into an NMD file for NMWiz
  # for this we will need to parse the protein structure as well
  prot = parsePDB('protein.pdb')
  prot_ca = prot.select('calpha')
  writeNMD('md_eda.nmd', eda[:3], prot_ca)

  # Let's print fraction of variance for top raking 4 essential modes
  for mode in eda[:4]:
      print mode.getFractOfVariance()
