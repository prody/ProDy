.. currentmodule:: prody.dynamics

.. _eda-md:

*******************************************************************************
EDA of MD Trajectories
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform essential dynamics analysis of molecular
dynamics (MD) trajectories.

Input
-------------------------------------------------------------------------------


User needs to provide coordinate sets from an MD trajectory as a 
:class:`numpy.ndarray`. 

In this example, we show how to extract coordinate sets
from a trajectory in :term:`DCD` file format using MDAnalysis package.
In this case, user needs to provide  :term:`PDB`, :term:`PSF`, and :term:`DCD` 
files from a completed simulation, but other formats will also do.  
   
Output
-------------------------------------------------------------------------------

A :class:`EDA` instance that stores covariance matrix and principal modes
that describes the essential dynamics of the system observed in the simulation. 
:class:`EDA` instance and principal modes (:class:`Mode`) can be used as input 
to functions in :mod:`~prody.dynamics` module for further analysis.

Notes
-------------------------------------------------------------------------------

MDAnalysis Python package is required for parsing coordinates from DCD
trajectory files. It can be used to parse trajectories in several other file 
formats as well. For more information interested user should consult the wiki 
pages of the package. See |mdanalysis| for more information.


ProDy Code
===============================================================================

We start by importing everything from ProDy and MDAnalysis packages::
  
  from prody import *
  import MDAnalysis

Prepare ensemble
-------------------------------------------------------------------------------

Instantiate a Universe for the simulated system::

  universe = MDAnalysis.Universe('protein.psf', 'protein.dcd')
  # Select atoms of interest
  universe_ca = universe.selectAtoms('name CA')
  # Get coordinates of CA atoms
  ca_coords = universe.dcd.timeseries(universe_ca, format='fac')


Instantiate a ProDy :class:`~prody.ensemble.Ensemble`::

  ensemble = Ensemble('MD snapshots')
  # Add all coordinate sets to ensemble
  ensemble.addCoordset(ca_coords)
  # Set reference coordinates 
  ensemble.setCoordinates(ca_coords[0])
  # Perform iterative sueprimposition
  ensemble.iterpose()

EDA calculations
-------------------------------------------------------------------------------

Instantiate :class:`EDA` and perform calculations::

  eda = EDA('EDA')
  eda.buildCovariance(ensemble)
  eda.calcModes()


Write NMD file
-------------------------------------------------------------------------------

We can write essential modes into an :term:`NMD` file for NMWiz.
For this we will need to parse the protein structure as well::

  prot = parsePDB('protein.pdb')
  prot_ca = prot.select('calpha')
  writeNMD('md_eda.nmd', eda[:3], prot_ca)

Print data
-------------------------------------------------------------------------------

Let's print fraction of variance for top raking 4 essential modes::

  for mode in eda[:4]:
      print mode.getFractOfVariance()


See Also
===============================================================================
   
User is referred to other examples in :ref:`pca-xray` for illustration of 
comparative analysis of theoretical and computational data.


|questions|

|suggestions|
