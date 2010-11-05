.. module:: prody.dynamics

.. _egeda:

*******************************************************************************
EDA of MD Trajectories
*******************************************************************************

This example shows how to perform essential dynamics analysis of a molecular
dynamics trajector.

Following ProDy classes and functions are used in the example:

Classes:
  * :class:`PCA` (aliased as :class:`EDA`)
  * :class:`Ensemble`
Functions:
  * :func:`prody.proteins.parsePDB`

Also required:

  * PSF, PDB, and DCD files from a simulation
  * MDAnalysis Python package is required for parsing coordinates from DCD
    trajectory files.
   
.. note::
   |mdanalysis| is a very useful Python package for analyzing MD trajectories.
   It can be used to parse trajectories in several different file formats. 
   For more information interested user should consult the wiki pages of the 
   code.
  
.. literalinclude:: eda_md.py
