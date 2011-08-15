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

User needs to provide a trajectory in DCD file format and PDB file of the 
system.

Example input: 

* :download:`MDM2 reference structure </doctest/mdm2.pdb>` 
* :download:`MDM2 trajectory </doctest/mdm2.dcd>` 

Output
-------------------------------------------------------------------------------

A :class:`EDA` instance that stores covariance matrix and principal modes
that describes the essential dynamics of the system observed in the simulation. 
:class:`EDA` instance and principal modes (:class:`Mode`) can be used as input 
to functions in :mod:`~prody.dynamics` module for further analysis.

Notes
-------------------------------------------------------------------------------

Analysis of trajectories in some other formats can be performed with the help
of MDAnalysis Python package. Interested user should consult |mdanalysis| for 
more information.


ProDy Code
===============================================================================

We start by importing everything from ProDy:
  
>>> from prody import *
   
Parse reference structure
-------------------------------------------------------------------------------

The PDB file provided with this example contains and X-ray structure which will 
be useful in a number of places, so let's start with parsing this file first:

>>> structure = parsePDB('mdm2.pdb')
>>> structure
<AtomGroup: mdm2 (1449 atoms; 1 coordinate sets, active set index: 0)>

This function returned a :class:`~prody.atomic.AtomGroup` instance that
stores all atomic data parsed from the PDB file.

EDA calculations
-------------------------------------------------------------------------------

Essential dynamics analysis (EDA or PCA) of a trajectory can be performed in 
two ways. 

**Small trajectory files**

If you are analyzing a small trajectory, you can use an 
:class:`Ensemble` instance obtained by parsing the trajectory at once using
:func:`parseDCD`:

>>> ensemble = parseDCD('mdm2.dcd')
>>> ensemble.setAtomGroup( structure )
>>> ensemble.select('calpha')
<Selection: "calpha" from mdm2 (85 atoms; 1 coordinate sets, active set index: 0)>
>>> #ensemble.superpose()
>>> ensemble.iterpose()
>>> eda_ensemble = EDA('MDM2 Ensemble')
>>> eda_ensemble.buildCovariance( ensemble )
>>> eda_ensemble.calcModes()
>>> eda_ensemble
<EDA: MDM2 Ensemble (20 modes, 85 atoms)>

**Large trajectory files**

If you are analyzing a large trajectory, you can pass the trajectory instance
to the :meth:`~dynamics.PCA.buildCovariance` method as follows:

>>> dcd = DCDFile('mdm2.dcd')
>>> dcd.setAtomGroup( structure )
>>> dcd.select('calpha')
<Selection: "calpha" from mdm2 (85 atoms; 1 coordinate sets, active set index: 0)>
>>> eda_trajectory = EDA('MDM2 Trajectory')
>>> eda_trajectory.buildCovariance( dcd )
>>> eda_trajectory.calcModes()
>>> eda_trajectory
<EDA: MDM2 Trajectory (20 modes, 85 atoms)>

**Compare two methods**

>>> printOverlapTable(eda_ensemble[:3], eda_trajectory[:3])
Overlap Table
                       EDA MDM2 Trajectory
                         #1     #2     #3
EDA MDM2 Ensemble #1   +1.00  +0.00  -0.00
EDA MDM2 Ensemble #2   -0.00  +1.00  +0.00
EDA MDM2 Ensemble #3   +0.00  -0.00  +1.00
<BLANKLINE>

Overlap values of +1 along the diagonal of the table shows that top ranking
3 essential (principal) modes are the same.


Write NMD file
-------------------------------------------------------------------------------

We can write essential modes into an :term:`NMD` file for NMWiz.
For this we will need to parse the protein structure as well::

>>> writeNMD('mdm2_eda.nmd', eda_trajectory[:3], structure.select('calpha'))
'mdm2_eda.nmd'

Print data
-------------------------------------------------------------------------------

Let's print fraction of variance for top raking 4 essential modes::

>>> for mode in eda_trajectory[:4]:
...     print mode.getFractOfVariance().round(2)
0.26
0.11
0.08
0.06

See Also
===============================================================================
   
See other examples in :ref:`pca-xray` for illustration of 
comparative analysis of theoretical and computational data.

See also :ref:`trajectory` for more analysis examples. 


|questions|

|suggestions|
