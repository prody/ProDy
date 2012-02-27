.. _generate-conformers:

*******************************************************************************
Generate conformers along normal modes
*******************************************************************************

Synopsis
===============================================================================

This example shows how to generate conformers along normal modes.  Normal modes
may be from :class:`~.ANM`, :class:`~.PCA`, or any other user provided 
:class:`~.NMA` models. 

Input
-------------------------------------------------------------------------------

A PDB structure will be used to calculate ANM modes. 

Output
-------------------------------------------------------------------------------

Output is conformers along selected normal modes. Conformers can be save in 
PDB format.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Normal mode analysis
-------------------------------------------------------------------------------

Conformers can be generated along any set of normal modes. In this example,
we will calculate normal modes for unbound structure of p38 MAP kinase and
generate backbone trace conformations. 

>>> p38ca = parsePDB('1p38', subset='ca')
>>> anm = ANM('1p38')
>>> anm.buildHessian(p38ca)
>>> anm.calcModes()

Traverse a mode
-------------------------------------------------------------------------------

:func:`~.traverseMode` function generates conformations along a single normal
mode. Conformations are generated in both directions along the given mode.
*rmsd* argument is used to set the RMSD distance to the farthest conformation.

Let's generate 10 conformations along ANM mode 1:

>>> trajectory = traverseMode(anm[0], p38ca, n_steps=5, rmsd=2.0)
>>> trajectory 
<Ensemble: Conformations along Mode 1 from ANM 1p38 (11 conformations; 351 atoms)>
>>> calcRMSD(trajectory).round(2)
array([ 2. ,  1.6,  1.2,  0.8,  0.4,  0. ,  0.4,  0.8,  1.2,  1.6,  2. ])

Note that the resulting trajectory contains 11 coordinate set including
the initial coordinates. 

**Write PDB file**

The newly generated Cα trajectory can be written in PDB format as 
follows:

>>> p38traj = p38ca.copy()
>>> p38traj.delCoordset(0)
>>> p38traj.addCoordset( trajectory )
>>> writePDB('p38_mode1_trajectory.pdb', p38traj)
'p38_mode1_trajectory.pdb'

Sample along modes
-------------------------------------------------------------------------------

:func:`~.sampleModes` function generates conformations using random 
combinations of given modes. The nice thing about this function is that 
user can preset the average RMSD of the generated ensemble to the initial 
coordinates. 

Let's generate 40 p38 MAP kinase Cα conformers using slowest 3 ANM modes  
with average RMSD 1.5 Å to the initial coordinates:

>>> ensemble = sampleModes(anm[:3], p38ca, n_confs=40, rmsd=1.5)
>>> ensemble
<Ensemble: Conformations along 3 modes from ANM 1p38 (40 conformations; 351 atoms)>
>>> round( calcRMSD(ensemble).mean(), 2)
1.5

**Write PDB file**

The newly generated Cα ensemble can be written in PDB format as follows:

>>> p38ens = p38ca.copy()
>>> p38ens.delCoordset(0)
>>> p38ens.addCoordset( ensemble )
>>> writePDB('p38_mode123_ensemble.pdb', p38ens)
'p38_mode123_ensemble.pdb'


|questions|

|suggestions|
