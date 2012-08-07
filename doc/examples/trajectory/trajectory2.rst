.. _trajectory2:

*******************************************************************************
Trajectory analysis II
*******************************************************************************

Synopsis
===============================================================================

This example shows how to perform a more elaborate calculations simultaneously. 
Radius of gyration, distance, psi angle calculated will be calculated 
using trajectory frames.

 
Input files
===============================================================================

Two DCD trajectory files and a PDB structure file is provided for this example:

* :download:`MDM2 files </doctest/mdm2.tar.gz>` 

Parse structure 
===============================================================================

We start by importing everything from the ProDy package and also the NumPy
package:

>>> from prody import *
>>> import numpy as np

The PDB file provided with this example contains and X-ray structure which will 
be useful in a number of places, so let's start with parsing this file first:

>>> structure = parsePDB('mdm2.pdb')
>>> structure
<AtomGroup: mdm2 (1449 atoms)>

This function returned a :class:`.AtomGroup` instance that stores all atomic 
data parsed from the PDB file.

Handling multiple files
===============================================================================

:class:`.Trajectory` is designed for handling multiple trajectory files:

>>> traj = Trajectory('mdm2.dcd')
>>> traj
<Trajectory: mdm2 (next 0 of 500 frames; 1449 atoms)>
>>> traj.addFile('mdm2sim2.dcd')
>>> traj 
<Trajectory: mdm2 (2 files; next 0 of 1000 frames; 1449 atoms)>

Link trajectory to atoms
===============================================================================

Atoms can be linked to the trajectory as follows:

>>> traj.link(structure)
>>> traj.setCoords(structure)

When an atom group is linked to a trajectory, frame coordinates parsed from
trajectory files will overwrite coordinates of the atom group. By making
atom selections, you can calculate and analyze different properties. 


Setup for calculations
===============================================================================

Let's make atom selections for different types of calculations:

**End-to-end calculation**

We select atoms from terminal residues and make an empty array whose length
equal to the number of frames: 

>>> nter = structure.select('name CA and resnum 25')
>>> cter = structure.select('name CA and resnum 109')
>>> e2e = np.zeros(traj.numFrames())

**Protein radius of gyration**

We select atoms protein atoms this calculation and make an empty array: 


>>> protein = structure.select('noh and protein') 
>>> rgyr = np.zeros(traj.numFrames())

**Psi angle of a residue**

We select a residue an make an empty array:

>>> res30 = structure['PPP', 'P', 30]
>>> res30
<Residue: PRO 30 from Chain P from mdm2 (14 atoms)>
>>> res30psi = np.zeros(traj.numFrames())

Perform calculations
===============================================================================

We perform all calculations simultaneously as follows:

>>> for i, frame in enumerate(traj):
...     e2e[i] = calcDistance(nter, cter)
...     res30psi[i] = calcPsi(res30)
...     rgyr[i] = calcGyradius(protein)

Let's print results:

>>> print e2e.round(2) # doctest: +ELLIPSIS
[ 11.79  14.13  15.66  14.52  16.46  17.21  16.45  14.29  11.6   12.66
  ...
  12.59  11.2   11.26  11.89  11.36]
>>> print rgyr.round(2) # doctest: +ELLIPSIS
[ 12.86  12.98  12.83  12.92  12.87  12.92  12.76  12.86  12.82  12.76
  ...
  12.91  12.88  12.73  12.85  12.88  12.86  12.9   12.99  12.8   12.84
  12.87  12.84]

>>> print res30psi.round(2) # doctest: +ELLIPSIS
[ 149.81  170.66  139.94  156.37  139.49  151.11  147.68  151.82  143.42
  ...
  159.33  126.08  125.54  139.35  133.5   129.46  132.58  147.61  145.03
  151.92]

See Also
===============================================================================

See :ref:`trajectory`, :ref:`outputtraj`, and :ref:`frame` for more usage
examples and :ref:`eda` for essential dynamics analysis example. 

|questions|

|suggestions|
