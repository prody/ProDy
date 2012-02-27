.. _outputtraj:

*******************************************************************************
Trajectory output
*******************************************************************************

Synopsis
===============================================================================

This example shows how to output trajectories.


Input files
===============================================================================

Currently, ProDy supports only DCD format files. Two DCD trajectory files and 
corresponding PDB structure file is needed for this example.

Example input:
 
* :download:`MDM2 files </doctest/mdm2.tar.gz>` 


Load structure
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

The PDB file provided with this example contains an X-ray structure:

>>> mdm2 = parsePDB('mdm2.pdb')
>>> mdm2
<AtomGroup: mdm2 (1449 atoms)>

This function returned a :class:`~.AtomGroup` instance that stores all atomic 
data parsed from the PDB file.

Open trajectories
===============================================================================

:class:`~.Trajectory` is designed for handling multiple trajectory files:

>>> traj = Trajectory('mdm2.dcd')
>>> traj
<Trajectory: mdm2 (1 files; next 0 of 500 frames; 1449 atoms)>
>>> traj.addFile('mdm2sim2.dcd')
>>> traj 
<Trajectory: mdm2 (2 files; next 0 of 1000 frames; 1449 atoms)>

Now we link the trajectory (*traj*) with the atom group (*mdm2*): 

>>> traj.setAtoms(mdm2) 

.. note::
   Note that when a frame (coordinate set) is parsed from the trajectory file,
   coordinates of the atom group will be updated.


Output selected atoms
===============================================================================

You can write a trajectory in DCD format using :func:`~.writeDCD` function.
Let's select non-hydrogen protein atoms and write a merged trajectory for
MDM2:

>>> traj.select('noh')
<Selection: "noh" from mdm2 (706 atoms)>
>>> writeDCD('mdm2_merged_noh.dcd', traj)
'mdm2_merged_noh.dcd'

Parsing this file returns:

>>> DCDFile('mdm2_merged_noh.dcd')
<DCDFile: mdm2_merged_noh (next 0 of 1000 frames; 706 atoms)>


Output aligned frames
===============================================================================

You can write a trajectory in DCD format after aligning the frames.
Let's return to the first frame by resetting the trajectory:

>>> traj.reset()
>>> traj
<Trajectory: mdm2 (2 files; next 0 of 1000 frames; selected 706 of 1449 atoms)>

It is possible to write multiple DCD files at the same time.  We open two DCD 
files in write mode, one for all atoms, and another for backbone atoms:

>>> out = DCDFile('mdm2_aligned.dcd', 'w')
>>> out_bb = DCDFile('mdm2_bb_aligned.dcd', 'w')
>>> mdm2_bb = mdm2.backbone

Let's align and write frames one by one: 

>>> for frame in traj:
...     frame.superpose()
...     out.write(mdm2)
...     out_bb.write(mdm2_bb)

Let's open these files to show number of atoms in each:

>>> DCDFile('mdm2_aligned.dcd')
<DCDFile: mdm2_aligned (next 0 of 1000 frames; 1449 atoms)>
>>> DCDFile('mdm2_bb_aligned.dcd')
<DCDFile: mdm2_bb_aligned (next 0 of 1000 frames; 339 atoms)>

See Also
===============================================================================

See :ref:`trajectory`, :ref:`trajectory2`, and :ref:`atomsframes` for more 
usage examples and :ref:`eda` for essential dynamics analysis example. 

|questions|

|suggestions|
