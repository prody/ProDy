.. _outputtraj:

Trajectory output
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example shows how to output trajectories.


Input files
-------------------------------------------------------------------------------

Currently, ProDy supports only DCD format files. Two DCD trajectory files and
corresponding PDB structure file is needed for this example.

Example input:

* :download:`MDM2 structure <trajectory_analysis_files/mdm2.pdb>`
* :download:`MDM2 trajectory I <trajectory_analysis_files/mdm2.dcd>`
* :download:`MDM2 trajectory II <trajectory_analysis_files/mdm2sim2.dcd>`


Setup environment
-------------------------------------------------------------------------------

We start by importing everything from ProDy:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Load structure
-------------------------------------------------------------------------------

The PDB file provided with this example contains an X-ray structure:

.. ipython:: python

   mdm2 = parsePDB('trajectory_analysis_files/mdm2.pdb')
   repr(mdm2)

This function returned a :class:`.AtomGroup` instance that stores all atomic
data parsed from the PDB file.


Open trajectories
-------------------------------------------------------------------------------

:class:`.Trajectory` is designed for handling multiple trajectory files:

.. ipython:: python

   traj = Trajectory('trajectory_analysis_files/mdm2.dcd')
   traj
   traj.addFile('trajectory_analysis_files/mdm2sim2.dcd')
   traj

Now we link the trajectory (*traj*) with the atom group (*mdm2*):

.. ipython:: python

   traj.link(mdm2)

.. note::
   Note that when a frame (coordinate set) is parsed from the trajectory file,
   coordinates of the atom group will be updated.


Output selected atoms
-------------------------------------------------------------------------------

You can write a trajectory in DCD format using :func:`.writeDCD` function.
Let's select non-hydrogen protein atoms and write a merged trajectory for
MDM2:

.. ipython:: python

   traj.setAtoms(mdm2.noh)
   traj
   writeDCD('mdm2_merged_noh.dcd', traj)

Parsing this file returns:

.. ipython:: python

   DCDFile('mdm2_merged_noh.dcd')


Output aligned frames
-------------------------------------------------------------------------------

You can write a trajectory in DCD format after aligning the frames.
Let's return to the first frame by resetting the trajectory:

.. ipython:: python

   traj.reset()
   traj

It is possible to write multiple DCD files at the same time.  We open two DCD
files in write mode, one for all atoms, and another for backbone atoms:

.. ipython:: python

   out = DCDFile('mdm2_aligned.dcd', 'w')
   out_bb = DCDFile('mdm2_bb_aligned.dcd', 'w')
   mdm2_bb = mdm2.backbone

Let's align and write frames one by one:

.. ipython:: python

   for frame in traj:
       frame.superpose()
       out.write(mdm2)
       out_bb.write(mdm2_bb)

Let's open these files to show number of atoms in each:

.. ipython:: python

   DCDFile('mdm2_aligned.dcd')
   DCDFile('mdm2_bb_aligned.dcd')