.. _trajectory:


Trajectory analysis
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example shows how to analyze a trajectory in DCD format. RMSD, RMSF,
radius of gyration, and distance will be calculated from trajectory frames.


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

Parse structure
-------------------------------------------------------------------------------

The PDB file provided with this example contains an X-ray structure which will
be useful in a number of places, so let's start with parsing this file first:

.. ipython:: python

   structure = parsePDB('trajectory_analysis_files/mdm2.pdb')
   repr(structure)

This function returned a :class:`.AtomGroup` instance that
stores all atomic data parsed from the PDB file.

Parse all frames
-------------------------------------------------------------------------------

Using :func:`.parseDCD` function all coordinate data in the DCD file can
be parsed at once. This function returns an :class:`.Ensemble` instance:

.. ipython:: python

   ensemble = parseDCD('trajectory_analysis_files/mdm2.dcd')
   repr(ensemble)

.. note:: When parsing large DCD files at once memory may become an issue.
   If the size of the DCD file is larger than half of the RAM in your machine,
   consider parsing DCD files frame-by-frame. See the following subsection for
   details.

Let's associate this ensemble with the *structure* we parsed from the PDB file:

.. ipython:: python

   ensemble.setAtoms(structure)
   ensemble.setCoords(structure)

This operation set the coordinates of the *structure* as the reference
coordinates of the *ensemble*. Now we can :meth:`.Ensemble.superpose`
the *ensemble* onto the coordinates of the *structure*.

.. ipython:: python

   ensemble.superpose()

Now, we can get calculate RMSDs and RMSFs as follows:

.. ipython:: python

   rmsd = ensemble.getRMSDs()
   rmsd[:10]
   rmsf = ensemble.getRMSFs()
   rmsf

Preceding calculations used all atoms in the structure. When we are interested
in a subset of atoms, let's say Cα atoms, we can make a selection before
performing calculations:

.. ipython:: python

   ensemble.setAtoms(structure.calpha)
   repr(ensemble)
   ensemble.superpose()

In this case, superposition was based on Cα atom coordinates.

.. ipython:: python

   rmsd = ensemble.getRMSDs()
   rmsd[:10]
   rmsf = ensemble.getRMSFs()
   rmsf


The :class:`.Ensemble` instance can also be used in :class:`.PCA`
calculations. See the examples in :ref:`pca` for more information.

Parse frames one-by-one
-------------------------------------------------------------------------------

.. ipython:: python

   dcd = DCDFile('trajectory_analysis_files/mdm2.dcd')
   repr(dcd)

.. ipython:: python

   structure = parsePDB('trajectory_analysis_files/mdm2.pdb')
   dcd.setCoords(structure)
   dcd.link(structure)

   dcd.nextIndex()
   frame = dcd.next()
   repr(frame)
   dcd.nextIndex()

.. ipython:: python

   frame.getRMSD()
   frame.superpose()
   frame.getRMSD()

   calcGyradius(frame)

We can perform these calculations for all frames in a for loop. Let's reset
*dcd* to return to the 0th frame:

.. ipython:: python

   dcd.reset()
   rgyr = zeros(len(dcd))
   rmsd = zeros(len(dcd))
   for i, frame in enumerate(dcd):
       rgyr[i] = calcGyradius(frame)
       frame.superpose()
       rmsd[i] = frame.getRMSD()
   rmsd[:10]
   rgyr[:10]

Handling multiple files
-------------------------------------------------------------------------------

:class:`.Trajectory` is designed for handling multiple trajectory files:

.. ipython:: python

   traj = Trajectory('trajectory_analysis_files/mdm2.dcd')
   repr(traj)
   traj.addFile('trajectory_analysis_files/mdm2sim2.dcd')
   repr(traj)

Instances of this class are also suitable for previous calculations:

.. ipython:: python

   structure = parsePDB('trajectory_analysis_files/mdm2.pdb')
   traj.link(structure)
   traj.setCoords(structure)
   rgyr = zeros(len(traj))
   rmsd = zeros(len(traj))
   for i, frame in enumerate(traj):
       rgyr[i] = calcGyradius( frame )
       frame.superpose()
       rmsd[i] = frame.getRMSD()
   rmsd[:10]
   rgyr[:10]
