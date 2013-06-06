.. _trajectory2:

Trajectory Analysis II
===============================================================================

This example shows how to perform a more elaborate calculations simultaneously.
Radius of gyration, distance, psi angle calculated will be calculated
using trajectory frames.


Input files
-------------------------------------------------------------------------------

Two DCD trajectory files and a PDB structure file is provided for this example:

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

The PDB file provided with this example contains and X-ray structure which will
be useful in a number of places, so let's start with parsing this file first:

.. ipython:: python

   structure = parsePDB('trajectory_analysis_files/mdm2.pdb')
   structure

This function returned a :class:`.AtomGroup` instance that stores all atomic
data parsed from the PDB file.


Handling multiple files
-------------------------------------------------------------------------------

:class:`.Trajectory` is designed for handling multiple trajectory files:

.. ipython:: python

   traj = Trajectory('trajectory_analysis_files/mdm2.dcd')
   traj
   traj.addFile('trajectory_analysis_files/mdm2sim2.dcd')
   traj


Link trajectory to atoms
-------------------------------------------------------------------------------

Atoms can be linked to the trajectory as follows:

.. ipython:: python

   traj.link(structure)
   traj.setCoords(structure)

When an atom group is linked to a trajectory, frame coordinates parsed from
trajectory files will overwrite coordinates of the atom group. By making
atom selections, you can calculate and analyze different properties.


Setup for calculations
-------------------------------------------------------------------------------

Let's make atom selections for different types of calculations:

End-to-end distance
^^^^^^^^^^^^^^^^^^^

We select atoms from terminal residues and make an empty array whose length
equal to the number of frames:

.. ipython:: python

   nter = structure.select('name CA and resnum 25')
   cter = structure.select('name CA and resnum 109')
   e2e = zeros(traj.numFrames())

Radius of gyration
^^^^^^^^^^^^^^^^^^

We select atoms protein atoms this calculation and make an empty array:


.. ipython:: python

   protein = structure.select('noh and protein')
   rgyr = zeros(traj.numFrames())

A psi angle
^^^^^^^^^^^

We select a residue an make an empty array:

.. ipython:: python

   res30 = structure['PPP', 'P', 30]
   res30
   res30psi = zeros(traj.numFrames())


Perform calculations
-------------------------------------------------------------------------------

We perform all calculations simultaneously as follows:

.. ipython:: python

   for i, frame in enumerate(traj):
       e2e[i] = calcDistance(nter, cter)
       res30psi[i] = calcPsi(res30)
       rgyr[i] = calcGyradius(protein)

Let's print part of results:

.. ipython:: python

   e2e[:10]
   rgyr[:10]
   res30psi[:10]


Plot results
-------------------------------------------------------------------------------

End-to-end distance
^^^^^^^^^^^^^^^^^^^
.. ipython:: python

   plot(e2e);
   xlabel('Frame index');
   @savefig trajectory_analysis_end2end.png width=4in
   ylabel('End-to-end distance (A)');

Radius of gyration
^^^^^^^^^^^^^^^^^^

.. ipython:: python

   plot(rgyr);
   xlabel('Frame index');
   @savefig trajectory_analysis_gyradius.png width=4in
   ylabel('Radius of gyration (A)');

A psi angle
^^^^^^^^^^^

.. ipython:: python

   plot(res30psi);
   xlabel('Frame index');
   @savefig trajectory_analysis_res30psi.png width=4in
   ylabel('Residue 30 psi angle');
