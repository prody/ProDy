.. _eda:


Essential Dynamics Analysis
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example shows how to perform essential dynamics analysis of molecular
dynamics (MD) trajectories.  A :class:`.EDA` instance that stores covariance
matrix and principal modes that describes the essential dynamics of the system
observed in the simulation will be built.  :class:`.EDA` and principal modes
(:class:`.Mode`) can be used as input to functions in :mod:`.dynamics` module
for further analysis.


User needs to provide trajectory in DCD file format and PDB file of the system.

Example input:

  * :download:`MDM2 files (ZIP) <trajectory_analysis_files.zip>`
  * :download:`MDM2 files (TGZ) <trajectory_analysis_files.tgz>`


Setup environment
-------------------------------------------------------------------------------

We start by importing everything from ProDy:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Parse reference structure
-------------------------------------------------------------------------------


The PDB file provided with this example contains and X-ray structure which will
be useful in a number of places, so let's start with parsing this file first:

.. ipython:: python

   structure = parsePDB('trajectory_analysis_files/mdm2.pdb')
   structure

This function returned a :class:`.AtomGroup` instance that
stores all atomic data parsed from the PDB file.

EDA calculations
-------------------------------------------------------------------------------

Essential dynamics analysis (EDA or PCA) of a trajectory can be performed in
two ways.

Small files
^^^^^^^^^^^

If you are analyzing a small trajectory, you can use an :class:`.Ensemble`
instance obtained by parsing the trajectory at once using :func:`.parseDCD`:

.. ipython:: python

   ensemble = parseDCD('trajectory_analysis_files/mdm2.dcd')
   ensemble.setCoords(structure)
   ensemble.setAtoms(structure.calpha)
   ensemble
   ensemble.superpose()
   eda_ensemble = EDA('MDM2 Ensemble')
   eda_ensemble.buildCovariance( ensemble )
   eda_ensemble.calcModes()
   eda_ensemble

Large files
^^^^^^^^^^^

If you are analyzing a large trajectory, you can pass the trajectory instance
to the :meth:`.PCA.buildCovariance` method as follows:

.. ipython:: python

   dcd = DCDFile('trajectory_analysis_files/mdm2.dcd')
   dcd.link(structure)
   dcd.setAtoms(structure.calpha)
   dcd

   eda_trajectory = EDA('MDM2 Trajectory')
   eda_trajectory.buildCovariance( dcd )
   eda_trajectory.calcModes()
   eda_trajectory

Comparison
^^^^^^^^^^

.. ipython:: python

   printOverlapTable(eda_ensemble[:3], eda_trajectory[:3])

Overlap values of +1 along the diagonal of the table shows that top ranking
3 essential (principal) modes are the same.

Multiple files
-------------------------------------------------------------------------------

It is also possible to analyze multiple trajectory files without concatenating
them. In this case we will use data from two independent simulations

.. ipython:: python

   trajectory = Trajectory('trajectory_analysis_files/mdm2.dcd')
   trajectory.addFile('trajectory_analysis_files/mdm2sim2.dcd')
   trajectory

   trajectory.link(structure)
   trajectory.setCoords(structure)
   trajectory.setAtoms(structure.calpha)
   trajectory

   eda = EDA('mdm2')
   eda.buildCovariance( trajectory )
   eda.calcModes()
   eda

Save your work
^^^^^^^^^^^^^^

You can save your work using ProDy function :func:`.saveModel`. This will
allow you to avoid repeating calculations when you return to your work later:

.. ipython:: python

   saveModel(eda)

:func:`.loadModel` function can be used to load this object without any loss.

Analysis
-------------------------------------------------------------------------------

Let's print fraction of variance for top raking 4 essential modes:

.. ipython:: python

   for mode in eda_trajectory[:4]:
       print calcFractVariance(mode).round(2)

You can find more analysis functions in :ref:`dynamics`.

Plotting
-------------------------------------------------------------------------------

Now, let's project the trajectories onto top three essential modes:

.. ipython:: python

   mdm2ca_sim1 = trajectory[:500]
   mdm2ca_sim1.superpose()
   mdm2ca_sim2 = trajectory[500:]
   mdm2ca_sim2.superpose()

   # We project independent trajectories in different color
   showProjection(mdm2ca_sim1, eda[:3], color='red', marker='.');
   showProjection(mdm2ca_sim2, eda[:3], color='blue', marker='.');
   # Now let's mark the beginning of the trajectory with a circle
   showProjection(mdm2ca_sim1[0], eda[:3], color='red', marker='o', ms=12);
   showProjection(mdm2ca_sim2[0], eda[:3], color='blue', marker='o', ms=12);
   # Now let's mark the end of the trajectory with a square
   showProjection(mdm2ca_sim1[-1], eda[:3], color='red', marker='s', ms=12);
   @savefig trajectory_analysis_eda_projection.png width=4in
   showProjection(mdm2ca_sim2[-1], eda[:3], color='blue', marker='s', ms=12);

You can find more plotting functions in :ref:`dynamics` and :ref:`measure`
modules.

Visualization
-------------------------------------------------------------------------------

The above projection is shown for illustration. Interpreting the essential
modes and projection of snapshots onto them is case dependent. One should know
what kind of motion the top essential modes describe. You can use :ref:`nmwiz`
for visualizing essential mode shapes and fluctuations along these modes.

We can write essential modes in :ref:`nmd-format` for NMWiz as follows:

.. ipython:: python

   writeNMD('mdm2_eda.nmd', eda[:3], structure.select('calpha'))

