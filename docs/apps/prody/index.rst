.. _prody-apps:

ProDy Applications
==================

ProDy applications are command line programs that automates structure
processing and structural dynamics analysis:

.. toctree::
   :maxdepth: 1
   :glob:

   prody_*

Running :command:`prody` command will provide a description of applications::

  $ prody

.. literalinclude:: prody.txt



Detailed information on a specific application can be obtained
by typing the command and application names as :command:`prody anm -h`.

Running :command:`prody anm` application as follows will perform ANM
calculations for the p38 MAP kinase structure, and will write
eigenvalues/vectors in plain text and :ref:`nmd-format`::

  $ prody anm 1p38

In the above example, the default parameters (``cutoff=15.`` and ``gamma=1.``)
and all of the Cα atoms of the protein structure 1p38 are used.

In the example below, the *cutoff* distance is changed to 14 Å,
and the Cα atoms of residues with numbers smaller than 340 are used,
the output files are prefixed with :file:`p38_anm`::

  $ prody anm -c 14 -s "calpha resnum < 340" -p p38_anm 1p38

The output file :file:`p38_anm.nmd` can be visualized using NMWiz_.
