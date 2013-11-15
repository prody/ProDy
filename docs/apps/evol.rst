.. _evol-apps:

Evol Applications
=================

Evol applications are command line programs that automate retrieval,
refinement, and analysis of multiple sequence alignments:

.. toctree::
   :maxdepth: 1
   :glob:

   evol_*

Running :command:`evol` command will provide a description of applications::

  $ evol

.. literalinclude:: evol.txt


Detailed information on a specific application can be obtained
by typing the command and application names as :command:`evol search -h`.

Running :command:`prody search` application as follows will search Pfam
database for protein families that match the proteins in PDB structure 2w5i::

  $ evol search 2w5i
