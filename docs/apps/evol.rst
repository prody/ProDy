.. _evol-apps:

Evol Applications
=================

Evol applications are command line programs that automate retrieval,
refinement, and analysis of multiple sequence alignments.  See

Running :command:`evol` command will provide a description of applications::

  $ evol

.. literalinclude:: evol.txt

See usage details of these applications below:

.. toctree::
   :maxdepth: 1
   :glob:

   evol_*


Application setup
-----------------

Evol applications come with ProDy package.  On Linux, when installing ProDy
from source, application scripts are placed into a default folder that is
included in :envvar:`PATH` environment variable, e.g. :file:`/usr/local/bin/`.
On Windows, installer places the scripts into the :file:`Scripts` folder under
Python distribution folder, e.g. :file:`C:\\Python27\\Scripts`.  You may need
to add this path to :envvar:`PATH` environment variable yourself.


Usage example
-------------

Detailed information on a specific application can be obtained
by typing the command and application names as :command:`evol search -h`.

Running :command:`prody search` application as follows will search Pfam
database for protein families that match the proteins in PDB structure 2w5i::

  $ evol search 2w5i
