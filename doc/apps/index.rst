.. _prody-apps:

*******************************************************************************
ProDy Applications
*******************************************************************************

ProDy applications are command line programs that automates structural data
analysis, theoretical calculations, and plotting resulting data.  

Running :command:`prody` command will provide a description of applications::
  
  $ prody
  
.. literalinclude:: prody.txt


See usage details of these applications below:

.. toctree::
   :maxdepth: 1
   :glob:

   prody_*


Application Setup
===============================================================================

ProDy comes with several handy applications.  On Linux, when installing ProDy 
from source, application scripts are placed into a default folder that is 
included in :envvar:`PATH` environment variable, e.g. :file:`/usr/local/bin/`.  
On Windows, installer places the scripts into the :file:`Scripts` folder under
Python distribution folder, e.g. :file:`C:\\Python27\\Scripts`.  You may need 
to add this path to :envvar:`PATH` environment variable yourself. 


Usage Example
===============================================================================

Detailed information on a specific application can be obtained 
by typing the command and application names as :command:`prody anm -h`.

Running :command:`prody anm` application as follows will perform ANM 
calculations for the p38 MAP kinase structure, and will write 
eigenvalues/vectors in plain text and :term:`NMD` formats::

  $ prody anm 1p38
  
In the above example, the default parameters (``cutoff=15.`` and ``gamma=1.``)
and all of the Cα atoms of the protein structure 1p38 are used.

In the example below, the *cutoff* distance is changed to 14 Å, 
and the Cα atoms of residues with numbers smaller than 340 are used, 
the output files are prefixed with :file:`p38_anm`::

  $ prody anm -c 14 -s "calpha resnum < 340" -p p38_anm 1p38

The output file :file:`p38_anm.nmd` can be visualized using NMWiz (|nmwiz|). 
