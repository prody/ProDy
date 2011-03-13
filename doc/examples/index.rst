.. _examples:

.. currentmodule:: prody

*******************************************************************************
Examples
*******************************************************************************

How to Use
===============================================================================

Examples in this section can be used in the following ways:

* in an interactive Python session (start a Python interpreter and insert
  code lines one by one).
* in a Python script to perform all calculations at once (copy the code
  into a Python file and run it). 

Notes
-------------------------------------------------------------------------------

* Some examples need internet connectivity for blast searching the PDB and 
  retrieving files from the PDB FTP server.

* Some examples attempt to download well over 100 structure 
  files, which may take several minutes depending on connection speed.  

Access PDB Data 
===============================================================================

.. toctree::
   :maxdepth: 1

   proteins/fetchpdb
   proteins/parsepdb
   proteins/writepdb
   proteins/blastpdb
   
Protein Structure 
===============================================================================

.. toctree::
   :maxdepth: 1

   select/selections
   select/contacts
   atomic/hierview
   measure/alignment
   proteins/biomolt
   compare/compare
   measure/deformation

Dynamics Analysis
===============================================================================

.. toctree::
   :maxdepth: 2
   :titlesonly:

   dynamics/enm/index
   dynamics/pca/index
   dynamics/generate_conformers
   dynamics/extrapolate

More on :class:`~atomic.AtomGroup` 
===============================================================================

.. toctree::
   :maxdepth: 1
   
   atomic/water

