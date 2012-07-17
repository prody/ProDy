.. _examples:

.. currentmodule:: prody

*******************************************************************************
Usage Examples
*******************************************************************************

How to Use
===============================================================================

Examples in this section can be used in the following ways:

* in an interactive Python session (start a Python interpreter and insert
  code lines one by one)
* in a Python script to perform all calculations at once (copy the code
  into a Python file and run it)

Notes
-------------------------------------------------------------------------------

* Some examples need internet connection for blast searching the PDB and 
  retrieving files from the `wwPDB <http://www.wwpdb.org/>`_ FTP servers.

* Some examples attempt to download well over 100 structure files, which 
  may take several minutes depending on connection speed.  

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

   measure/contacts
   measure/alignment
   proteins/biomolt
   proteins/compare
   proteins/extract_ligands
   measure/deformation

Dynamics Analysis
===============================================================================

.. toctree::
   :maxdepth: 2
   :titlesonly:

   dynamics/enm/index
   dynamics/pca/index
   dynamics/sampling
   dynamics/extendmodel
   dynamics/editmodel

Trajectory Analysis
===============================================================================

.. toctree::
   :maxdepth: 2
   :titlesonly:
   :glob:

   trajectory/*


Atomic classes and data 
===============================================================================

.. toctree::
   :maxdepth: 1
   
   atomic/userdata   

