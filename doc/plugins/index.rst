.. _nmwiz:

*******************************************************************************
NMWiz
*******************************************************************************

Normal Mode Wizard (NMWiz) is a `VMD <www.ks.uiuc.edu/Research/vmd/>`_ 
plugin designed for visual comparative analysis of normal mode data, 
i.e. modes may come from principal component, essential dynamics, normal 
mode analysis or may be any vector describing a molecular motion. 

NMWiz can be used for:

  * drawing normal modes arrows
  * making animations (conformations along a normal mode)
  * plotting square-fluctuations (labeling and highlighting residues)
  * comparing two structures and drawing deformation arrows
  
Following molecular representations are prepared using NMWiz:
  
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
|                                                  | Example figures                                |                                                |
+==================================================+================================================+================================================+
| .. image:: /_static/gallery/p38_modes_123_sm.png | .. image:: /_static/gallery/p38_anm_pca_sm.png | .. image:: /_static/gallery/p38_network_sm.png |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
| ANM modes 1-3 for p38 MAPK                       | ANM and PCA modes for p38                      | p38 network model                              |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+

NMWiz can also be used to generate trajectories on the fly.  The movie shows 
normal mode representation and animation generated using NMWiz.  Anisotropic 
network model modes were calculated using ProDy.  Movie was generated using 
`VMD Movie Plugin <http://www.ks.uiuc.edu/Research/vmd/plugins/vmdmovie/>`_.

.. only:: html

   .. youtube:: 1OUzdzm68YY
      :width: 400

See `documentation <http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/>`_  
NMWiz usage for details. NMWiz recognizes :ref:`nmd-format`, which is 
a simple data format described below.

Availability
===============================================================================

You can obtain a copy of NMWiz with *VMD 1.9.1*. Please see `VMD downloads 
<http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD>`_.


.. 
    NMWiz also works with VMD version 1.8.7 or higher. You can 
    NMWiz works on all platforms that VMD is available for, including Linux, 
    Mac OS, and Windows.


..
    Manual Installation
    -------------------------------------------------------------------------------

    Following instructions apply to all computer architectures and operating 
    systems that VMD runs on, but may require root (or administrator) access.

    #. Extract tarball/zip (:file:`NMWiz1.{x}.tar.gz`) contents.

    #. Copy :file:`nmwiz1.{x}` folder into VMD plugins directory 
       (:file:`$VMDDIR/plugins/noarch/tcl/`).

    #. Insert the following line to :file:`$VMDDIR/scripts/vmd/loadplugins.tcl` 
       file in VMD directory at line number 143 (or wherever you like)::

        vmd_install_extension nmwiz nmwiz_tk "Analysis/Normal Mode Wizard"


    If you are not sure where VMD directory is located, run :program:`vmd`, and 
    type the following command line in the VMD console::

        global env; puts $env(VMDDIR)

    Once you perform these steps, NMWiz GUI will show up in 
    :menuselection:`Extensions --> Analysis` menu of VMD main window. 
    It is also possible to make it appear in another :menuselection:`Extensions` 
    submenu by replacing *Analysis* in step 3 with another submenu name.

..      
    Using Installer Script
    -------------------------------------------------------------------------------

    Alternatively, you can use :file:`install_nmwiz.py` script. This script
    will locate VMD plugins directory, copy the files, remove older versions if 
    found, and update the :file:`loadplugins.tcl` file. Again, this script
    also requires write access to the VMD folders. On Linux, following command
    should work:: 
     
      $ sudo python install_nmwiz.py
      
    This installer script works in Linux and Windows. Mac OS users, please
    follow the manual installation instructions.

..    
    Updates
    -------------------------------------------------------------------------------

    To install a newer version, you need to run the installer script again or
    delete the existing version and repeat the manual instructions.


.. _nmd-format:

NMD format
===============================================================================

Description
-------------------------------------------------------------------------------

NMD files (extension :file:`.nmd`) are plain text files that contain at 
least normal mode and system coordinate data.

NMD files can be visualized using :ref:`nmwiz`. 
ProDy functions :func:`~.writeNMD` and :func:`~.parseNMD` can be used to read 
and write NMD files. 

Data fields
-------------------------------------------------------------------------------

Data fields in bold face are required. All data arrays and lists must be in a 
single line and items must be separated by one or more space characters.

**coordinates**: system coordinates as a list of decimal numbers
  Coordinate array is the most important line in an NMD file. All mode array 
  lengths must match the length of the coordinate array. Also, number of atoms
  in the system is deduced from the length of the coordinate array.

::

  coordinates 27.552 4.354 23.629 24.179 4.807 21.907 ...

**mode**: normal mode array as a list of decimal numbers
  Optionally, mode index and a scaling factor may be provided
  in the same line as a mode array. Both of these must precede the mode array.
  Providing a scaling factor enables relative scaling of the mode arrows and
  the amplitude of the fluctuations in animations. For NMA, scaling factors
  may be chosen to be the square-root of the inverse-eigenvalue associated
  with the mode. Analogously, for PCA data, scaling factor would be the 
  square-root of the eigenvalue.
  
  If a mode line contains numbers preceding the mode array, they are evaluated 
  based on their type. If an integer is encountered, it is considered the mode 
  index. If a decimal number is encountered, it is considered the scaling 
  factor. Scaling factor may be the square-root of the inverse eigenvalue
  if data is from an elastic network model, or the square-root of the 
  eigenvalue if data is from an essential dynamics (or principal component) 
  analysis.
  
  For example, all of the following lines are valid. The first line contains
  mode index and scaling factor. Second and third lines contain mode index or
  scaling factor. Last line contains only the mode array.

::

  mode 1 2.37    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 1    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 2.37    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 0.039 0.009 0.058 0.038 -0.011 0.052 0.043  ...
  
*name*: name of the model

The length of all following data fields must be equal to the number of atoms in
the system. NMWiz uses such data when writing a temporary PDB files for
loading coordinate data into VMD.

*atomnames*: list of atom names
  If not provided, all atom names are set to "CA".
  
*resnames*: list of residue names
  If not provided, all residue names are set to "GLY".
  
*chainids*: list of chain identifiers
  If not provided, all chain identifiers are set to "A".

*resids*: list of residue numbers
  If not provided, residue numbers are started from 1 and incremented by one 
  for each atom.

*bfactors*: list of experimental beta-factors
  If not provided, all beta-factors are set to zero. 
  Beta-factors can be used to color the protein representation.
  
NMD files may contain additional lines. Only lines that start with one of the 
above field names are evaluated by NMWiz.

Examples
-------------------------------------------------------------------------------


**Example 1**

File: :download:`p38anm.nmd <p38anm.nmd.gz>`

This example contains normal modes from Anisotropic Network Model calculations
for p38 MAP kinase (PDB identifier 1P38). This example
contains all required and optional data fields. Mode arrays are preceded by 
mode indices and square-root of inverse eigenvalues.

::

  name 1p38.anm
  atomnames CA CA CA ...
  resnames GLU ARG PRO ...
  chainids A A A ...
  resids 4 5 6 ...
  bfactors 69.99 59.83 47.29 ...
  coordinates 27.552 4.354 23.629 24.179 4.807 21.907 ...
  mode 1 2.37 0.039 0.009 0.058 0.038 -0.011 0.052 ...
  mode 2 1.73 -0.045 -0.096 -0.009 -0.040 -0.076 -0.010 ...
  mode 3 1.70 0.007 -0.044 0.080 0.015 -0.037 0.062 0.012 ...
  mode 4 1.12 0.010 0.024 0.003 0.007 0.017 0.004 0.010 ...
  mode 5 1.03 0.006 0.010 0.025 0.007 0.003 0.017 0.007 ...
  mode 6 0.99 -0.063 -0.066 0.060 -0.054 -0.045 0.049 ...

**Example 2** 

File: :download:`xyzeros.nmd <xyzeros.nmd.gz>`

This example contains minimal amount of data sufficient for visualizing modes.
Mode data comes from *Example 1*.  Coordinates are set to zero. NMWiz Logo on 
the upper left corner of this documentation is generated using this NMD file.

::

  coordinates 0 0 0 0 0 0 ...
  mode 0.039 0.009 0.058 0.038 -0.011 0.052 ...
  mode -0.045 -0.096 -0.009 -0.040 -0.076 -0.010 ...
  mode 0.007 -0.044 0.080 0.015 -0.037 0.062 ...


Autoload Trick
-------------------------------------------------------------------------------

By adding a special line in an NMD file, file content can be automatically 
loaded into VMD at startup. The first line calls a NMWiz function to load the 
file itself (:file:`xyzeros.nmd`).

::

  nmwiz_load xyzeros.nmd
  coordinates 0 0 0 0 0 0  ...
  mode 0.039 0.009 0.058 0.038 -0.011 0.052 ...
  mode -0.045 -0.096 -0.009 -0.040 -0.076 -0.010 ...
  mode 0.007 -0.044 0.080 0.015 -0.037 0.062 ...


In this case, VMD must be started from the command line by typing 
:program:`vmd -e xyzeros.nmd`.


.. toctree::
   :glob:
   :maxdepth: 2
   :hidden:

