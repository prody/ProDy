.. _nmwiz:

*******************************************************************************
NMWiz
*******************************************************************************

Normal Mode Wizard (NMWiz) is a `VMD <www.ks.uiuc.edu/Research/vmd/>`_ 
plugin designed for visual and comparative analysis of normal mode data, 
i.e. modes from principal component, essential dynamics, or normal mode 
analysis. 

NMWiz can be used for:

  * drawing normal modes arrows
  * making animations (conformations along a normal mode)
  * plotting square-fluctuations (labeling and highlighting residues)
  * comparing two structures and drawing deformation arrows
  
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
|                                                  | Example figures                                |                                                |
+==================================================+================================================+================================================+
| .. image:: /_static/gallery/p38_modes_123_sm.png | .. image:: /_static/gallery/p38_anm_pca_sm.png | .. image:: /_static/gallery/p38_network_sm.png |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
| ANM modes 1-3 for p38 MAPK                       | ANM and PCA modes for p38                      | p38 network model                              |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+


Downloads
===============================================================================

.. note:: NMWiz is incorporated into VMD and will be available in VMD version
   1.9.1. Until this version is released, NMWiz files can be obtained from here.

NMWiz is written in `TCL <http://tcl.tk/>`_. To be able to use it, 
you need to have VMD version 1.8.7 or higher installed on your computer.
NMWiz works on all platforms that VMD is available for, including Linux, 
Mac OS, and Windows.

See |vmd| for obtaining VMD.

NMWiz files are available in the following:
 
  * :download:`NMWiz1.0.tar.gz <NMWiz1.0.tar.gz>`
  * :download:`NMWiz1.0.zip <NMWiz1.0.zip>`

See `documentation <nmwiz/index.html>`_  for details.


Installation
===============================================================================

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
  
Updates
-------------------------------------------------------------------------------

To install a newer version, you need to run the installer script again or
delete the existing version and repeat the manual instructions.


.. toctree::
   :glob:
   :maxdepth: 2
   :hidden:

