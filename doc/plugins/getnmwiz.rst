.. _getnmwiz:

*******************************************************************************
NMWiz Installation
*******************************************************************************

NMWiz is written in `TCL <http://tcl.tk/>`_. To be able to use it, 
you need to have VMD version 1.8.7 or higher installed on your computer.

See |vmd| for obtaining VMD.

Downloads
===============================================================================

Plug-in files:
 
  * Windows: :download:`NMWiz-0.6.zip <nmwiz/NMWiz-0.7.zip>`
  * Linux: :download:`NMWiz-0.6.tar.gz <nmwiz/NMWiz-0.7.tar.gz>`


Installation
===============================================================================

Manual Installation
-------------------------------------------------------------------------------

Following instructions apply to all computer architectures and operating 
systems that VMD runs on, but may require root (or administrator) access.

#. Extract tarball/zip (:file:`NMWiz-0.{x}.tar.gz`) contents.

#. Copy :file:`nmwiz0.{x}` folder into VMD plug-ins directory 
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
will locate VMD plug-ins directory, copy the files, remove older versions if 
found, and update the :file:`loadplugins.tcl` file. Again, this script
also requires write access to the VMD folders. On Linux, following command
should work:: 
 
  $ sudo python install_nmwiz.py
  
Updates
===============================================================================

To install a newer version, you will need to delete the existing plug-in 
directory and copy the new release to the same location.
