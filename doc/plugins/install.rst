.. _getnmwiz:

*******************************************************************************
Plug-in Installation
*******************************************************************************

NMWiz is a plugin written in `TCL <http://tcl.tk/>`_. To be able to use it, 
you need to have VMD version 1.8.7 or higher installed on your computer.

If you have VMD installed when installing ProDy, ProDy will try to install
NMWiz. Alternatively, you can follow these instructions to install it yourself.

See |vmd| for obtaining VMD.

Downloads
===============================================================================

Plug-in installation file: `NMWiz 0.6 <http://pypi/NMWiz/NMWiz0.6.tar.gz>`_



Installation
===============================================================================

Following instructions apply to all computer architectures and operating 
systems that VMD runs on, but may require root access.

#. Extract tarball (:file:`NMWiz0.x.y.tar.gz`) contents.

#. Copy :file:`nmwiz0.x` folder into VMD plugin directory 
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
 
Updates
===============================================================================

To install a newer version, you will need to delete the existing plug-in 
directory and copy the new release to the same location.
