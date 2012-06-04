.. _nmwiz:

*******************************************************************************
Getting NMWiz
*******************************************************************************

Normal Mode Wizard (NMWiz) is a VMD plugin designed for visual comparative 
analysis of normal mode data.  NMWiz is available with `VMD`_ 1.9.1 or later.

.. _VMD: http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

See :ref:`nmwiz-tutorial` or `NMWiz documentation`_ for usage details.  NMWiz 
recognizes :ref:`nmd-format`.

.. _NMWiz documentation: http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/


Manual Updates
-------------------------------------------------------------------------------

NMWiz plugin distributed with VMD is updated when changes and improvements 
are available (see :ref:`changes`).  In the case that updates are not yet 
available in the latest version of VMD, you can follow below instructions 
to update your NMWiz copy.  These instructions apply to all computer 
architectures and operating systems that VMD runs on, but may require 
root (or administrator) access.

#. Retrieve plugin files from 
   https://github.com/abakan/ProDy/tree/master/plugins/nmwiz
  
  * :file:`nmwiz.tcl`
  * :file:`pkgIndex.tcl`

#. Check the version number in :file:`pkgIndex.tcl`.

#. Copy files into :file:`nmwiz{x.y}` folder in the VMD plugins directory 
   (:file:`$VMDDIR/plugins/noarch/tcl/`).  You may need to make the directory
   first.

#. If this is not an update, insert the following line to 
   :file:`$VMDDIR/scripts/vmd/loadplugins.tcl` file in VMD directory at 
   line number 143 (or a suitable place you like)::

    vmd_install_extension nmwiz nmwiz_tk "Analysis/Normal Mode Wizard"


If you are not sure where VMD directory is located, run :program:`vmd`, and 
type the following command line in the VMD console::

    global env; puts $env(VMDDIR)

Once you perform these steps, NMWiz GUI will show up in 
:menuselection:`Extensions --> Analysis` menu of VMD main window. 
It is also possible to make it appear in another :menuselection:`Extensions` 
submenu by replacing *Analysis* in step 3 with another submenu name.

Alternatively, instructions for installing `3rd-party`_ plugins may be helpful
too.

.. _3rd-party: http://physiology.med.cornell.edu/faculty/hweinstein/vmdplugins/installation.html
