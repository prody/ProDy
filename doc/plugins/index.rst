.. _nmwiz:

*******************************************************************************
NMWiz
*******************************************************************************


Getting NMWiz
===============================================================================

Normal Mode Wizard (NMWiz) is a VMD plugin for comparative visual analysis of 
protein dynamics modeled using theory or inferred from experimental structural
ensembles.  NMWiz is available with `VMD`_ 1.9.1.  For updates follow below 
instructions.

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/

See :ref:`nmwiz-tutorial`.  NMWiz recognizes :ref:`nmd-format`.

.. _NMWiz documentation: http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/


Manual Updates
-------------------------------------------------------------------------------

There has been several improvements in NMWiz with ProDy release v1.4.
Improved version of NMWiz and other related plotting plugins can be 
found in the following files:

#. Download plugin files
  
  * :download:`nmwiz_multiplot_heatmapper.zip`
  * :download:`nmwiz_multiplot_heatmapper.tgz`

#. Update files in VMD plugins directory :file:`$VMDDIR/plugins/noarch/tcl/`.  

If you are not sure where VMD directory is located, run :program:`vmd`, and 
type the following command line in the VMD console::

    global env; puts $env(VMDDIR)

Alternatively, instructions for installing `3rd-party`_ plugins may be helpful
too.

.. _3rd-party: http://physiology.med.cornell.edu/faculty/hweinstein/vmdplugins/installation.html
