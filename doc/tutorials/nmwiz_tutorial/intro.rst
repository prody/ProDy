.. _nmwiz:

Normal Mode Wizard
===============================================================================

Normal Mode Wizard (NMWiz) is a VMD plugin for depiction, animation, and
comparative analysis of normal modes.  Normal modes may come from principal
component of structural ensembles, essential dynamics analysis of simulation
trajectories, or normal mode analysis of protein structures.  In addition,
NMWiz can be used to depict any vector that describes a molecular motion.


Required Programs
-------------------------------------------------------------------------------

Latest version of `ProDy`_ and `VMD`_ are required.

.. _ProDy: http://csb.pitt.edu/ProDy/index.html#downloads

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/


Manual Updates
-------------------------------------------------------------------------------

There has been several improvements in NMWiz with ProDy release v1.4.
Improved version of NMWiz and other related plotting plugins can be
found in the following files:

1. Download plugin files:

  * :download:`nmwiz_multiplot_heatmapper.zip`
  * :download:`nmwiz_multiplot_heatmapper.tgz`

2. Update files in VMD plugins directory :file:`$VMDDIR/plugins/noarch/tcl/`.

If you are not sure where VMD directory is located, run :program:`vmd`, and
type the following command line in the VMD console::

    global env; puts $env(VMDDIR)


Getting Started
-------------------------------------------------------------------------------

To follow this tutorial, you will need the following files which can be
downloaded from :ref:`tutorials`.

.. literalinclude:: files.txt
