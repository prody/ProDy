.. _nmwiz:

*******************************************************************************
NMWiz
*******************************************************************************

Normal Mode Wizard (NMWiz) is a `VMD <www.ks.uiuc.edu/Research/vmd/>`_ 
plug-in designed for visual and comparative analysis of normal mode data, 
i.e. modes from principal component, essential dynamics, or normal mode 
analysis. 

NMWiz can be used for:

  * visualizing normal modes with arrows
  * generating alternate conformations along a normal mode (animation)
  * plotting squared-fluctuations (allows labeling and highlighting residues)
  * comparative analysis of normal modes from different sources/methods
  * preparing publication quality images of normal modes with the power of VMD
  
See examples figures on the side panel. 
For **downloads** and **installation** instructions see :ref:`getnmwiz`.


.. contents:: NMWiz User Guide
   :local:
   :backlinks: none

Input for NMWiz
===============================================================================

NMWiz recognizes :ref:`nmd-format` (:file:`.nmd`). NMD files can be 
generated using ProDy or obtained from the 
`ANM server <http:/www.csb.pitt.edu/ANM/>`_.  Alternatively, you can prepare 
your results in NMD format for analysis with NMWiz. NMD is a simple plain text 
format and can be easily prepared using a text editor. See 
:ref:`nmd-format` for details.

Example Input
===============================================================================

Examples :term:`nmd` files (*right click* and select *save as* option):

  * :download:`p38 MAP kinase </_downloads/p38_MAPK_1p38_anm_modes.nmd>`
  * :download:`p38 modes with zero coordinates </_downloads/xyzeros.nmd>` 
    that was used to make :download:`NMWiz logo </_static/nm.png>`.


NMWiz GUIs
===============================================================================

This section describes NMWiz Interfaces briefly. Help on functionality can be
obtained by simply clicking on the question marks (:guilabel:`?`) on the 
graphical user interface. In addition, a tutorial that was developed for an 
earlier version of the plug-in can be obtained from here:
:download:`Tutorial </_downloads/NMWiz_tutorial.pdf>`. 

Main
-------------------------------------------------------------------------------

*NMWiz - Main* enables loading :term:`nmd` files and submitting
ProDy and ANM server jobs. Settings and online documentation (at ProDy website)
can also be accessed from this window.

.. image:: /_static/NMWiz-Main.png

Additionally, for each loaded nmd file two buttons will appear in the main 
window. :guilabel:`Show GUI` button recovers a closed interface for a dataset.
:guilabel:`Remove` button completely removes the dataset from
the VMD session.

.. note::
   
   NMWiz does not allow for loading the same file twice. To get around this, 
   file may be renamed or the dataset may be removed from from the session.


NMWiz
-------------------------------------------------------------------------------

Below interface is generated for all datasets loaded into NMWiz. These 
interfaces are independent of each other, hence allows comparative 
analysis of distinct data sets.

NMWiz interfaces allow visualizing, animating, and plotting modes. 
At startup, various options are hidden from user. 
Arrow size, animation length, graphics quality etc. can be 
adjusted using option frames after they are switched on.

*Active Mode and Color*: The top panel enables choosing the active mode and 
color. When the active mode is changed, previously drawn mode arrows will 
be undisplayed. This behavior and other graphics options can be
changed from the "Arrow Graphics Option" panel.

.. image:: /_static/NMWiz-GUI.png


Mobility Plots
-------------------------------------------------------------------------------

User can plot squared-fluctuations along the active normal mode by clicking on
the :guilabel:`Plot` button. Plots will be generated using a modified version of 
`MultiPlot <http://www.ks.uiuc.edu/Research/vmd/plugins/multiplot/>`_ plug-in.
Clicking on the plot will label and highlight the residue (or atom) in the VMD
display.

.. image:: /_static/NMWiz-Plot.png


Settings
-------------------------------------------------------------------------------

Below window allows users to specify the path to ProDy scripts and select
the default color for displaying arrows.  

.. image:: /_static/NMWiz-Settings.png


ProDy Interface
-------------------------------------------------------------------------------

Below interface allows users to submit ProDy ANM and PCA jobs for proteins
loaded in VMD. Upon completion of the calculations, NMWiz automatically
loads the results.

.. image:: /_static/NMWiz-ProDy.png


ANM Server
-------------------------------------------------------------------------------

Finally, ANM jobs can be submitted to ANM server using the below interface 
(Linux only). User needs to provide PDB and chain identifiers. Cutoff and 
distance weight parameters can also be adjusted by the user.

.. image:: /_static/NMWiz-ANMServer.png

After ANM Server completes calculations, the user needs to download normal
mode data in an NMD file (see *download files* link in results page) 
and load it into VMD.


.. toctree::
   :glob:
   :maxdepth: 2
   :hidden:

   *
