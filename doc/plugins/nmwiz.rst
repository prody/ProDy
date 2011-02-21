.. _nmwiz:

*******************************************************************************
NMWiz Usage
*******************************************************************************

This section describes NMWiz Interfaces briefly. Help on functionality can be
obtained by simply clicking on the question marks (**?**) on the graphical 
user interface.

Examples
===============================================================================

Examples :term:`nmd` files:

  * :download:`p38 MAP kinase <p38_MAPK_1p38_anm_modes.nmd>`
  * :download:`p38 modes with zero coordinates <xyzeros.nmd>` 
    that was used to make NMWiz logo.

Tutorial: 
  
  * :download:`Tutorial <NMWiz_tutorial.pdf>`


Main
===============================================================================

NMWiz Main enables loading data from a :ref:`nmdformat` file or retreiving data 
from `ANM server <http://ignmtest.ccbb.pitt.edu/cgi-bin/anm/anm1.cgi>`_. 
License information and documentation (at NMWiz website)
can also be accessed from the :guilabel:`Main` window.

.. image:: /_static/gui_main.png

Additionally, for each loaded data set two buttons will appear in the main 
window. :guilabel:`Show GUI` button recovers a closed interface for a dataset.
:guilabel:`Remove Data & GUI` button completely removes the dataset from
the VMD session.

ANM Server
===============================================================================


Normal mode data can be retreived from ANM server using this interface. User
needs to provide PDB and chain identifiers. Cutoff and distance weight 
parameters can also be adjusted by the user.

.. image:: /_static/gui_anm.png

After ANM Server completes calculations, the user needs to download normal
mode data in an NMD file and load it to VMD.

NMWiz
===============================================================================


This interface is generated for all datasets loaded into NMWiz. These interfaces
are independent of each other, hence allows comparative analysis of different
data sets.

NMWiz interaces allow visualizing, animating, and plotting modes. 
At startup, variuos options are hidden from user. 
Arrow size, amination length, graphics quality etc. can be 
adjusted using option frames after they are switched on.

**Active Mode and Color**

The top panel enables choosing the active mode and display color of the active
mode. Initially, when the active mode is changed, previously drawn mode arrows
will be removed. This can be changed from the "Arrow Graphics Option" panel.
Other arrow graphics parameters can also be adjusted from this panel.

.. image:: /_static/gui_nmwiz.png


Plot
===============================================================================


User can plot squared-fluctuations along the active normal mode by clicking on
the :guilabel:`Plot` button. Plots will
be generated using a modified version of 
`MultiPlot <http://www.ks.uiuc.edu/Research/vmd/plugins/multiplot/>`_.
Clicking on the plot will label and highlight the residue (or atom) in the VMD
display.

.. image:: /_static/gui_plot.png
