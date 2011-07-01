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
  
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
|                                                  | Example figures                                |                                                |
+==================================================+================================================+================================================+
| .. image:: /_static/gallery/p38_modes_123_sm.png | .. image:: /_static/gallery/p38_anm_pca_sm.png | .. image:: /_static/gallery/p38_network_sm.png |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
| ANM modes 1-3 for p38 MAPK                       | ANM and PCA modes for p38                      | p38 network model                              |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+

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

|new| NMWiz can retrieve coordinate and normal mode data from molecules
loaded in VMD. See the :guilabel:`From Molecule` interface.

Example Input
===============================================================================

Examples :term:`nmd` files (*right click* and select *save as* option):

  * :download:`p38 MAP kinase <nmwiz/p38_MAPK_1p38_anm_modes.nmd.zip>`
  * :download:`p38 modes with zero coordinates <nmwiz/xyzeros.nmd.zip>` 
    that was used to make :download:`NMWiz logo </_static/nm.png>`.


NMWiz GUI
===============================================================================

This section describes NMWiz Interfaces briefly. Help on functionality can be
obtained by simply clicking on the question marks (:guilabel:`?`) on the 
graphical user interface. In addition, a tutorial that was developed for an 
earlier version of the plug-in can be obtained from here:
:download:`Tutorial <nmwiz/NMWiz_tutorial.pdf>`. 

Main Window
-------------------------------------------------------------------------------

|

========================= ====================================================================================================================================
.. image:: nmwiz/main.png Main window enables loading :term:`nmd` files and submitting ProDy and ANM server jobs. 
                          Settings and online documentation (at ProDy website) can also be accessed from this window.
                                     
                          Additionally, for each loaded nmd file two buttons will appear in the main window. 
                          :guilabel:`Show GUI` button recovers a closed interface for a dataset.
                          :guilabel:`Remove` button completely removes the dataset from the VMD session.
                                     
                          |bulb| Normal mode data can be saved in NMD format from this window, which then can be parsed with ProDy for further analysis. 
========================= ====================================================================================================================================

|

From Molecule Window
-------------------------------------------------------------------------------

|

============================ ====================================================================================================================================
.. image:: nmwiz/frommol.png |new| Use this interface to retrieve data from a molecule which contains normal modes as frames.
                             In the following two examples mode data are provided in Gromacs TRR or PDB formats. 
                             
                             |example| Lysozyme example contains 10 modes from all-atom NMA calculations (courtesy of Guang Hu). Coordinate data is in GRO format
                             and mode data is in TRR format.
                             
                             * :download:`Lysozyme dataset <nmwiz/lysozyme.zip>`
                             
                             |example| Ubiquitin example contains 10 modes from all-atom PCA calculations for the structure 2K39. Both coordinate and mode data 
                             are in PDB format. The zip archive also contains the same data in NMD format. 
                             
                             * :download:`Ubiquitin dataset <nmwiz/ubiquitin.zip>`
          
                             |bulb| When loading data from a large molecular system, you can choose to get data for select atoms, e.g. 
                             ``name CA and protein`` will obtain parts of normal modes matching carbon alpha atoms.    
                             
============================ ====================================================================================================================================

|

NMWiz Window
-------------------------------------------------------------------------------

|

======================== ====================================================================================================================================
.. image:: nmwiz/gui.png NMWiz interface is generated for all datasets loaded into NMWiz. These
                         interfaces are independent of each other, hence allows comparison of different data sets.
                                     
                         NMWiz interfaces allow visualizing, animating, and plotting modes.
                         At startup, settings panels are hidden from user.
                         Arrow size, animation length, graphics quality etc. can be
                         adjusted using option frames after they are switched on.
                                     
                         *Active Mode and Color*
                         
                         The top panel enables choosing the active mode and
                         color. When the active mode is changed, previously drawn mode arrows will
                         be undisplayed. This behavior and other graphics options can be
                         changed from the "Arrow Graphics Option" panel.
                         
                         |bulb| Are you visualizing a large system? Keep the resolution low to draw arrows faster.
                         
                         |bulb| Are the arrow graphics too crowded? Draw arrows for a subset of residues that are evenly spaced, e.g try 
                         the selection string ``residue % 4 == 0``, which will draw an arrow for every fourth residue.
======================== ====================================================================================================================================

|


Plot Window
-------------------------------------------------------------------------------

|

========================= ====================================================================================================================================
.. image:: nmwiz/plot.png User can plot squared-fluctuations along the active normal mode by clicking on
                          the :guilabel:`Plot` button. Plots will be generated using a modified version of
                          `MultiPlot <http://www.ks.uiuc.edu/Research/vmd/plugins/multiplot/>`_ plug-in.
                          
                          |bulb| Clicking on the plot will label and highlight the residue (or atom) in the VMD
                          display.
========================= ====================================================================================================================================

Settings Window
-------------------------------------------------------------------------------

|

============================= ====================================================================================================================================
.. image:: nmwiz/settings.png Settings window allows users to specify the path to ProDy scripts and select the default color for displaying arrows.
============================= ====================================================================================================================================

|

ProDy Interface
-------------------------------------------------------------------------------

|

========================== ====================================================================================================================================
.. image:: nmwiz/prody.png ProDy interface allows users to submit ProDy ANM and PCA jobs for proteins 
                           loaded in VMD. Upon completion of the calculations, NMWiz automatically loads the results.
========================== ====================================================================================================================================

.. toctree::
   :glob:
   :maxdepth: 2
   :hidden:

   *
