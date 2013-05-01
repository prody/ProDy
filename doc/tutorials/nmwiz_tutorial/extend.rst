Extending a Model
===============================================================================

In previous calculations, we used Cα atoms and the results retrieved from
ProDy contained only a trace of the structure.  VMD requires more information
(at least a complete backbone) for displaying cartoon and ribbon 
representation of proteins which are suitable for publications.  In this
part, we will use :guilabel:`Extend model to` option for extending the
model to backbone atoms of the protein.

.. figure:: /_static/nmwiz_1dlo_ANM1.png
   :align: right
   :scale: 50 %
   
   ANM mode 1 for HIV Reverse Transcriptase
   
ANM Calculation
-------------------------------------------------------------------------------

Let's fetch an X-ray structure of the protein HIV reverse transcriptase (RT)
and load into VMD::

  $ prody fetch 1dlo
  $ vmd 1dlo.pdb
  
In the :guilabel:`ProDy Interface`, we select :guilabel:`ANM Calculation`,
check :guilabel:`backbone` option, and click :guilabel:`Submit Job`.  
Model will be calculated for 971 selected Cα atoms, but the normal modes will 
be extended to all backbone atoms.


Visualization
-------------------------------------------------------------------------------

When the results are loaded, you will see four arrows per residue (or node).
Change the :guilabel:`Selection` string to read ``name CA`` and click 
:guilabel:`Redraw`.  This will draw only one arrow per mode.

RT is a large structure and updating the display with every little change you
make might be time consuming.  You can uncheck :guilabel:`auto update graphics`
option in :guilabel:`Mode Graphics Options` panel.  

To get the view displayed in the figure, you will need to hide arrows that
are shorter than a given length using :guilabel:`Draw if longer than` option
and draw an arrow for every forth residue using the selection
``name CA and resid % 4 == 0``. The protein representation is *NewCartoon*.

Animation
-------------------------------------------------------------------------------

You can generate a trajectory along the selected mode by clicking 
:guilabel:`Make` in :guilabel:`Animation` row. For large proteins,
keeping the :guilabel:`Graphics resolution` low (10) will make
the animation run smoother.
