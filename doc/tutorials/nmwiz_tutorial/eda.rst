Essential Dynamics Analysis
===============================================================================

In this part, we will perform essential dynamics analysis calculations
for a p38 MAP kinase trajectory and compare results with normal modes from
ANM calculations.  p38 files can be downloaded from :download:`p38 trajectory
<nmwiz_tutorial_files/p38_100frames.dcd>`.  The :file:`p38_100frames.dcd` is
from a 40 ns long simulation of p38.  Detailed analysis of this trajectory can
be found in [AB11]_.

Load the contents of this file into VMD as follows::

  $ tar -xzf p38_trajectory.tar.gz
  $ vmd p38.pdb p38_100frames.dcd


Click on :guilabel:`ProDy Interface` for performing ANM and EDA jobs.

.. [AB11] Bakan A, Bahar I. Computational generation of inhibitor-bound
   conformers of p38 MAP kinase and comparison with experiments. *Pacific
   Symposium on Biocomputing* **2011** 16 181-192.

EDA Calculation
-------------------------------------------------------------------------------

Select :guilabel:`PCA Calculation` in the :guilabel:`ProDy Job Settings` panel.
Set :guilabel:`First frame` 1 to exclude the X-ray coordinates from the
calculation.  You may also check :guilabel:`aligned` to make ProDy skip
alignment step in PCA/EDA calculations when you already have done the alignment
in VMD. In this case, the frames are already aligned.
Click :guilabel:`Submit Job` and results will be loaded automatically when
:ref:`prody-pca` command finishes the calculations.

ANM Calculation
-------------------------------------------------------------------------------

We will perform ANM calculations for all Cα atoms and keep the rest of the
parameters unchanged.  Click :guilabel:`Submit Job` and results obtained from
:ref:`prody-anm` command will load automatically.

.. figure:: /_static/nmwiz_p38_EDA1vsANM1.png
   :align: right
   :scale: 50 %

   EDA 1 (orange) vs. ANM mode 2 (lime green)

Comparison
-------------------------------------------------------------------------------

For each dataset you load into or generate via NMWiz, a GUI will pop up with
independent controls for normal mode display, animation, and plotting. Select
PC 2 and ANM mode 2 and try to get the view in the image in VMD display.


Suggestions
-------------------------------------------------------------------------------

NMWiz writes a DCD or PDB file for PCA/EDA calculations.  For large systems
and long trajectories you may try one or more of the following for speedier
calculations:

  * select a subset of atoms, e.g. Cα atoms
  * select a subset of frames, e.g. set :guilabel:`Skip frame` a value
    greater than 0
  * use :guilabel:`DCD file` for faster IO operations and less disk usage
  * alternatively, if you have trajectories in DCD format, use :ref:`prody-pca`
    directly to obtain results in :ref:`nmd-format`
