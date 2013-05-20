Introduction
===============================================================================

*ProDy* is an application programming interface (API) designed for
structure-based analysis of protein dynamics, in particular for inferring
protein dynamics from large heterogeneous structural ensembles.  It comes with
several command line applications (:ref:`prody-apps`) and graphical user
interface for visualization (GUI, :ref:`nmwiz`).  This tutorial shows core
features of *ProDy* and some basic analysis tasks.  You can find links to more
detailed and advanced tutorials below.


Structural Ensemble Analysis
-------------------------------------------------------------------------------

*ProDy* is primarily designed for analysis of *large* heterogeneous structural
datasets that may be composed of sequence homologs, mutants, or ligand bound
forms of a protein that have with missing loops or terminal residues.  Dominant
patterns in structural variability are extracted by principal component
analysis (PCA) of the ensemble.  Helper functions allow for comparison of
dynamics inferred from experiments with theoretically models and simulation
data.  For detailed usage examples see :ref:`ensemble-analysis`.


Elastic Network Models
-------------------------------------------------------------------------------

*ProDy* can be used for normal mode analysis (NMA) of protein dynamics based
on elastic network models (ENMs).  Flexible classes allow for developing and
using customized gamma functions in ENMs and numerous helper functions allow
for comparative analysis of experimental and theoretical datasets.  See
:ref:`enm-analysis` for detailed usage examples.


Trajectory Analysis
-------------------------------------------------------------------------------

In addition to analysis of experimental data and theoretical models, *ProDy*
can be used to analyze trajectories from molecular dynamics simulations, such
as for performing essential dynamics analysis (EDA).  *ProDy* supports DCD
file format, but trajectories in other formats can be parsed using other
Python packages and analyzed using ProDy.  See :ref:`trajectory-analysis`
for detailed usage examples.


Visualization
-------------------------------------------------------------------------------

Finally, results from *ProDy* calculations can be visualized using NMWiz,
which is a `VMD`_ plugin GUI. NMWiz can also be used for submitting *ProDy*
calculations for molecules in VMD.  See :ref:`nmwiz-tutorial` for analysis
of various types of datasets and visualization of protein dynamics.

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/