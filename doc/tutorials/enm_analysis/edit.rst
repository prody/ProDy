.. _reduce-slice:


Editing a Model
===============================================================================

This example shows how to analyze the normal modes corresponding to
a system of interest.  In this example, ANM calculations will be performed for
HIV-1 reverse transcriptase (RT) subunits p66 and p51. Analysis will be made
for subunit p66.  Output is a reduced/sliced model that can be used as input
to analysis and plotting functions.

ANM calculations
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()


We start with parsing the Cα atoms of the RT structure 1DLO and performing ANM
calculations for them:

.. ipython:: python

   rt = parsePDB('1dlo', subset="ca")
   anm, sel = calcANM(rt)
   anm
   saveModel(anm, 'rt_anm')
   anm[:5].getEigvals().round(3)
   (anm[0].getArray() ** 2).sum() ** 0.5


Analysis
^^^^^^^^

We can plot the cross-correlations and square fluctuations for the full model
as follows:

Cross-correlations
""""""""""""""""""

.. ipython:: python


   @savefig enm_analysis_edit_crosscorr.png width=4in
   showCrossCorr(anm);


Square fluctuations
"""""""""""""""""""

.. ipython:: python

   @savefig enm_analysis_edit_sqflucts.png width=4in
   showSqFlucts(anm[0]);


Slicing a model
-------------------------------------------------------------------------------

Slicing a model is analogous to slicing a list, i.e.:

.. ipython:: python

   numbers = list(range(10))
   numbers
   slice_first_half = numbers[:10]
   slice_first_half

In this case, we want to slice normal modes, so that we will handle mode
data corresponding to subunit p66, which is chain A in the structure.
We use :func:`.sliceModel` function:

.. ipython:: python

   anm_slc_p66, sel_p66 = sliceModel(anm, rt, 'chain A')
   anm_slc_p66

You see that now the sliced model contains 556 nodes out of the
971 nodes in the original model.

.. ipython:: python

   saveModel(anm_slc_p66, 'rt_anm_sliced')
   anm_slc_p66[:5].getEigvals().round(3)
   '%.3f' % (anm_slc_p66[0].getArray() ** 2).sum() ** 0.5

Note that slicing does not change anything in the model apart from taking parts
of the modes matching the selection. The sliced model contains fewer nodes,
has the same eigenvalues, and modes in the model are not normalized.

Analysis
^^^^^^^^

We plot the cross-correlations and square fluctuations for the sliced model
in the same way. Note that the plots contain the selected part of the model
without any change:

Cross-correlations
""""""""""""""""""

.. ipython:: python

   showCrossCorr(anm_slc_p66);

   @savefig enm_analysis_edit_slice_cc.png width=4in
   title('Cross-correlations for ANM slice');



Square fluctuations
"""""""""""""""""""

.. ipython:: python

   @savefig enm_analysis_edit_slice_sqf.png width=4in
   showSqFlucts(anm_slc_p66[0]);


Reducing a model
-------------------------------------------------------------------------------

We reduce the ANM model to subunit p66 using :func:`.reduceModel` function.
This function implements the method described in 2000 paper of Hinsen et al.
[KH00]_

.. ipython:: python

   anm_red_p66, sel_p66 = reduceModel(anm, rt, 'chain A')
   anm_red_p66.calcModes()
   anm_red_p66
   saveModel(anm_red_p66, 'rt_anm_reduced')
   anm_red_p66[:5].getEigvals().round(3)
   '%.3f' % (anm_red_p66[0].getArray() ** 2).sum() ** 0.5


Analysis
^^^^^^^^

We plot the cross-correlations and square fluctuations for the reduced model
in the same way. Note that in this case the plots are not identical to the
full model:

Cross-correlations
""""""""""""""""""

.. ipython:: python

   @savefig enm_analysis_edit_reduce_cc.png width=4in
   showCrossCorr(anm_red_p66);

Square fluctuations
"""""""""""""""""""

.. ipython:: python

   @savefig enm_analysis_edit_reduce_sqf.png width=4in
   showSqFlucts(anm_red_p66[0]);


Compare reduced and sliced models
-------------------------------------------------------------------------------

We can compare the sliced and reduced models by plotting the overlap table
between modes:

.. ipython:: python

   @savefig enm_analysis_edit_overlap.png width=4in
   showOverlapTable(anm_slc_p66, anm_red_p66);


The sliced and reduced models are not the same. While the purpose of slicing is
simply enabling easy plotting/analysis of properties of a part of the system,
reducing has other uses as in [WZ05]_.

.. [WZ05] Zheng W, Brooks BR. Probing the Local Dynamics of Nucleotide-Binding
   Pocket Coupled to the Global Dynamics: Myosin versus Kinesin.
   *Biophysical Journal*  **2005** 89:167–178.
