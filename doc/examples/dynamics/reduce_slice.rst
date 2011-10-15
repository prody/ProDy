.. currentmodule:: prody.dynamics

.. _reduce-slice:

*******************************************************************************
Reducing/slicing a model
*******************************************************************************

Synopsis
=============================================================================

This example shows how to analyze part of a system that is to interest. 
In this example, ANM calculations will be performed for HIV-1 reverse
transcriptase (RT) subunits p66 and p51. Analysis will be made for 
subunit p66.

Input
-------------------------------------------------------------------------------

A PDB structure.

Output
-------------------------------------------------------------------------------

Output is a reduced/sliced model that can be used as input to analysis and
plotting functions.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

ANM calculations
-------------------------------------------------------------------------------

We start with parsing the Cα atoms of the RT structure 1DLO and performing ANM
calculations for them:

>>> rt = parsePDB('1dlo', subset="ca")
>>> anm, sel = calcANM(rt)
>>> anm
<ANM: 1dlo_ca (20 modes, 971 nodes)>
>>> saveModel(anm, 'rt_anm')
'rt_anm.anm.npz'
>>> anm[:5].getEigenvalues()
array([ 0.039,  0.063,  0.126,  0.181,  0.221])
>>> '%.3f' % (anm[0].getArray() ** 2).sum() ** 0.5
'1.000'


.. plot::
   :context:
   :nofigs:

   from prody import *
   rt = parsePDB('1dlo', subset="ca")
   anm = loadModel('rt_anm.anm.npz')
   anm_slc_p66 = loadModel('rt_anm_sliced.anm.npz')
   anm_red_p66 = loadModel('rt_anm_reduced.anm.npz')
   import matplotlib.pyplot as plt
   plt.close('all')
   
Analysis of full model
-------------------------------------------------------------------------------

We can plot the cross-correlations and square fluctuations for the full model
as follows:

**Cross-correlations**

.. plot::
   :context:
   :include-source:
   
   import matplotlib.pyplot as plt
   plt.figure(figsize=(5,4))
   showCrossCorrelations(anm)

   
.. plot::
   :context:
   :nofigs:

   plt.close('all')


**Square fluctuations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showSqFlucts(anm[0])
   
.. plot::
   :context:
   :nofigs:

   plt.close('all') 


Slicing a model
-------------------------------------------------------------------------------

We take the slice of the ANM model corresponding to subunit p66, which is 
chain A in the structure, using :func:`sliceModel` function:

>>> anm_slc_p66, sel_p66 = sliceModel(anm, rt, 'chain A')
>>> anm_slc_p66
<ANM: 1dlo_ca slice "chain A" (20 modes, 556 nodes)>
>>> saveModel(anm_slc_p66, 'rt_anm_sliced')
'rt_anm_sliced.anm.npz'
>>> anm_slc_p66[:5].getEigenvalues()
array([ 0.039,  0.063,  0.126,  0.181,  0.221])
>>> '%.3f' % (anm_slc_p66[0].getArray() ** 2).sum() ** 0.5
'0.895'

Slicing do not change anything in the model apart from taking parts of the 
modes matching the selection. Note that the sliced model contains fewer nodes, 
has the same eigenvalues, and that the sliced modes are not normal.

Analysis of the slice
-------------------------------------------------------------------------------

We plot the cross-correlations and square fluctuations for the sliced model
in the same way. Note that the plots contain selected part of the model
without any change:

**Cross-correlations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showCrossCorrelations(anm_slc_p66) 
   
   # If title does not fit the figure, can be changed as follows
   plt.title('Cross-correlations for ANM slice')
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')


**Square fluctuations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showSqFlucts(anm_slc_p66[0])
   
.. plot::
   :context:
   :nofigs:

   plt.close('all') 

Reducing a model
-------------------------------------------------------------------------------

We reduce the ANM model to subunit p66 using :func:`reduceModel` function. This
function implements the method described in 2000 paper of Hinsen et al. [KH00]_

>>> anm_red_p66, sel_p66 = reduceModel(anm, rt, 'chain A')
>>> anm_red_p66.calcModes()
>>> anm_red_p66
<ANM: 1dlo_ca reduced (20 modes, 556 nodes)>
>>> saveModel(anm_red_p66, 'rt_anm_reduced')
'rt_anm_reduced.anm.npz'
>>> anm_red_p66[:5].getEigenvalues()
array([ 0.05 ,  0.098,  0.214,  0.289,  0.423])
>>> '%.3f' % (anm_red_p66[0].getArray() ** 2).sum() ** 0.5
'1.000'


Analysis of the slice
-------------------------------------------------------------------------------

We plot the cross-correlations and square fluctuations for the reduced model
in the same way. Note that in this case the plots are not identical to the
full model:

**Cross-correlations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showCrossCorrelations(anm_red_p66) 
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')


**Square fluctuations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showSqFlucts(anm_red_p66[0])
   
.. plot::
   :context:
   :nofigs:

   plt.close('all') 

Compare reduced and sliced models
-------------------------------------------------------------------------------

We can compare the sliced and reduced models by plotting the overlap table
between modes:

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showOverlapTable(anm_slc_p66, anm_red_p66)
   
.. plot::
   :context:
   :nofigs:

   plt.close('all') 

The sliced and reduced models are not the same. While the purpose of slicing is 
simply enabling easy plotting/analysis of properties of a part of the system, 
reducing has other uses as in [WZ05]_.  

|questions|

|suggestions|
