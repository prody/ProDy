.. currentmodule:: prody.dynamics

.. _pca-xray-visualization:

*******************************************************************************
PCA of X-ray structures: Visualization
*******************************************************************************

Synopsis
===============================================================================

This example is continued from :ref:`pca-xray-plotting`.
The aim of this part is to generate molecular graphics that was presented
in our paper [AB09]_.

Input
-------------------------------------------------------------------------------

PCA and ANM data calculated in :ref:`pca-xray-calculations`.

Output
-------------------------------------------------------------------------------

Qualitative molecular graphics comparing experimental and theoretical modes.

Notes
-------------------------------------------------------------------------------

To make a comparative visual analysis of PCA and ANM modes
that were calculated in the previous parts, NMWiz needs to be installed.
NMWiz is a VMD plugin designed to complement ProDy.
You can get it from |nmwiz|. 


ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Then we load data saved in :ref:`pca-xray-calculations`:

>>> pca = loadModel('p38_xray.pca.npz')
>>> anm = loadModel('1p38.anm.npz')
>>> ref_chain = parsePDB('p38_ref_chain.pdb')

Write NMD files
-------------------------------------------------------------------------------

We will save PCA and ANM data in NMD format. 
NMWiz can read and visualize multiple NMD files at once. Interested
user is referred to NMWiz documentation for more information. NMD files
are saved as follows using :func:`writeNMD` functions:

>>> writeNMD('p38_pca.nmd', pca[:3], ref_chain)
'p38_pca.nmd'
>>> writeNMD('p38_anm.nmd', anm[:3], ref_chain)
'p38_anm.nmd'
   

It is also possible to load VMD to visualize normal mode data 
from within an interactive Python session. For this to work, you need
VMD and NMWiz plugin installed. Check if VMD path is correct using :func:`getVMDpath`:
    
>>> getVMDpath()
'/usr/local/bin/vmd' 

If this is not the correct path to your VMD executable you can change it 
using :func:`setVMDpath`.
    
::

    viewNMDinVMD( '1p38_pca.nmd' )


|questions|

|suggestions|
