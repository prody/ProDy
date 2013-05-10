.. _pca-xray-visualization:


Visualization
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example is continued from :ref:`pca-xray-plotting`.  The aim of this part
is visual comparison of experimental and theoretical modes.
We will generate molecular graphics that was presented in our paper [AB09]_.

Notes
^^^^^

To make a comparative visual analysis of PCA and ANM modes that were calculated
in the previous parts, NMWiz needs to be installed.  NMWiz is a VMD plugin
designed to complement ProDy.  You can get it from |nmwiz|.


Load data
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Then we load data saved in :ref:`pca-xray-calculations`:

.. ipython:: python

   pca = loadModel('p38_xray.pca.npz')
   anm = loadModel('1p38.anm.npz')
   ensemble = loadEnsemble('p38_X-ray.ens.npz')
   ref_chain = parsePDB('p38_ref_chain.pdb')

Write NMD files
-------------------------------------------------------------------------------

We will save PCA and ANM data in NMD format.
NMWiz can read and visualize multiple NMD files at once. Interested
user is referred to NMWiz documentation for more information. NMD files
are saved as follows using :func:`.writeNMD` functions:

.. ipython:: python

   writeNMD('p38_pca.nmd', pca[:3], ref_chain)
   writeNMD('p38_anm.nmd', anm[:3], ref_chain)


It is also possible to load VMD to visualize normal mode data
from within an interactive Python session. For this to work, you need
VMD and NMWiz plugin installed. Check if VMD path is correct using
:func:`.pathVMD`:

.. ipython:: python

   pathVMD()

If this is not the correct path to your VMD executable you can change it
using the same function.

.. ipython:: python
   :verbatim:

   viewNMDinVMD('1p38_pca.nmd')

