.. currentmodule:: prody.dynamics

.. _gnm:

*******************************************************************************
Gaussian Network Model analysis
*******************************************************************************

This example shows how to perform GNM calculations using a ubiquitin NMR
model.  

Prepare protein
===============================================================================

.. plot::
   :context:
   :nofigs:
   :include-source:
    
   from prody import *
    
   # parse pdb file by identifier
   # if file is not found in current working directory, it will be downloaded
   ubi = parsePDB('2k39', model=1)
   print ubi # print p38 to see number of atoms and number of coordinate sets

   # we only want to use alpha carbons, so we select them
   calphas = ubi.select('calpha')
   # "calpha" is a keyword and is equivalent to "protein and name CA"
   print calphas # this will show the number of atoms in the selection 

   # Note that, atoms other than alpha carbons can be selected and used in GNM 
   # calculations

Perform calculations
===============================================================================

**Build Kirchoff**
 
.. plot::
   :context:
   :nofigs:
   :include-source:
    
   # instantiate GNM instance
   gnm = GNM('Ubiquitin')

   # build Kirchhoff matrix using selected atoms (351 alpha carbons)
   gnm.buildKirchhoff(calphas)

   # get a copy of the Kirchhoff matrix
   gnm.getKirchhoff()

**How to change cutoff and gamma parameters?**


We didn't pass any parameters, but :meth:`~GNM.buildKirchhoff` method accepts two of them
by default ``cutoff=10.0`` and ``gamma=1.0`` is passed
that is, ``buildKirchhoff(calphas, cutoff=10., gamma=1.)`` 
 
.. plot::
   :context:
   :nofigs:
   :include-source:

   gnm.getCutoff()
   # this prints 15.0
   gnm.getGamma()
   # this prints 1.0

Note that it is also possible to use an externally calculated Kirchhoff 
matrix. Just pass it to the GNM instance using ``setKirchhoff`` method.

**Calculate normal modes**

.. plot::
   :context:
   :nofigs:
   :include-source:

   # calculate modes (by default slowest 20 will be calculated)
   gnm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate different number
of non-zero modes or to keep zero modes, try ``gnm.calcModes(50, zeros=True)``

Access calculated data
===============================================================================

::

   # get covariance matrix
   # note that covariance calculated using 20 modes when you call this method
   gnm.getCovariance()

   gnm.getEigenvalues()
   gnm.getEigenvectors()
   
.. note::
   Covariance matrices are calculated using available modes in the model. 
   If user calculated M slowest modes, only they will be used in the 
   calculation of covariance.

Inspect individual modes
===============================================================================

::

   slowest_mode = gnm[0]
   print slowest_mode.getEigenvalue()
   print slowest_mode.getEigenvector()

Normal mode indices start from 0, so slowest mode has index 0. By default,
modes with 0 eigenvalue are excluded. If they were retained, slowest 
non-trivial mode would have index 6.

Ploting results
===============================================================================

ProDy plotting functions are prefixed with ``show``. Let's use some of them
to plot data:

**Contact Map**

.. plot::
   :context:
   :include-source:
   
   # We import plotting library
   import matplotlib.pyplot as plt
   
   plt.figure(figsize=(5,4))
   showContactMap(gnm)
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  

**Cross-correlations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showCrossCorrelations(gnm)
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  

**Slow mode shape**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(5,4))
   showMode(gnm[0])
   plt.grid()
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  

**Square fluctuations**

.. plot::
   :context:
   :include-source:
   
   plt.figure(figsize=(7,4))
   showSqFlucts(gnm[0])
   
.. plot::
   :context:
   :nofigs:

   plt.close('all')  
   import os
   os.remove('2k39.pdb.gz')
