# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""Protein dynamics
===============================================================================

This module defines classes and functions for protein dynamics analysis. 

Modeling and analysis 
-------------------------------------------------------------------------------

Following classes are designed for modeling and analysis of protein dynamics:

==============  ===============================================================
Model           Description
==============  ===============================================================
:class:`~.ANM`  Anisotropic network model, for coarse-grained NMA
:class:`~.GNM`  Gaussian network model, for coarse-grained dynamics analysis
:class:`~.PCA`  Principal component analysis of conformation ensembles
:class:`~.EDA`  Essential dynamics analysis of dynamics trajectories
:class:`~.NMA`  Normal mode analysis, for analyzing data from external programs
==============  ===============================================================

Following classes are for supporting above models by allowing analysis of 
individual modes or subsets of modes:

==================  ===========================================================
Class               Description
==================  ===========================================================
:class:`~.Mode`     analyze individual normal/principal/essential modes
:class:`~.ModeSet`  analyze subset of modes from one of aforementioned models
:class:`~.Vector`   analyze modified modes or deformation vectors
==================  ===========================================================

Customized ENMs 
-------------------------------------------------------------------------------

Following classes allow for using structure or distance based, or other custom 
force constants and cutoff distances in :class:`~.ANM` and :class:`~.GNM` 
calculations: 

==============================  ===============================================
Class                           Description
==============================  ===============================================
:class:`~.Gamma`                base class for developing property custom force 
                                constant calculation methods
:class:`~.GammaStructureBased`  secondary structure based force constants
:class:`~.GammaVariableCutoff`  atom type based variable cutoff function
==============================  ===============================================

Function library
===============================================================================

Dynamics of the functions described below accept a *modes* argument (may also 
appear in different names), which may refer to one or more of the following:

  * a dynamics model, :class:`~.ANM`, :class:`~.GNM`, :class:`~.NMA`, 
    :class:`~.PCA`, or :class:`~.EDA` 
  * a :class:`~.Mode` obtained by indexing an NMA model, e.g. ``anm[0]``
  * a :class:`~.ModeSet` obtained by slicing an NMA model, e.g. ``anm[0:10]``

Some of these functions may also accept :class:`Vector` instances as *mode* 
argument.  These are noted in function documentations. 


Analyze models
-------------------------------------------------------------------------------

Following functions are for calculating atomic properties from normal modes:

=============================  ================================================
Function                       Calculated data or property
=============================  ================================================
:func:`~.calcCollectivity`     degree of collectivity of a mode
:func:`~.calcCovariance`       covariance matrix for given modes
:func:`~.calcCrossCorr`        cross-correlations of fluctuations
:func:`~.calcPerturbResponse`  response to perturbations in positions
:func:`~.calcProjection`       projection of conformations onto modes
:func:`~.calcSqFlucts`         square-fluctuations
:func:`~.calcTempFactors`      temperature factors fitted to exp. data
=============================  ================================================

Compare models
-------------------------------------------------------------------------------

Following functions are for comparing normal modes or dynamics models:

=============================  ================================================
Function                       Output
=============================  ================================================
:func:`~.calcOverlap`          overlap (correlation) between modes
:func:`~.calcCumOverlap`       cumulative overlap between modes
:func:`~.calcCumOverlapArray`  incremental cumulative overlap
:func:`~.calcSubspaceOverlap`  overlap between normal mode subspaces 
:func:`~.calcCovOverlap`       covariance overlap between models
:func:`~.printOverlapTable`    formatted overlap table printed on screen
=============================  ================================================

Generate conformers
-------------------------------------------------------------------------------

Following functions can be used to generate conformers along normal modes:

======================  =======================================================
Function                Method
======================  =======================================================
:func:`~.deformAtoms`   deform atoms along a mode
:func:`~.sampleModes`   deform along random combination of a set of modes 
:func:`~.traverseMode`  traverse a mode along both directions
======================  =======================================================

Editing models
-------------------------------------------------------------------------------

Following functions can be used to reduce, slice, or extrapolate models:

==========================  ===================================================
Function                    Description
==========================  ===================================================
:func:`~.sliceMode`         take a slice of the normal mode      
:func:`~.sliceModel`        take a slice of a model
:func:`~.sliceVector`       take a slice of a vector
:func:`~.reduceModel`       reduce a model to a subset of atoms
:func:`~.extrapolateModel`  extrapolate a coarse-grained model to all-atoms  
==========================  ===================================================

Parse/write data
-------------------------------------------------------------------------------

Following functions are parsing or writing normal mode data:

===========================  ==================================================
Function                     Input/output
===========================  ==================================================
:func:`~.parseArray`         numeric arrays, e.g. coordinates, eigenvectors
:func:`~.parseModes`         normal modes
:func:`~.parseNMD`           normal mode, coordinate, and atomic data for NMWiz
:func:`~.parseSparseMatrix`  matrix data in sparse coordinate list format
:func:`~.writeArray`         numeric arrays, e.g. coordinates, eigenvectors
:func:`~.writeModes`         normal modes
:func:`~.writeNMD`           normal mode, coordinate, and atomic data
:func:`~.writeOverlapTable`  overlap between modes in a formatted table
===========================  ==================================================

Save/load models 
-------------------------------------------------------------------------------

Dynamics objects can be efficiently saved and loaded in later Python sessions 
using the following functions:

| :func:`~.loadModel`, :func:`~.saveModel` - load/save dynamics models
| :func:`~.loadVector`, :func:`~.saveVector` - load/save modes or vectors
  

Short-hand functions
-------------------------------------------------------------------------------

Following allow for performing some dynamics calculations in one function call:

| :func:`~.calcANM` - perform ANM calculations
| :func:`~.calcGNM` - perform GNM calculations

Plotting functions
-------------------------------------------------------------------------------

Plotting functions are called by the name of the plotted data/property 
and are prefixed with "show".  Function documentations refers to the 
:mod:`matplotlib.pyplot` function utilized for actual plotting. 
Arguments and keyword arguments are passed to the Matplotlib functions.  


=============================  ================================================
Function                       Plotted data
=============================  ================================================
:func:`~.showMode`             mode shape
:func:`~.showOverlap`          overlap between modes
:func:`~.showSqFlucts`         square-fluctuations
:func:`~.showEllipsoid`        depict projection of a normal mode space on 
                               another 
:func:`~.showContactMap`       contact map based on a Kirchhoff matrix
:func:`~.showProjection`       projection of conformations onto normal modes
:func:`~.showOverlapTable`     overlaps between two models
:func:`~.showScaledSqFlucts`   square-fluctuations fitted to experimental data
:func:`~.showNormedSqFlucts`   normalized square-fluctuations
:func:`~.showCrossProjection`  project conformations onto modes from different 
                               models
:func:`~.showFractOfVar`       fraction of variances
:func:`~.showCrossCorr`        cross-correlations between fluctuations
                               in atomic positions
:func:`~.showCumOverlap`       cumulative overlap of a mode with multiple modes 
                               from another model
:func:`~.showCumFractOfVar`    cumulative fraction of variances 
:func:`~.resetTicks`           change ticks in a plot
=============================  ================================================

Visualize modes
-------------------------------------------------------------------------------

Finally, normal modes can be visualized and animated using VMD plugin 
:ref:`nmwiz`. Following functions allow for running NMWiz from within Python: 
 
| :func:`~.viewNMDinVMD` - run VMD and load normal mode data
| :func:`~.getVMDpath`, :func:`~.setVMDpath` - get/set path to VMD executable


Examples
-------------------------------------------------------------------------------

Results from the example :ref:`pca-xray-calculations` will be used to
illustrate class methods and functions in the module.
  

>>> from prody import *
>>> import matplotlib.pyplot as plt
>>> import numpy as np

>>> p38_pca = loadModel('p38_xray.pca.npz')
>>> p38_anm = loadModel('1p38.anm.npz') 
>>> p38_ensemble = loadEnsemble('p38_X-ray.ens.npz')
>>> p38_structure = parsePDB('p38_ref_chain.pdb')


.. plot::
   :nofigs: 
   :context: 
    
   from prody import *
   import matplotlib.pyplot as plt
   import numpy as np
   
   plt.close('all')
   
   p38_pca = loadModel('p38_xray.pca.npz')
   p38_anm = loadModel('1p38.anm.npz') 
   p38_ensemble = loadEnsemble('p38_X-ray.ens.npz')
   p38_structure = parsePDB('p38_ref_chain.pdb')

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = []

import analysis
from analysis import *
__all__.extend(analysis.__all__)

import compare
from compare import *
__all__.extend(compare.__all__)

import editing
from editing import *
__all__.extend(editing.__all__)

import functions
from functions import *
__all__.extend(functions.__all__)

import plotting
from plotting import *
__all__.extend(plotting.__all__)

import sampling
from sampling import *
__all__.extend(sampling.__all__)

import nma
from nma import *
__all__.extend(nma.__all__)

import mode
from mode import *
__all__.extend(mode.__all__)

mode.calcCollectivity = calcCollectivity

import modeset
from modeset import *
__all__.extend(modeset.__all__)

modeset.NMA = NMA

import gamma
from gamma import *
__all__.extend(gamma.__all__)

import anm
from anm import *
__all__.extend(anm.__all__)

import gnm
from gnm import *
__all__.extend(gnm.__all__)

import pca
from pca import *
__all__.extend(pca.__all__)

