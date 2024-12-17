# -*- coding: utf-8 -*-
"""This module defines classes and functions for protein dynamics analysis.

Dynamics Models
===============

The following classes are designed for modeling and analysis of protein dynamics:

  * :class:`.ANM` - Anisotropic network model, for coarse-grained NMA
  * :class:`.GNM` - Gaussian network model, for coarse-grained dynamics
    analysis
  * :class:`.PCA` - Principal component analysis of conformation ensembles
  * :class:`.EDA` - Essential dynamics analysis of dynamics trajectories
  * :class:`.NMA` - Normal mode analysis, for analyzing data from external
    programs
  * :class:`.RTB` - Rotations and Translation of Blocks method

Usage of these classes are shown in :ref:`anm`, :ref:`gnm`, :ref:`pca`, and
:ref:`eda` examples.

The following classes are for analysis of individual modes or subsets of modes:

  * :class:`.Mode` - analyze individual normal/principal/essential modes
  * :class:`.ModeSet` - analyze subset of modes from a dynamics model
  * :class:`.Vector` - analyze modified modes or deformation vectors

Customize ENMs
==============

The following classes allow for using structure or distance based, or other custom
force constants and cutoff distances in :class:`~.ANM` and :class:`~.GNM`
calculations:

  * :class:`.Gamma` - base class for developing property custom force constant
    calculation methods
  * :class:`.GammaStructureBased` - secondary structure based force constants
  * :class:`.GammaVariableCutoff` - atom type based variable cutoff function

Membrane Models
===============

The following classes are designed for modeling and analysis of protein dynamics in membranes:

  * :class:`.exANM` - perform explicit membrane ANM, creating an explicit membrane lattice
  * :class:`.exGNM` - perform explicit membrane GNM, creating an explicit membrane lattice
  * :func:`.imANM` - perform implicit membrane GNM, creating anisotropic constraints in the membrane domain using RTB

Usage of these classes are shown in :ref:`exanm` and :ref:`imanm` examples.

Signature Dynamics (SignDy)
===============

The following classes are designed for signature dynamics analysis of protein/domain families, 
together with those in the database module:

  * :func:`.calcEnsembleENMs` - perform NMA on a protein family ensemble using ENMs
  * :class:`.ModeEnsemble` - handle outputs of ensemble NMA, an ensemble of normal modes
  * :class:`.sdarray` - handle signature dynamics data in an array based on a numpy array

There are many other functions starting `showSignature` or `calcSignature` for plotting and analysis.
There are also load and save functions for mode ensembles and signature arrays.

Usage of these classes are shown in :ref:`signdy-overview`, :ref:`signdy-core`, and :ref:`signdy-class` examples.


Function library
================

Dynamics of the functions in this library accept a *modes* argument (may also
appear in different names), which may refer to one or more of the following:

  * a dynamics model, :class:`.ANM`, :class:`SM`, :class:`.GNM`, :class:`.NMA`,
    :class:`.PCA`, or :class:`.EDA`
  * a :class:`.Mode` obtained by indexing an NMA model, e.g. ``anm[0]``
  * a :class:`.ModeSet` obtained by slicing an NMA model, e.g. ``anm[0:10]``

Some of these functions may also accept :class:`Vector` instances as *mode*
argument.  These are noted in function documentations.


Analyze models
==============

The following functions are for calculating atomic properties from normal modes:

  * :func:`.calcCollectivity` - degree of collectivity of a mode
  * :func:`.calcCovariance` - covariance matrix for given modes
  * :func:`.calcCrossCorr` - cross-correlations of fluctuations
  * :func:`.calcFractVariance` - fraction of variance explained by a mode
  * :func:`.calcPerturbResponse` - response to perturbations in positions
  * :func:`.calcProjection` - projection of conformations onto modes
  * :func:`.calcSqFlucts` - square-fluctuations
  * :func:`.calcTempFactors` - temperature factors fitted to exp. data

Compare models
==============

The following functions are for comparing normal modes or dynamics models:

  * :func:`.calcOverlap` - overlap (correlation) between modes
  * :func:`.calcCumulOverlap` - cumulative overlap between modes
  * :func:`.calcSubspaceOverlap` - overlap between normal mode subspaces
  * :func:`.calcCovOverlap` - covariance overlap between models
  * :func:`.printOverlapTable` - formatted overlap table printed on screen

Generate conformers
===================

The following functions can be used to generate conformers along normal modes:

  * :func:`.deformAtoms` - deform atoms along a mode
  * :func:`.sampleModes` - deform along random combination of a set of modes
  * :func:`.traverseMode` - traverse a mode along both directions

Adaptive ANM
===================

The following class and its functions can be used to generate conformers using adaptive ANM:

  * :class:`.AdaptiveANM` - generate transitions between two conformers using best overlapping modes 

ENM-MD hybrid methods
===================

The following classes and functions can be used to generate conformers using ENM-MD hybrid methods:

  * :class:`.ClustENM` - generate conformers by exploring combinations of modes and clustering, using minimisation and optionally MD
  * :function:`.ANMD` - generate conformers by traversing a number of individual modes, applying minimisation to each conformer

Essential Site Scanning Analysis (ESSA)
========================================

The following class and its functions can be used to perform Essential Site Scanning Analysis:

  * :class:`.ESSA`

Editing models
==============

The following functions can be used to reduce, slice, or extrapolate models:

  * :func:`.sliceMode` - take a slice of the normal mode
  * :func:`.extendMode` - extend a coarse-grained mode to all-atoms
  * :func:`.sliceModel` - take a slice of a model
  * :func:`.extendModel` - extend a coarse-grained model to all-atoms
  * :func:`.reduceModel` - reduce a model to a subset of atoms
  * :func:`.sliceVector` - take a slice of a vector
  * :func:`.extendVector` - extend a coarse-grained vector to all-atoms

Parse/write data
================

The following functions are parsing or writing normal mode data:

  * :func:`.parseArray` - numeric arrays, e.g. coordinates, eigenvectors
  * :func:`.parseModes` - normal modes
  * :func:`.parseNMD` - normal mode, coordinate, and atomic data for NMWiz
  * :func:`.parseSparseMatrix` - matrix data in sparse coordinate list format
  * :func:`.writeArray` - numeric arrays, e.g. coordinates, eigenvectors
  * :func:`.writeModes` - normal modes
  * :func:`.writeNMD` - normal mode, coordinate, and atomic data
  * :func:`.writeOverlapTable` - overlap between modes in a formatted table
  * :func:`.writeBILD` - normal mode and coordinate data for ChimeraX

Save/load models
================

Dynamics objects can be efficiently saved and loaded in later Python sessions
using the following functions:

  * :func:`.loadModel`, :func:`.saveModel` - load/save dynamics models
  * :func:`.loadVector`, :func:`.saveVector` - load/save modes or vectors


Short-hand functions
====================

Following allow for performing some dynamics calculations in one function call:

  * :func:`.calcANM` - perform ANM calculations
  * :func:`.calcGNM` - perform GNM calculations
  * :func:`.calcSM` - perform GNM calculations

Plotting functions
==================

Plotting functions are called by the name of the plotted data/property
and are prefixed with "show".  Function documentations refers to the
:mod:`matplotlib.pyplot` function utilized for actual plotting.
Arguments and keyword arguments are passed to the Matplotlib functions.


  * :func:`.showMode` - mode shape
  * :func:`.showOverlap` - overlap between modes
  * :func:`.showSqFlucts` - square-fluctuations
  * :func:`.showEllipsoid` - depict projection of a normal mode space on
    another
  * :func:`.showContactMap` - contact map based on a Kirchhoff matrix
  * :func:`.showProjection` - projection of conformations onto normal modes
  * :func:`.showOverlapTable` - overlaps between two models
  * :func:`.showScaledSqFlucts` - square-fluctuations fitted to experimental
    data
  * :func:`.showNormedSqFlucts` - normalized square-fluctuations
  * :func:`.showCrossProjection` - project conformations onto modes from
    different models
  * :func:`.showCrossCorr` - cross-correlations between fluctuations
    in atomic positions
  * :func:`.showCumulOverlap` - cumulative overlap of a mode with multiple
    modes from another model
  * :func:`.showFractVars` - fraction of variances
  * :func:`.showCumulFractVars` - cumulative fraction of variances
  * :func:`.resetTicks` - change ticks in a plot


Heat Mapper support
===================

The following functions can be used to read, write, and plot VMD plugin
`Heat Mapper`_ files.

  * :func:`.showHeatmap`
  * :func:`.parseHeatmap`
  * :func:`.writeHeatmap`

.. _Heat Mapper: http://www.ks.uiuc.edu/Research/vmd/plugins/heatmapper/


Visualize modes
===============

Finally, normal modes can be visualized and animated using VMD plugin
:ref:`nmwiz`. The following functions allow for running NMWiz from within Python:

  * :func:`.viewNMDinVMD` - run VMD and load normal mode data
  * :func:`.pathVMD` - get/set path to VMD executable
  
"""

__all__ = []

from . import entropy
from .entropy import *
__all__.extend(entropy.__all__)


from . import analysis
from .analysis import *
__all__.extend(analysis.__all__)


from . import compare
from .compare import *
__all__.extend(compare.__all__)

from . import editing
from .editing import *
__all__.extend(editing.__all__)

from . import functions
from .functions import *
__all__.extend(functions.__all__)

from . import perturb
from .perturb import *
__all__.extend(perturb.__all__)

from . import plotting
from .plotting import *
__all__.extend(plotting.__all__)

from . import sampling
from .sampling import *
__all__.extend(sampling.__all__)

from . import nma
from .nma import *
__all__.extend(nma.__all__)

from . import mode
from .mode import *
__all__.extend(mode.__all__)

mode.calcCollectivity = calcCollectivity

from . import modeset
from .modeset import *
__all__.extend(modeset.__all__)

nma.ModeSet = ModeSet

from . import gamma
from .gamma import *
__all__.extend(gamma.__all__)

from . import anm
from .anm import *
__all__.extend(anm.__all__)

from . import rtb
from .rtb import *
__all__.extend(rtb.__all__)

from . import gnm
from .gnm import *
__all__.extend(gnm.__all__)

from . import pca
from .pca import *
__all__.extend(pca.__all__)

from . import heatmapper
from .heatmapper import *
__all__.extend(heatmapper.__all__)

from . import mechstiff
from .mechstiff import *
__all__.extend(mechstiff.__all__)

from . import nmdfile
from .nmdfile import *
__all__.extend(nmdfile.__all__)

from . import vmdfile
from .vmdfile import *
__all__.extend(vmdfile.__all__)

from . import bildfile
from .bildfile import *
__all__.extend(bildfile.__all__)

#from . import bbenm
#from .bbenm import *
#__all__.extend(bbenm.__all__)

from . import exanm
from .exanm import *
__all__.extend(exanm.__all__)

from . import imanm
from .imanm import *
__all__.extend(imanm.__all__)

from . import exgnm
from .exgnm import *
__all__.extend(exgnm.__all__)

#from . import saxs
#from .saxs import *
#__all__.extend(saxs.__all__)

from . import signature
from .signature import *
__all__.extend(signature.__all__)

from . import adaptive
from .adaptive import *
__all__.extend(adaptive.__all__)

from . import clustenm
from .clustenm import *
__all__.extend(clustenm.__all__)

from . import essa
from .essa import *
__all__.extend(essa.__all__)

# workaround for circular dependency to accommodate original design style 
from prody.ensemble import functions
functions.ClustENM = ClustENM

from . import lda
from .lda import *
__all__.extend(lda.__all__)

from . import logistic
from .logistic import *
__all__.extend(logistic.__all__)

from . import anmd
from .anmd import *
__all__.extend(anmd.__all__)