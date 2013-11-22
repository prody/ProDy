# -*- coding: utf-8 -*-
"""This module defines classes measuring quantities, transforming coordinates,
and identifying contacts.


Identify contacts
=================

Following class and functions are for contact identifications:

  * :class:`.Contacts` - identify intermolecular contacts
  * :func:`.findNeighbors` - identify interacting atom pairs
  * :func:`.iterNeighbors` - identify interacting atom pairs

Measure quantities
==================

Following functions are for measuring simple quantities:

  * :func:`.calcDistance` - calculate distance(s)
  * :func:`.calcAngle` - calculate bond angle
  * :func:`.calcDihedral` - calculate dihedral angle
  * :func:`.calcOmega` - calculate omega (ω) angle
  * :func:`.calcPhi` - calculate phi (φ) angle
  * :func:`.calcPsi` - calculate psi (ψ) angle
  * :func:`.calcGyradius` - calculate radius of gyration
  * :func:`.calcCenter` - calculate geometric (or mass) center
  * :func:`.calcDeformVector` - calculate deformation vector


Anisotropic factors
===================

Following functions handle anisotropic displacement parameter (ADP) present
in some X-ray structures.

  * :func:`.buildADPMatrix` - build ADP matrix
  * :func:`.calcADPAxes` - calculate ADP axes
  * :func:`.calcADPs` - calculate ADPs

Transformations
===============

Following class and functions are for handling coordinate transformations:

  * :class:`.Transformation` - store transformation matrix
  * :func:`.alignCoordsets` - align multiple coordinate sets
  * :func:`.applyTransformation` - apply a transformation
  * :func:`.calcTransformation` - calculate a transformation
  * :func:`.calcRMSD` - calculate root-mean-square distance
  * :func:`.superpose` - superpose atoms or coordinate sets
  * :func:`.moveAtoms` - move atoms by given offset
"""

__all__ = []

from . import measure
from .measure import *
__all__.extend(measure.__all__)

from .measure import getDihedral, getPhiAtoms
from .measure import getAngle, getCenter, getCentral

from . import contacts
from .contacts import *
__all__.extend(contacts.__all__)

from . import transform
from .transform import *
__all__.extend(transform.__all__)

from .transform import getRMSD, getTransformation
