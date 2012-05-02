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

"""This module defines classes measuring quantities, transforming coordinates,
and identifying contacts.


Identify contacts
=================

Following class and functions are for contact identifications:

  * :class:`~.Contacts` - identify intermolecular contacts
  * :func:`~.iterNeighbors` - yield interacting atom pairs

Measure quantities
==================

Following functions are for measuring simple quantities:

  * :func:`~.calcDistance` - calculate distance(s)
  * :func:`~.calcAngle` - calculate bond angle
  * :func:`~.calcDihedral` - calculate dihedral angle
  * :func:`~.calcOmega` - calculate omega (ω) angle
  * :func:`~.calcPhi` - calculate phi (φ) angle
  * :func:`~.calcPsi` - calculate psi (ψ) angle
  * :func:`~.calcGyradius` - calculate radius of gyration
  * :func:`~.calcCenter` - calculate geometric (or mass) center
  * :func:`~.calcDeformVector` - calculate deformation vector


Anisotropic factors
===================

Following functions handle anisotropic displacement parameter (ADP) present
in some X-ray structures.

  * :func:`~.buildADPMatrix` - build ADP matrix 
  * :func:`~.calcADPAxes` - calculate ADP axes
  * :func:`~.calcADPs` - calculate ADPs

Transformations
===============

Following class and functions are for handling coordinate transformations:
      
  * :class:`~.Transformation` - store transformation matrix
  * :func:`~.alignCoordsets` - align multiple coordinate sets
  * :func:`~.applyTransformation` - apply a transformation
  * :func:`~.calcTransformation` - calculate a transformation
  * :func:`~.calcRMSD` - calculate root-mean-square distance
  * :func:`~.superpose` - superpose atoms or coordinate sets
  * :func:`~.moveAtoms` - move atoms by given offset
"""

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = []

import measure
from measure import *
__all__.extend(measure.__all__)

from measure import getDihedral, getPhiAtoms
from measure import getAngle, getCenter, getCentral

import contacts
from contacts import *
__all__.extend(contacts.__all__)

import transform
from transform import *
__all__.extend(transform.__all__)

from transform import getRMSD, getTransformation
