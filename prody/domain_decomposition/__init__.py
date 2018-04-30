# -*- coding: utf-8 -*-
"""This module defines functions for computing structural/dynamical
domain decompositions, and related properties, from either ANM modes 
or analysis of structural ensembles.
"""

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = []

from . import rtb_domain
__all__.extend(rtb_domain.__all__)
from .rtb_domain import *

from . import spectrus 
__all__.extend(spectrus.__all__)
from .spectrus import *
