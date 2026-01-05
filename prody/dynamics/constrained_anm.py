# -*- coding: utf-8 -*-
"""Backward compatibility shim for constrained_anm module.

This module has been renamed to generalized_anm. The old names are kept for
backward compatibility but will issue a DeprecationWarning.

Please update your code to use:
  - generalized_anm.genANM instead of constrained_anm.cANM
  - generalized_anm.calcGenANM instead of constrained_anm.calcCANM
"""

import warnings

# Issue deprecation warning
warnings.warn(
    "The 'constrained_anm' module has been renamed to 'generalized_anm'. "
    "Please update your imports:\n"
    "  - Use 'from prody.dynamics.generalized_anm import genANM' instead of 'from prody.dynamics.constrained_anm import cANM'\n"
    "  - Use 'calcGenANM' instead of 'calcCANM'\n"
    "The old names will be removed in a future version.",
    DeprecationWarning,
    stacklevel=2
)

# Import and re-export with old names
from .generalized_anm import genANM as cANM
from .generalized_anm import calcGenANM as calcCANM
from .generalized_anm import build_generalized_hessian as build_constrained_hessian

__all__ = ['cANM', 'calcCANM', 'build_constrained_hessian']
