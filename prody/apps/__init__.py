"""This module contains ProDy applications.

Dynamics analysis
===============================================================================

  * :func:`.prody_anm`
  * :func:`.prody_gnm`
  * :func:`.prody_pca`

Structure analysis
===============================================================================

  * :func:`.prody_align`
  * :func:`.prody_biomol`
  * :func:`.prody_blast`
  * :func:`.prody_catdcd`
  * :func:`.prody_contacts`
  * :func:`.prody_fetch`
  * :func:`.prody_select`

Sequence analysis
===============================================================================

  * :func:`.evol_search`
  * :func:`.evol_fetch`
  * :func:`.evol_filter`
  * :func:`.evol_refine`
  * :func:`.evol_merge`
  * :func:`.evol_conserv`
  * :func:`.evol_coevol`
  * :func:`.evol_occupancy`
  * :func:`.evol_rankorder`
"""

from .evol_apps import evol_main, evol_parser, EVOL_APPS
from .prody_apps import prody_main, prody_parser, PRODY_APPS
