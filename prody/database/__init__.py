# -*- coding: utf-8 -*-
"""This module contains features for accessing databases containing protein
related data.

Pfam
====

Following functions can be used to search and retrieve Pfam_ data:

  * :func:`.fetchPfamMSA` - download MSA files
  * :func:`.searchPfam` - search families of a protein


.. _Pfam: http://pfam.sanger.ac.uk/"""

__all__ = []

from . import pfam
from .pfam import *
__all__.extend(pfam.__all__)

from . import uniprot
from .uniprot import *
__all__.extend(uniprot.__all__)
