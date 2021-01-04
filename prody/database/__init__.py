# -*- coding: utf-8 -*-
"""This module contains features for accessing databases containing protein
related data.

Pfam
====

The following functions can be used to search and retrieve Pfam_ data:

  * :func:`.fetchPfamMSA` - download MSA files
  * :func:`.searchPfam` - search for domain families of a protein


.. _Pfam: http://pfam.sanger.ac.uk/

UniProt
========

The following functions and class can be used to search and retrieve UniProt_ data:

  * :func:`.queryUniprot` - query UniProt and parse the results as a dictionary
  * :class:`.UniprotRecord` - a wrapper from UniProt data with functions including parsing PDBs
  * :func:`.searchUniprot` - search UniProt and return a UniprotRecord


.. _UniProt: https://www.uniprot.org/

CATH
====

The following class and its functions can be used to search and retrieve CATH_ data:

  * :class:`.CATHDB` - parse, handle and navigate the tree-like structure of the CATH database

.. _CATH: http://download.cathdb.info

DALI
====

The following class and functions can be used to search and retrieve data using the DALI_ structure alignment server:

  * :func:`.searchDali` - search for similar structures using DALI
  * :class:`.DaliRecord` - fetch and handle outputs from DALI searches
  * :func:`.daliFilterMultimers` - filter DALI results to obtain multimers of a particular size

.. _DALI: http://ekhidna2.biocenter.helsinki.fi/dali/

QuartataWeb
============

The following classes and functions can be used to search and retrieve data using the QuartataWeb_ structure alignment server:

  * :class:`.QuartataWebBrowser` - class based on the Splinter web browser package to search QuartataWeb
  * :class:`.QuartataChemicalRecord` - class to handle the outputs of QuartataWeb searches
  * :func:`.searchQuartataWeb` - perform QuartataWeb searches and return the output in a QuartataChemicalRecord

.. _QuartataWeb: http://quartata.csb.pitt.edu

Gene Ontology Annotation (GOA)
================================

The following classes and functions can be used to search and retrieve data from the EBI GOA_ database:

  * :func:`.queryGOA` - query GOA using a PDB ID
  * :class:`.GOADictList` - class to handle data from GOA queries
  * :func:`.parseOBO` - parse an OBO file containing the Gene Ontology.
  * :func:`.parseGAF` - parse a Gene Association File (GAF)
  * :func:`.showGoLineage` - visualize GO tree
  * :func:`.calcGoOverlap` - Calculate overlap between GO terms from their distance in the graph

.. _GOA: https://www.ebi.ac.uk/GOA/


"""

__all__ = []

from . import pfam
from .pfam import *
__all__.extend(pfam.__all__)

from . import uniprot
from .uniprot import *
__all__.extend(uniprot.__all__)

from . import cath
from .cath import *
__all__.extend(cath.__all__)

from . import dali
from .dali import *
__all__.extend(dali.__all__)

from . import goa
from .goa import *
__all__.extend(goa.__all__)

from . import quartataweb
from .quartataweb import *
__all__.extend(quartataweb.__all__)
