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

"""This module contains features for analyzing protein sequences.

Classes
=======

  * :class:`.MSA` - store MSA data indexed by label
  * :class:`.Sequence` - store sequence data


MSA IO
======

  * :class:`.MSAFile` - read/write MSA files in FASTA/SELEX/Stockholm formats
  * :func:`.parseMSA` - parse MSA files
  * :func:`.writeMSA` - parse MSA files

Editing
========

  * :func:`.mergeMSA` - merge MSA data for multi-domain proteins
  * :func:`.refineMSA` - refine MSA by removing gapped columns and/or sequences

Analysis
========

  * :func:`.calcMSAOccupancy` - calculate row (sequence) or column occupancy
  * :func:`.calcShannonEntropy` - calculate Shannon entropy
  * :func:`.buildMutinfoMatrix` - build mutual information matrix
  * :func:`.buildOMESMatrix` - build mutual observed minus expected squared
    covariance matrix
  * :func:`.buildSCAMatrix`- build statistical coupling analysis matrix
  * :func:`.buildSeqidMatrix`- build sequence identity matrix
  * :func:`.uniqueSequences` - select unique sequences
  * :func:`.applyMutinfoCorr` - apply correction to mutual information matrix
  * :func:`.applyMutinfoNorm` - apply normalization to mutual information
    matrix
  * :func:`.calcRankorder` - rank order scores


Plotting
========

  * :func:`.showShannonEntropy` - plot Shannon entropy
  * :func:`.showMSAOccupancy` - plot row (sequence) or column occupancy
  * :func:`.showMutinfoMatrix` - show mutual information matrix
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

__all__ = []

from . import msa
from .msa import *
__all__.extend(msa.__all__)

from . import msafile
from .msafile import *
__all__.extend(msafile.__all__)

from . import analysis
from .analysis import *
__all__.extend(analysis.__all__)

from . import plotting
from .plotting import *
__all__.extend(plotting.__all__)

from . import sequence
from .sequence import *
__all__.extend(sequence.__all__)

