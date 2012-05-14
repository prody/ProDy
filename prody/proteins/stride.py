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

"""This module defines functions for executing STRIDE program and parsing
its output."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

from prody.atomic import ATOMIC_FIELDS
from prody.atomic import AtomGroup
from prody.utilities import gunzip, which, PLATFORM

from pdbfile import parsePDB
from wwpdbftp import fetchPDB

__all__ = ['execSTRIDE', 'parseSTRIDE', 'performSTRIDE']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

def execSTRIDE(pdb, outputname=None, outputdir=None):
    """Execute STRIDE program for given *pdb*.  *pdb* can be an identifier or 
    a PDB file path.  If *pdb* is a compressed file, it will be decompressed 
    using Python :mod:`gzip` library.  When no *outputname* is given, output 
    name will be :file:`pdb.stride`.  :file:`.stride` extension will be 
    appended automatically to *outputname*.  If :file:`outputdir` is given, 
    STRIDE output and uncompressed PDB file will be written into this folder.
    Upon successful execution of :command:`stride pdb > out` command, output
    filename is returned. 
    
    For more information on STRIDE see http://webclu.bio.wzw.tum.de/stride/.
    If you benefited from STRIDE, please consider citing [DF95]_."""
    
    stride = which('stride')
    if stride is None:
        raise EnvironmentError('command not found: stride executable is not '
                               'found in one of system paths')
    assert outputname is None or isinstance(outputname, str),\
        'outputname must be a string'
    assert outputdir is None or isinstance(outputdir, str),\
        'outputdir must be a string'
    if not os.path.isfile(pdb):
        pdb = fetchPDB(pdb, compressed=False)
    if pdb is None:
        raise ValueError('pdb is not a valid PDB identifier or filename')
    if os.path.splitext(pdb)[1] == '.gz':
        if outputdir is None:
            pdb = gunzip(pdb, os.path.splitext(pdb)[0])
        else:
            pdb = gunzip(pdb, os.path.join(outputdir, 
                                os.path.split(os.path.splitext(pdb)[0])[1]))
    if outputdir is None:
        outputdir = '.'
    if outputname is None:
        out = os.path.join(outputdir,
                        os.path.splitext(os.path.split(pdb)[1])[0] + '.stride')
    else:
        out = os.path.join(outputdir, outputname + '.stride')
        
    status = os.system('{0:s} {1:s} > {2:s}'.format(stride, pdb, out))
    if status == 0:
        return out
    
def parseSTRIDE(stride, ag):
    """Parse STRIDE output from file *stride* into :class:`~.AtomGroup` 
    instance *ag*.  STRIDE output file must be in the new format used 
    from July 1995 and onwards.  When *stride* file is parsed, following 
    attributes are added to *ag*:
        
    * *stride_resnum*: STRIDE's sequential residue number, starting at the 
      first residue actually in the data set.
    
    * *stride_phi*, *stride_psi*: peptide backbone torsion angles phi and psi
    
    * *stride_area*: residue solvent accessible area"""
    
    if not os.path.isfile(stride):
        raise IOError('{0:s} is not a valid file path'.format(stride))
    if not isinstance(ag, AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')
        
    stride = open(stride)
    
    n_atoms = ag.numAtoms()
    NUMBER = np.zeros(n_atoms, int)
    AREA = np.zeros(n_atoms, float)
    PHI = np.zeros(n_atoms, float)
    PSI = np.zeros(n_atoms, float)

    ag.setSecstrs(np.zeros(n_atoms), dtype=ATOMIC_FIELDS['secondary'].dtype)
    for line in stride:
        if not line.startswith('ASG '):
            continue
        res = ag[(line[9], int(line[10:15]), line[15].strip())]
        if res is None:
            continue
        indices = res.getIndices()
        res.setSecstrs(line[24].strip())
        NUMBER[indices] = int(line[16:20])
        PHI[indices] = float(line[42:49])
        PSI[indices] = float(line[52:59])
        AREA[indices] = float(line[64:69])
    ag.setData('stride_resnum', NUMBER)
    ag.setData('stride_phi', PHI)
    ag.setData('stride_psi', PSI)
    ag.setData('stride_area', AREA)
    return ag

def performSTRIDE(pdb):
    """Perform STRIDE calculations and parse results.  STRIDE data is 
    returned in an :class:`~.AtomGroup` instance.  See also 
    :func:`execSTRIDE` and :func:`parseSTRIDE`."""
    
    pdb = fetchPDB(pdb, compressed=False)
    return parseSTRIDE(execSTRIDE(pdb), parsePDB(pdb))

