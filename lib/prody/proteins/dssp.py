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

"""This module defines functions for executing DSSP program and parsing
its output."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

from prody import LOGGER
from prody.atomic import ATOMIC_FIELDS
from prody.atomic import AtomGroup
from prody.utilities import gunzip, which, PLATFORM

from .pdbfile import parsePDB
from .localpdb import fetchPDB

__all__ = ['execDSSP', 'parseDSSP', 'performDSSP']

def execDSSP(pdb, outputname=None, outputdir=None, stderr=True):
    """Execute DSSP for given *pdb*.  *pdb* can be a PDB identifier or a PDB 
    file path.  If *pdb* is a compressed file, it will be decompressed using
    Python :mod:`gzip` library.  When no *outputname* is given, output name 
    will be :file:`pdb.dssp`.  :file:`.dssp` extension will be appended 
    automatically to *outputname*.  If :file:`outputdir` is given, DSSP 
    output and uncompressed PDB file will be written into this folder.
    Upon successful execution of :command:`dssp pdb > out` command, output
    filename is returned.  On Linux platforms, when *stderr* is false, 
    standard error messages are suppressed, i.e.
    ``dssp pdb > outputname 2> /dev/null``.
    
    For more information on DSSP see http://swift.cmbi.ru.nl/gv/dssp/.
    If you benefited from DSSP, please consider citing [WK83]_."""
    
    dssp = which('dssp')
    if dssp is None:
        raise EnvironmentError('command not found: dssp executable is not '
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
                        os.path.splitext(os.path.split(pdb)[1])[0] + '.dssp')
    else:
        out = os.path.join(outputdir, outputname + '.dssp')
        
    if not stderr and PLATFORM != 'Windows':
        status = os.system('{0} {1} > {2} 2> /dev/null'.format(
                            dssp, pdb, out))
    else:
        status = os.system('{0} {1} > {2}'.format(dssp, pdb, out))

    if status == 0:
        return out
    
def parseDSSP(dssp, ag, parseall=False):
    """Parse DSSP data from file *dssp* into :class:`~.AtomGroup` instance 
    *ag*.  DSSP output file must be in the new format used from July 1995 
    and onwards.  When *dssp* file is parsed, following attributes are added 
    to *ag*:
        
    * *dssp_resnum*: DSSP's sequential residue number, starting at the first 
      residue actually in the data set and including chain breaks; this number 
      is used to refer to residues throughout.
    
    * *dssp_acc*: number of water molecules in contact with this residue \*10. 
      or residue water exposed surface in Angstrom^2.
      
    * *dssp_kappa*: virtual bond angle (bend angle) defined by the three Cα 
      atoms of residues I-2,I,I+2.  Used to define bend (structure code 'S').
        
    * *dssp_alpha*: virtual torsion angle (dihedral angle) defined by the four 
      Cα atoms of residues I-1,I,I+1,I+2.Used to define chirality (structure 
      code '+' or '-').
        
    * *dssp_phi* and *dssp_psi*: IUPAC peptide backbone torsion angles 

    The following attributes are parsed when ``parseall=True`` is passed: 

    * *dssp_bp1*, *dssp_bp2*, and *dssp_sheet_label*: residue number of first 
      and second bridge partner followed by one letter sheet label
      
    * *dssp_tco*: cosine of angle between C=O of residue I and C=O of residue 
      I-1.  For α-helices, TCO is near +1, for β-sheets TCO is near -1.  Not 
      used for structure definition.
      
    * *dssp_NH_O_1_index*, *dssp_NH_O_1_energy*, etc.: hydrogen bonds; e.g. 
      -3,-1.4 means: if this residue is residue i then N-H of I is h-bonded to 
      C=O of I-3 with an electrostatic H-bond energy of -1.4 kcal/mol.  There 
      are two columns for each type of H-bond, to allow for bifurcated H-bonds.
        
    See http://swift.cmbi.ru.nl/gv/dssp/DSSP_3.html for details."""
    
    if not os.path.isfile(dssp):
        raise IOError('{0} is not a valid file path'.format(dssp))
    if not isinstance(ag, AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')
        
    dssp = open(dssp)
    
    n_atoms = ag.numAtoms()
    NUMBER = np.zeros(n_atoms, int)
    SHEETLABEL = np.zeros(n_atoms, np.array(['a']).dtype.char + '1')
    ACC = np.zeros(n_atoms, float)
    KAPPA = np.zeros(n_atoms, float)
    ALPHA = np.zeros(n_atoms, float)
    PHI = np.zeros(n_atoms, float)
    PSI = np.zeros(n_atoms, float)

    if parseall:
        BP1 = np.zeros(n_atoms, int)
        BP2 = np.zeros(n_atoms, int)
        NH_O_1 = np.zeros(n_atoms, int)
        NH_O_1_nrg = np.zeros(n_atoms, float)
        O_HN_1 = np.zeros(n_atoms, int)
        O_HN_1_nrg = np.zeros(n_atoms, float)
        NH_O_2 = np.zeros(n_atoms, int)
        NH_O_2_nrg = np.zeros(n_atoms, float)
        O_HN_2 = np.zeros(n_atoms, int)
        O_HN_2_nrg = np.zeros(n_atoms, float)
        TCO = np.zeros(n_atoms, float)

    ag.setSecstrs(np.zeros(n_atoms, dtype=ATOMIC_FIELDS['secondary'].dtype))
    for line in dssp:
        if line.startswith('  #  RESIDUE'):
            break
    for line in dssp:
        if line[13] == '!':
            continue
        res = ag[(line[11], int(line[5:10]), line[10].strip())]
        if res is None:
            continue
        indices = res.getIndices()
        res.setSecstrs(line[16].strip())
        NUMBER[indices] = int(line[:5])
        SHEETLABEL[indices] = line[33].strip()
        ACC[indices] = int(line[35:38])
        KAPPA[indices] = float(line[91:97])
        ALPHA[indices] = float(line[97:103])
        PHI[indices] = float(line[103:109])
        PSI[indices] = float(line[109:115])

        if parseall:
            BP1[indices] = int(line[25:29])
            BP2[indices] = int(line[29:33])
            NH_O_1[indices] = int(line[38:45])
            NH_O_1_nrg[indices] = float(line[46:50]) 
            O_HN_1[indices] = int(line[50:56])
            O_HN_1_nrg[indices] = float(line[57:61])
            NH_O_2[indices] = int(line[61:67])
            NH_O_2_nrg[indices] = float(line[68:72])
            O_HN_2[indices] = int(line[72:78])
            O_HN_2_nrg[indices] = float(line[79:83])
            TCO[indices] = float(line[85:91])
    
    ag.setData('dssp_resnum', NUMBER)
    ag.setData('dssp_sheet_label', SHEETLABEL)
    ag.setData('dssp_acc', ACC)
    ag.setData('dssp_kappa', KAPPA)
    ag.setData('dssp_alpha', ALPHA)
    ag.setData('dssp_phi', PHI)
    ag.setData('dssp_psi', PSI)

    if parseall:
        ag.setData('dssp_bp1', BP1)
        ag.setData('dssp_bp2', BP2)
        ag.setData('dssp_NH_O_1_index', NH_O_1)
        ag.setData('dssp_NH_O_1_energy', NH_O_1_nrg)
        ag.setData('dssp_O_NH_1_index', O_HN_1)
        ag.setData('dssp_O_NH_1_energy', O_HN_1_nrg)    
        ag.setData('dssp_NH_O_2_index', NH_O_2)
        ag.setData('dssp_NH_O_2_energy', NH_O_2_nrg)
        ag.setData('dssp_O_NH_2_index', O_HN_2)
        ag.setData('dssp_O_NH_2_energy', O_HN_2_nrg)
        ag.setData('dssp_tco', TCO)
    return ag

def performDSSP(pdb, parseall=False, stderr=True):
    """Perform DSSP calculations and parse results.  DSSP data is returned 
    in an :class:`~.AtomGroup` instance.  See also :func:`execDSSP` 
    and :func:`parseDSSP`."""
    
    pdb = fetchPDB(pdb, compressed=False)
    return parseDSSP(execDSSP(pdb, stderr=stderr), parsePDB(pdb), parseall)
