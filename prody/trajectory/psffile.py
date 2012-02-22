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

"""This module defines a function for parsing protein structure files in PSF 
format."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

from prody.atomic import ATOMIC_DATA_FIELDS
from prody.atomic import AtomGroup
from prody.tools import openFile

__all__ = ['parsePSF']

def parsePSF(filename, title=None, ag=None):
    """Return an :class:`~.AtomGroup` instance storing data 
    parsed from X-PLOR format PSF file *filename*.  Atom and bond information
    is parsed from the file.  If *title* is not given, *filename* will be set 
    as the title of the :class:`AtomGroup` instance.  An :class:`AtomGroup` 
    instance may be provided as *ag* argument.  When provided, *ag* must have 
    the same number of atoms in the same order as the file.  Data from PSF 
    file will be added to the *ag*.  This may overwrite present data if it 
    overlaps with PSF file content. Note that this function does not evaluate 
    angles, dihedrals, and impropers sections."""
    
    if ag is not None:
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance') 
    
    psf = openFile(filename)
    line = psf.readline()
    i_line = 1
    while line:
        line = line.strip()
        if line.endswith('!NATOM'):
            n_atoms = int(line.split('!')[0])
            break
        line = psf.readline()
        i_line += 1
    if title is None:
        title = os.path.splitext(os.path.split(filename)[1])[0]
    else:
        title = str(title)
    if ag is None:
        ag = AtomGroup(title)
    else:
        if n_atoms != ag.numAtoms():
            raise ValueError('ag and PSF file must have same number of atoms')
        
    serials = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['serial'].dtype)
    segnames = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['segment'].dtype)
    resnums = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['resnum'].dtype)
    resnames = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['resname'].dtype)
    atomnames = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['name'].dtype)
    atomtypes = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['type'].dtype)
    charges = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['charge'].dtype)
    masses = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['mass'].dtype)
    
    lines = psf.readlines(71 * (n_atoms + 5))
    if len(lines) < n_atoms:
        raise IOError('number of lines in PSF is less than the number of '
                      'atoms')
    
    for i, line in enumerate(lines):
        if i == n_atoms:
            break
        i_line += 1
        if len(line) <= 71:
            serials[i] = line[:8]
            segnames[i] = line[9:13].strip()
            resnums[i] = line[14:19]
            resnames[i] = line[19:23].strip()
            atomnames[i] = line[24:28].strip()
            atomtypes[i] = line[29:35].strip()
            charges[i] = line[35:44]
            masses[i] = line[50:60]
        else:
            items = line.split()
            serials[i] = items[0]
            segnames[i] = items[1]
            resnums[i] = items[2]
            resnames[i] = items[3]
            atomnames[i] = items[4]
            atomtypes[i] = items[5]
            charges[i] = items[6]
            masses[i] = items[7]
    
    i = n_atoms
    while 1:
        line = lines[i].split()
        if len(line) >= 2 and line[1] == '!NBOND:':
             n_bonds = int(line[0])
             break
        i += 1
    lines = ''.join(lines[i+1:]) + psf.read(n_bonds/4 * 71)
    array = np.fromstring(lines, count=n_bonds*2, dtype=int, sep=' ')
    if len(array) != n_bonds*2:
        raise IOError('number of bonds expected and parsed do not match')

    psf.close()
    ag.setSerials(serials)
    ag.setSegnames(segnames)
    ag.setResnums(resnums)
    ag.setResnames(resnames)
    ag.setNames(atomnames)
    ag.setTypes(atomtypes)
    ag.setCharges(charges)
    ag.setMasses(masses)

    array = np.add(array, -1, array)
    ag.setBonds(array.reshape((n_bonds, 2)))

    return ag
