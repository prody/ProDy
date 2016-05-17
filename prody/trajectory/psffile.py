# -*- coding: utf-8 -*-
"""This module defines a function for parsing protein structure files in
`PSF format`_.

.. _PSF format:
   http://www.ks.uiuc.edu/Training/Tutorials/namd/
   namd-tutorial-unix-html/node21.html"""

import os.path

from numpy import fromstring, zeros, ones, array, add

from prody import PY2K
from prody.atomic import ATOMIC_FIELDS, AtomGroup
from prody.utilities import openFile

if PY2K:
    range = xrange

__all__ = ['parsePSF', 'writePSF']

def parsePSF(filename, title=None, ag=None):
    """Return an :class:`.AtomGroup` instance storing data parsed from X-PLOR
    format PSF file *filename*.  Atom and bond information is parsed from the
    file.  If *title* is not given, *filename* will be set as the title of the
    :class:`.AtomGroup` instance.  An :class:`.AtomGroup` instance may be
    provided as *ag* argument.  When provided, *ag* must have the same number
    of atoms in the same order as the file.  Data from PSF file will be added
    to the *ag*.  This may overwrite present data if it overlaps with PSF file
    content.  Note that this function does not evaluate angles, dihedrals, and
    impropers sections."""

    if ag is not None:
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')

    psf = openFile(filename, 'rb')
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

    serials = zeros(n_atoms, ATOMIC_FIELDS['serial'].dtype)
    segnames = zeros(n_atoms, ATOMIC_FIELDS['segment'].dtype)
    resnums = zeros(n_atoms, ATOMIC_FIELDS['resnum'].dtype)
    resnames = zeros(n_atoms, ATOMIC_FIELDS['resname'].dtype)
    atomnames = zeros(n_atoms, ATOMIC_FIELDS['name'].dtype)
    atomtypes = zeros(n_atoms, ATOMIC_FIELDS['type'].dtype)
    charges = zeros(n_atoms, ATOMIC_FIELDS['charge'].dtype)
    masses = zeros(n_atoms, ATOMIC_FIELDS['mass'].dtype)
    
    #lines = psf.readlines(71 * (n_atoms + 5))
    i_line = 0
    n_bonds = 0
    for i, line in enumerate(psf):
        if line.strip() == '':
            continue
        if '!NBOND:' in line.upper():
            items = line.split()
            n_bonds = int(items[0])
            break
        if i_line + 1 > n_atoms:
            continue
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
    
    if i_line < n_atoms:
        raise IOError('number of lines in PSF is less than the number of '
                      'atoms')
                      
#    i = n_atoms
#    while 1:
#        line = lines[i].split()
#        if len(line) >= 2 and line[1] == '!NBOND:':
#             n_bonds = int(line[0])
#             break
#        i += 1
#    lines = ''.join(lines[i+1:]) + psf.read(n_bonds/4 * 71)
    lines = []
    for i, line in enumerate(psf):
        if line.strip() == '':
            continue
        if '!' in line:
            break
        lines.append(line)
    
    lines = ''.join(lines)
    array = fromstring(lines, count=n_bonds*2, dtype=int, sep=' ')
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

    array = add(array, -1, array)
    ag.setBonds(array.reshape((n_bonds, 2)))

    return ag

PSFLINE = ('%8d %-4s %-4d %-4s %-4s %-4s %10.6f %13.4f %11d\n')

def writePSF(filename, atoms):
    """Write atoms in X-PLOR format PSF file with name *filename* and return
    *filename*.  This function will write available atom and bond information
    only."""

    try:
        n_atoms, segments, rnums, rnames, names, types, charges, masses = (
        atoms.numAtoms(), atoms._getSegnames(), atoms._getResnums(),
        atoms._getResnames(), atoms._getNames(), atoms._getTypes(),
        atoms._getCharges(), atoms._getMasses())

    except AttributeError:
        raise TypeError('atoms must be an Atomic instance')

    if segments is None:
        segments = atoms._getChids()
        if segments is None:
            segments = ['UNK'] * n_atoms

    if rnums is None:
        rnums = ones(n_atoms, int)

    if rnames is None:
        rnames = ['UNK'] * n_atoms

    if names is None:
        raise ValueError('atom names are not set')

    if types is None:
        atomtypes = zeros(n_atoms, array(['a']).dtype.char + '1')

    long_fields = array([len(tp) for tp in types]).max() > 4

    out = openFile(filename, 'w')
    write = out.write
    write('PSF{0}\n'.format( ' NAMD' if long_fields else ''))
    write('\n')
    write('{0:8d} !NTITLE\n'.format(1))
    write(' REMARKS {0}\n'.format(str(atoms)))
    write('\n')
    write('{0:8d} !NATOM\n'.format(n_atoms))

    for i in range(n_atoms):
        write(PSFLINE % (i + 1, segments[i], rnums[i], rnames[i], names[i],
                        types[i], charges[i], masses[i], 0))
    bonds = list(atoms._iterBonds())
    if bonds:
        bonds = array(bonds, int) + 1
        write('\n')
        write('{0:8d} !NBOND: bonds\n'.format(len(bonds)))
        for i, bond in enumerate(bonds):
            write('%8s%8s' % (bond[0], bond[1]))
            if i % 4 == 3:
                write('\n')
        if i % 4 != 3:
            write('\n')
    out.close()
    return filename
