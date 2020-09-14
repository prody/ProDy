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
    """Returns an :class:`.AtomGroup` instance storing data parsed from X-PLOR
    format PSF file *filename*.  Atom and bond information is parsed from the
    file.  If *title* is not given, *filename* will be set as the title of the
    :class:`.AtomGroup` instance.  An :class:`.AtomGroup` instance may be
    provided as *ag* argument.  When provided, *ag* must have the same number
    of atoms in the same order as the file.  Data from PSF file will be added
    to the *ag*.  This may overwrite present data if it overlaps with PSF file
    content.  
    
    This function now includes the angles, dihedrals, and impropers sections
    as well as donors, acceptors and crossterms!"""

    if ag is not None:
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')

    psf = openFile(filename, 'rb')
    line = psf.readline()
    while line:
        line = line.strip()
        if line.endswith(b'!NATOM'):
            n_atoms = int(line.split(b'!')[0])
            break
        line = psf.readline()
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
    
    n = 0
    n_bonds = 0
    for line in psf:
        if line.strip() == b'':
            continue
        if b'!NBOND:' in line.upper():
            items = line.split()
            n_bonds = int(items[0])
            break
        if n + 1 > n_atoms:
            continue

        if len(line) <= 71:
            serials[n] = line[:8]
            segnames[n] = line[9:13].strip()
            resnums[n] = line[14:19]
            resnames[n] = line[19:23].strip()
            atomnames[n] = line[24:28].strip()
            atomtypes[n] = line[29:35].strip()
            charges[n] = line[35:44]
            masses[n] = line[50:60]
        else:
            items = line.split()
            serials[n] = items[0]
            segnames[n] = items[1]
            resnums[n] = items[2]
            resnames[n] = items[3]
            atomnames[n] = items[4]
            atomtypes[n] = items[5]
            charges[n] = items[6]
            masses[n] = items[7]
        n += 1
    
    if n < n_atoms:
        raise IOError('number of lines in PSF atoms block is less than the number of '
                      'atoms')

    n_angles = 0
    lines = []
    for i, line in enumerate(psf):
        if line.strip() == b'':
            continue
        if b'!NTHETA' in line:
            items = line.split()
            n_angles = int(items[0])
            break
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    b_array = fromstring(lines, count=n_bonds*2, dtype=int, sep=' ')
    if len(b_array) != n_bonds*2:
        raise IOError('number of bonds expected and parsed do not match')

    n_dihedrals = 0
    lines = []
    for i, line in enumerate(psf):
        if line.strip() == b'':
            continue
        if b'!NPHI' in line:
            items = line.split()
            n_dihedrals = int(items[0])
            break
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    a_array = fromstring(lines, count=n_angles*3, dtype=int, sep=' ')
    if len(a_array) != n_angles*3:
        raise IOError('number of angles expected and parsed do not match')

    n_impropers = 0
    lines = []
    for i, line in enumerate(psf):
        if line.strip() == b'':
            continue
        if b'!NIMPHI' in line:
            items = line.split()
            n_impropers = int(items[0])
            break
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    d_array = fromstring(lines, count=n_dihedrals*4, dtype=int, sep=' ')
    if len(d_array) != n_dihedrals*4:
        raise IOError('number of dihedrals expected and parsed do not match')

    n_donors = 0
    lines = []
    for i, line in enumerate(psf):
        if line.strip() == b'':
            continue
        if b'!NDON' in line:
            items = line.split()
            n_donors = int(items[0])
            break
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    i_array = fromstring(lines, count=n_impropers*4, dtype=int, sep=' ')
    if len(i_array) != n_impropers*4:
        raise IOError('number of impropers expected and parsed do not match')

    n_acceptors = 0
    lines = []
    for i, line in enumerate(psf):
        if line.strip() == b'':
            continue
        if b'!NACC' in line:
            items = line.split()
            n_acceptors = int(items[0])
            break
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    do_array = fromstring(lines, count=n_donors*2, dtype=int, sep=' ')
    if len(do_array) != n_donors*2:
        raise IOError('number of donors expected and parsed do not match')

    n_exclusions = 0
    lines = []
    for i, line in enumerate(psf):
        if line.strip() == b'':
            continue
        if b'!NNB' in line:
            items = line.split()
            n_exclusions = int(items[0])
            break
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    ac_array = fromstring(lines, count=n_acceptors*2, dtype=int, sep=' ')
    if len(ac_array) != n_acceptors*2:
        raise IOError('number of acceptors expected and parsed do not match')

    lines = []
    for i, line in enumerate(psf):
        if line.strip() == b'':
            continue
        if b'!' in line:
            break
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    nbe_array = fromstring(lines, count=n_exclusions*2, dtype=int, sep=' ')
    if len(nbe_array) != n_exclusions*2:
        raise IOError('number of nonbonded exclusions expected and parsed do not match')

    n_crossterms = 0
    for i, line in enumerate(psf):
        if b'!NCRTERM' in line:
            items = line.split()
            n_crossterms = int(items[0])
            break

    lines = []
    for i, line in enumerate(psf):
        lines.append(line.decode(encoding='UTF-8'))
    
    lines = ''.join(lines)
    c_array = fromstring(lines, count=n_crossterms*4, dtype=int, sep=' ')
    if len(c_array) != n_crossterms*4:
        raise IOError('number of crossterms expected and parsed do not match')

    psf.close()
    ag.setSerials(serials)
    ag.setSegnames(segnames)
    ag.setResnums(resnums)
    ag.setResnames(resnames)
    ag.setNames(atomnames)
    ag.setTypes(atomtypes)
    ag.setCharges(charges)
    ag.setMasses(masses)

    if n_bonds > 0:
        b_array = add(b_array, -1, b_array)
        ag.setBonds(b_array.reshape((n_bonds, 2)))

    if n_angles > 0:
        a_array = add(a_array, -1, a_array)
        ag.setAngles(a_array.reshape((n_angles, 3)))

    if n_dihedrals > 0:
        d_array = add(d_array, -1, d_array)
        ag.setDihedrals(d_array.reshape((n_dihedrals, 4)))

    if n_impropers > 0:
        i_array = add(i_array, -1, i_array)
        ag.setImpropers(i_array.reshape((n_impropers, 4)))

    if n_donors > 0:
        do_array = add(do_array, -1, do_array)
        ag.setDonors(do_array.reshape((n_donors, 2)))

    if n_acceptors > 0:
        ac_array = add(ac_array, -1, ac_array)
        ag.setAcceptors(ac_array.reshape((n_acceptors, 2)))

    if n_exclusions > 0:
        nbe_array = add(nbe_array, -1, nbe_array)
        ag.setNBExclusions(nbe_array.reshape((n_exclusions, 2)))

    if n_crossterms > 0:
        c_array = add(c_array, -1, c_array)
        ag.setCrossterms(c_array.reshape((n_crossterms, 4)))

    return ag


PSFLINE = ('%8d %-4s %-4d %-4s %-4s %-4s %10.6f %13.4f %11d\n')

def writePSF(filename, atoms):
    """Write atoms in X-PLOR format PSF file with name *filename* and return
    *filename*.  This function will write available atom and bond information
    only."""

    if not filename.endswith('.psf'):
        filename = filename + '.psf'

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
        types = zeros(n_atoms, ATOMIC_FIELDS['type'].dtype)

    if charges is None:
        charges = zeros(n_atoms, ATOMIC_FIELDS['charge'].dtype)

    if masses is None:
        masses = zeros(n_atoms, ATOMIC_FIELDS['mass'].dtype)

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
    if len(bonds) > 0:
        bonds = array(bonds, int) + 1
        write('\n')
        write('{0:8d} !NBOND: bonds\n'.format(len(bonds)))
        for i, bond in enumerate(bonds):
            write('%8s%8s' % (bond[0], bond[1]))
            if i % 4 == 3:
                write('\n')
        if i % 4 != 3:
            write('\n')

    angles = list(atoms._iterAngles())
    if len(angles) > 0:
        angles = array(angles, int) + 1
        write('\n')
        write('{0:8d} !NTHETA: angles\n'.format(len(angles)))
        for i, angle in enumerate(angles):
            write('%8s%8s%8s' % (angle[0], angle[1], angle[2]))
            if i % 3 == 2:
                write('\n')
        if i % 3 != 2:
            write('\n')

    dihedrals = list(atoms._iterDihedrals())
    if len(dihedrals) > 0:
        dihedrals = array(dihedrals, int) + 1
        write('\n')
        write('{0:8d} !NPHI: dihedrals\n'.format(len(dihedrals)))
        for i, dihedral in enumerate(dihedrals):
            write('%8s%8s%8s%8s' % (dihedral[0], dihedral[1], dihedral[2], dihedral[3]))
            if i % 4 == 3:
                write('\n')
        if i % 4 != 3:
            write('\n')

    impropers = list(atoms._iterImpropers())
    if len(impropers) > 0:
        impropers = array(impropers, int) + 1
        write('\n')
        write('{0:8d} !NIMPHI: impropers\n'.format(len(impropers)))
        for i, improper in enumerate(impropers):
            write('%8s%8s%8s%8s' % (improper[0], improper[1], improper[2], improper[3]))
            if i % 2 == 1:
                write('\n')
        if i % 2 != 1:
            write('\n')

    write('\n')

    donors = list(atoms._iterDonors())
    if len(donors) > 0:
        donors = array(donors, int) + 1
        write('\n')
        write('{0:8d} !NDON: donors\n'.format(len(donors)))
        for i, donor in enumerate(donors):
            write('%8s%8s' % (donor[0], donor[1]))
            if i % 4 == 3:
                write('\n')
        if i % 4 != 3:
            write('\n')
    else:
        write('{0:8d} !NDON: donors\n'.format(0))
        write('\n')
    
    write('\n')

    acceptors = list(atoms._iterAcceptors())
    if len(acceptors) > 0:
        acceptors = array(acceptors, int) + 1
        write('\n')
        write('{0:8d} !NACC: acceptors\n'.format(len(acceptors)))
        for i, acceptor in enumerate(acceptors):
            write('%8s%8s' % (acceptor[0], acceptor[1]))
            if i % 4 == 3:
                write('\n')
        if i % 4 != 3:
            write('\n')
    else:
        write('{0:8d} !NACC: acceptors\n'.format(0))
        write('\n')

    nbexclusions = list(atoms._iterNBExclusions())
    if len(nbexclusions) > 0:
        nbexclusions = array(nbexclusions, int) + 1
        write('\n')
        write('{0:8d} !NNB\n'.format(len(nbexclusions)))
        for i, nbexclusion in enumerate(nbexclusions):
            write('%8s%8s' % (nbexclusion[0], nbexclusion[1]))
            if i % 4 == 3:
                write('\n')
        if i % 4 != 3:
            write('\n')
    else:
        write('{0:8d} !NNB\n'.format(0))
        write('\n')

    crossterms = list(atoms._iterCrossterms())
    if len(crossterms) > 0:
        crossterms = array(crossterms, int) + 1
        write('\n')
        write('{0:8d} !NCRTERM: crossterms\n'.format(len(crossterms)))
        for i, crossterm in enumerate(crossterms):
            write('%8s%8s%8s%8s' % (crossterm[0], crossterm[1], crossterm[2], crossterm[3]))
            if i % 2 == 1:
                write('\n')
        if i % 2 != 1:
            write('\n')

    write('\n')
    out.close()
    return filename
