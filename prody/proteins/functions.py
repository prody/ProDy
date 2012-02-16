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

"""This module defines miscellaneous functions dealing with protein data."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

from prody.atomic import Atomic, Atom, AtomGroup, Selection, HierView
from prody.atomic import ATOMIC_DATA_FIELDS
from prody.tools import openFile

__all__ = ['parsePSF', 'showProtein', 'writePQR', ]

def writePQR(filename, atoms):
    """Write *atoms* in PQR format to a file with name *filename*.  Only 
    current coordinate set is written.  Returns *filename* upon success.  If 
    *filename* ends with :file:`.gz`, a compressed file will be written.."""
    
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms does not have a valid type')
    if isinstance(atoms, Atom):
        atoms = Selection(atoms.getAtomGroup(), [atoms.getIndex()], 
                                atoms.getACSIndex(), 
                                'index ' + str(atoms.getIndex()))
    stream = openFile(filename, 'w')
    n_atoms = atoms.numAtoms()
    atomnames = atoms.getNames()
    if atomnames is None:
        raise RuntimeError('atom names are not set')
    for i, an in enumerate(atomnames):
        lenan = len(an)
        if lenan < 4:
            atomnames[i] = ' ' + an
        elif lenan > 4:
            atomnames[i] = an[:4]
    resnames = atoms._getResnames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms
    resnums = atoms._getResnums()
    if resnums is None:
        resnums = np.ones(n_atoms, int)
    chainids = atoms._getChids()
    if chainids is None: 
        chainids = np.zeros(n_atoms, '|S1')
    charges = atoms._getCharges()
    if charges is None:
        charges = np.zeros(n_atoms, float)
    radii = atoms._getRadii()
    if radii is None:
        radii = np.zeros(n_atoms, float)
    icodes = atoms._getIcodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, '|S1')
    hetero = ['ATOM'] * n_atoms 
    heteroflags = atoms._getHeteros()
    if heteroflags is not None:
        hetero = np.array(hetero, '|S6')
        hetero[heteroflags] = 'HETATM'
    altlocs = atoms._getAltlocs()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, '|S1')
    
    format = ('{0:6s}{1:5d} {2:4s}{3:1s}' +
              '{4:4s}{5:1s}{6:4d}{7:1s}   ' + 
              '{8:8.3f}{9:8.3f}{10:8.3f}' +
              '{11:8.4f}{12:7.4f}\n').format
    coords = atoms._getCoords()
    write = stream.write
    for i, xyz in enumerate(coords):
        write(format(hetero[i], i+1, atomnames[i], altlocs[i], 
                     resnames[i], chainids[i], int(resnums[i]), 
                     icodes[i], xyz[0], xyz[1], xyz[2], charges[i], radii[i]))
    write('TER\nEND')
    stream.close()
    return filename

def showProtein(*atoms, **kwargs):
    """Show protein representation using :meth:`~mpl_toolkits.mplot3d.Axes3D`.
    This function is designed for generating a quick view of the contents of a 
    :class:`~.AtomGroup` or :class:`~.Selection`.
           
    Protein atoms matching ``"calpha"`` selection are displayed using solid 
    lines by picking a random and unique color per chain.  Line with can 
    be adjusted using *lw* argument, e.g. ``lw=12``. Default width is 4.  
    Chain colors can be overwritten using chain identifier as in ``A='green'``.  
      
    Water molecule oxygen atoms are represented by red colored circles.  Color 
    can be changed using *water* keyword argument, e.g. ``water='aqua'``.
    Water marker and size can be changed using *wmarker* and *wsize* keywords, 
    defaults values are ``wmarker='.', wsize=6``.
    
    Hetero atoms matching ``"hetero and noh"`` selection are represented by 
    circles and unique colors are picked at random on a per residue basis.  
    Colors can be customized using residue name as in ``NAH='purple'``.  Note 
    that this will color all distinct residues with the same name in the same 
    color.  Hetero atom marker and size can be changed using *hmarker* and 
    *hsize* keywords, default values are ``hmarker='o', hsize=6``. 

    ProDy will set the size of axis so the representation is not distorted when
    the figure window is close to a square.  Colors are picked at random,
    except for water oxygens which will always be colored red.
    

    .. plot::
       :nofigs: 
       :context: 
        
       from prody import *
       import matplotlib.pyplot as plt
       import numpy as np
        
       plt.close('all')    

    .. plot::
       :context:
       :include-source:
       
       p38 = parsePDB('1p38')
       p38inh = parsePDB('1zz2')
       matchAlign(p38inh, p38)
       plt.figure(figsize=(5,4))
       showProtein(p38, p38inh)
       plt.legend(prop={'size': 10})
        
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    """
    
    alist = atoms
    for atoms in alist:    
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    cf = plt.gcf()
    show = None
    for child in cf.get_children():
        if isinstance(child, Axes3D):
            show = child
            break 
    if show is None:
        show = Axes3D(cf)
    from matplotlib import colors
    cnames = dict(colors.cnames)
    wcolor = kwargs.get('water', 'red').lower()
    avoid = np.array(colors.hex2color(cnames.pop(wcolor, cnames.pop('red'))))
    for cn, val in cnames.items():
        clr = np.array(colors.hex2color(val))
        if clr.sum() > 2.4:
            cnames.pop(cn)
        elif np.abs(avoid - clr).sum() <= 0.6:
            cnames.pop(cn)
    cnames = cnames.keys()
    import random
    random.shuffle(cnames)
    min_ = list()
    max_ = list()
    for atoms in alist:
        if isinstance(atoms, AtomGroup):
            title = atoms.getTitle()
        else:
            title = atoms.getAtomGroup().getTitle()
        calpha = atoms.select('calpha')
        if calpha:
            for ch in HierView(calpha, chain=True):
                xyz = ch._getCoords()
                chid = ch.getChid()
                show.plot(xyz[:,0], xyz[:,1], xyz[:,2], 
                          label=title + '_' + chid,
                          color=kwargs.get(chid, cnames.pop()).lower(),
                          lw=kwargs.get('lw', 4))
        water = atoms.select('water and noh')
        if water: 
            xyz = atoms.select('water')._getCoords()
            show.plot(xyz[:,0], xyz[:,1], xyz[:,2], label=title + '_water',
                      color=wcolor, 
                      ls='None', marker=kwargs.get('wmarker', '.'), 
                      ms=kwargs.get('wsize', 6))
        hetero = atoms.select('not protein and not nucleic and not water')
        if hetero:
            for res in HierView(hetero).iterResidues():
                xyz = res._getCoords()
                resname = res.getResname()
                resnum = str(res.getResnum())
                chid = res.getChid()
                show.plot(xyz[:,0], xyz[:,1], xyz[:,2], ls='None',
                          color=kwargs.get(resname, cnames.pop()).lower(), 
                          label=title + '_' + chid + '_' + resname + resnum, 
                          marker=kwargs.get('hmarker', 'o'), 
                          ms=kwargs.get('hsize', 6))
        xyz = atoms._getCoords()
        min_.append(xyz.min(0))
        max_.append(xyz.max(0))

    show.set_xlabel('x')
    show.set_ylabel('y')
    show.set_zlabel('z')
    min_ = np.array(min_).min(0)
    max_ = np.array(max_).max(0)
    center = (max_ + min_) / 2
    half = (max_ - min_).max() / 2
    show.set_xlim3d(center[0]-half, center[0]+half)
    show.set_ylim3d(center[1]-half, center[1]+half)
    show.set_zlim3d(center[2]-half, center[2]+half)
    if kwargs.get('legend', False):
        show.legend(prop={'size': 10})
    return show

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
