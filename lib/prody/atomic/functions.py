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

"""This module defines some functions for handling atomic classes and data."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from textwrap import wrap

from numpy import load, savez, zeros

from prody.utilities import openFile, rangeString
from prody import LOGGER

from . import flags
from . import select

from .atomic import Atomic
from .atomgroup import AtomGroup
from .atommap import AtomMap
from .bond import trimBonds, evalBonds
from .fields import ATOMIC_FIELDS
from .selection import Selection

__all__ = ['iterFragments', 'findFragments', 'loadAtoms', 'saveAtoms',
           'isReserved', 'listReservedWords', 'getReservedWords', 'sortAtoms']


SAVE_SKIP_ATOMGROUP = set(['numbonds', 'fragindex'])
SAVE_SKIP_POINTER = set(['numbonds', 'fragindex', 'segindex', 'chindex', 
                         'resindex'])

def saveAtoms(atoms, filename=None, **kwargs):
    """Save *atoms* in ProDy internal format.  All :class:`.Atomic` classes are
    accepted as *atoms* argument.  This function saves user set atomic data as 
    well.  Note that title of the :class:`.AtomGroup` instance is used as the 
    filename when *atoms* is not an :class:`.AtomGroup`.  To avoid overwriting
    an existing file with the same name, specify a *filename*."""
    
    try:
        atoms.getACSIndex()
    except AttributeError:
        raise TypeError('atoms must be Atomic instance, not {0:s}'
                        .format(type(atoms)))

    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        ag = atoms
        title = ag.getTitle()
        SKIP = SAVE_SKIP_ATOMGROUP
    else:
        SKIP = SAVE_SKIP_POINTER
        title = str(atoms)
    
    if filename is None:
        filename = ag.getTitle().replace(' ', '_')
    if '.ag.npz' not in filename:
        filename += '.ag.npz'
        
    attr_dict = {'title': title}
    attr_dict['n_atoms'] = atoms.numAtoms()
    attr_dict['n_csets'] = atoms.numCoordsets()
    attr_dict['cslabels'] = atoms.getCSLabels()
    attr_dict['flagsts'] = ag._flagsts
    coords = atoms._getCoordsets()
    if coords is not None:
        attr_dict['coordinates'] = coords
    bonds = ag._bonds
    bmap = ag._bmap
    if bonds is not None and bmap is not None:
        if atoms == ag:
            attr_dict['bonds'] = bonds
            attr_dict['bmap'] = bmap
            attr_dict['numbonds'] = ag._data['numbonds']
            frags = ag._data['fragindex']
            if frags is not None:
                attr_dict['fragindex'] = frags
        else:
            bonds = trimBonds(bonds, atoms._getIndices())
            attr_dict['bonds'] = bonds
            attr_dict['bmap'], attr_dict['numbonds'] = \
                evalBonds(bonds, len(atoms))
    
    for label in atoms.getDataLabels():
        if label in SKIP:
            continue
        attr_dict[label] = atoms._getData(label)
    for label in atoms.getFlagLabels():
        if label in SKIP:
            continue
        attr_dict[label] = atoms._getFlags(label)
    
    ostream = openFile(filename, 'wb', **kwargs)
    savez(ostream, **attr_dict)
    ostream.close()
    return filename


SKIPLOAD = set(['title', 'n_atoms', 'n_csets', 'bonds', 'bmap',
                'coordinates', 'cslabels', 'numbonds', 'flagsts',
                'segindex', 'chindex', 'resindex'])


def loadAtoms(filename):
    """Return :class:`.AtomGroup` instance loaded from *filename* using 
    :func:`numpy.load` function.  See also :func:`saveAtoms`."""
    
    LOGGER.timeit('_prody_loadatoms')
    attr_dict = load(filename)
    files = set(attr_dict.files)

    if not 'n_atoms' in files:
        raise ValueError('{0:s} is not a valid atomic data file'
                         .format(repr(filename)))
    title = str(attr_dict['title'])
    
    if 'coordinates' in files:
        coords = attr_dict['coordinates']
        ag = AtomGroup(title)
        ag._n_csets = int(attr_dict['n_csets'])
        ag._coords = coords
    ag._n_atoms = int(attr_dict['n_atoms'])
    ag._setTimeStamp()
    if 'flagsts' in files:
        ag._flagsts = int(attr_dict['flagsts'])
    
    if 'bonds' in files and 'bmap' in files and 'numbonds' in files:
        ag._bonds = attr_dict['bonds']
        ag._bmap = attr_dict['bmap']
        ag._data['numbonds'] = attr_dict['numbonds']
    
    skip_flags = set()
    
    for label, data in attr_dict.items():
        if label in SKIPLOAD:
            continue
        if data.ndim == 1 and data.dtype == bool:
            if label in skip_flags:
                continue
            else:
                ag._setFlags(label, data)
                skip_flags.update(flags.ALIASES.get(label, [label]))
        else:
            ag.setData(label, data)
    
    for label in ['segindex', 'chindex', 'resindex']:
        if label in attr_dict:
            ag._data[label] = attr_dict[label]
    
    if ag.numCoordsets() > 0:
        ag._acsi = 0
    
    if 'cslabels' in files:
        ag.setCSLabels(list(attr_dict['cslabels']))
    
    LOGGER.report('Atom group was loaded in %.2fs.', '_prody_loadatoms')
    return ag


def iterFragments(atoms):
    """Yield fragments, connected subsets in *atoms*, as :class:`.Selection` 
    instances."""
    
    try:
        return atoms.iterFragments()
    except AttributeError:
        pass
    
    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        raise TypeError('atoms must be an Atomic instance')
        
    bonds = atoms._iterBonds()
    return _iterFragments(atoms, ag, bonds)
    
def _iterFragments(atoms, ag, bonds):
    
    fids = zeros((len(ag)), int)
    fdict = {}
    c = 0
    for a, b in bonds:
        af = fids[a]
        bf = fids[b]
        if af and bf:
            if af != bf:
                frag = fdict[af]
                temp = fdict[bf]
                fids[temp] = af
                frag.extend(temp)
                fdict.pop(bf)
        elif af:
            fdict[af].append(b)
            fids[b] = af
        elif bf:
            fdict[bf].append(a)
            fids[a] = bf
        else:
            c += 1
            fdict[c] = [a, b]
            fids[a] = fids[b] = c
    fragments = []
    append = fragments.append
    fidset = set()
    indices = atoms._getIndices()
    for i, fid in zip(indices, fids[indices]):
        if fid in fidset:
            continue
        elif fid:
            fidset.add(fid)
            indices = fdict[fid]
            indices.sort()
            append(indices)
        else:
            # these are non-bonded atoms, e.g. ions
            append([i])

    acsi = atoms.getACSIndex()
    for indices in fragments:
        yield Selection(ag, indices, 'index ' + rangeString(indices), acsi, 
                        unique=True)


def findFragments(atoms):
    """Return list of fragments, connected subsets in *atoms*.  See also 
    :func:`iterFragments`."""
    
    return list(iterFragments(atoms))


RESERVED = set(ATOMIC_FIELDS)
RESERVED.update(['and', 'or', 'not', 'within', 'of', 'exwithin', 'same', 'as',
                 'bonded', 'exbonded', 'to', 'all', 'none',
                 'index', 'sequence', 'x', 'y', 'z']) 
RESERVED.update(flags.PLANTERS)
RESERVED.update(select.FUNCTIONS) 
RESERVED.update(select.FIELDS_SYNONYMS) 
RESERVED.update(['n_atoms', 'n_csets', 'cslabels', 'title', 'coordinates',
                 'bonds', 'bmap'])


def isReserved(word):
    """Return **True** if *word* is reserved for internal data labeling or atom
    selections.  See :func:`listReservedWords` for a list of reserved words."""
    
    return word in RESERVED
        
        
def getReservedWords():
    """Deprecated for removal in v1.4, use :func:`listReservedWords` instead.
    """
    
    from prody import deprecate
    deprecate('getReservedWords', 'listReservedWords')
    
    return listReservedWords()
        
        
def listReservedWords():
    """Return list of words that are reserved for atom selections and internal 
    variables. These words are: """

    words = list(RESERVED)
    words.sort()
    return words

_ = listReservedWords.__doc__ + '*' + '*, *'.join(listReservedWords()) + '*.'
listReservedWords.__doc__ = '\n'.join(wrap(_, 79))


def sortAtoms(atoms, label, reverse=False):
    """Return an :class:`.AtomMap` pointing to *atoms* sorted in ascending 
    data *label* order, or optionally in *reverse* order."""
    
    try:
        data, acsi = atoms.getData(label), atoms.getACSIndex()
    except AttributeError:
        raise TypeError('atoms must be an Atomic instance')
    else:
        if data is None:
            raise ValueError('{0:s} data is not set for {1:s}'
                                .format(repr(label), atoms))
    sort = data.argsort()
    if reverse:
        sort = sort[::-1]
    try:
        indices = atoms.getIndices()
    except AttributeError:
        ag = atoms
    else:
        ag = atoms.getAtomGroup()
        sort = indices[sort]
    return AtomMap(ag, sort, acsi)
