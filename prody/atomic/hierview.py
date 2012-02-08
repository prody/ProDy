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
"""This module contains a class for building hierarchical views of atom groups.

.. currentmodule:: prody.atomic"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from atomic import MultiCoordset
from atomgroup import AtomGroup
from selection import Selection
from chain import Chain, MCChain
from residue import Residue, MCResidue
from segment import Segment, MCSegment

__all__ = ['HierView']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER 
SETTINGS = pkg.SETTINGS


class HierView(object):
    
    """Hierarchical views can be generated for :class:`AtomGroup` and 
    :class:`Selection` instances.  Indexing a :class:`HierView` instance 
    returns a :class:`Chain` instance.
    
    >>> from prody import *
    >>> pdb = parsePDB('1p38')
    >>> hv = pdb.getHierView()
    >>> chA = hv['A']
    >>> chA
    <Chain: A from 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>
    >>> print hv['B'] # Chain B does not exist in 1p38
    None
    
    Note that is the atom group instance have distinct segments, they will
    be considered when building the hierarchical view.  A :class:`Segment`
    instance will be generated for each distinct segment name.  Then,
    for each segment chains and residues will be evaluated.  Having
    segments in the structure will not change most behaviors of this class,
    except indexing.  For example, when indexing a hierarchical view for 
    chain P in segment PROT needs to be indexed as ``hv['PROT', 'P']``."""
    
    def __init__(self, atoms, **kwargs):
        """Build hierarchical view for *atoms*."""
        
        
        if not isinstance(atoms, (AtomGroup, Selection)):
            raise TypeError('atoms must be an AtomGroup or Selection instance')
        self._atoms = atoms
        self.update(**kwargs)
        
    def __repr__(self):
        
        if self._segments:
            return ('<HierView: {0:s} ({1:d} segments, {2:d} chains, {3:d} '
                    'residues)>').format(str(self._atoms), self.numSegments(),
                                         self.numChains(), self.numResidues())
        else:
            return ('<HierView: {0:s} ({1:d} chains, {2:d} residues)>'
               ).format(str(self._atoms), self.numChains(), self.numResidues())
    
    def __str__(self):
        return 'HierView of {0:s}'.format(str(self._atoms))
    
    def __iter__(self):
        """Iterate over chains."""
        
        return self._chains.__iter__()
    
    def __len__(self):
        return len(self._chains)
    
    def __getitem__(self, key):
        
        _dict = self._dict
        if isinstance(key, str):
            return _dict.get(key, self._dict.get((None, key)))
        
        elif isinstance(key, tuple):
            length = len(key)
        
            if length == 1:
                return self.__getitem__(key[0])

            if length == 2:
                return _dict.get(key,
                        _dict.get((None, None, key[0], key[1] or None), 
                        _dict.get((None, key[0], key[1], None))))
        
            if length == 3:
                return _dict.get((None, key[0] or None, 
                                       key[1], key[2] or None), 
                          _dict.get((key[0] or None, key[1] or None, 
                                     key[2], None)))
        
            if length == 4:
                return _dict.get((key[0] or None, key[1] or None, 
                                  key[2], key[3] or None))
        
        elif isinstance(key, int):
            return _dict.get((None, None, key, None))

    def getAtoms(self):
        """Return atoms for which the hierarchical view was built.
        
        .. versionadded:: 0.6.2"""
        
        return self._atoms
    
    def update(self, **kwargs):
        """Rebuild hierarchical view of atoms.  This method is called at 
        instantiation, but can be used to rebuild the hierarchical view 
        when attributes of atoms change."""
        
        array = np.array
        atoms = self._atoms
        if isinstance(atoms, AtomGroup):
            ag = atoms
            _indices = np.arange(ag._n_atoms)
            selstr = False
        else:
            ag = atoms.getAtomGroup()
            _indices = atoms._indices
            selstr = atoms.getSelstr()
        
        if isinstance(ag, MultiCoordset):
            acsi = self._atoms.getACSIndex()
            Segment_ = MCSegment
            Chain_ = MCChain
            Residue_ = MCResidue
        else:
            acsi = None
            Segment_ = Segment
            Chain_ = Chain
            Residue_ = Residue
        
        n_atoms = len(ag)
        self._dict = _dict = dict()
        self._chains = _chains = list()
        self._residues = _residues = list()
        self._segments = _segments = list()

        segindex = -1
        segindices = np.zeros(n_atoms, int)
        chindex = -1
        chindices = np.zeros(n_atoms, int)
        resindex = -1
        resindices = np.zeros(n_atoms, int)
        
        sgnms = ag._getSegnames()
        if sgnms is None:
            _segments = None
        else:
            if selstr:
                sgnms = sgnms[_indices]
            unique = np.unique(sgnms)
            s = sgnms[0]
            if len(unique) == 1:
                if  s != '':
                    segment = Segment_(ag, _indices, acsi=acsi, unique=True, 
                                      selstr=selstr)
                    _dict[s] = segment
                    _segments.append(segment)
                    LOGGER.info('Hierarchical view contains segments.')
                else: 
                    _segments = None
            else:
                ps = None
                for i, s in enumerate(sgnms):
                    if s == ps or s in _dict:
                        continue
                    ps = s
                    segindex += 1
                    idx = _indices[i:][sgnms[i:] == s]
                    segment = Segment_(ag, idx, acsi=acsi, unique=True, 
                                       selstr=selstr)
                    segindices[idx] = segindex
                    _dict[s] = segment
                    _segments.append(segment)
                LOGGER.info('Hierarchical view contains segments.')

        chids = ag._getChids()
        if chids is None:
            _chains = None
        else:
            if selstr:
                chids = chids[_indices]
            if _segments is None:
                if len(np.unique(chids)) == 1:
                    chain = Chain_(ag, _indices, acsi=acsi, unique=True)
                    _dict[(None, chids[0] or None)] = chain
                    _chains.append(chain)
                else:
                    pc = None
                    for i, c in enumerate(chids):
                        if c == pc or (None, c) in _dict:
                            continue
                        pc = c
                        chindex += 1
                        idx = _indices[i:][chids[i:] == c]
                        chain = Chain_(ag, idx, acsi=acsi, unique=True)
                        chindices[idx] = chindex
                        _dict[(None, c)] = chain
                        _chains.append(chain)
            else:
                pc = chids[0]
                ps = sgnms[0]
                _i = 0
                for i, c in enumerate(chids):
                    s = sgnms[i]
                    if c == pc and s == ps:
                        continue
                    s_c = (ps, pc or None)
                    chain = _dict.get(s_c)
                    if chain is None:
                        segment = _dict[ps]
                        chindex += 1
                        idx = _indices[_i:i]
                        chain = Chain_(ag, idx, acsi=acsi, segment=segment, 
                                       unique=True)
                        chindices[idx] = chindex
                        _dict[s_c] = chain
                        segment._dict[pc] = chain
                        _chains.append(chain)
                    else:
                        idx = _indices[_i:i]
                        chindices[idx] = chain._indices[0]
                        chain._indices = np.concatenate((chain._indices, idx))
                    pc = c
                    ps = s
                    _i = i
                s_c = (ps, pc or None)
                chain = _dict.get(s_c)
                idx = _indices[_i:]
                if chain is None:
                    segment = _dict[ps]
                    chindex += 1
                    chindices[idx] = chindex
                    chain = Chain_(ag, idx, acsi=acsi, segment=segment, 
                                   unique=True)
                    _dict[s_c] = chain
                    segment._dict[pc] = chain
                    _chains.append(chain)
                else:
                    chindices[idx] = chain._indices[0]
                    chain._indices = np.concatenate((chain._indices, idx)) 
        
        if kwargs.get('chain') == True:
            return
    
        rnums = ag._getResnums()
        if rnums is None:
            raise ValueError('resnums are not set')
        if selstr:
            rnums = rnums[_indices]
        nones = None
        if _segments is None:
            if nones is None:
                nones = [None] * len(rnums)
            sgnms = nones
        if _chains is None:
            if nones is None:
                nones = [None] * len(rnums)
            chids = nones
        icods = ag._getIcodes()
        if icods is None:
            if nones is None:
                nones = [None] * len(rnums)
            icods = nones
        elif selstr:
                icods = icods[_indices]

        pr = rnums[0]
        pi = icods[0] or None          
        pc = chids[0]
        ps = sgnms[0] or None
        _j = 0
        for j, r in enumerate(rnums):
            i = icods[j] or None
            c = chids[j] or None
            s = sgnms[j]
            if r != pr or i != pi or c != pc or s != ps:
                s_c_r_i = (ps, pc, pr, pi)
                res = _dict.get(s_c_r_i)
                if res is None:
                    chain = _dict.get((ps, pc))
                    resindex += 1
                    idx = _indices[_j:j]
                    res = Residue_(ag, idx, acsi=acsi, chain=chain, 
                                   unique=True, selstr=selstr)
                    resindices[idx] = resindex
                    if chain is not None:
                        chain._dict[(pr, pi)] = len(chain._list)
                        chain._list.append(res)
                    _residues.append(res)
                    _dict[s_c_r_i] = res
                else:
                    res._indices = np.concatenate((res._indices, 
                                                   _indices[_j:j]))
                ps = s
                pc = c
                pr = r
                pi = i
                _j = j 
        s_c_r_i = (ps, pc, pr, pi)
        res = _dict.get(s_c_r_i)
        idx = _indices[_j:]
        if res is None:
            chain = _dict.get((ps, pc))
            resindex += 1
            res = Residue_(ag, idx, acsi=acsi, chain=chain, unique=True, 
                           selstr=selstr)
            resindices[idx] = resindex
            if chain is not None:
                chain._dict[(pr, pi)] = len(chain._list)
                chain._list.append(res)
            _residues.append(res)
            _dict[s_c_r_i] = res
        else:
            resindices[idx] = res._indices[0]
            res._indices = np.concatenate((res._indices, idx))
        
        ag._data['segindices'] = segindices
        ag._data['chindices'] = chindices
        ag._data['resindices'] = resindices

    def getResidue(self, chid, resnum, icode=None, segname=None):
        """Return residue with number *resnum* and insertion code *icode* from 
        the chain with identifier *chid* in segment with name *segname*."""
        
        return self._dict.get((segname or None, chid or None, 
                               resnum, icode or None))

    def numResidues(self):
        """Return number of residues."""
        
        return len(self._residues)    

    def iterResidues(self):
        """Iterate over residues."""
        
        return self._residues.__iter__()
                
    def getChain(self, chid, segname=None):
        """Return chain with identifier *chid*, if it exists."""
        
        return self._dict.get((segname or None, chid or None))

    def iterChains(self):
        """Iterate over chains."""

        return self._chains.__iter__()
    
    def numChains(self):
        """Return number of chains."""
        
        return len(self._chains)

    def getSegment(self, segname):
        """Return segment with name *segname*, if it exists."""
        
        return self._dict.get(segname or None)

    def numSegments(self):
        """Return number of chains."""
        
        return len(self._segments)

    def iterSegments(self):
        """Iterate over segments."""

        return self._segments.__iter__()
