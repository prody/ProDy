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

"""This module defines :class:`HierView` class that builds a hierarchical 
views of atom groups.

Below example shows how to use hierarchical views.  We start by importing 
everything from the ProDy package:

>>> from prody import *

Then we parses a structure to get an :class:`.AtomGroup` instance which has a 
plain view of atoms: 

>>> structure = parsePDB('3mkb')
>>> structure
<AtomGroup: 3mkb (4776 atoms)>

.. _hierview:

Hierarchical view
===============================================================================

A hierarchical view of the structure can be simply get by calling the
:meth:`.AtomGroup.getHierView` method:

>>> hv = structure.getHierView()
>>> hv
<HierView: AtomGroup 3mkb (4 chains, 946 residues)>

Indexing
-------------------------------------------------------------------------------

Indexing :class:`HierView` instances return :class:`.Chain`:

>>> hv['A']
<Chain: A from 3mkb (254 residues, 1198 atoms)>
>>> hv['B']
<Chain: B from 3mkb (216 residues, 1193 atoms)>
>>> hv['Z'] # This will return None, which means chain Z does not exist

The length of the *hv* variable gives the number of chains in the structure:

>>> len(hv)
4
>>> hv.numChains()
4

It is also possible to get a :class:`.Residue` by 
directly indexing the :class:`HierView` instance:

>>> hv['A', 100]
<Residue: MET 100 from Chain A from 3mkb (8 atoms)>

Insertion codes can also be passed:

>>> hv['A', 100, 'B']

But this does not return anything, since residue 100B does not exist.

Iterations
-------------------------------------------------------------------------------

One can iterate over :class:`HierView` instances to get chains:

>>> for chain in hv:
...     chain # doctest: +ELLIPSIS
<Chain: A from 3mkb (254 residues, 1198 atoms)>
...
<Chain: D from 3mkb (231 residues, 1196 atoms)>
    
It is also possible to get a :func:`list` of chains simply as follows:

>>> chains = list( hv )
>>> chains # doctest: +SKIP
[<Chain: A from 3mkb (1198 atoms)>, 
 <Chain: B from 3mkb (1193 atoms)>, 
 <Chain: C from 3mkb (1189 atoms)>, 
 <Chain: D from 3mkb (1196 atoms)>]

Residues
-------------------------------------------------------------------------------

In addition, one can also iterate over all residues:

>>> for residue in hv.iterResidues():
...     residue # doctest: +ELLIPSIS
<Residue: ALA 1 from Chain A from 3mkb (5 atoms)>
<Residue: PHE 2 from Chain A from 3mkb (11 atoms)>
...
<Residue: HOH 475 from Chain D from 3mkb (1 atoms)>
<Residue: HOH 492 from Chain D from 3mkb (1 atoms)>

Chains
===============================================================================


>>> chA = hv['A']
>>> chA
<Chain: A from 3mkb (254 residues, 1198 atoms)>

Length of the chain equals to the number of residues in it:

>>> len(chA)
254
>>> chA.numResidues()
254

Indexing
-------------------------------------------------------------------------------

Indexing a :class:`.Chain` instance returns a :class:`.Residue` instance.

>>> chA[1]
<Residue: ALA 1 from Chain A from 3mkb (5 atoms)>

If a residue does not exist, ``None`` is returned:

>>> chA[1000]
>>> chA[1, 'A'] # Residue 1 with insertion code A also does not exist  

If residue with given integer number does not exist, ``None`` is returned. 

Iterations
-------------------------------------------------------------------------------

Iterating over a chain yields residues:

>>> for residue in chA:
...     residue # doctest: +ELLIPSIS
<Residue: ALA 1 from Chain A from 3mkb (5 atoms)>
<Residue: PHE 2 from Chain A from 3mkb (11 atoms)>
...
<Residue: HOH 490 from Chain A from 3mkb (1 atoms)>
<Residue: HOH 493 from Chain A from 3mkb (1 atoms)>

Note that water atoms, each constituting a residue, are also part of a chain
if they are labeled with that chain's identifier.

This enables getting a :func:`list` of residues simply as follows:

>>> chA_residues = list(chA)
>>> chA_residues # doctest: +SKIP
[<Residue: ALA 1 from Chain A from 3mkb (5 atoms)>,
 ...,
 <Residue: HOH 493 from Chain A from 3mkb (1 atoms)>]

Get data
-------------------------------------------------------------------------------

All methods defined for :class:`.AtomGroup` class are also defined for 
:class:`.Chain` and :class:`.Residue` classes:

>>> print( chA.getCoords() ) # doctest: +ELLIPSIS
[[ -2.139  17.026 -13.287]
 [ -1.769  15.572 -13.111]
 [ -0.296  15.257 -13.467]
 ...
 [ -5.843  17.181 -16.86 ]
 [-13.199  -9.21  -49.692]
 [ -0.459   0.378 -46.156]]
>>> print( chA.getBetas() )
[ 59.35  59.14  58.5  ...,  57.79  47.77  40.77]

Selections
-------------------------------------------------------------------------------

Finally, you can select atoms from a :class:`.Chain` instance:

>>> chA_backbone = chA.select('backbone')
>>> chA_backbone
<Selection: '(backbone) and (chain A)' from 3mkb (560 atoms)>
>>> chA_backbone.getSelstr()
'(backbone) and (chain A)'

As you see, the selection string passed by the user is augmented with 
"chain" keyword and identifier automatically to provide internal
consistency:

>>> structure.select( chA_backbone.getSelstr() )
<Selection: '(backbone) and (chain A)' from 3mkb (560 atoms)>
 

Residues
===============================================================================

>>> chA_res1 = chA[1]
>>> chA_res1
<Residue: ALA 1 from Chain A from 3mkb (5 atoms)>

Indexing
-------------------------------------------------------------------------------

:class:`.Residue` instances can be indexed to get individual atoms:

>>> chA_res1['CA']
<Atom: CA from 3mkb (index 1)>
>>> chA_res1['CB']
<Atom: CB from 3mkb (index 4)>
>>> chA_res1['X'] # if atom does not exist, None is returned

Iterations
-------------------------------------------------------------------------------

Iterating over a residue instance yields :class:`Atom` instances:

>>> for atom in chA_res1:
...     atom # doctest: +ELLIPSIS
<Atom: N from 3mkb (index 0)>
...
<Atom: CB from 3mkb (index 4)>

This makes it easy to get a :func:`list` of atoms:

>>> list( chA_res1 ) # doctest: +SKIP
[<Atom: N from 3mkb (index 0)>,
 <Atom: CA from 3mkb (index 1)>,
 <Atom: C from 3mkb (index 2)>,
 <Atom: O from 3mkb (index 3)>,
 <Atom: CB from 3mkb (index 4)>]

Get data
-------------------------------------------------------------------------------

All methods defined for :class:`.AtomGroup` class are also defined for 
:class:`.Residue` class:

>>> print( chA_res1.getCoords() )
[[ -2.139  17.026 -13.287]
 [ -1.769  15.572 -13.111]
 [ -0.296  15.257 -13.467]
 [  0.199  14.155 -13.155]
 [ -2.752  14.639 -13.898]]
>>> print( chA_res1.getBetas() )
[ 59.35  59.14  58.5   59.13  59.02]

Selections
-------------------------------------------------------------------------------

Finally, you can select atoms from a :class:`.Residue` instance:

>>> chA_res1_bb = chA_res1.select('backbone')
>>> chA_res1_bb
<Selection: '(backbone) and ... and (chain A))' from 3mkb (4 atoms)>
>>> chA_res1_bb.getSelstr()
'(backbone) and (resnum 1 and (chain A))'

Again, the selection string is augmented with the chain identifier and 
residue number ("resnum").

Atoms
===============================================================================

The lowest level of the hierarchical view contains :class:`Atom` instances.

>>> chA_res1_CA = chA_res1['CA']
>>> chA_res1_CA
<Atom: CA from 3mkb (index 1)>

*Get atomic data*

All methods defined for :class:`.AtomGroup` class are also defined for 
:class:`.Atom` class with the difference that method names are singular 
(except for coordinates):

>>> print( chA_res1_CA.getCoords() )
[ -1.769  15.572 -13.111]
>>> print( chA_res1_CA.getBeta() )
59.14

State Changes 
===============================================================================

A :class:`HierView` instance represents the state of an :class:`.AtomGroup` 
instance at the time it is built.  When chain identifiers or residue numbers 
change, the state that hierarchical view represents may not match the current 
state of the atom group:

>>> chA.setChid('X')
>>> chA
<Chain: X from 3mkb (254 residues, 1198 atoms)>
>>> hv['X'] # returns None, since hierarchical view is not updated
>>> hv.update() # this updates hierarchical view
>>> hv['X']
<Chain: X from 3mkb (254 residues, 1198 atoms)>

When this is the case, :meth:`HierView.update` method can be used to update 
hierarchical view."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import unique, zeros, arange, concatenate

from atomgroup import AtomGroup
from selection import Selection
from chain import Chain
from residue import Residue
from segment import Segment

from prody import LOGGER, SETTINGS

__all__ = ['HierView']


class HierView(object):
    
    """Hierarchical views can be generated for :class:`.AtomGroup` and 
    :class:`.Selection` instances.  Indexing a :class:`HierView` instance 
    returns a :class:`.Chain`  instance.
    
    Some :class:`object` methods are customized as follows:
    
    * :func:`len` returns the number of atoms, i.e. :meth:`numChains`
    * :func:`iter` yields :class:`.Chain` instances
    * indexing by:
         - *segment name* (:func:`str`), e.g. ``"PROT"``, returns
           a :class:`.Segment` 
         - *chain identifier* (:func:`str`), e.g. ``"A"``, returns 
           a :class:`.Chain`
         - *[segment name,] chain identifier, residue number[, insertion code]* 
           (:func:`tuple`), e.g. ``"A", 10`` or  ``"A", 10, "B"`` or
           ``"PROT", "A", 10, "B"``, returns a :class:`.Residue`
        
    Note that when an :class:`.AtomGroup` instance have distinct segments, 
    they will be considered when building the hierarchical view.  
    A :class:`.Segment` instance will be generated for each distinct segment 
    name.  Then, for each segment chains and residues will be evaluated.  
    Having segments in the structure will not change most behaviors of this 
    class, except indexing.  For example, when indexing a hierarchical view 
    for chain P in segment PROT needs to be indexed as ``hv['PROT', 'P']``."""
    
    def __init__(self, atoms, **kwargs):
        
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
    
    def __len__(self):
        
        return len(self._chains)
    
    def __getitem__(self, key):
        
        if isinstance(key, str):
            return self.getSegment(key) or self.getChain(key)
        
        elif isinstance(key, tuple):
            length = len(key)
        
            if length == 1:
                return self.__getitem__(key[0])

            if length == 2:
                return (self.getChain(key[1], key[0]) or 
                         self.getResidue(*key) or
                         self.getResidue(None, key[0], key[1]))
        
            if length == 3:
                return self.getResidue(*key) or self.getResidue(key[1],
                                                key[2], None, key[0]) 
        
            if length == 4:
                return getResidue(key[1], key[2], key[3], key[0])
        
        elif isinstance(key, int):
            return self.getResidue(None, key)

    def _getSegname(self):
        """Return name of the segment when there is only one segment."""
        
        if self.numSegments() == 1:
            return self._ag._getSegnames()[0]

    def _getChid(self):
        """Return identifier of the chain when there is only one chain."""
        
        if self.numChains() == 1:
            return self._ag._getChids()[0]

    def _getResidue(self, index):
        
        try:
            residue = self._residues[index]
        except IndexError:
            pass
        else:
            if residue is not None:
                try:
                    residue.getAtomGroup()
                except AttributeError:
                    residue = self._residues[index] = Residue(self._ag, 
                                            residue, self, self._acsi, 
                                            unique=True, selstr=self._selstr)
            return residue

    def _getChain(self, index):
        
        try:
            chain = self._chains[index]
        except IndexError:
            pass
        else:
            if chain is not None:
                try:
                    chain.getAtomGroup()
                except AttributeError:
                    chain = self._chains[index] = Chain(self._ag, 
                                            chain, self, self._acsi, 
                                            unique=True, selstr=self._selstr)
            return chain

    def _getSegment(self, index):
        
        try:
            segment = self._segments[index]
        except IndexError:
            pass
        else:
            if segment is not None:
                try:
                    segment.getAtomGroup()
                except AttributeError:
                    segment = self._segments[index] = Segment(self._ag, 
                                            segment, self, self._acsi, 
                                            unique=True, selstr=self._selstr)
            return segment

    def getAtoms(self):
        """Return atoms for which the hierarchical view was built."""
        
        return self._atoms
    
    def update(self, **kwargs):
        """Update (or build) hierarchical view of atoms.  This method is called 
        at instantiation, but can be used to rebuild the hierarchical view when
        attributes of atoms change."""

        self._acsi = self._atoms.getACSIndex()
        try:
            self._ag = self._atoms.getAtomGroup()
        except AttributeError:
            self._selstr = None
            self._update(**kwargs)
        else:
            self._selhv(**kwargs)

    def _selhv(self, **kwargs):

        atoms = self._atoms
        ag = self._ag
        indices = atoms._getIndices()
        self._selstr = atoms.getSelstr()
        
        self._dict = ag.getHierView()._dict

        self._segments = _segments = [None] * ag.numSegments()
        self._residues = _residues = [None] * ag.numResidues() 
        self._chains = _chains = [None] * ag.numChains()

        for hvidx, _list in [(atoms._getSegindices(), _segments),
                             (atoms._getChindices(), _chains),
                             (atoms._getResindices(), _residues),]: 
            if not _list: continue
            pidx = hvidx[0]
            pi = 0
            for i, idx in enumerate(hvidx):
                if pidx == idx: continue
                subset = _list[pidx]
                if subset is None:
                    _list[pidx] = indices[pi:i]
                else:
                    _list[pidx] = concatenate((subset, indices[pi:i]))
                pidx, pi = idx, i
            subset = _list[pidx]
            if subset is None:
                _list[pidx] = indices[pi:]
            else:
                _list[pidx] = concatenate((subset, indices[pi:]))
            
    def _update(self, **kwargs):
            
        atoms = ag = self._ag = self._atoms
        n_atoms = len(ag)
        _indices = arange(n_atoms)

        self._dict = _dict = dict()
        self._residues = _residues = list()
        self._segments = _segments = list()
        self._chains = _chains = list()

        # identify segments
        segindex = -1
        segindices = zeros(n_atoms, int)

        sgnms = ag._getSegnames()
        if sgnms is None:
            _segments = None
        else:
            s = sgnms[0]
            if len(unique(sgnms)) == 1:
                # 1 segment
                if s:
                    _segments.append(_indices)
                    _dict[s] = 0
                else: 
                    _segments = None
            else:
                ps = None       # previous segment name 
                for i, s in enumerate(sgnms):
                    if s == ps or s in _dict:
                        continue
                    ps = s
                    idx = _indices[i:][sgnms[i:] == s]
                    segindex += 1
                    segindices[idx] = segindex
                    _dict[s] = segindex
                    _segments.append(idx)

        ag._data['segindex'] = segindices

        # identify chains
        chindex = -1
        chindices = zeros(n_atoms, int)

        chids = ag._getChids()
        if chids is None:
            _chains = None
        else:
            if _segments is None:
                if len(unique(chids)) == 1:
                    _dict[(None, chids[0] or None)] = 0
                    _chains.append(_indices)
                else:
                    pc = None
                    for i, c in enumerate(chids):
                        if c == pc or (None, c) in _dict:
                            continue
                        pc = c
                        idx = _indices[i:][chids[i:] == c]
                        chindex += 1
                        chindices[idx] = chindex
                        _dict[(None, c)] = chindex
                        _chains.append(idx)
            else:
                pc = chids[0]
                ps = sgnms[0]
                _i = 0
                for i, c in enumerate(chids):
                    s = sgnms[i]
                    if c == pc and s == ps:
                        continue
                    s_c = (ps, pc or None)
                    cid = _dict.get(s_c)
                    idx = _indices[_i:i]
                    if cid is None:
                        segment = _dict[ps]
                        chindex += 1
                        chindices[idx] = chindex
                        _dict[s_c] = chindex
                        _chains.append(idx)
                    else:
                        chain = _chains[cid]
                        chindices[idx] = chain._indices[0]
                        _chains[cid] = concatenate((chain, idx))
                    pc = c
                    ps = s
                    _i = i
                s_c = (ps, pc or None)
                cid = _dict.get(s_c)
                idx = _indices[_i:]
                if cid is None:
                    segment = _dict[ps]
                    chindex += 1
                    chindices[idx] = chindex
                    _dict[s_c] = chindex
                    _chains.append(idx)
                else:
                    chain = _chains[cid]
                    chindices[idx] = chain._indices[0]
                    _chains[cid] = concatenate((chain, idx))

        ag._data['chindex'] = chindices

        if kwargs.get('chain') == True:
            return
    
        # identify residues
        resindex = -1
        resindices = zeros(n_atoms, int)
    
        rnums = ag._getResnums()
        if rnums is None:
            raise ValueError('resnums are not set')
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

        pr = rnums[0]
        pi = icods[0] or None          
        pc = chids[0]
        ps = sgnms[0]
        _j = 0
        _append = _residues.append
        _get = _dict.get
        _set = _dict.__setitem__
        for j, r in enumerate(rnums):
            i = icods[j] or None
            c = chids[j]
            s = sgnms[j]
            if r != pr or i != pi or c != pc or s != ps:
                s_c_r_i = (ps, pc, pr, pi)
                rid = _get(s_c_r_i)
                idx = _indices[_j:j]
                if rid is None:
                    resindex += 1
                    resindices[idx] = resindex
                    _set(s_c_r_i, resindex)
                    _append(idx)
                else:
                    residue = _residues[rid]
                    resindices[idx] = residue[0]
                    _residues[rid] = concatenate((residue, idx))
                ps = s
                pc = c
                pr = r
                pi = i
                _j = j 
        s_c_r_i = (ps, pc, pr, pi)
        rid = _get(s_c_r_i)
        idx = _indices[_j:]
        if rid is None:
            resindex += 1
            resindices[idx] = resindex
            _append(idx)
            _set(s_c_r_i, resindex)
        else:
            residue = _residues[rid]
            resindices[idx] = residue[0]
            _residues[rid] = concatenate((residue, idx))
        
        ag._data['resindex'] = resindices

    def getResidue(self, chid, resnum, icode=None, segname=None):
        """Return residue with number *resnum* and insertion code *icode* from 
        the chain with identifier *chid* in segment with name *segname*."""
        
        try:
            index = self._dict[(segname or self._getSegname(), 
                                chid or self._getChid(), 
                                resnum, icode or None)]
        except KeyError:
            pass
        else:
            return self._getResidue(index)

    def numResidues(self):
        """Return number of residues."""
        
        return (len(self._residues) if self._ag is self._atoms else
                len(self._residues) - self._residues.count(None))

    def iterResidues(self):
        """Yield residues."""
        
        alist = self._residues
        ag = self._ag
        acsi = self._acsi
        selstr = self._selstr
        for i, item in enumerate(alist):
            if item is None:
                continue
            try:
                item.dtype
            except AttributeError:
                pass
            else:
                item = alist[i] = Residue(ag, item, self, acsi, selstr=selstr,
                                          unique=True)
            yield item
                
    def getChain(self, chid, segname=None):
        """Return chain with identifier *chid*, if it is present."""
        
        try:
            index = self._dict[(segname or self._getSegname(), 
                                chid or None)]
        except KeyError:
            pass
        else:
            return self._getChain(index)

    def iterChains(self):
        """Yield chains."""

        alist = self._chains
        ag = self._ag
        acsi = self._acsi
        selstr = self._selstr
        for i, item in enumerate(alist):
            if item is None:
                continue
            try:
                item.dtype
            except AttributeError:
                pass
            else:
                item = alist[i] = Chain(ag, item, self, acsi, selstr=selstr,
                                        unique=True)
            yield item
    
    __iter__ = iterChains
    
    def numChains(self):
        """Return number of chains."""
        
        return (len(self._chains) if self._ag is self._atoms else
                len(self._chains) - self._chains.count(None))

    def getSegment(self, segname):
        """Return segment with name *segname*, if it is present."""
        
        try:
            index = self._dict[segname or None]
        except KeyError:
            pass
        else:
            return self._getSegment(index)

    def numSegments(self):
        """Return number of chains."""
        
        return (len(self._segments) if self._ag is self._atoms else
                len(self._segments) - self._segments.count(None))

    def iterSegments(self):
        """Yield segments."""

        alist = self._segments
        ag = self._ag
        acsi = self._acsi
        selstr = self._selstr
        for i, item in enumerate(alist):
            if item is None:
                continue
            try:
                item.dtype
            except AttributeError:
                pass
            else:
                item = alist[i] = Segment(ag, item, self, acsi, selstr=selstr,
                                          unique=True)
            yield item
