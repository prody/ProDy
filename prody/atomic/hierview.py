# -*- coding: utf-8 -*-
"""This module defines :class:`HierView` class that builds a hierarchical
views of atom groups."""

from numpy import unique, zeros, arange, concatenate
from prody.utilities.misctools import count

from .atomgroup import AtomGroup
from .selection import Selection
from .chain import Chain
from .residue import Residue
from .segment import Segment

__all__ = ['HierView']


class HierView(object):

    """Hierarchical views can be generated for :class:`.AtomGroup`,
    :class:`.Selection`, and :class:`.Chain` instances. Indexing a 
    :class:`HierView` instance returns a :class:`.Chain` instance.

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

        if not isinstance(atoms, (AtomGroup, Selection, Chain, Segment)):
            raise TypeError('atoms must be an AtomGroup, Selection, Chain, or Segment instance')

        self._atoms = atoms
        self.update(**kwargs)

    def __repr__(self):

        if self._segments:
            return ('<HierView: {0} ({1} segments, {2} chains, {3} '
                    'residues)>').format(str(self._atoms), self.numSegments(),
                                         self.numChains(), self.numResidues())
        else:
            return ('<HierView: {0} ({1} chains, {2} residues)>'
               ).format(str(self._atoms), self.numChains(), self.numResidues())

    def __str__(self):

        return 'HierView of {0}'.format(str(self._atoms))

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
                return self.getResidue(key[1], key[2], key[3], key[0])

        elif isinstance(key, int):
            return self.getResidue(None, key)

    def _getSegname(self):
        """Returns name of the segment when there is only one segment."""

        if self.numSegments() == 1:
            return self._ag._getSegnames()[0]

    def _getChid(self):
        """Returns identifier of the chain when there is only one chain."""

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
        """Returns atoms for which the hierarchical view was built."""

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
        """Build hierarchical view for :class:`.Selection` instances."""

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
        """Build hierarchical view for :class:`.AtomGroup` instances."""

        ag = self._ag = self._atoms
        n_atoms = len(ag)
        _indices = arange(n_atoms)

        self._dict = _dict = {}
        self._residues = _residues = []
        self._segments = _segments = []
        self._chains = _chains = []

        nones = None
        getnones = lambda: [None] * n_atoms if nones is None else nones
        termini = ag.getFlags('pdbter')
        if termini is None:
            termini = getnones()

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
                        #segment = _dict[ps]
                        chindex += 1
                        chindices[idx] = chindex
                        _dict[s_c] = chindex
                        _chains.append(idx)
                    else:
                        chain = _chains[cid]
                        chindices[idx] = cid
                        _chains[cid] = concatenate((chain, idx))
                    pc = c
                    ps = s
                    _i = i
                s_c = (ps, pc or None)
                cid = _dict.get(s_c)
                idx = _indices[_i:]
                if cid is None:
                    #segment = _dict[ps]
                    chindex += 1
                    chindices[idx] = chindex
                    _dict[s_c] = chindex
                    _chains.append(idx)
                else:
                    chain = _chains[cid]
                    chindices[idx] = cid
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
        if _segments is None:
            sgnms = nones = getnones()
        if _chains is None:
            chids = nones = getnones()
        icods = ag._getIcodes()
        if icods is None:
            icods = nones = getnones()
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
            if r != pr or i != pi or c != pc or s != ps or (j and termini[j-1]):
                s_c_r_i = (ps, pc, pr, pi)
                rid = _get(s_c_r_i)
                idx = _indices[_j:j]
                if (rid is None or isinstance(rid, list) or
                    termini[_residues[rid][-1]]):
                    resindex += 1
                    resindices[idx] = resindex
                    _append(idx)
                    if rid is None:
                        _set(s_c_r_i, resindex)
                    elif isinstance(rid, list):
                        rid.append(resindex)
                    else:
                        _set(s_c_r_i, [rid, resindex])
                else:
                    residue = _residues[rid]
                    resindices[idx] = rid
                    _residues[rid] = concatenate((residue, idx))
                ps = s
                pc = c
                pr = r
                pi = i
                _j = j
        s_c_r_i = (ps, pc, pr, pi)
        rid = _get(s_c_r_i)
        idx = _indices[_j:]
        if rid is None or isinstance(rid, list) or termini[_residues[rid][-1]]:
            resindex += 1
            resindices[idx] = resindex
            _append(idx)
            if rid is None:
                _set(s_c_r_i, resindex)
            elif isinstance(rid, list):
                rid.append(resindex)
            else:
                _set(s_c_r_i, [rid, resindex])
        else:
            residue = _residues[rid]
            resindices[idx] = rid
            _residues[rid] = concatenate((residue, idx))

        ag._data['resindex'] = resindices

    def getResidue(self, chid, resnum, icode=None, segname=None):
        """Returns residue with number *resnum* and insertion code *icode* from
        the chain with identifier *chid* in segment with name *segname*."""

        try:
            index = self._dict[(segname or self._getSegname(),
                                chid or self._getChid(),
                                resnum, icode or None)]
        except KeyError:
            pass
        else:
            if isinstance(index, list):
                return [r for r in [self._getResidue(i) for i in index]
                        if r is not None]
            else:
                return self._getResidue(index)

    def numResidues(self):
        """Returns number of residues."""

        return (len(self._residues) if self._ag is self._atoms else
                len(self._residues) - count(self._residues, None))

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
        """Returns chain with identifier *chid*, if it is present."""

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
        """Returns number of chains."""

        return (len(self._chains) if self._ag is self._atoms else
                len(self._chains) - count(self._chains, None))

    def getSegment(self, segname):
        """Returns segment with name *segname*, if it is present."""

        try:
            index = self._dict[segname or None]
        except KeyError:
            pass
        else:
            return self._getSegment(index)

    def numSegments(self):
        """Returns number of chains."""

        return (len(self._segments) if self._ag is self._atoms else
                len(self._segments) - count(self._segments, None))

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
