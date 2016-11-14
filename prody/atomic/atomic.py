# -*- coding: utf-8 -*-
"""This module defines base class :class:`Atomic` that all other
:mod:`~prody.atomic` classes are derived from."""

from numpy import all, arange

from prody import LOGGER

from . import flags
from .bond import trimBonds
from .fields import READONLY


__all__ = ['Atomic']

SELECT = None
isSelectionMacro = None
NOTALLNONE = set(['not', 'all', 'none', 'index', 'sequence', 'x', 'y', 'z'])


class Atomic(object):

    """Base class for all atomic classes that can be used for type checking."""

    __slots__ = []

    def __getattribute__(self, name):

        try:
            return object.__getattribute__(self, name)

        except AttributeError:
            if name.startswith('is') and self.isFlagLabel(name[2:]):
                return all(self._getFlags(name[2:]))
            else:
                if name == 'all':
                    try:
                        ag = self.getAtomGroup()
                    except AttributeError:
                        ag = self
                        selstr = name
                        return Selection(ag, arange(self.numAtoms()), 'all',
                                         self._acsi, unique=True)
                    else:
                        try:
                            dummies = self.numDummies()
                        except AttributeError:
                            return Selection(ag, self.getIndices(),
                                             self.getSelstr(),
                                             self._acsi, unique=True)
                        else:
                            return AtomMap(ag, self.getIndices(), self._acsi,
                                           intarrays=True, dummies=dummies,
                                           title=self.getTitle())
                elif name == 'none':
                    return None
                elif self.isFlagLabel(name):
                    try:
                        ag = self.getAtomGroup()
                    except AttributeError:
                        ag = self
                        selstr = name
                    else:
                        selstr = '({0}) and ({1})'.format(name,
                                                          self.getSelstr())
                    try:
                        dummies = self.numDummies()
                    except AttributeError:
                        indices = self._getSubset(name)
                        if len(indices):
                            return Selection(ag, indices, selstr,
                                             self._acsi, unique=True)
                        else:
                            return None
                    else:
                        indices = self._getSubset(name)
                        if len(indices):
                            return AtomMap(ag, indices, self._acsi,
                                           intarrays=True, dummies=dummies,
                                           title='Selection ' + repr(name) +
                                                 ' from ' + str(self))
                        else:
                            return None
                else:
                    selstr = name
                    items = name.split('_')
                    word = items[0]
                    if (self.isFlagLabel(word) or self.isDataLabel(word) or
                       word in NOTALLNONE or isSelectionMacro(word)):
                        selstr = ' '.join(items)
                        return SELECT.select(self, selstr)

        raise AttributeError('{0} object has no attribute `{1}` and {2} '
                             'is not a valid selection string'
                             .format(self.__class__.__name__, name,
                                     repr(selstr)))

    def __getstate__(self):

        return dict([(slot, getattr(self, slot))
                     for slot in self.__class__.__slots__])

    def __setstate__(self, state):

        for slot in self.__class__.__slots__:
            try:
                value = state[slot]
            except KeyError:
                pass
            else:
                setattr(self, slot, value)

    def copy(self):
        """Returns a copy of atoms (and atomic data) in an :class:`.AtomGroup`
        instance."""

        dummies = None
        indices = None
        readonly = False
        try:
            ag = self.getAtomGroup()
        except AttributeError:
            ag = self
            readonly = True
            new = AtomGroup(ag.getTitle())
        else:
            indices = self.getIndices()
            new = AtomGroup(ag.getTitle() + ' ' + str(self))
            try:
                dummies = self.numDummies()
            except AttributeError:
                pass
            else:
                if dummies:
                    dummy = self.getFlags('dummy')
                    mapped = self.getFlags('mapped')

        try:
            self.getIndex()
        except AttributeError:
            this = self
        else:
            this = self.all

        if self.numCoordsets():
            new.setCoords(this.getCoordsets(), label=ag.getCSLabels())

        for label in ag.getDataLabels():
            if label in READONLY:
                if readonly:
                    new._data[label] = this.getData(label)
            else:
                new.setData(label, this.getData(label))

        #if readonly:
        #    for label in READONLY:
        #        data = this.getData(label)
        #        if data is not None:
        #            new._data[label] = data

        skip_flags = set()
        for label in ag.getFlagLabels():
            if label in skip_flags:
                continue
            else:
                new._setFlags(label, this.getFlags(label))
                skip_flags.update(flags.ALIASES.get(label, [label]))

        if dummies:
            new._setFlags('dummy', dummy)
            new._setFlags('mapped', mapped)

        bonds = ag._bonds
        bmap = ag._bmap
        if bonds is not None and bmap is not None:
            if indices is None:
                new._bonds = bonds.copy()
                new._bmap = bmap.copy()
                new._data['numbonds'] = ag._data['numbonds'].copy()
            elif dummies:
                if dummies:
                    indices = indices[self._getMapping()]
                if len(set(indices)) == len(indices):
                    new.setBonds(trimBonds(bonds, indices))
                else:
                    LOGGER.warn('Duplicate atoms in mapping, bonds are '
                                'not copied.')
            else:
                bonds = trimBonds(bonds, indices)
                if bonds is not None:
                    new.setBonds(bonds)
        return new

    __copy__ = copy

    def select(self, selstr, **kwargs):
        """Returns atoms matching *selstr* criteria.  See :mod:`~.select` module
        documentation for details and usage examples."""

        return SELECT.select(self, selstr, **kwargs)
