# -*- coding: utf-8 -*-
"""This module defines :class:`Selection` class for handling arbitrary subsets
of atom."""

from .subset import AtomSubset

__all__ = ['Selection']

SELECT = None

ellipsis = lambda s: s[:15] + '...' + s[-15:] if len(s) > 33 else s


class Selection(AtomSubset):

    """A class for accessing and manipulating attributes of selection of atoms
    in an :class:`.AtomGroup` instance.  Instances can be generated using
    :meth:`~.AtomGroup.select` method.  Following built-in functions are
    customized for this class:

    * :func:`len` returns the number of selected atoms
    * :func:`iter` yields :class:`.Atom` instances"""

    __slots__ = ['_ag', '_indices', '_acsi', '_selstr']

    def __init__(self, ag, indices, selstr, acsi=None, **kwargs):

        kwargs['selstr'] = selstr
        AtomSubset.__init__(self, ag, indices, acsi, **kwargs)

    def __repr__(self):

        n_csets = self._ag.numCoordsets()
        selstr = ellipsis(self._selstr)
        if n_csets:
            if n_csets == 1:
                return ('<Selection: {0} from {1} ({2} atoms)>'
                        ).format(repr(selstr), self._ag.getTitle(),
                                 len(self), n_csets)
            else:
                return ('<Selection: {0} from {1} ({2} atoms; '
                        'active #{3} of {4} coordsets)>'
                        ).format(repr(selstr), self._ag.getTitle(), len(self),
                                 self.getACSIndex(), n_csets)
        else:
            return ('<Selection: {0} from {1} ({2} atoms; no '
                    'coordinates)>').format(repr(selstr), self._ag.getTitle(),
                                            len(self))

    def __str__(self):

        return 'Selection {0}'.format(repr(ellipsis(self._selstr)))

    def getSelstr(self):
        """Returns selection string that selects this atom subset."""

        return self._selstr

    def getHierView(self, **kwargs):
        """Returns a hierarchical view of the atom selection."""

        return HierView(self, **kwargs)

    def update(self):
        """Update selection."""

        self._indices = SELECT.getIndices(self._ag, self._selstr)
