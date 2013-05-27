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

"""This module defines :class:`Selection` class for handling arbitrary subsets
of atom.

.. _selection-operations:

Operations on Selections
===============================================================================

Before reading this section, familiarity with :ref:`selections` can be helpful.

Let's import all classes and functions from ProDy and parse coordinates from
the PDB structure 1p38:

.. ipython:: python

   from prody import *
   prot = parsePDB('1p38')

We will use :meth:`~.AtomGroup.select` method to get :class:`Selection`
instances as follows:

.. ipython:: python

   water = prot.select('water')

Union
-------------------------------------------------------------------------------

Let's select β-carbon atoms for non-GLY amino acid residues, and
α-carbons for GLYs in two steps:

.. ipython:: python

   betas = prot.select('name CB and protein')
   len(betas)
   gly_alphas = prot.select('name CA and resname GLY')
   len(gly_alphas)

The above shows that the p38 structure contains 15 GLY residues.

These two selections can be combined as follows:

.. ipython:: python

   betas_gly_alphas = betas | gly_alphas
   betas_gly_alphas
   len(betas_gly_alphas)

The selection string for the union of selections becomes:

.. ipython:: python

   betas_gly_alphas.getSelstr()

Note that it is also possible to yield the same selection using selection
string ``(name CB and protein) or (name CA and resname GLY)``.


Intersection
-------------------------------------------------------------------------------

It is as easy to get the intersection of two selections. Let's find
charged and medium size residues in a protein:

.. ipython:: python

   charged = prot.select('charged')
   charged
   medium = prot.select('medium')
   medium

.. ipython:: python

   medium_charged = medium & charged
   medium_charged
   medium_charged.getSelstr()

Let's see which amino acids are considered charged and medium:

.. ipython:: python

   set(medium_charged.getResnames())

What about amino acids that are medium or charged:

.. ipython:: python

   set((medium | charged).getResnames())


Inversion
-------------------------------------------------------------------------------

It is also possible to invert a selection:

.. ipython:: python

   only_protein = prot.select('protein')
   only_protein
   only_non_protein = ~only_protein
   only_non_protein
   water = prot.select('water')
   water

The above shows that 1p38 does not contain any non-water
hetero atoms.

Addition
-------------------------------------------------------------------------------

Another operation defined on the :class:`~.Select` object is addition
(also on other :class:`~.AtomPointer` derived classes).

This may be useful if you want to yield atoms in an :class:`~.AtomGroup` in a
specific order.
Let's think of a simple case, where we want to output atoms in 1p38 in a
specific order:

.. ipython:: python

   protein = prot.select('protein')
   water = prot.select('water')
   water_protein = water + protein
   writePDB('1p38_water_protein.pdb', water_protein)

In the resulting file, the water atoms will precedes the
protein atoms.


Membership
-------------------------------------------------------------------------------

Selections also allows membership test operations:

.. ipython:: python

   backbone = prot.select('protein')
   calpha = prot.select('calpha')

Is :term:`calpha` a subset of :term:`backbone`?

.. ipython:: python

   calpha in backbone

Or, is water in protein selection?

.. ipython:: python

   water in protein

Other tests include:

.. ipython:: python

   protein in prot
   backbone in prot
   prot in prot
   calpha in calpha


Equality
-------------------------------------------------------------------------------

You can also check the equality of selections. Comparison will return
``True`` if both selections refer to the same atoms.

.. ipython:: python

   calpha = prot.select('protein and name CA')
   calpha2 = prot.select('calpha')
   calpha == calpha2"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

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
        """Return selection string that selects this atom subset."""

        return self._selstr

    def getHierView(self):
        """Return a hierarchical view of the atom selection."""

        return HierView(self)

    def update(self):
        """Update selection."""

        self._indices = SELECT.getIndices(self._ag, self._selstr)
