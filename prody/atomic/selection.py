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

>>> from prody import *
>>> prot = parsePDB('1p38')

We will use :meth:`~.AtomGroup.select` method to get :class:`Selection` 
instances as follows:
    
>>> water = prot.select('water') 

Union
-------------------------------------------------------------------------------

Let's select β-carbon atoms for non-GLY amino acid residues, and 
α-carbons for GLYs in two steps:

>>> betas = prot.select('name CB and protein')
>>> print( len(betas) )
336
>>> gly_alphas = prot.select('name CA and resname GLY')
>>> print( len(gly_alphas) )
15

The above shows that the p38 structure contains 15 GLY residues.

These two selections can be combined as follows:

>>> betas_gly_alphas = betas | gly_alphas
>>> print( betas_gly_alphas )
Selection '(name CB and pr...nd resname GLY)' from 1p38
>>> print( len(betas_gly_alphas) )
351

The selection string for the union of selections becomes:

>>> print( betas_gly_alphas.getSelstr() )
(name CB and protein) or (name CA and resname GLY)

Note that it is also possible to yield the same selection using selection 
string ``(name CB and protein) or (name CA and resname GLY)``.


Intersection
-------------------------------------------------------------------------------

It is as easy to get the intersection of two selections. Let's find 
charged and medium size residues in a protein:

>>> charged = prot.select('charged')
>>> print( charged )
Selection 'charged' from 1p38
>>> medium = prot.select('medium')
>>> print( medium )
Selection 'medium' from 1p38

>>> medium_charged = medium & charged
>>> print( medium_charged )
Selection '(medium) and (charged)' from 1p38
>>> print( medium_charged.getSelstr() )
(medium) and (charged)

Let's see which amino acids are considered charged and medium:

>>> print( set(medium_charged.getResnames()) )
set(['ASP'])

What about amino acids that are medium or charged:

>>> print( set((medium | charged).getResnames()) )
set(['CYS', 'ASP', 'VAL', 'LYS', 'PRO', 'THR', 'GLU', 'HIS', 'ARG', 'ASN'])


Inversion
-------------------------------------------------------------------------------

It is also possible to invert a selection:

>>> only_protein = prot.select('protein')
>>> print( only_protein )
Selection 'protein' from 1p38
>>> only_non_protein = ~only_protein
>>> print( only_non_protein )
Selection 'not (protein) ' from 1p38

>>> water = prot.select('water')
>>> print( water )
Selection 'water' from 1p38

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

>>> protein = prot.select('protein')
>>> water = prot.select('water')
>>> water_protein = water + protein
>>> writePDB('1p38_water_protein.pdb', water_protein)
'1p38_water_protein.pdb'

In the resulting file, the water atoms will precedes the 
protein atoms.


Membership
-------------------------------------------------------------------------------

Selections also allows membership test operations:

>>> backbone = prot.select('protein') 
>>> calpha = prot.select('calpha')

Is ``calpha`` a subset of ``backbone``?

>>> calpha in backbone
True

Or, is water in protein selection?

>>> water in protein
False

Other tests include:

>>> protein in prot
True
>>> backbone in prot
True
>>> prot in prot
True
>>> calpha in calpha
True

Equality
-------------------------------------------------------------------------------

You can also check the equality of selections. Comparison will return
``True`` if both selections refer to the same atoms.

>>> calpha = prot.select('protein and name CA') 
>>> calpha2 = prot.select('calpha')
>>> calpha == calpha2
True


"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from atom import Atom
from subset import AtomSubset

__all__ = ['Selection']

SELECT = None

ellipsis = lambda s: s[:15] + '...' + s[-15:] if len(s) > 33 else s

class Selection(AtomSubset):
    
    """A class for accessing and manipulating attributes of selection of atoms 
    in an :class:`~.AtomGroup` instance.  Instances can be generated using
    :meth:`~.AtomGroup.select` method.  Following built-in functions are 
    customized for this class:
    
    * :func:`len` returns the number of selected atoms
    * :func:`iter` yields :class:`~.Atom` instances"""
    
    __slots__ = ['_ag', '_indices', '_acsi', '_selstr']
    
    def __init__(self, ag, indices, selstr, acsi=None, **kwargs):
        
        kwargs['selstr'] = selstr
        AtomSubset.__init__(self, ag, indices, acsi, **kwargs)

    def __repr__(self):
        
        n_csets = self._ag.numCoordsets()
        selstr = ellipsis(self._selstr)
        if n_csets:
            if n_csets == 1:
                return ('<Selection: {0:s} from {1:s} ({2:d} atoms)>').format(
                    repr(selstr), self._ag.getTitle(), len(self), n_csets)
            else:
                return ('<Selection: {0:s} from {1:s} ({2:d} atoms; '
                        'active of {3:d} coordsets)>').format(repr(selstr), 
                        self._ag.getTitle(), len(self), self.getACSIndex(), 
                        n_csets)
        else:
            return ('<Selection: {0:s} from {1:s} ({2:d} atoms; no '
                    'coordinates)>').format(repr(selstr), self._ag.getTitle(), 
                    len(self))

    def __str__(self):

        return 'Selection {0:s} from {1:s}'.format(
                repr(ellipsis(self._selstr)), self._ag.getTitle())
    
    def getSelstr(self):
        """Return selection string that selects this atom subset."""
        
        return self._selstr

    def getHierView(self):
        """Return a hierarchical view of the atom selection."""
        
        return HierView(self)

    def update(self):    
        """Update selection."""
        
        self._indices = SELECT.getIndices(self._ag, self._selstr)
