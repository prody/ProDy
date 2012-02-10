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

"""This module defines base class :class:`Atomic` that all other 
:mod:`~prody.atomic` classes are derived from."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

__all__ = ['Atomic']

isMacro = lambda none: None
isKeyword = lambda none: None
SELECT = None

class Atomic(object):
    
    """Base class for all atomic classes.  This class can be used for type
    checking:
    
    >>> from prody import *
    >>> ag = parsePDB('1aar')
    >>> isinstance(ag, Atomic)
    True
    >>> prot = ag.select('protein')
    >>> isinstance(prot, Atomic)
    True"""
    
    __slots__ = []
    
    def __getattribute__(self, name):
        
        try:
            return object.__getattribute__(self, name)

        except AttributeError:
            selstr = name
            items = name.split('_')
            word = items[0]
            if (isKeyword(word) or items == 'not' or isMacro(word)):
                selstr = ' '.join(items)
                return SELECT.select(self, selstr)

        raise AttributeError("'{0:s}' object has no attribute '{1:s}' "
                             "and '{2:s}' is not a valid selection string"
                             .format(self.__class__.__name__, name, selstr))

    def select(self, selstr, **kwargs):
        """Return atoms matching *selstr* criteria.
        
        .. seealso:: :mod:`~prody.atomic.select` module documentation for 
        details and usage examples."""
        
        return SELECT.select(self, selstr, **kwargs)
