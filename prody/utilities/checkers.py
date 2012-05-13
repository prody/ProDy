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

"""This module defines functions for type, value, and/or attribute checking."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import array

__all__ = ['checkCoords', 'checkTypes']


def checkCoords(coords, name='coords', cset=False, n_atoms=None, reshape=None, 
                dtype=float):
    """Return a copy of *coords* as a :class:`numpy.ndarray` instance after its
    dimensionality and shape are checked for fitness to be a coordinate set.
    
    :arg coords: coordinate data, any array like object 
    
    :arg name: name of *coords* argument
    
    :arg cset: allow multiple coordinate sets, default is **False**
    
    :arg n_atoms: number of atoms
    
    :arg reshape: reshape *coords* to multiple coordinate set  
    
    :arg dtype: allowed data types, default is :func:`float`
    
    :returns: coordinate array
    
    .. note::  Note that the array that this function returns is not the array 
        that is given as input.  This helps ProDy classes to internalize the
        data.
    """

    if not isinstance(dtype, tuple):
        dtype = (dtype, )

    coords = array(coords, dtype=dtype[0])
    ndim = coords.ndim
    shape = coords.shape
    
    if cset and ndim not in (2,3): 
        raise ValueError(str(name) + '.ndim must be 2 or 3')
        
    elif not cset and ndim != 2:
        raise ValueError(str(name) + '.ndim must be 2')
        
    elif shape[-1] != 3:
        raise ValueError(str(name) + '.shape[-1] must be 3')
        
    if n_atoms and shape[-2] != n_atoms:
        raise ValueError(str(name) + ' size do not match number of atoms')
        
    if cset and reshape and ndim == 2:
        coords = coords.reshape([1, shape[0], 3])
        
    return coords


def checkTypes(args, **types):
    """Return **True** if types of *args* match those given in *types*, 
    otherwise raise a :class:`TypeError`.
    
    ::
        
        def incr(n, i):
            '''Return sum of *n* and *i*.'''
            
            checkTypes(locals(), n=(float, int), i=(float, int))
            return n + i"""
        
    for arg, allowed in types.iteritems():
        if arg in args and not isinstance(args[arg], types[arg]):

            val = args[arg]
            if isinstance(allowed, (list, tuple)):
                if len(allowed) > 1:
                    tstr = ', '.join([repr(tp.__name__) for tp in allowed[:-1]]
                                     ) + ', or ' + repr(allowed[-1].__name__)
                
                else:
                    tstr = repr(allowed[0].__name__)
            
            else:
                tstr = repr(allowed.__name__)
            
            raise TypeError('{0:s} must be an instance of {1:s}, not {2:s}'
                            .format(repr(arg), tstr, repr(type(val).__name__)))

    return True
