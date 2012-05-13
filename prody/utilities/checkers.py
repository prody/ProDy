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

from numpy import float32

__all__ = ['checkCoords', 'checkTypes']

COORDS_NDIM = set([2])
CSETS_NDIMS = set([2, 3])

def checkCoords(coords, csets=False, natoms=None, dtype=(float, float32), 
                name='coords'):
    """Return **True** if shape, dimensionality, and data type of *coords*
    array are as expected.
    
    :arg coords: coordinate array 
    
    :arg csets: whether multiple coordinate sets (i.e. ``.ndim in (2, 3)``) are 
        allowed, default is **False**
    
    :arg natoms: number of atoms, if **None** number of atoms is not checked
    
    :arg dtype: allowed data type(s), default is ``(float, numpy.float32)``, 
        if **None** data type is not checked
    
    :arg name: name of the coordinate argument
    
    :raises: :exc:`TypeError` when *coords* is not an instance of 
        :class:`numpy.ndarray`
        
    :raises: :exc:`ValueError` when wrong shape, dimensionality, or data type
        is encountered"""

    try:
        ndim = coords.ndim
        shape = coords.shape
    except AttributeError:
        raise TypeError('coords must be a numpy.ndarray instance')

    ndims = CSETS_NDIMS if csets else COORDS_NDIM    
    if ndim not in ndims: 
        raise ValueError(str(name) + '.ndim must be ' + 
                         ' or '.join([str(d) for d in ndims]))

    elif shape[-1] != 3:
        raise ValueError(str(name) + '.shape[-1] must be 3')
        
    elif natoms and shape[-2] != natoms:
        raise ValueError(str(name) + '.shape[-2] must match number of atoms')
        
    elif dtype: 
        if isinstance(dtype, type) and coords.dtype != dtype:
            raise ValueError(str(name) + '.dtype must be ' + dtype.__name__)
        elif coords.dtype not in dtype:
            if len(dtype) > 1:
                msg = ', '.join([repr(dt.__name__) for dt in dtype[:-1]]
                                ) + ', or ' + repr(dtype[-1].__name__)
            else:
                msg = dtype[0].__name__
            raise ValueError(str(name) + '.dtype must be ' + msg)
  
    return True

def checkTypes(args, **types):
    """Return **True** if types of all *args* match those given in *types*.
   
    :raises: :exc:`TypeError` when type of an argument is not one of allowed 
        types
        
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
