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

from numpy import ndarray

__all__ = ['checkCoords', 'checkTypes']


def checkCoords(array, name='array', cset=False, n_atoms=None, reshape=None, 
                dtype=(float,)):
    """Check whether *array* is a coordinate set.  An exception is raised if 
    the array is not suitable.
    
    :arg array: coordinate array
    :type array: :class:`numpy.ndarray`
    
    :arg name: argument name used in the function where array is being checked
    :type name: str
    
    :arg cset: whether the arra
    :type cset: bool
    
    :arg n_atoms: number of atoms, when available is used to check the shape of
        the array
    :type n_atoms: int
    
    :arg reshape: when *cset* is **True**, reshape array
    :type reshape: bool
    
    :arg dtype: acceptable data types for array elements, default is 
        :func:`float`
    :type dtype: type
    
    :returns: coordinate array is returned, after its shape or element data 
        type is changed when needed"""

    assert isinstance(name, str), 'name must be a string'
    assert isinstance(cset, bool), 'cset must be a boolean'
    assert n_atoms is None or isinstance(n_atoms, int) and n_atoms >= 0, \
        'n_atoms must be a positive integer'
    assert reshape is None or isinstance(reshape, bool), \
        'reshape must be a boolean'
    if not isinstance(dtype, tuple):
        dtype = (dtype, )

    if not isinstance(array, ndarray):
        raise TypeError(name + ' must be a Numpy array')
        
    elif cset and array.ndim not in (2,3): 
        raise ValueError(name + '.ndim must be 2 or 3')
        
    elif not cset and array.ndim != 2:
        raise ValueError(name + '.ndim must be 2')
        
    elif array.shape[-1] != 3:
        raise ValueError(name + '.shape[-1] of 3, i.e. ([n_csets,]n_atoms,3)')
        
    if n_atoms is not None and n_atoms != 0 and array.shape[-2] != n_atoms:
        raise ValueError(name + ' size do not match number of atoms')
        
    if array.dtype not in dtype:
        try:
            array = array.astype(dtype[0])
        except ValueError:
            raise ValueError(name + '.astype(' + str(dtype[0]) + ') fails, '
                            'float type could not be assigned')
                            
    if cset and reshape and array.ndim == 2:
        array = array.reshape([1, array.shape[0], 3])
        
    return array


def checkTypes(**types):
    """Decorate a function with a type checker.  In this decorator design, a 
    single line needs to be added inside the function to make the call for
    type checking.  
    
    >>> from prody.utilities import checkTypes
    >>> @checkTypes(n=int, i=int)
    ... def incr(n, i):
    ...     '''Increment *n* by *i*, both arguments must be integers.'''
    ...     
    ...     globals()['incr'].checktypes(locals())
    ...     return n + i
    >>> incr(10, 1)
    11
    >>> incr(10., 1)
    Traceback (most recent call last):
      File "checkers.py", line 81, in <module>
        incr(10., 1)
      File "checkers.py", line 78, in incr
        globals()['incr'].func_typecheck(locals())
      File "checkers.py", line 65, in func_typecheck
        .format(arg, tstr, type(val).__name__))
    TypeError: `n` must be int, not float"""

    def decorate(func):
    
        def checktypes(local, types=types):
            
            for arg, val in local.iteritems():
                
                if arg in types and not isinstance(val, types[arg]):
                    
                    tps = types[arg]
                    if isinstance(tps, (list, tuple)):
                        if len(tps) > 1:
                            tstr = 'one of ' + ', '.join([tp.__name__ 
                                                for tp in tps[:-1]]
                                                ) + ', or ' + tps[-1].__name__
                        
                        else:
                            tstr = tps[0].__name__
                    
                    else:
                        tstr = tps.__name__
                    
                    raise TypeError('`{0:s}` must be an instance of {1:s}, '
                            'not {2:s}'.format(arg, tstr, type(val).__name__))
                
        func.checktypes = checktypes
        return func
    return decorate 
    
    
if __name__ == '__main__':

    @checkTypes(n=int)
    def incr(n, i):
        '''Increment *n* by *i*, both arguments must be integers.'''
        
        globals()['incr'].checktypes(locals())
        return n + i
    incr(10, 1)
    incr(10., 1)
