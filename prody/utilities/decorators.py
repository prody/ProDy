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

"""This module defines some decorator functions."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

__all__ = ['checkTypes']

def checkTypes(**types):
    """Decorate a function with a type checker.  In this decorator design, a 
    single line needs to be added inside the function to make the call for
    type checking.  
    
    >>> @checkTypes(n=int, i=int)
    >>> def incr(n, i):
    ...     '''Increment *n* by *i*, both arguments must be integers.'''
    ...     
    ...     globals()['incr'].func_typecheck(locals())
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
    
        def func_typecheck(local, types=types):
            
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
                
        func.func_typecheck = func_typecheck
        return func
    return decorate 
    
    
if __name__ == '__main__':

    @checkTypes(n=int)
    def incr(n, i):
        '''Increment *n* by *i*, both arguments must be integers.'''
        
        globals()['incr'].func_typecheck(locals())
        return n + i
    incr(10, 1)
    incr(10., 1)
