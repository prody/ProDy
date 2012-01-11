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

"""This module defines classes and methods for parsing volumetric data.

.. versionadded:: 0.6.1

Classes
-------
  
  * :class:`Volume`
    
Functions
---------
    
  * :func:`parseDX`
  * :func:`writeDX`


"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

import prody
LOGGER = prody.LOGGER

__all__ = ['Volume', 'parseDX', 'writeDX']

if prody._PY3K:
    range = xrange

class Volume(object):
    
    """
    
    .. versionadded:: 0.6.1
    
    """
    
    def __init__(self, name, array=None):
        self._name = str(name)
        self._array = None
        if array is not None:
            self.setArray(array)
        self._origin = None
        self._spacing = None
        self._scratch = {} # A dictionary for parsers to keep arbitrary data

    def __repr__(self):
        return '<{0:s}>'.format(str(self))
    
    def __str__(self):
        return ('Volume: {0:s} with shape ({1[0]:d},{1[1]:d},{1[2]:d}) at '
                '({2[0]:.2f},{2[1]:.2f},{2[2]:.2f})'.format(self._name,
                self.getShape(), self._origin))
    
    
    def __add__(self, other):
        """Compare grid attributes to determine if grids can be added.
        
        File type specific Grid derivatives need to implement method for 
        adding grid data array.
        
        Before adding data arrays of the grids, four attributes are compared. 
        These are :attr:`state`, :attr:`shape`, :attr:`spacing`, and 
        :attr:`offset`. All of these attributes must have the same values.  
        
        """
        if self._shape != other._shape:
            raise ValueError('both grids must have the same shape')
        elif all(self._spacing != other._spacing):
            raise ValueError('both grids must have the same spacing')
        elif all(self._offset != other._offset):
            raise ValueError('both grids must have the same offset')
        else:
            return True

    
    def getArray(self):
        """Return a copy of 3-dimensional NumPy volume array."""
        
        return self._array.copy()
    
    def setArray(self, array):
        """Set the volume data array.
        
        If origin and spacing of the volume data is not set yet, they are
        set to (0, 0, 0) and (1, 1, 1), respectively.
        
        """
        
        if not isinstance(array, np.ndarray): 
            raise TypeError('array must be an ndarray instance')
        elif array.ndim != 3:
            raise ValueError('array must a 3-dimensional array')
        if array.dtype != float:
            try:
                array = array.astype(float)
            except:
                raise ValueError('array must have data type float')
        self._array = array
        if self._origin is None:
            self._origin = np.zeros(3, float)
        if self._spacing is None:
            self._spacing = np.ones(3, float)
    
    def getOrigin(self):
        """Return the coordinates of the grid lower corner."""
    
        return self._origin.copy()

    def setOrigin(self, origin):
        """Set the coordinates of the grid lower corner."""
        if not isinstance(origin, np.ndarray): 
            try:
                origin = np.array(origin, float)
            except:
                raise TypeError('origin must be an ndarray, a list, or a '
                                'tuple instance')
        if origin.ndim != 1:
            raise ValueError('origin must a 1-dimensional array')
        if origin.dtype != float:
            try:
                origin = origin.astype(float)
            except:
                raise ValueError('origin must have data type float')
        self._origin = origin

    def getShape(self):
        """Return the shape of the volume array."""
        
        return self._array.shape
    
    def getSpacing(self):
        """Return the spacing (resolution) between array elements."""
        
        return self._spacing.copy()

    def setSpacing(self, spacing):
        """Set the spacing (resolution) between array elements."""
        
        if not isinstance(spacing, np.ndarray): 
            try:
                spacing = np.array(spacing, float)
            except:
                raise TypeError('spacing must be an ndarray, a list, or a '
                                'tuple instance')
        if spacing.ndim != 1:
            raise ValueError('spacing must a 1-dimensional array')
        if spacing.dtype != float:
            try:
                spacing = spacing.astype(float)
            except:
                raise ValueError('spacing must have data type float')

        self._spacing = spacing


    def getMetadata(self, key):
        """Return metadata associated with given *key*. 
        
        If data associated with *key* is not found, returns ``None``.
         
        """
        
        return self._scratch.get(key, None)
    
    def setMetadata(self, key, metadata):
        """Store *data* for later access using *key*. 
        
        If metadata associated with *key* exists, it is overwritten.
         
        """
        
        self._scratch[key] = metadata

def parseDX(filename):
    """Parse volume data from OpenDX files.
    
    .. versionadded:: 0.6.1
    
    """

    if not os.path.isfile(filename):
        raise IOError('{0:s} not found'.format(filename))

    volume = Volume(os.path.splitext(os.path.split(filename)[1])[0])

    opendx_file = open(filename)
    opendx = opendx_file.read()
    opendx_file.close()
    
    lindex = opendx.index('data follows') + len('data follows') + 1
    lines = opendx[:lindex].split('\n')
    comments = []
    spacing = np.zeros(3, float)
    for line in lines:
        if line.startswith('#'):
            comments.append(line)
        elif line.startswith('object 1'):
            items = line.strip().split()
            shape = np.array(items[-3:], int)
        elif line.startswith('origin'):
            items = line.strip().split()
            origin = np.array(items[1:], float)
        elif line.startswith('delta'):
            items = line.strip().split()
            spacing += np.array(items[1:], float)
    rindex = -1
    index = -1
    count = 0
    while count < 100:
        if opendx[index].isalpha():
            count = 0
            rindex = index - 1
        else:
            count += 1
        index -= 1
    epilog = opendx[rindex:].strip()
    #volume._offset = origin - volume._spacing / 2
    try:
        array = np.fromstring(opendx[lindex:rindex], dtype=float, 
                              sep=' ').reshape(tuple(shape))
    except:
        raise IOError('There was a problem when parsing the volume data. '
                      'Check file formatting and integrity.')
    volume.setArray(array)
    volume.setOrigin(origin)
    volume.setSpacing(spacing)
    volume.setMetadata('comments', comments)
    volume.setMetadata('epilog', epilog)
    return volume

def writeDX(filename, volume):
    """Write volume data in Open DX format.
    
    .. versionadded:: 0.6.1
    
    """
    
    opendx = open(filename, 'w')
    spacing = volume.getSpacing()
    shape = volume.getShape()
    origin = volume.getOrigin()
    comments = volume.getMetadata('comments')
    if comments is not None and isinstance(comments, list):
        for line in comments:
            line = line.strip()
            if not line.startswith('#'):
                line = '# ' + line
            opendx.write('{0:s}\n'.format(line))
    opendx.write('object 1 class gridpositions counts {0[0]:d} {0[1]:d} '
                 '{0[2]:d}\n'.format(shape))
    opendx.write('origin {0[0]:.9g} {0[1]:.9g} {0[2]:.9g}\n'
                 .format(origin))
    opendx.write('delta {0:.9g} 0 0\n'.format(spacing[0]))
    opendx.write('delta 0 {0:.9g} 0\n'.format(spacing[1]))
    opendx.write('delta 0 0 {0:.9g}\n'.format(spacing[2]))
    opendx.write('object 2 class gridconnections counts {0[0]:d} {0[1]:d} '
                 '{0[2]:d}\n'.format(shape))
    length = np.prod(shape)
    opendx.write('object 3 class array type double rank 0 items {0:d} data'
                 ' follows\n'.format(length))
    
    array = volume.getArray().flatten()
    string = ''
    times = length / 9
    for i in range(times):
        string += ('{0[0]:.9g} {0[1]:.9g} {0[2]:.9g}\n'
                   '{0[3]:.9g} {0[4]:.9g} {0[5]:.9g}\n'
                   '{0[6]:.9g} {0[7]:.9g} {0[8]:.9g}\n'
                  ).format(
                   array[i*9:i*9+9])
    length = length - times * 9
    times = length / 3
    for i in range(times):
        string += '{0[0]:.9g} {0[1]:.9g} {0[2]:.9g}\n'.format(
                    array[i*3:i*3+3])
    length = length - times * 3
    if length == 2:
        string += '{0[0]:.9g} {0[1]:.9g}\n'.format(array[-2:])            
    elif length == 1:
        string += '{0:.9g}\n'.format(array[-1])
    opendx.write(string)
    epilog = volume.getMetadata('epilog')
    if isinstance(epilog, str):
        opendx.write('\n{0:s}'.format(epilog))
    opendx.close()
    return filename

def parseXPLOR(filename):
    pass
    
def writeXPLOR(filename, volume):
    pass    
