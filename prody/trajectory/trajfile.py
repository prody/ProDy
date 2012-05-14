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

"""This module defines a base class for format specific trajectory classes."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

from trajbase import TrajBase

from prody.utilities import relpath

__all__ = ['TrajFile']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

class TrajFile(TrajBase):
    
    """A base class for providing a consistent interface over trajectories in 
    different formats.  Derived classes are:
    
      * :class:`~.DCDFile`"""
    
    
    def __init__(self, filename, mode='r'):
        """Open *filename* for reading (default, ``mode="r"``), writing 
        (``mode="w"``), or appending (``mode="r+"`` or ``mode="a"``)."""

        if not isinstance(filename, str):
            raise TypeError("filename argument must be a string")
        if not isinstance(mode, str): 
            TypeError('mode argument must be string')
        if not mode in ('r', 'w', 'a', 'r+'): 
            ValueError("mode string must begin with one of 'r', 'w', 'r+', or "
                       "'a'")
        if mode == 'r' and not os.path.isfile(filename):
            raise IOError("[Errno 2] No such file or directory: '{0:s}'"
                          .format(filename))
        self._filename = filename
        if mode in ('a', 'r+'):
            self._file = open(filename, 'r+b')
            self._file.seek(0)
            mode = 'a'
        else:
            self._file = open(filename, mode+'b')
        self._mode = mode
        name = os.path.splitext(os.path.split(filename)[1])[0]
        TrajBase.__init__(self, name)
        self._bytes_per_frame = None
        self._first_byte = None
        self._dtype = np.float32
        
        self._timestep = 1
        self._first_ts = 0
        self._framefreq = 1
        self._n_fixed = 0
        
    def __del__(self):
        
        self._file.close()
    
    def __repr__(self):

        if self._closed:
            return ('<{0:s}: {1:s} (closed)>').format(
                        self.__class__.__name__, self._title)
        if self._mode == 'r':
            if self._sel is None:            
                return ('<{0:s}: {1:s} (next {2:d} of {3:d} frames; '
                        '{4:d} atoms)>').format(
                        self.__class__.__name__, self._title, 
                        self._nfi, self._n_csets, self._n_atoms)
            else:
                return ('<{0:s}: {1:s} (next {2:d} of {3:d} frames; '
                        'selected {4:d} of {5:d} atoms)>').format(
                        self.__class__.__name__, self._title, 
                        self._nfi, self._n_csets, self.numSelected(),
                        self._n_atoms)
        else:
            return ('<{0:s}: {1:s} ({2:d} atoms, {3:d} frames written)>'
                    ).format(
                    self.__class__.__name__, self._title, 
                    self._n_atoms, self._n_csets)
    
    def getFilename(self, absolute=False):
        """Return relative path to the current file. For absolute path,
        pass ``absolute=True`` argument."""
        
        if absolute:
            return os.path.abspath(self._filename)    
        return relpath(self._filename)
    
    def getFrame(self, index):
        """Return frame at given *index*."""
        
        if self._closed:
            raise ValueError('I/O operation on closed file')
        if not isinstance(index, (int, long)):
            raise IndexError('index must be an integer')
        if not 0 <= index < self._n_csets:
            raise IndexError('index must be greater or equal to 0 and less '
                             'than number of frames')
        nfi = self._nfi
        if index > nfi:
            self.skip(index - nfi)
        elif index < nfi:
            self.reset()
            if index > 0:
                self.skip(index)
        return self.next()
    
    getFrame.__doc__ = TrajBase.getFrame.__doc__
                
    def getCoordsets(self, indices=None):
        
        if self._closed: 
            raise ValueError('I/O operation on closed file')
        if indices is None:
            indices = np.arange(self._n_csets)
        elif isinstance(indices, (int, long)):
            indices = np.array([indices])
        elif isinstance(indices, slice):
            indices = np.arange(*indices.indices(self._n_csets))
            indices.sort()
        elif isinstance(indices, (list, np.ndarray)):
            indices = np.unique(indices)
        else:
            raise TypeError('indices must be an integer or a list of integers')

        nfi = self._nfi
        self.reset()

        n_atoms = self.numSelected() 
        coords = np.zeros((len(indices), n_atoms, 3), self._dtype)

        prev = 0
        next = self.nextCoordset
        for i, index in enumerate(indices):
            diff = index - prev
            if diff > 1:
                self.skip(diff-1)
            xyz = next()
            if xyz is None:         
                LOGGER.warning('Expected {0:d} frames, but parsed {1:d}.'
                               .format(len(indices), i))
                self.goto(nfi)
                return coords[:i]
            coords[i] = xyz 
            prev = index

        self.goto(nfi)
        return coords
    
    getCoordsets.__doc__ = TrajBase.getCoordsets.__doc__
    
    def skip(self, n):
     
        if self._closed: 
            raise ValueError('I/O operation on closed file')
        if not isinstance(n, (int, long)):
            raise ValueError('n must be an integer')
        if n > 0:
            left = self._n_csets - self._nfi
            if n > left:
                n = left
            self._file.seek(n * self._bytes_per_frame, 1)                
            self._nfi += n
            
    skip.__doc__ = TrajBase.skip.__doc__
    
    def goto(self, n):
        
        if self._closed:
            raise ValueError('I/O operation on closed file')
        if not isinstance(n, (int, long)):
            raise ValueError('n must be an integer')
        n_csets = self._n_csets
        if n == 0:
            self.reset()
        else:
            if n < 0:
                n = n_csets + n
            if n < 0:
                n = 0
            elif n > n_csets: 
                n = n_csets
            self._file.seek(self._first_byte + n * self._bytes_per_frame)
            self._nfi = n
            
    goto.__doc__ = TrajBase.goto.__doc__  
    
    def reset(self):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        self._file.seek(self._first_byte)
        self._nfi = 0
            
    reset.__doc__ = TrajBase.reset.__doc__  
    
    def close(self):
       
        self._file.close()
        self._nfi = 0
        self._closed = True
    
    close.__doc__ = TrajBase.close.__doc__
    
    def getTimestep(self):
        """Return timestep size."""
        
        return self._timestep
    
    def getFirstTimestep(self):
        """Return first timestep value."""
        
        return self._first_ts
    
    def getFrameFreq(self):
        """Return timesteps between frames."""
        
        return self._framefreq
    
    def numFixed(self):
        """Return number of fixed atoms."""
        
        return self._n_fixed
