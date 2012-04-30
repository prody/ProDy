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

"""This module defines :class:`KDTree` class that can handle periodic boundary
conditions."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'
 
from numpy import array, ndarray

from _CKDTree import KDTree as CKDTree 

__all__ = ['KDTree']

class KDTree(object):
    
    """This is a Python interface to Thomas Hamelryck's KDTree distributed
    with Biopython."""
    
    def __init__(self, coords, unitcell=None, bucket_size=1):
        """
        :arg coords: coordinate array with shape ``(N, 3)``, where N is number 
            of atoms
        :type coords: :class:`numpy.ndarray`, :class:`.Atomic`
        
        :arg unitcell: unitcell array with shape ``(3,)``
        :type unitcell: :class:`numpy.ndarray`
        
        :arg bucket_size: number of points per tree node
        :type bucket_size: int"""
        
        if not isinstance(coords, ndarray):
            try:
                # using getCoords() because coords will be stored internally
                # and reused when needed, this will avoid unexpected results
                # due to changes made to coordinates externally 
                coords = coords.getCoords()
            except AttributeError: 
                raise TypeError('coords must be a Numpy array or must have '
                                'getCoords attribute')
        else:
            coords = coords.copy()
        
        if coords.ndim != 2:
                raise Exception('coords.ndim must be 2')
        if coords.min() <= -1e6 or coords.max() >= 1e6:
                raise Exception('coords must be between -1e6 and 1e6')
            
        self._point = None
        self._coords = None
        self._neighbors = None
        if unitcell is None:
            self._dim = coords.shape[-1]
            self._kdtree = CKDTree(self._dim, bucket_size)
            self._kdtree.set_data(coords)
        else:
            if not isinstance(self._unitcell, ndarray):
                raise TypeError('unitcell must be a Numpy array')
            if unitcell.shape != (3,):
                raise ValueError('unitcell.shape must be (3,)')
            self._coords = coords
            self._unitcell = unitcell            
        
    def __call__(self, radius, point=None):
        """Shorthand method for searching and retrieving results."""
        
        self.search(radius, point)       
        return self.getIndices(), self.getDistances()

    def search(self, radius, point=None):
        """Search pairs within *radius* of each other or within *radius* of
        *point*.
        
        :arg radius: distance (Å)
        :type radius: float

        :arg point: a point in Cartesian coordinate system
        :type point: :class:`numpy.ndarray`"""
        
        self._point = point
        if point is None:
            self._neighbors = self._kdtree.neighbor_search(radius)
        else:
            if not isinstance(point, ndarray): 
                raise TypeError('point must be a Numpy array instance')
            if point.shape != (self._dim,):
                raise ValueError('point.shape must be ' + repr((self._dim,)))
            self._kdtree.search_center_radius(point, radius)
    
    def getIndices(self):
        """Return array of indices or list of pairs of indices, depending on
        the type of the most recent search."""
        
        if self._point is not None:
            return self._kdtree.get_indices()
        else:
            return array([(n.index1, n.index2) for n in self._neighbors], int)
    
    def getDistances(self):
        """Return array of distances."""
        
        if self._point is not None:
            return self._kdtree.get_radii()
        else:
            return array([n.radius for n in self._neighbors])
    
    def getCount(self:
        """Return number of neighbors."""
        
        if self._point is not None:
            return self._kdtree.get_count()
        else:
            return self._kdtree.get_neighbor_count()

