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

"""This module defines :class:`KDTree` class for dealing with atomic coordinate
sets and handling periodic boundary conditions."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'
 
from numpy import array, ndarray, concatenate, mod

from prody import LOGGER

from _CKDTree import KDTree as CKDTree 

__all__ = ['KDTree']

_ = array([-1., 0., 1.])
REPLICATE = array([[x, y, z] for x in _ for y in _ for z in _])

class KDTree(object):
    
    """An interface to Thomas Hamelryck's C KDTree module that can handle 
    periodic boundary conditions.  Both point and pair search are performed 
    using the single :meth:`search` method and results are retrieved using 
    :meth:`getIndices` and :meth:`getDistances`.
    
    **Periodic Boundary Conditions**
    
    *Point search*

    A point search around a *center*, indicated with a question mark (``?``) 
    below, involves making images of the point in cells sharing a wall or an 
    edge with the unitcell that contains the system.  The search is performed 
    for all images of the *center* (27 in 3-dimensional space) and unique 
    indices with the minimum distance from them to the *center* are returned.
    ::

          _____________________________
         |        1|        2|        3|
         |       ? |       ? |      ?  |
         |_________|_________|_________|
         |        4|o  h h  5|        6| ? and H interact in periodic image 4
         |       ?H| h  o  ? |      ?  | but not in the original unitcell (0)   
         |_________|_________|_________| 
         |        7|        8|        9|
         |       ? |       ? |      ?  |
         |_________|_________|_________|
    
    There are two requirements for this approach to work: (i) the *center* must 
    be in the original unitcell, and (ii) the system must be in the original 
    unitcell with parts in its immediate periodic images. 
    
    *Pair search*
    
    A pair search involves making 26 (or 8 in 2-d) replicas of the system 
    coordinates.  A KDTree is built for the system (``O`` and ``H``) and all 
    its replicas (``o`` and ``h``).  After pair search is performed, unique
    pairs of indices and minimum distance between them are returned.
    ::

          _____________________________
         |o  h h  1|o  h h  2|o  h h  3|
        h| h  o   h| h  o   h| h  o    |
         |_________|_________|_________|
         |o  h h  4|O  H H  5|o  h h  6|
        h| h  o   H| H  O   h| h  o    |   
         |_________|_________|_________| 
         |o  h h  7|o  h h  8|o  h h  9|
        h| h  o   h| h  o   h| h  o    |
         |_________|_________|_________|
    
    Only requirement for this approach to work is that the system must be
    in the original unitcell with parts in its immediate periodic images.      
    

    """
    
    def __init__(self, coords, **kwargs):
        """
        :arg coords: coordinate array with shape ``(N, 3)``, where N is number 
            of atoms
        :type coords: :class:`numpy.ndarray`, :class:`.Atomic`, :class:`.Frame`
        
        :arg unitcell: orthorhombic unitcell dimension array with shape 
            ``(3,)``
        :type unitcell: :class:`numpy.ndarray`
        
        :arg bucketsize: number of points per tree node, default is 10
        :type bucketsize: int"""
        
        unitcell = kwargs.get('unitcell')
        if not isinstance(coords, ndarray):
            if unitcell is None:
                try:
                    unitcell = coords.getUnitcell()
                except AttributeError:
                    pass
                else:
                    if unitcell is not None:
                        LOGGER.info('Unitcell information from {0:s} will be '
                                    'used.'.format(str(coords)))
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
        if coords.shape[-1] != 3:
                raise Exception('coords.shape must be (N,3)')
        if coords.min() <= -1e6 or coords.max() >= 1e6:
                raise Exception('coords must be between -1e6 and 1e6')

        self._bucketsize = kwargs.get('bucketsize', 10)
        
        if not isinstance(self._bucketsize, int):
            raise TypeError('bucketsize must be an integer')
        if self._bucketsize < 1:
            raise ValueError('bucketsize must be a positive integer')
            
        self._coords = None
        self._unitcell = None
        self._neighbors = None
        if unitcell is None:
            self._kdtree = CKDTree(3, self._bucketsize)
            self._kdtree.set_data(coords)
        else:
            if not isinstance(unitcell, ndarray):
                raise TypeError('unitcell must be a Numpy array')
            if unitcell.shape != (3,):
                raise ValueError('unitcell.shape must be (3,)')
            self._kdtree = CKDTree(3, self._bucketsize)
            self._kdtree.set_data(coords)
            self._coords = coords
            self._unitcell = unitcell
            self._replicate = REPLICATE * unitcell
            self._kdtree2 = None
            self._pbcdict = {}
            self._pbckeys = []   
            self._n_atoms = coords.shape[0]         
        
    def __call__(self, radius, center=None):
        """Shorthand method for searching and retrieving results."""
        
        self.search(radius, center)       
        return self.getIndices(), self.getDistances()

    def search(self, radius, center=None):
        """Search pairs within *radius* of each other or points within *radius*
        of *center*.
        
        :arg radius: distance (Å)
        :type radius: float

        :arg center: a point in Cartesian coordinate system
        :type center: :class:`numpy.ndarray`"""
        
        if not isinstance(radius, (float, int)):
            raise TypeError('radius must be a number')
        if radius <= 0:
            raise TypeError('radius must be a positive number')
        
        if center is not None:
            if not isinstance(center, ndarray): 
                raise TypeError('center must be a Numpy array instance')
            if center.shape != (3,):
                raise ValueError('center.shape must be (3,)')
            
            if self._unitcell is None:
                self._kdtree.search_center_radius(center, radius)
                self._neighbors = None
                
            else:
                kdtree = self._kdtree
                search = kdtree.search_center_radius
                get_radii = kdtree.get_radii
                get_indices = kdtree.get_indices
                get_count = kdtree.get_count

                _dict = {}
                _dict_get = _dict.get
                _dict_set = _dict.__setitem__
                for center in center + self._replicate:
                    search(center, radius)
                    if get_count():
                        [_dict_set(i, min(r, _dict_get(i, 1e6)))
                         for i, r in zip(get_indices(), get_radii())]
                self._pbcdict = _dict
                self._pdbkeys = _dict.keys() 
                
        else:
            if self._unitcell is None:
                self._neighbors = self._kdtree.neighbor_search(radius)
            else:
                kdtree = self._kdtree2
                if kdtree is None:
                    coords = self._coords
                    coords = concatenate([coords + rep 
                                          for rep in self._replicate])
                    kdtree = CKDTree(3, self._bucketsize)
                    kdtree.set_data(coords)
                    self._kdtree2 = kdtree
                n_atoms = len(self._coords)
                _dict = {}
                neighbors = kdtree.neighbor_search(radius)
                if kdtree.neighbor_get_count():
                    _get = _dict.get
                    _set = _dict.__setitem__
                    
                    for nb in neighbors:
                        i = nb.index1 % n_atoms
                        j = nb.index2 % n_atoms
                        if i < j: 
                            _set((i, j), min(nb.radius, _get((i, j), 1e6)))
                        elif j < i:
                            _set((j, i), min(nb.radius, _get((j, i), 1e6)))
                self._pbcdict = _dict
                self._pdbkeys = _dict.keys() 
                    
    
    def getIndices(self):
        """Return array of indices for points or pairs, depending on the type 
        of the most recent search."""

        if self.getCount():
            if self._unitcell is None:        
                if self._neighbors is None:
                    return self._kdtree.get_indices()
                else:
                    return array([(n.index1, n.index2) 
                                  for n in self._neighbors], int)
            else:
                return array(self._pdbkeys)
            
    
    def getDistances(self):
        """Return array of distances."""

        if self.getCount():        
            if self._unitcell is None:        
                if self._neighbors is None:
                    return self._kdtree.get_radii()
                else:
                    return array([n.radius for n in self._neighbors])
            else:
                _dict = self._pbcdict
                return array([_dict[i] for i in self._pdbkeys])
    
    def getCount(self):
        """Return number of points or pairs."""

        if self._unitcell is None:        
            if self._neighbors is None:
                return self._kdtree.get_count()
            else:
                return self._kdtree.neighbor_get_count()
        else:
            return len(self._pbcdict)
