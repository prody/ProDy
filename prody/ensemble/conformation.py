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

"""This module defines classes handling individual conformations.

  * :class:`Conformation`
  * :class:`PDBConformation`
  
Examples
--------

Results from the example :ref:`pca-xray-calculations` will be used to
illustrate class methods and functions in the module.

>>> from prody import *
>>> ens = loadEnsemble('p38_X-ray.ens.npz')

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from prody.measure import measure

__all__ = ['Conformation', 'PDBConformation']

class Conformation(object):

    """A class to provide methods on a conformation in an ensemble.  
    Instances of this class do not keep coordinate and weights data."""

    __slots__ = ['_ensemble', '_index', '_sel', '_indices']

    def __init__(self, ensemble, index):
        self._ensemble = ensemble
        self._index = index
        self._sel = ensemble.getSelection()
        self._indices = ensemble._getSelIndices()
        
    def __repr__(self):
        return ('<{0:s}: {1:d} from {2:s} (selected {3:d} of {4:d} atoms)>'
               ).format(self.__class__.__name__, self._index, 
                        self._ensemble.getTitle(), self.numSelected(), 
                        self._ensemble.numAtoms())

    def __str__(self):
        return '{0:s} {1:d} from {2:s}'.format(self.__class__.__name__,
                    self._index, self._ensemble.getTitle())

    def numAtoms(self):
        """Return number of atoms."""
        
        return self._ensemble.numAtoms()
    
    def numSelected(self):
        """Return number of selected atoms."""
        
        if self._sel is None:
            return self._ensemble.numAtoms()
        else:
            return len(self._indices)
    
    def getIndex(self):
        """Return index."""
        
        return self._index
    
    def getWeights(self):
        """Return coordinate weights for selected atoms."""
        
        if self._sel is None:
            return self._ensemble.getWeights()
        else:
            return self._ensemble.getWeights()[self._indices]

    def _getWeights(self):
        
        if self._sel is None:
            return self._ensemble._getWeights()
        else:
            return self._ensemble._getWeights()[self._indices]


    def getEnsemble(self):
        """Return the ensemble that this conformation belongs to."""
        
        return self._ensemble
    
    def getCoords(self):
        """Return a copy of the coordinates of the conformation. If a subset
        of atoms are selected in the ensemble, coordinates for selected
        atoms will be returned."""
        
        if self._ensemble._confs is None:
            return None
        indices = self._indices
        if indices is None:
            return self._ensemble._confs[self._index].copy()
        else:
            return self._ensemble._confs[self._index, indices].copy()
    
    def _getCoords(self):

        if self._ensemble._confs is None:
            return None
        indices = self._indices
        if indices is None:
            return self._ensemble._confs[self._index]
        else:
            return self._ensemble._confs[self._index, indices]

    
    def getDeviations(self):
        """Return deviations from the ensemble reference coordinates. 
        Deviations are calculated for selected atoms."""
        
        ensemble = self._ensemble 
        if ensemble._confs is None:
            return None
        indices = ensemble._indices
        if indices is None:
            return ensemble._confs[self._index] - ensemble._coords 
        else:
            return (ensemble._confs[self._index, indices].copy() - 
                    ensemble._coords[indices])

    def getRMSD(self):
        """Return RMSD from the ensemble reference coordinates. RMSD is
        calculated for selected atoms."""
        
        ensemble = self._ensemble 
        indices = self._indices
        if indices is None:
            return measure._calcRMSD(ensemble._coords,
                                     ensemble._confs[self._index], 
                                     self.getWeights())[0]
        else:
            return measure._calcRMSD(ensemble._coords[indices],
                                     ensemble._confs[self._index, indices], 
                                     self.getWeights())[0]
    
class PDBConformation(Conformation):
    
    """This class is the same as :class:`Conformation`, except that the 
    conformation has a name (or identifier), e.g. PDB identifier.
    
    .. versionadded:: 0.8

    >>> conf = ens[0] 
    >>> conf
    <PDB Conformation: 1a9u_ca from p38 X-ray (index: 0; selected 321 of 321 atoms)>

    """
    
    def __repr__(self):
        return ('<PDB Conformation: {0:s} from {1:s} (index: {2:d}; '
                'selected {3:d} of {4:d} atoms)>').format(
                    self._ensemble._labels[self._index], 
                    self._ensemble.getTitle(), self._index, 
                    self.numSelected(), self.numAtoms())
    
    def __str__(self):
        return 'PDB Conformation {0:s} from {1:s}'.format(
                    self._ensemble._labels[self._index], 
                    self._ensemble.getTitle())
    
    def getLabel(self):
        """Return the label of the conformation.
        
        >>> print( conf.getLabel() )
        1a9u_ca"""
        
        return self._ensemble._labels[self._index]
    
    def setLabel(self, label):
        """Set the label of the conformation.
        
        >>> conf.setLabel('1a9u')
        >>> print( conf.getLabel() )
        1a9u"""
        
        self._ensemble._labels[self._index] = str(label)
        
    def getCoords(self):
        """Return a copy of the coordinates of the conformation. If a subset
        of atoms are selected in the ensemble, coordinates for selected
        atoms will be returned.
        
        .. warning:: When there are atoms with weights equal to zero (0),
           their coordinates will be replaced with the coordinates of the
           ensemble reference coordinate set.
        """
        
        ensemble = self._ensemble
        if ensemble._confs is None:
            return None
        indices = ensemble._indices
        index = self._index
        if indices is None:
            coords = ensemble._confs[index].copy()
            which = ensemble._weights[index].flatten()==0
            coords[which] = ensemble._coords[which]
        else: 
            coords = ensemble._confs[index, indices].copy()
            which = ensemble._weights[index, indices].flatten()==0
            coords[which] = ensemble._coords[indices][which]
        return coords
    
    _getCoords = getCoords
    
    def getDeviations(self):
        """Return deviations from the ensemble reference coordinates. 
        Deviations are calculated for selected atoms."""
        
        indices = self._indices
        if indices is None:
            return self.getCoords() - self._ensemble._coords
        else:
            return self.getCoords() - self._ensemble._coords[indices]
    
    def getRMSD(self):
        """Return RMSD from the ensemble reference coordinates. RMSD is
        calculated for selected atoms.
        
        >>> print( conf.getRMSD().round(2) )
        0.74
        """
        
        ensemble = self._ensemble
        index = self._index
        indices = ensemble._indices
        if indices is None:
            return measure._calcRMSD(ensemble._coords,
                                     ensemble._confs[index], 
                                     ensemble._weights[index])
        else:
            return measure._calcRMSD(ensemble._coords[indices],
                                     ensemble._confs[index, indices], 
                                     ensemble._weights[index, indices])
    
    def getWeights(self):
        """Return coordinate weights for selected atoms."""
        
        ensemble = self._ensemble
        indices = ensemble._indices
        if indices is None:
            return ensemble._weights[self._index].copy()
        else:
            return ensemble._weights[self._index, indices]
    
    def _getWeights(self):
        
        ensemble = self._ensemble
        indices = ensemble._indices
        if indices is None:
            return ensemble._weights[self._index]
        else:
            return ensemble._weights[self._index, indices]


