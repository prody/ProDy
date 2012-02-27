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

"""This module defines classes handling individual conformations."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from prody.measure import getRMSD

__all__ = ['Conformation', 'PDBConformation']

class Conformation(object):

    """A class to provide methods on a conformation in an ensemble.  
    Instances of this class do not keep coordinate and weights data."""

    __slots__ = ['_ensemble', '_index']

    def __init__(self, ensemble, index):
        
        self._ensemble = ensemble
        self._index = index
        
    def __repr__(self):
        
        if self.getSelection() is None:
            return ('<Conformation: {0:d} from {1:s} ({2:d} atoms)>').format(
                    self._index, self._ensemble.getTitle(), self.numAtoms())
        else:
            return ('<Conformation: {0:d} from {1:s} (selected {2:d} of {3:d} '
                    'atoms)>').format(self._index, self._ensemble.getTitle(), 
                    self.numSelected(), self.numAtoms())

    def __str__(self):
        
        return 'Conformation {0:d} from {1:s}'.format(
                self._index, self._ensemble.getTitle())

    def numAtoms(self):
        """Return number of atoms."""
        
        return self._ensemble.numAtoms()
    
    def numSelected(self):
        """Return number of selected atoms."""
        
        return self._ensemble.numSelected()
    
    def getAtoms(self):
        """Return associated atom group."""
        
        return self._ensemble.getAtoms()        
    
    def getSelection(self):
        """Return the current selection. If ``None`` is returned, it means
        that all atoms are selected."""
        
        return self._ensemble.getSelection()        

    def getIndex(self):
        """Return index."""
        
        return self._index
    
    def getWeights(self):
        """Return coordinate weights for (selected) atoms."""
        
        ensemble = self._ensemble
        indices = ensemble._indices
        if indices is None:
            return ensemble.getWeights()
        else:
            return ensemble.getWeights()[indices]

    def _getWeights(self):
        
        ensemble = self._ensemble
        indices = ensemble._indices
        if indices is None:
            return ensemble._getWeights()
        else:
            return ensemble._getWeights()[indices]

    def getEnsemble(self):
        """Return the ensemble that this conformation belongs to."""
        
        return self._ensemble
    
    def getCoords(self):
        """Return a copy of the coordinates of the conformation. If a subset
        of atoms are selected in the ensemble, coordinates for selected
        atoms will be returned."""
        
        ensemble = self._ensemble
        if ensemble._confs is None:
            return None
        indices = ensemble._indices
        if indices is None:
            return ensemble._confs[self._index].copy()
        else:
            return ensemble._confs[self._index, indices].copy()
    
    def _getCoords(self):

        ensemble = self._ensemble
        if ensemble._confs is None:
            return None
        indices = ensemble._indices
        if indices is None:
            return ensemble._confs[self._index]
        else:
            return ensemble._confs[self._index, indices]
    
    def getDeviations(self):
        """Return deviations from the ensemble reference coordinates. 
        Deviations are calculated for (selected) atoms."""
        
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
        calculated for (selected) atoms."""
        
        ensemble = self._ensemble 
        indices = ensemble._indices
        if indices is None:
            return getRMSD(ensemble._coords, 
                           ensemble._confs[self._index], 
                           self.getWeights())[0]
        else:
            return getRMSD(ensemble._coords[indices],
                           ensemble._confs[self._index, indices], 
                           self.getWeights())[0]
    
class PDBConformation(Conformation):
    
    """This class is the same as :class:`Conformation`, except that the 
    conformation has a name (or identifier), e.g. PDB identifier."""
    
    def __repr__(self):
        
        if self.getSelection() is None:
            return ('<PDBConformation: {0:s} from {1:s} (index: {2:d}; '
                    '{3:d} atoms)>').format(
                    self._ensemble._labels[self._index], 
                    self._ensemble.getTitle(), self._index, self.numAtoms())
        else:
            return ('<PDBConformation: {0:s} from {1:s} (index: {2:d}; '
                    'selected {3:d} of {4:d} atoms)>').format(
                        self._ensemble._labels[self._index], 
                        self._ensemble.getTitle(), self._index, 
                        self.numSelected(), self.numAtoms())
    
    def __str__(self):
    
        return 'PDBConformation {0:s} from {1:s}'.format(
                    self._ensemble._labels[self._index], 
                    self._ensemble.getTitle())
    
    def getLabel(self):
        """Return the label of the conformation."""
        
        return self._ensemble._labels[self._index]
    
    def setLabel(self, label):
        """Set the label of the conformation."""
        
        self._ensemble._labels[self._index] = str(label)
        
    def getCoords(self):
        """Return a copy of the coordinates of the conformation. If a subset
        of atoms are selected in the ensemble, coordinates for selected
        atoms will be returned.
        
        .. warning:: When there are atoms with weights equal to zero (0),
           their coordinates will be replaced with the coordinates of the
           ensemble reference coordinate set."""
        
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
        Deviations are calculated for (selected) atoms."""
        
        ensemble = self._ensemble
        indices = ensemble._indices
        if indices is None:
            return self.getCoords() - ensemble._coords
        else:
            return self.getCoords() - ensemble._coords[indices]
    
    def getRMSD(self):
        """Return RMSD from the ensemble reference coordinates. RMSD is
        calculated for (selected) atoms."""
        
        ensemble = self._ensemble
        index = self._index
        indices = ensemble._indices
        if indices is None:
            return getRMSD(ensemble._coords,
                           ensemble._confs[index], 
                           ensemble._weights[index])
        else:
            return getRMSD(ensemble._coords[indices],
                           ensemble._confs[index, indices], 
                           ensemble._weights[index, indices])
    
    def getWeights(self):
        """Return coordinate weights for (selected) atoms."""
        
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


