"""This module defines classes handling individual conformations."""

from prody.measure import getRMSD
from prody.measure import Transformation

__all__ = ['Conformation', 'PDBConformation']

class Conformation(object):

    """A class to provide methods on a conformation in an ensemble.
    Instances of this class do not keep coordinate and weights data."""

    __slots__ = ['_ensemble', '_index']

    def __init__(self, ensemble, index):

        self._ensemble = ensemble
        self._index = index

    def __repr__(self):

        if self._ensemble._indices is None:
            return ('<Conformation: {0} from {1} ({2} atoms)>').format(
                    self._index, self._ensemble.getTitle(), self.numAtoms())
        else:
            return ('<Conformation: {0} from {1} (selected {2} of {3} '
                    'atoms)>').format(self._index, self._ensemble.getTitle(),
                    self.numSelected(), self.numAtoms())

    def __str__(self):

        return 'Conformation {0} from {1}'.format(
                self._index, self._ensemble.getTitle())

    def numAtoms(self):
        """Returns number of atoms."""

        return self._ensemble.numAtoms()

    def numSelected(self):
        """Returns number of selected atoms."""

        return self._ensemble.numSelected()

    def getAtoms(self):
        """Returns associated atom group."""

        return self._ensemble.getAtoms()

    def getIndex(self):
        """Returns conformation index."""

        return self._index

    def getWeights(self):
        """Returns coordinate weights for (selected) atoms."""

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
        """Returns the ensemble that this conformation belongs to."""

        return self._ensemble

    def getCoords(self):
        """Returns a copy of the coordinates of the conformation. If a subset
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
        """Returns deviations from the ensemble reference coordinates.
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
        """Returns RMSD from the ensemble reference coordinates. RMSD is
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

        if self.numSelected() == self.numAtoms():
            return ('<PDBConformation: {0} from {1} (index: {2}; '
                    '{3} atoms)>').format(
                    self._ensemble._labels[self._index],
                    self._ensemble.getTitle(), self._index, self.numAtoms())
        else:
            return ('<PDBConformation: {0} from {1} (index: {2}; '
                    'selected {3} of {4} atoms)>').format(
                        self._ensemble._labels[self._index],
                        self._ensemble.getTitle(), self._index,
                        self.numSelected(), self.numAtoms())

    def __str__(self):

        return 'PDBConformation {0} from {1}'.format(
                    self._ensemble._labels[self._index],
                    self._ensemble.getTitle())

    def getLabel(self):
        """Returns the label of the conformation."""

        return self._ensemble._labels[self._index]

    def setLabel(self, label):
        """Set the label of the conformation."""

        self._ensemble._labels[self._index] = str(label)

    def getCoords(self):
        """Returns a copy of the coordinates of the conformation. If a subset
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
        """Returns deviations from the ensemble reference coordinates.
        Deviations are calculated for (selected) atoms."""

        ensemble = self._ensemble
        indices = ensemble._indices
        if indices is None:
            return self.getCoords() - ensemble._coords
        else:
            return self.getCoords() - ensemble._coords[indices]

    def getRMSD(self):
        """Returns RMSD from the ensemble reference coordinates. RMSD is
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
        """Returns coordinate weights for (selected) atoms."""

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

    def getTransformation(self):
        """Returns the :class:`~.Transformation` used to superpose this
        conformation onto reference coordinates.  The transformation can
        be used to superpose original PDB file onto the reference PDB file."""

        if self._ensemble._trans is not None:
            return Transformation(self._ensemble._trans[self._index])

