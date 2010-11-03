# ProDy: A Python Package for Protein Structural Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
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

""":mod:`atommap` module defines a class to map atom data with missing atoms.

Classes:
    
  * :class:`AtomMap`
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np


from . import DTYPES, Atom, AtomGroup, AtomSubset
from .select import ProDyAtomSelect as SELECT

__all__ = ['AtomMap']

class AtomMap(object):
    """A class for mapping atomic data.
    
    This class stores a reference to an :class:`AtomGroup` instance, a set of 
    atom indices, active coordinate set index, mapping for indices, and
    indices of unmapped atoms.
    
    """
    
    __slots__ = ['_ag', '_indices', '_acsi', '_name', '_mapping', '_unmapped',
                 '_len']
    
    def __init__(self, atomgroup, indices, mapping, unmapped, name='Unnamed', acsi=None):
        """Instantiate with an AtomMap with following arguments:        
        
        :arg atomgroup: the atomgroup instance from which atoms are mapped
        :arg indices: indices of mapped atoms
        :arg mapping: mapping of the atoms as a list of indices
        :arg unmapped: list of indices for unmapped atoms
        :arg acsi: active coordinate set index, if ``None`` defaults to that of *atomgrup*
        :arg name: name of the AtomMap instance
        
        Length of *mapping* must be equal to length of *indices*. Number of 
        atoms (including unmapped dummy atoms) are determined from the 
        sum of lengths of *mapping* and *unmapped* arrays.         
        
        """
        if not isinstance(atomgroup, AtomGroup):
            raise TypeError('atomgroup must be AtomGroup, not {0:s}'
                            .format(type(atomgroup)))
            
        self._ag = atomgroup

        if not isinstance(indices, np.ndarray):
            self._indices = np.array(indices, np.int64)
        elif not indices.dtype == np.int64:
            self._indices = indices.astype(np.int64)
        else:
            self._indices = indices

        if not isinstance(mapping, np.ndarray):
            self._mapping = np.array(mapping, np.int64)
        elif not mapping.dtype == np.int64:
            self._mapping = mapping.astype(np.int64)
        else:
            self._mapping = mapping

        if not isinstance(unmapped, np.ndarray):
            self._unmapped = np.array(unmapped, np.int64)
        elif not unmapped.dtype == np.int64:
            self._unmapped = unmapped.astype(np.int64)
        else:
            self._unmapped = unmapped
        
        self._name = str(name)
        
        if acsi is None:
            self._acsi = atomgroup.getActiveCoordsetIndex()
        else:
            self._acsi = int(acsi)

        self._len = len(self._unmapped) + len(self._mapping)

        
    def __repr__(self):
        return ('<AtomMap: {0:s} (from {1:s}; {2:d} atoms; '
                '{3:d} mapped; {4:d} unmapped; {5:d} coordinate sets, '
                'active set index: {6:d})>').format(self._name,
                self._ag.getName(), self._len, len(self._mapping), 
                len(self._unmapped), self.getNumOfCoordsets(), self._acsi)
    
    def __str__(self):
        return 'AtomMap {0:s}'.format(self._name)
    
    def __iter__(self):
        indices = np.zeros(self._len, np.int64)
        indices[self._unmapped] = -1
        indices[self._mapping] = self._indices
        ag = self._ag
        acsi = self._acsi
        for index in indices:
            if index > -1:
                yield Atom(ag, index, acsi)
            else:
                yield None
    
    def __len__(self):
        return self._len
    
    def __add__(self, other):
        if not isinstance(other, AtomMap):
            raise TypeError('other must be an AtomMap instance')
        if self._ag != other._ag:
            raise ValueError('both AtomMaps must be from the same AtomGroup')
        acsi = self._acsi
        if acsi != other._acsi:
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero')
            acsi = 0
        indices = np.concatenate((self._indices, other._indices))
        #if isinstance(other, AtomMap): 
        name = '({0:s}) + ({1:s})'.format(self._name, other._name)
        mapping = np.concatenate((self._mapping, other._mapping + self._len))
        unmapped = np.concatenate((self._unmapped, other._unmapped + self._len))
        #else:
        #    name = '({0:s}) + ({1:s})'.format(self._name, other.getSelectionString())
        #    mapping = np.concatenate((self._mapping, np.arange(len(other)) + self._len))
        #    unmapped = self._unmapped.copy()
            
        return AtomMap(self._ag, indices, mapping, unmapped, name, acsi)
    
    def getName(self):
        return self._name
    
    def setName(self, name):
        self._name = name


    def getAtomGroup(self):
        """Return the atom group from which the atoms are mapped."""
        return self._ag
    
    def getNumOfAtoms(self):
        """Return number of mapped atoms."""
        return self._len

    def getNumOfUnmapped(self):
        """Return number of unmapped atoms."""
        return len(self._unmapped)

    def getNumOfMapped(self):
        """Return number of mapped atoms."""
        return len(self._mapping)

    def getIndices(self):
        """Return indices of mapped atoms."""
        return self._indices.copy()

    def getMapping(self):
        """Return mapping of indices."""
        return self._mapping.copy()


    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each coordinate set."""
        for i in range(self._ag._n_coordsets):
            coordinates = np.zeros((self._len, 3), DTYPES['coordinates'])
            coordinates[self._mapping] = self._ag._coordinates[i, self._indices] 
            yield coordinates

    def getNumOfCoordsets(self):
        """Return number of coordinate sets."""
        return self._ag._n_coordsets
    
    def getActiveCoordsetIndex(self):
        """Return the index of the active coordinate set."""
        return self._acsi
    
    def setActiveCoordsetIndex(self, index):
        """Set the index of the active coordinate set."""
        if self._ag._coordinates is None:
            return
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if self._ag._n_coordsets <= index or \
           self._ag._n_coordsets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += self._ag._n_coordsets
        self._acsi = index
    
    def getCoordinates(self):
        """Return coordinates from the active coordinate set."""
        coordinates = np.zeros((self._len, 3), DTYPES['coordinates'])
        coordinates[self._mapping] = self._ag._coordinates[self._acsi, self._indices] 
        return coordinates
    
    def setCoordinates(self, coordinates):
        """Set coordinates in the active coordinate set."""
        self._ag._coordinates[self._acsi, self._indices] = coordinates
    
    def getCoordsets(self, indices):
        """Return coordinate sets at given indices.
        
        *indices* may be an integer or a list of integers.
        
        """
        if self._ag._coordinates is None:
            return None
        try: 
            return self._ag._coordinates[indices, self._indices].copy()
        except IndexError:
            raise IndexError('indices may be an integer or a list of integers')

    def getAtomNames(self):
        """Return atom names of the atoms."""
        if self._ag._atomnames is None:
            raise AttributeError('atom names are not set')
        atomnames = np.zeros(self._len, DTYPES['atomnames'])
        atomnames[self._mapping] = self._ag._atomnames[self._indices]
        atomnames[self._unmapped] = 'DUMMY'
        return atomnames
    
    def getAlternateLocationIndicators(self):
        """Return alternate location indicators of the atoms."""
        if self._ag._altlocs is None:
            return None
        altlocs = np.zeros(self._len, DTYPES['altlocs'])
        altlocs[self._mapping] = self._ag._altlocs[self._indices]
        return altlocs
    
    def getAnisotropicTemperatureFactors(self):
        """Return anisotropic temperature factors of the atoms."""
        if self._ag._anisou is None:
            return None
        anisou = np.zeros((self._len, 3), DTYPES['anisou'])
        anisou[self._mapping] = self._ag._anisou[self._indices]
        return anisou
    
    def getChainIdentifiers(self):
        """Return chain identifiers of the atoms."""
        if self._ag._chainids is None:
            raise AttributeError('chain identifiers are not set')
        chainids = np.zeros(self._len, DTYPES['chainids'])
        chainids[self._mapping] = self._ag._chainids[self._indices]
        chainids[self._unmapped] = ' '
        return chainids
    
    def getElementSymbols(self):
        """Return element symbols of the atoms."""
        if self._ag._elements is None:
            return None
        elements = np.zeros(self._len, DTYPES['elements'])
        elements[self._mapping] = self._ag._elements[self._indices]
        elements[self._unmapped] = 'Du'
        return elements
    
    def getHeteroFlags(self):
        """Return hetero flags of the atoms."""
        if self._ag._hetero is None:
            return None
        hetero = np.zeros(self._len, DTYPES['hetero'])
        hetero[self._mapping] = self._ag._hetero[self._indices]
        return hetero
    
    def getOccupancies(self):
        """Return occupancies of the atoms."""
        if self._ag._occupancies is None:
            return None
        occupancies = np.zeros(self._len, DTYPES['occupancies'])
        occupancies[self._mapping] = self._ag._occupancies[self._indices]
        return occupancies
    
    def getResidueNames(self):
        """Return residue names of the atoms."""
        if self._ag._resnames is None:
            return None
        resnames = np.zeros(self._len, DTYPES['resnames'])
        resnames[self._mapping] = self._ag._resnames[self._indices]
        resnames[self._unmapped] = 'NONE'
        return resnames
    
    def getResidueNumbers(self):
        """Return residue numbers of the atoms."""
        if self._ag._resnums is None:
            return None
        resnums = np.zeros(self._len, DTYPES['resnums'])
        resnums[self._mapping] = self._ag._resnums[self._indices]
        return resnums
    
    def getSecondaryStructureAssignments(self):
        """Return secondary structure assignments of the atoms."""
        if self._ag._secondary is None:
            return None
        secondary = np.zeros(self._len, DTYPES['secondary'])
        secondary[self._mapping] = self._ag._secondary[self._indices]
        return secondary
    
    def getSegmentNames(self):
        """Return segment names of the atoms."""
        if self._ag._segnames is None:
            return None
        segnames = np.zeros((self._len, 3), DTYPES['segnames'])
        segnames[self._mapping] = self._ag._segnames[self._indices]
        return segnames
    
    def getAnisotropicStandardDeviations(self):
        """Return standard deviations for the anisotropic temperature
        factors of the atoms."""
        if self._ag._siguij is None:
            return None
        siguij = np.zeros((self._len, 3), DTYPES['siguij'])
        siguij[self._mapping] = self._ag._siguij[self._indices]
        return siguij
    
    def getTemperatureFactors(self):
        """Return temperature (B) factors of the atoms."""
        if self._ag._bfactors is None:
            return None
        bfactors = np.zeros(self._len, DTYPES['bfactors'])
        bfactors[self._mapping] = self._ag._bfactors[self._indices]
        return bfactors
    
    def getRadii(self):
        """Return atomic radii of the atoms."""
        if self._ag._radii is None:
            return None
        radii = np.zeros(self._len, DTYPES['radii'])
        radii[self._mapping] = self._ag._radii[self._indices]
        return radii
    
    def getMasses(self):
        """Return atomic masses of the atoms."""
        if self._ag._masses is None:
            return None
        masses = np.zeros(self._len, DTYPES['masses'])
        masses[self._mapping] = self._ag._masses[self._indices]
        return masses
    
    def getCharges(self):
        """Return atomic partial charges of the atoms."""
        if self._ag._charges is None:
            return None
        charges = np.zeros(self._len, DTYPES['charges'])
        charges[self._mapping] = self._ag._charges[self._indices]
        return charges
    
    def getAtomTypes(self):
        """Return atom types of the atoms."""
        if self._ag._atomtypes is None:
            return None
        atomtypes = np.zeros(self._len, DTYPES['atomtypes'])
        atomtypes[self._mapping] = self._ag._atomtypes[self._indices]
        return atomtypes

    def getInsertionCodes(self):
        """Return insertion codes of the atoms."""
        if self._ag._icodes is None:
            return None
        icodes = np.zeros(self._len, DTYPES['icodes'])
        icodes[self._mapping] = self._ag._icodes[self._indices]
        return icodes
    
    def getUnmappedFlags(self):
        """Return an array with 1s for unmapped atoms."""
        flags = np.zeros(self._len)
        if len(self._unmapped):
            flags[self._unmapped] = 1
        return flags
    
    def getMappedFlags(self):
        """Return an array with 1s for mapped atoms."""
        flags = np.ones(self._len)
        if len(self._unmapped):
            flags[self._unmapped] = 0
        return flags

    def select(self, selstr):
        """Return a atom map matching the criteria given by *selstr*.
        
        Note that this is a special case for making atom selections. Unmapped
        atoms will not be included in the returned :class:`AtomMap` instance.
        The order of atoms will be preserved.
        
        """
        return SELECT.select(self, selstr)
