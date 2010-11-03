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

""":mod:`atom` module defines a class for accessing individual atom data.

Classes:

  * :class:`Atom`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import prody

__all__ = ['Atom']

class Atom(object):
    """A class for accessing and manipulating attributes of an atom 
    in a :class:`AtomGroup` instance."""
    
    __slots__ = ('_ag', '_index', '_acsi')
    
    def __init__(self, atomgroup, index, acsi=None):
        if not isinstance(atomgroup, prody.AtomGroup):
            raise TypeError('atomgroup must be AtomGroup, not {0:s}'
                            .format(type(atomgroup)))
        self._ag = atomgroup
        self._index = int(index)
        if acsi is None:
            self._acsi = atomgroup.getActiveCoordsetIndex()
        else: 
            self._acsi = int(acsi)
        
    def __repr__(self):
        return ('<Atom: {0:s} from {1:s} (index {2:d}; {3:d} '
                'coordinate sets, active set index: {4:d})>').format(
                self.getAtomName(), self._ag._name, self._index,  
                self._ag._n_coordsets, self._acsi)

    def __str__(self):
        return ('Atom {0:s} from {1:s} (index {2:d})').format(
                self.getAtomName(), self._ag.getName(), self._index)
        return ('{0:s} from {2:s} (index {1:d}; {3:d} '
                'coordinate sets, active set index: {4:d})').format(
                self.getAtomName(), self._index, self._ag._name, 
                self._ag._n_coordsets, self._acsi)

    def __len__(self):
        return 1
    
    def getAtomGroup(self):
        """Return associated atom group."""
        return self._ag
    
    def getIndex(self):
        """Return index of the atom."""
        return self._index
    
    def getNumOfCoordsets(self):
        """Return number of coordinate sets."""
        return self._ag._n_coordsets
    
    def getActiveCoordsetIndex(self):
        """Return the index of the active coordinate set for the atom."""
        return self._acsi
    
    def setActiveCoordsetIndex(self, index):
        """Set the index of the active coordinate set for the atom."""
        if self._ag._coordinates is None:
            raise AttributeError('coordinates are not set')
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if self._ag._n_coordsets <= index or \
           self._ag._n_coordsets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += self._ag._n_coordsets
        self._acsi = index
    
    def getCoordinates(self):
        """Return a copy of coordinates of the atom from the active coordinate set."""
        return self._ag._coordinates[self._acsi, self._index].copy()
    
    def setCoordinates(self, coordinates):
        """Set coordinates of the atom in the active coordinate set."""
        self._ag._coordinates[self._acsi, self._index] = coordinates
        
    def getCoordsets(self, indices):
        """Return a copy of coordinate sets at given indices.
        
        *indices* may be an integer or a list of integers.
        
        """
        if self._ag._coordinates is None:
            raise AttributeError('coordinates are not set')
        try: 
            return self._ag._coordinates[indices, self._index].copy()
        except IndexError:
            raise IndexError('indices may be an integer or a list of integers')

    def iterCoordsets(self):
        """Iterate over coordinate sets."""
        for i in range(self._ag._n_coordsets):
            yield self._ag._coordinates[i, self._index].copy()

    def getAtomName(self):
        """Return name of the atom."""
        if self._ag._atomnames is None:
            return None
        return self._ag._atomnames[self._index]
    
    getName = getAtomName
    
    def setAtomName(self, atom_name):
        """Set name of the atom."""
        if self._ag._atomnames is None:
            raise AttributeError('atom names are not set')
        self._ag._atomnames[self._index] = name
        
    setName = setAtomName
    
    def getAlternateLocationIndicator(self):
        """Return alternate location indicator of the atom."""
        if self._ag._altlocs is None:
            return None
        return self._ag._altlocs[self._index]
    
    def setAlternateLocationIndicator(self, alternate_location_indicator):
        """Set alternate location indicator of the atom."""
        if self._ag._altlocs is None:
            raise AttributeError('alternate location indicators are not set')
        self._ag._altlocs[self._index] = alternate_location_indicator
    
    def getAnisotropicTemperatureFactors(self):
        """Return anisotropic temperature factors of the atom."""
        if self._ag._anisou is None:
            return None
        return self._ag._anisou[self._index]
    
    def setAnisotropicTemperatureFactors(self, anisotropic_temperature_factors):
        """Set anisotropic temperature factors of the atom."""
        if self._ag._anisou is None:
            raise AttributeError('anisotropic temperature factors are not set')
        self._ag._anisou[self._index] = anisotropic_temperature_factors
        
    def getChainIdentifier(self):
        """Return chain identifier of the atom."""
        if self._ag._chainids is None:
            return None
        return self._ag._chainids[self._index]
    
    def setChainIdentifier(self, chain_identifier):
        """Set chain identifier of the atom."""
        if self._ag._chainids is None:
            raise AttributeError('chain identifiers are not set')
        self._ag._chainids[self._index] = chain_identifier
    
    def getElementSymbol(self):
        """Return element symbol of the atom."""
        if self._ag._elements is None:
            return None
        return self._ag._elements[self._index]
    
    def setElementSymbol(self, element_symbol):
        """Set element symbol of the atom."""
        if self._ag._elements is None:
            raise AttributeError('element symbols are not set')
        self._ag._elements[self._index] = element_symbol
    
    def getHeteroFlag(self):
        """Return hetero flag of the atom."""
        if self._ag._hetero is None:
            return None
        return self._ag._hetero[self._index]
    
    def setHeteroFlag(self, hetero_flag):
        """Set hetero flag of the atom."""
        if self._ag._hetero is None:
            raise AttributeError('hetero flags are not set')
        self._ag._hetero[self._index] = hetero_flag
    
    def getOccupancy(self):
        """Return occupancy of the atom."""
        if self._ag._occupancies is None:
            return None
        return self._ag._occupancies[self._index]
    
    def setOccupancy(self, occupancy):
        """Set occupancy of the atom."""
        if self._ag._occupancies is None:
            raise AttributeError('occupancies are not set')
        self._ag._occupancies[self._index] = occupancy
    
    def getResidueName(self):
        """Return residue name of the atom."""
        if self._ag._resnames is None:
            return None
        return self._ag._resnames[self._index]
    
    def setResidueName(self, residue_name):
        """Set residue name of the atom."""
        if self._ag._resnames is None:
            raise AttributeError('residue names are not set')
        self._ag._resnames[self._index] = residue_name
    
    def getResidueNumber(self):
        """Return residue number of the atom."""
        if self._ag._resnums is None:
            return None
        return self._ag._resnums[self._index]
    
    def setResidueNumber(self, residue_number):
        """Set residue number of the atom."""
        if self._ag._resnums is None:
            raise AttributeError('residue numbers are not set')
        self._ag._resnums[self._index] = residue_number
    
    def getSecondaryStructureAssignment(self):
        """Return secondary structure assignment of the atom."""
        if self._ag._secondary is None:
            return None
        return self._ag._secondary[self._index]
    
    def setSecondaryStructureAssignment(self, secondary_structure_assignment):
        """Set secondary structure assignment of the atom."""
        if self._ag._secondary is None:
            raise AttributeError('secondary structure assignments are not set')
        self._ag._secondary[self._index] = secondary_structure_assignment
    
    def getSegmentName(self):
        """Return segment name of the atom."""
        if self._ag._segnames is None:
            return None
        return self._ag._segnames[self._index]
    
    def setSegmentName(self, segment_name):
        """Set segment name of the atom."""
        if self._ag._segnames is None:
            raise AttributeError('segment names are not set')
        self._ag._segnames[self._index] = segment_name
    
    def getAnisotropicStandardDeviations(self):
        """Return standard deviations for the anisotropic temperature
        factors of the atom."""
        if self._ag._siguij is None:
            return None
        return self._ag._siguij[self._index]
    
    def setAnisotropicStandardDeviations(self, anisotropic_standard_deviations):
        """Set standard deviations for the anisotropic temperature factors of 
        the atom."""
        if self._ag._siguij is None:
            raise AttributeError('standard deviations for the anisotropic temperature factors are not set')
        self._ag._siguij[self._index] = anisotropic_standard_deviations
    
    def getTemperatureFactor(self):
        """Return temperature (B) factor of the atom."""
        if self._ag._bfactors is None:
            return None
        return self._ag._bfactors[self._index]
    
    def setTemperatureFactor(self, temperature_factor):
        """Set temperature (B) factor of the atom."""
        if self._ag._bfactors is None:
            raise AttributeError('temperature (B) factors are not set')
        self._ag._bfactors[self._index] = temperature_factor
    
    def getRadius(self):
        """Return radius of the atom."""
        if self._ag._radii is None:
            return None
        return self._ag._radii[self._index]
    
    def setRadius(self, radius):
        """Set atomic radii of the atom."""
        if self._ag._radii is None:
            raise AttributeError('atomic radii are not set')
        self._ag._radii[self._index] = radius
    
    def getMass(self):
        """Return mass of the atom."""
        if self._ag._masses is None:
            return None
        return self._ag._masses[self._index]
    
    def setMass(self, mass):
        """Set mass of the atom."""
        if self._ag._masses is None:
            raise AttributeError('atomic masses are not set')
        self._ag._masses[self._index] = mass
    
    def getCharge(self):
        """Return partial charges of the atom."""
        if self._ag._charges is None:
            return None
        return self._ag._charges[self._index]
    
    def setCharge(self, charge):
        """Set partial charges of the atom."""
        if self._ag._atomtypes is None:
            raise AttributeError('atomic partial charges are not set')
        self._ag._charges[self._index] = charge
    
    def getAtomType(self):
        """Return type of the atom."""
        if self._ag._atomtypes is None:
            return None
        return self._ag._atomtypes[self._index]
    
    def setAtomType(self, atom_type):
        """Set type of the atom."""
        if self._ag._atomtypes is None:
            raise AttributeError('atom types are not set')
        self._ag._atomtypes[self._index] = atom_type

    def getInsertionCode(self):
        """Return insertion code of the atom."""
        if self._ag._icodes is None:
            return None
        return self._ag._icodes[self._index]
    
    def setInsertionCode(self, insertion_code):
        """Set insertion code of the atom."""
        if self._ag._icodes is None:
            raise AttributeError('atom types are not set')
        self._ag._icodes[self._index] = insertion_code


    def getSelectionString(self):
        """Return selection string that will select this atom."""
        return 'index {0:d}'.format(self._index)

