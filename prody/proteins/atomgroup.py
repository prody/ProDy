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

""":mod:`atomgroup` module defines a class for mainly storing atomic data.

Classes:

    * :class:`AtomGroup`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np

import prody
from prody import ProDyLogger as LOGGER
from . import Atom, Selection
from .select import ProDyAtomSelect as SELECT
from . import DTYPES
from .hierview import HierView

__all__ = ['AtomGroup']


class AtomGroup(object):
    
    """A class for storing and accessing atomic data.
    
    
    The number of atoms of the atom group is inferred at the first set method
    call from the size of the data array. 

    **Atomic Data**
    
    All atomic data is stored in :class:`numpy.ndarray` instances.

    **Get and Set Methods** 
    
    *get_attribute()* methods return copies of the data arrays. 
    
    *set_attribute()* methods accepts data contained in :class:`list` or 
    :class:`numpy.ndarray` instances. The length of the list or array must 
    match the number of atoms in the atom group. Set method sets attributes of 
    all atoms at once.
    
    Atom groups with multiple coordinate sets may have one of these sets as 
    the active coordinate set. The active coordinate set may be changed using
    :meth:`setActiveCoordsetIndex()` method. :meth:`getCoordinates` returns
    coordinates from the active set.
    
    To access and modify data associated with a subset of atoms in an atom group, 
    :class:`Selection` instances may be used. A selection from an atom 
    group has initially the same coordinate set as the active coordinate set.
    
    User can iterate over atoms and coordinate sets in an atom group. To 
    iterate over residues and chains, get a hierarchical view of the atom 
    group by calling :meth:`getHierView()`. 

    """
    
    def __init__(self, name):
        """Instantiate an AtomGroup with a *name*."""
        self._name = str(name)
        self._n_atoms = 0
        self._coordinates = None
        self._acsi = 0                  # Active Coordinate Set Index
        self._n_coordsets = 0
        
        self._atomnames = None          # Attributes parsed from PDB files 
        self._altlocs = None
        self._anisou = None
        self._chainids = None
        self._elements = None
        self._occupancies = None
        self._resnums = None
        self._resnames = None
        self._hetero = None
        self._secondary = None
        self._segnames = None
        self._siguij = None
        self._bfactors = None
        self._icodes = None

        self._charges = None            # Some useful atomic attributes
        self._masses = None
        self._radii = None
        self._atomtypes = None          # Typically used in force fields
        
    def __repr__(self):
        return ('<AtomGroup: {0:s} ({1:d} atoms; {2:d} coordinate sets, active '
               'set index: {3:d})>').format(self._name, 
              self._n_atoms, self._n_coordsets, self._acsi)
        return ('<AtomGroup: {0:s}>').format(str(self))
        
    def __str__(self):
        return ('AtomGroup {0:s}').format(self._name)
        return ('{0:s} ({1:d} atoms; {2:d} coordinate sets, active '
               'set index: {3:d})').format(self._name, 
              self._n_atoms, self._n_coordsets, self._acsi)
    
    def __getitem__(self, indices):
        if isinstance(indices, int):
            if indices < 0:
                indices = self._n_atoms + indices
            return Atom(self, indices, self._acsi)
        elif isinstance(indices, slice):
            start, stop, step = indices.indices(self._n_atoms)
            if start is None:
                start = 0
            if step is None:
                step = 1
            selstr = 'index {0:d}:{1:d}:{2:d}'.format(start, stop, step)
            return Selection(self, 
                             np.arange(start, stop, step), 
                             selstr,
                             self._acsi)
        elif isinstance(indices, (list, np.ndarray)):
            return Selection(self, 
                             np.array(indices), 
                             'Some atoms', 
                             'index {0:s}'.format(
                                            ' '.join(np.array(indices, '|S'))),
                             self._acsi)
        else:
            raise IndexError('invalid index') 
    
    def __iter__(self):
        """Iterate over atoms in the atom group."""
        acsi = self._acsi
        for index in xrange(self._n_atoms):
            yield Atom(self, index, acsi)

    def __len__(self):
        return self._n_atoms
    
    def __add__(self, other):
        if not isinstance(other, AtomGroup):
            raise TypeError('type mismatch')
        if self == other:
            raise ValueError('an atom group cannot be added to itself')
        
        new = AtomGroup(self._name + ' + ' + other._name)
        n_coordsets = self._n_coordsets
        if n_coordsets != other._n_coordsets:
            LOGGER.warning('AtomGroups {0:s} and {1:s} do not have same number '
                'of coordinate sets. First from both AtomGroups will be merged.'
                .format(str(self._name), str(other._name), n_coordsets))
            n_coordsets = 1
        coordset_range = range(n_coordsets)
        new.setCoordinates(np.concatenate((self._coordinates[coordset_range],
                                        other._coordinates[coordset_range]), 1))
        
        if self._atomnames is not None and other._atomnames is not None:
            new._atomnames = np.concatenate((self._atomnames, other._atomnames))

        if self._altlocs is not None and other._altlocs is not None:
            new._altlocs = np.concatenate((self._altlocs, other._altlocs))

        if self._anisou is not None and other._anisou is not None:
            new._anisou = np.concatenate((self._anisou, other._anisou))

        if self._resnames is not None and other._resnames is not None:
            new._resnames = np.concatenate((self._resnames, other._resnames))
            
        if self._resnums is not None and other._resnums is not None:
            new._resnums = np.concatenate((self._resnums, other._resnums))
            
        if self._chainids is not None and other._chainids is not None:
            new._chainids = np.concatenate((self._chainids, other._chainids))
            
        if self._bfactors is not None and other._bfactors is not None:
            new._bfactors = np.concatenate((self._bfactors, other._bfactors))
            
        if self._occupancies is not None and other._occupancies is not None:
            new._occupancies = np.concatenate((self._occupancies, other._occupancies))
            
        if self._hetero is not None and other._hetero is not None:
            new._hetero = np.concatenate((self._hetero, other._hetero))
            
        if self._elements is not None and other._elements is not None:
            new._elements = np.concatenate((self._elements, other._elements))
            
        if self._segnames is not None and other._segnames is not None:
            new._segnames = np.concatenate((self._segnames, other._segnames))
            
        if self._secondary is not None and other._secondary is not None:
            new._secondary = np.concatenate((self._secondary, other._secondary))
            
        if self._siguij is not None and other._siguij is not None:
            new._siguij = np.concatenate((self._siguij, other._siguij))
            
        if self._charges is not None and other._charges is not None:
            new._charges = np.concatenate((self._charges, other._charges))

        if self._radii is not None and other._radii is not None:
            new._radii = np.concatenate((self._radii, other._radii))

        if self._masses is not None and other._masses is not None:
            new._masses = np.concatenate((self._masses, other._masses))

        if self._atomtypes is not None and other._atomtypes is not None:
            new._atomtypes = np.concatenate((self._atomtypes, other._atomtypes))

        return new

    def getName(self):
        """Return name of the atom group instance."""
        return self._name
    
    def setName(self, name):
        """Set name of the atom group instance."""
        self._name = str(name)
    
    def getNumOfAtoms(self):
        """Return number of atoms."""
        return self._n_atoms
    
    def getCoordinates(self): 
        """Return coordinates from active coordinate set."""
        if self._coordinates is None:
            return None
        return self._coordinates[self._acsi].copy()
    
    def setCoordinates(self, coordinates):
        """Set coordinates.
        
        Coordinates must be a NumPy ndarray instance.
        
        If the shape of the coordinates array is (n_coordsets,n_atoms,3),
        the given array will replace all coordinate sets. To avoid it,
        :meth:`addCoordset` may be used.
        
        If the shape of the coordinates array is (n_atoms,3) or (1,n_atoms,3),
        the coordinate set will replace the coordinates of the currently active 
        coordinate set.
        
        """
        if not isinstance(coordinates, np.ndarray):
            raise TypeError('coordinates must be an ndarray instance')

        if not coordinates.ndim in (2, 3):
            raise ValueError('coordinates must be a 2d or a 3d array')
            
        if coordinates.shape[-1] != 3:
            raise ValueError('shape of coordinates must be (n_atoms,3) or '
                             '(n_coordsets,n_atoms,3)')
        
        if coordinates.dtype != DTYPES['coordinates']:
            try:
                coordinates.astype(DTYPES['coordinates'])
            except ValueError:
                raise ValueError('coordinate array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['coordinates']))

        if self._n_atoms == 0:
            self._n_atoms = coordinates.shape[-2] 
        elif coordinates.shape[-2] != self._n_atoms:
            raise ValueError('length of coordinate array must match n_atoms')
        
        if self._coordinates is None:
            if coordinates.ndim == 2:
                coordinates = coordinates.reshape(
                                (1, coordinates.shape[0], coordinates.shape[1]))
            self._coordinates = coordinates
            self._n_coordsets = self._coordinates.shape[0]
            self._acsi = 0
        else:
            if coordinates.ndim == 2:
                self._coordinates[self._acsi] = coordinates
            elif coordinates.shape[0] == 1:
                self._coordinates[self._acsi] = coordinates[0]
            else:
                self._coordinates = coordinates
                self._n_coordsets = self._coordinates.shape[0]
                self._acsi = min(self._n_coordsets-1,
                                                    self._acsi)
    
    def addCoordset(self, coords):
        """Add a coordinate set to the atom group."""
        if self._coordinates is None:
            self.setCoordinates(coords)
        if not isinstance(coords, np.ndarray):
            raise TypeError('coords must be an ndarray instance')
        elif not coords.ndim in (2, 3):
            raise ValueError('coords must be a 2d or a 3d array')
        elif coords.shape[-2:] != self._coordinates.shape[1:]:
            raise ValueError('shape of coords must be ([n_coordsets,] n_atoms, 3)')
        elif coords.dtype != DTYPES['coordinates']:
            try:
                coords.astype(DTYPES['coordinates'])
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['coordinates']))
        if coords.ndim == 2:
            coords = coords.reshape((1, coords.shape[0], coords.shape[1]))
        
        self._coordinates = np.concatenate((self._coordinates, coords), axis=0)
        self._n_coordsets = self._coordinates.shape[0]

    def delCoordset(self, index):
        """Delete a coordinate set from the atom group."""
        which = np.ones(self._n_coordsets, np.bool)
        which[index] = False
        if which.sum() == self._n_coordsets:
            self._coordinates = None
            self._n_coordsets = 0
        else:
            self._coordinates = self._coordinates[which]
            self._n_coordsets = self._coordinates.shape[0]
        self._acsi = 0

    def getCoordsets(self, indices=None):
        """Return a copy of coordinate sets at given indices.
        
        *indices* may be an integer, a list of integers or ``None``. ``None``
        returns all coordinate sets. 
        
        """
        if indices is None:
            indices = slice(None)
        if self._coordinates is None:
            return None
        try: 
            return self._coordinates[indices].copy()
        except IndexError:
            raise IndexError('indices must be an integer, a list of integers, or None')

    def getNumOfCoordsets(self):
        """Return number of coordinate sets."""
        return self._n_coordsets
    
    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each 
        coordinate set."""
        for i in range(self._n_coordsets):
            yield self._coordinates[i].copy()
    
    def getActiveCoordsetIndex(self):
        """Return index of the active coordinate set."""
        return self._acsi
    
    def setActiveCoordsetIndex(self, index):
        """Set the index of the active coordinate set."""
        if self._n_coordsets == 0:
            self._acsi = 0
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if self._n_coordsets <= index or self._n_coordsets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += self._n_coordsets 
        #self._kdtree = None
        self._acsi = index
    
    def getAtomNames(self):
        """Return a copy of atom names."""
        if self._atomnames is None:
            return None
        return self._atomnames.copy()

    def setAtomNames(self, atom_names):
        """Set atom names."""
        if self._n_atoms == 0:
            self._n_atoms = len(atom_names)
        elif len(atom_names) != self._n_atoms:
            raise ValueError('length of atom_names must match n_atoms')

        if isinstance(atom_names, list):
            atom_names = np.array(atom_names, DTYPES['atomnames'])
        elif not isinstance(atom_names, np.ndarray):
            raise TypeError('atom_names must be an ndarray or a list')
        elif atom_names.dtype != DTYPES['atomnames']:
            try:
                atom_names.astype(DTYPES['atomnames'])
            except ValueError:
                raise ValueError('atom_names array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['atomnames']))
        self._atomnames = atom_names
        
    def getAlternateLocationIndicators(self):
        """Return a copy of alternate location indicators."""
        if self._altlocs is None:
            return None
        return self._altlocs.copy()

    def setAlternateLocationIndicators(self, alternate_location_indicators):
        """Set alternate location identifiers."""
        if self._n_atoms == 0:
            self._n_atoms = len(alternate_location_indicators)
        elif len(alternate_location_indicators) != self._n_atoms:
            raise ValueError('length of alternate_location_indicators must match n_atoms')
            
        if isinstance(alternate_location_indicators, list):
            alternate_location_indicators = np.array(alternate_location_indicators, DTYPES['altlocs'])
        elif not isinstance(alternate_location_indicators, np.ndarray):
            raise TypeError('alternate_location_indicators must be an array or a list')
        elif alternate_location_indicators.dtype != DTYPES['altlocs']:
            try:
                alternate_location_indicators.astype(DTYPES['altlocs'])
            except ValueError:
                raise ValueError('alternate_location_indicators array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['altlocs']))
        self._altlocs = alternate_location_indicators
        
    def getAnisotropicTemperatureFactors(self):
        """Return a copy of anisotropic temperature factors."""
        if self._anisou is None:
            return None
        return self._anisou.copy()
    
    def setAnisotropicTemperatureFactors(self, anisotropic_temperature_factors):
        """Set anisotropic temperature factors."""
        if self._n_atoms == 0:
            self._n_atoms = len(anisotropic_temperature_factors)
        elif len(anisotropic_temperature_factors) != self._n_atoms:
            raise ValueError('length of anisotropic_temperature_factors must match n_atoms')

        if isinstance(anisotropic_temperature_factors, list):
            anisotropic_temperature_factors = np.array(anisotropic_temperature_factors, DTYPES['anisou'])
        elif not isinstance(anisotropic_temperature_factors, np.ndarray):
            raise TypeError('anisotropic_temperature_factors must be an array or a list')
        elif anisotropic_temperature_factors.dtype != DTYPES['anisou']:
            try:
                anisotropic_temperature_factors.astype(DTYPES['anisou'])
            except ValueError:
                raise ValueError('anisotropic_temperature_factors array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['anisou']))
        self._anisou = anisotropic_temperature_factors
        
    def getChainIdentifiers(self):
        """Return a copy of chain identifiers."""
        if self._chainids is None:
            return None
        return self._chainids.copy()
    
    def setChainIdentifiers(self, chain_identifiers):
        """Set chain identifiers."""
        if self._n_atoms == 0:
            self._n_atoms = len(chain_identifiers)
        elif len(chain_identifiers) != self._n_atoms:
            raise ValueError('length of chain_identifiers must match n_atoms')
                
        if isinstance(chain_identifiers, list):
            chain_identifiers = np.array(chain_identifiers)
        elif not isinstance(chain_identifiers, np.ndarray):
            raise TypeError('chain_identifiers must be an array or a list')
        elif chain_identifiers.dtype != DTYPES['chainids']:
            try:
                chain_identifiers.astype(DTYPES['chainids'])
            except ValueError:
                raise ValueError('chain_identifiers array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['chainids']))
        self._chainids = chain_identifiers
        
    def getElementSymbols(self):
        """Return a copy of element symbols."""
        if self._elements is None:
            return None
        return self._elements.copy()

    def setElementSymbols(self, element_symbols):
        """Set element symbols."""
        if self._n_atoms == 0:
            self._n_atoms = len(element_symbols)
        elif len(element_symbols) != self._n_atoms:
            raise ValueError('length of element_symbols must match n_atoms')

        if isinstance(element_symbols, list):
            element_symbols = np.array(element_symbols, DTYPES['elements'])
        elif not isinstance(element_symbols, np.ndarray):
            raise TypeError('element_symbols must be an array or a list')
        elif element_symbols.dtype != DTYPES['elements']:
            try:
                element_symbols.astype(DTYPES['elements'])
            except ValueError:
                raise ValueError('element_symbols array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['elements']))
        self._elements = element_symbols
        
    def getHeteroFlags(self):
        """Return a copy of hetero flags."""
        if self._hetero is None:
            return None
        return self._hetero.copy()

    def setHeteroFlags(self, hetero_flags):
        """Set hetero flags."""
        if self._n_atoms == 0:
            self._n_atoms = len(hetero_flags)
        elif len(hetero_flags) != self._n_atoms:
            raise ValueError('length of hetero_flags must match n_atoms')
        
        if isinstance(hetero_flags, list):
            hetero_flags = np.array(hetero_flags, DTYPES['hetero'])
        elif not isinstance(hetero_flags, np.ndarray):
            raise TypeError('hetero_flags must be an array or a list')
        elif hetero_flags.dtype != DTYPES['hetero']:
            try:
                hetero_flags.astype(DTYPES['hetero'])
            except ValueError:
                raise ValueError('hetero_flags array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['hetero']))
        self._hetero = hetero_flags

    def getOccupancies(self):
        """Return a copy of occupancies."""
        if self._occupancies is None:
            return None
        return self._occupancies.copy()

    def setOccupancies(self, occupancies):
        """Set occupancies."""
        if self._n_atoms == 0:
            self._n_atoms = len(occupancies)
        elif len(occupancies) != self._n_atoms:
            raise ValueError('length of occupancies must match n_atoms')

        if isinstance(occupancies, list):
            occupancies = np.array(occupancies, DTYPES['occupancies'])
        elif not isinstance(occupancies, np.ndarray):
            raise TypeError('occupancies must be an array or a list')
        elif occupancies.dtype != DTYPES['occupancies']:
            try:
                occupancies.astype(DTYPES['occupancies'])
            except ValueError:
                raise ValueError('occupancies array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['occupancies']))
        self._occupancies = occupancies
        
    def getResidueNames(self):
        """Return a copy of residue names."""
        if self._resnames is None:
            return None
        return self._resnames.copy()
    
    def setResidueNames(self, residue_names):
        """Set residue names."""
        if len(residue_names) != self._n_atoms:
            raise ValueError('length of residue_names must match n_atoms')
        elif isinstance(residue_names, list):
            residue_names = np.array(residue_names, DTYPES['resnames'])
        elif not isinstance(residue_names, np.ndarray):
            raise TypeError('residue_names must be an array or a list')
        elif residue_names.dtype != DTYPES['resnames']:
            try:
                residue_names.astype(DTYPES['resnames'])
            except ValueError:
                raise ValueError('residue_names array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['resnames']))
        self._resnames = residue_names
        
    def getResidueNumbers(self):
        """Return a copy of residue numbers."""
        if self._resnums is None:
            return None
        return self._resnums.copy()

    def setResidueNumbers(self, residue_numbers):
        """Set residue numbers."""
        if self._n_atoms == 0:
            self._n_atoms = len(residue_numbers)
        elif len(residue_numbers) != self._n_atoms:
            raise ValueError('length of residue_numbers must match n_atoms')
        
        if isinstance(residue_numbers, list):
            residue_numbers = np.array(residue_numbers, DTYPES['resnums'])
        elif not isinstance(residue_numbers, np.ndarray):
            raise TypeError('residue_numbers must be an array or a list')
        elif residue_numbers.dtype != DTYPES['resnums']:
            try:
                residue_numbers.astype(DTYPES['resnums'])
            except ValueError:
                raise ValueError('residue_numbers array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['resnums']))
        self._resnums = residue_numbers
        
    def getSecondaryStructureAssignments(self):
        """Return a copy of secondary structure assignments."""
        if self._secondary is None:
            return None
        return self._secondary.copy()
    
    def setSecondaryStructureAssignments(self, secondary_structure_assignments):
        """Set secondary structure assignments."""
        if self._n_atoms == 0:
            self._n_atoms = len(secondary_structure_assignments)
        elif len(secondary_structure_assignments) != self._n_atoms:
            raise ValueError('length of secondary_structure_assignments must match n_atoms')
        
        if isinstance(secondary_structure_assignments, list):
            secondary_structure_assignments = np.array(secondary_structure_assignments, DTYPES['secondary'])
        elif not isinstance(secondary_structure_assignments, np.ndarray):
            raise TypeError('secondary_structure_assignments must be an array or a list')
        elif secondary_structure_assignments.dtype != DTYPES['secondary']:
            try:
                secondary_structure_assignments.astype(DTYPES['secondary'])
            except ValueError:
                raise ValueError('secondary_structure_assignments array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['secondary']))
        self._secondary = secondary_structure_assignments
        
    def getSegmentNames(self):
        """Return a copy of segment names."""
        if self._segnames is None:
            return None
        return self._segnames.copy()
    
    def setSegmentNames(self, segment_names):
        """Set segment names."""
        if self._n_atoms == 0:
            self._n_atoms = len(segment_names)
        elif len(segment_names) != self._n_atoms:
            raise ValueError('length of segment_names must match n_atoms')
        
        if isinstance(segment_names, list):
            segment_names = np.array(segment_names, DTYPES['segnames'])
        elif not isinstance(segment_names, np.ndarray):
            raise TypeError('segment_names must be an array or a list')
        elif segment_names.dtype != DTYPES['segnames']:
            try:
                segment_names.astype(DTYPES['segnames'])
            except ValueError:
                raise ValueError('segment_names array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['segnames']))
        self._segnames = segment_names
        
    def getAnisotropicStandardDeviations(self):
        """Return a copy of standard deviations for the anisotropic temperature factors."""
        if self._siguij is None:
            return None
        return self._siguij.copy()
    
    def setAnisotropicStandardDeviations(self, anisotropic_standard_deviations):
        """Set standard deviations for the anisotropic temperature factors."""
        if self._n_atoms == 0:
            self._n_atoms = len(anisotropic_standard_deviations)
        elif len(anisotropic_standard_deviations) != self._n_atoms:
            raise ValueError('length of siguij must match n_atoms')
        
        if isinstance(anisotropic_standard_deviations, list):
            anisotropic_standard_deviations = np.array(anisotropic_standard_deviations, DTYPES['siguij'])
        elif not isinstance(anisotropic_standard_deviations, np.ndarray):
            raise TypeError('anisotropic_standard_deviations must be an array or a list')
        elif anisotropic_standard_deviations.dtype != DTYPES['siguij']:
            try:
                anisotropic_standard_deviations.astype(DTYPES['siguij'])
            except ValueError:
                raise ValueError('anisotropic_standard_deviations array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['siguij']))
        self._siguij = anisotropic_standard_deviations
        
    def getTemperatureFactors(self):
        """Return a copy of temperature (B) factors."""
        if self._bfactors is None:
            return None
        return self._bfactors.copy()
                        
    def setTemperatureFactors(self, temperature_factors):
        """Set temperature (B) factors."""
        if self._n_atoms == 0:
            self._n_atoms = len(temperature_factors)
        elif len(temperature_factors) != self._n_atoms:
            raise ValueError('length of temperature_factors must match n_atoms')
        
        if isinstance(temperature_factors, list):
            temperature_factors = np.array(temperature_factors, DTYPES['bfactor'])
        elif not isinstance(temperature_factors, np.ndarray):
            raise TypeError('temperature_factors must be an array or a list')
        elif temperature_factors.dtype != DTYPES['bfactors']:
            try:
                temperature_factors.astype(DTYPES['bfactors'])
            except ValueError:
                raise ValueError('temperature_factors array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['bfactors']))
        self._bfactor = temperature_factors
        
    def getRadii(self):
        """Return a copy of atomic radii."""
        if self._radii is None:
            return None
        return self._radii.copy()
    
    def setRadii(self, radii):
        """Set atomic radii."""
        if self._n_atoms == 0:
            self._n_atoms = len(radii)
        elif len(radii) != self._n_atoms:
            raise ValueError('length of radii must match n_atoms')
        
        if isinstance(radii, list):
            radii = np.array(radii, DTYPES['radii'])
        elif not isinstance(radii, np.ndarray):
            raise TypeError('radii must be an array or a list')
        elif radii.dtype != DTYPES['radii']:
            try:
                radii.astype(DTYPES['radii'])
            except ValueError:
                raise ValueError('radii array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['radii']))
        self._radii = radii
        
    def getMasses(self):
        """Return a copy of atomic masses."""
        if self._masses is None:
            return None
        return self._masses.copy()
    
    def setMasses(self, masses):
        """Set atomic masses."""
        if self._n_atoms == 0:
            self._n_atoms = len(masses)
        elif len(masses) != self._n_atoms:
            raise ValueError('length of masses must match n_atoms')
        
        if isinstance(masses, list):
            masses = np.array(masses, DTYPES['masses'])
        elif not isinstance(masses, np.ndarray):
            raise TypeError('masses must be an array or a list')
        elif masses.dtype != DTYPES['masses']:
            try:
                masses.astype(DTYPES['masses'])
            except ValueError:
                raise ValueError('masses array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['masses']))
        self._masses = masses
        
    def getCharges(self):
        """Return a copy of atomic partial charges."""
        if self._charges is None:
            return None
        return self._charges.copy()
    
    def setCharges(self, charges):
        """Set atomic partial charges."""
        if self._n_atoms == 0:
            self._n_atoms = len(charges)
        elif len(charges) != self._n_atoms:
            raise ValueError('length of charges must match n_atoms')
        
        if isinstance(charges, list):
            charges = np.array(charges, DTYPES['charges'])
        elif not isinstance(charges, np.ndarray):
            raise TypeError('charges must be an array or a list')
        elif charges.dtype != DTYPES['charges']:
            try:
                charges.astype(DTYPES['charges'])
            except ValueError:
                raise ValueError('charges array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['charges']))
        self._charges = charges
        
    def getAtomTypes(self):
        """Return a copy of atom types."""
        if self._atomtypes is None:
            return None
        return self._atomtypes.copy()

    def setAtomTypes(self, atom_types):
        """Set atom types."""
        if self._n_atoms == 0:
            self._n_atoms = len(atom_types)
        elif len(atom_types) != self._n_atoms:
            raise ValueError('length of atom_types must match n_atoms')
        
        if isinstance(atom_types, list):
            atom_types = np.array(atom_types, DTYPES['atomtypes'])
        elif not isinstance(atom_types, np.ndarray):
            raise TypeError('atom_types must be an ndarray or a list')
        elif atom_types.dtype != DTYPES['atomtypes']:
            try:
                atom_types.astype(DTYPES['atomtypes'])
            except ValueError:
                raise ValueError('atom_types array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['atomtypes']))
        self._atomtypes = atom_types
        
    def getInsertionCodes(self):
        """Return a copy of insertion codes."""
        if self._icodes is None:
            return None
        return self._icodes.copy()

    def setInsertionCodes(self, insertion_codes):
        """Set atom types."""
        if self._n_atoms == 0:
            self._n_atoms = len(insertion_codes)
        elif len(insertion_codes) != self._n_atoms:
            raise ValueError('length of insertion_codes must match n_atoms')
        
        if isinstance(insertion_codes, list):
            insertion_codes = np.array(insertion_codes, DTYPES['icodes'])
        elif not isinstance(insertion_codes, np.ndarray):
            raise TypeError('insertion_codes must be an ndarray or a list')
        elif insertion_codes.dtype != DTYPES['icodes']:
            try:
                insertion_codes.astype(DTYPES['icodes'])
            except ValueError:
                raise ValueError('insertion_codes array cannot be assigned type '
                                 '{0:s}'.format(DTYPES['icodes']))
        self._icodes = insertion_codes


    def select(self, selstr):
        """Return a selection matching the criteria given by *selstr*."""
        return SELECT.select(self, selstr)
    
    
    def copy(self, which=None):
        """Return a copy of atoms indicated *which* as a new AtomGroup instance.
        
        *which* may be:
            
            * a Selection, Residue, Chain, or Atom instance
            * a list or an array of indices
            * a selection string
        
        """
        
        if which is None:
            indices = None
            newmol = AtomGroup('Copy of {0:s}'.format(self._name))
        elif isinstance(which, int):
            indices = [which]
            newmol = AtomGroup('Copy of {0:s} index {1:d}'.format(
                                 self._name, which))
        elif isinstance(which, str):
            indices = SELECT.select(self, which).getIndices()
            newmol = AtomGroup('Copy of {0:s} selection "{1:s}"'
                              .format(self._name, which))
        elif isinstance(which, (list, np.ndarray)):
            if isinstance(which, list):
                indices = np.array(which)
            else:
                indices = which
            newmol = AtomGroup('Copy of a {0:s} subset'
                              .format(self._name))
        elif isinstance(which, prody.Selection):
            indices = which.getIndices()
            newmol = AtomGroup('Copy of {0:s} selection "{1:s}"'
                              .format(self._name, which.getSelectionString()))
        elif isinstance(which, prody.Chain):
            indices = which.getIndices()
            newmol = AtomGroup('Copy of {0:s} chain {1:s}'
                              .format(self._name, which.getIdentifier()))
        elif isinstance(which, prody.Residue):
            indices = which.getIndices()
            newmol = AtomGroup('Copy of {0:s} residue {1:s}{2:d}'
                              .format(self._name, which.getName(), which.getNumber()))
        elif isinstance(which, prody.Atom):
            indices = [which.getIndex()]
            newmol = AtomGroup('Copy of {0:s} index {1:d}'.format(
                                 self._name, which.getIndex()))
        elif isinstance(which, prody.AtomMap):
            indices = which.getIndices()
            newmol = AtomGroup('Copy of {0:s} atom map {1:s}'.format(
                                 self._name, str(which)))
            

        if indices is None:
            newmol.setCoordinates(self._coordinates.copy())
            if self._atomnames is not None:
                newmol._atomnames = self._atomnames.copy()
            if self._resnames is not None:
                newmol._resnames = self._resnames.copy()
            if self._resnums is not None:
                newmol._resnums = self._resnums.copy()
            if self._chainids is not None:
                newmol._chainids = self._chainids.copy()
            if self._bfactors is not None:
                newmol._bfactors = self._bfactors.copy()
            if self._occupancies is not None:
                newmol._occupancies = self._occupancies.copy()
            if self._hetero is not None:
                newmol._hetero = self._hetero.copy()
            if self._altlocs is not None:
                newmol._altlocs = self._altlocs.copy()
            if self._elements is not None:
                newmol._elements = self._elements.copy()
            if self._segnames is not None:
                newmol._segnames = self._segnames.copy()
            if self._secondary is not None:
                newmol._secondary = self._secondary.copy()
            if self._anisou is not None:
                newmol._anisou = self._anisou.copy()
            if self._siguij is not None:
                newmol._siguij = self._siguij.copy()
            if self._charges is not None:
                newmol._charges = self._charges.copy()
            if self._masses is not None:
                newmol._masses = self._masses.copy()
            if self._radii is not None:
                newmol._radii = self._radii.copy()
            if self._atomtypes is not None:
                newmol._atomtypes = self._atomtypes.copy()
            if self._icodes is not None:
                newmol._icodes = self._icodes.copy()
        else:
            newmol.setCoordinates(
                    self._coordinates[:, indices].copy(
                        ).reshape((self._n_coordsets, len(indices), 3)))
            if self._atomnames is not None:
                newmol._atomnames = self._atomnames[indices].copy()
            if self._resnames is not None:
                newmol._resnames = self._resnames[indices].copy()
            if self._resnums is not None:
                newmol._resnums = self._resnums[indices].copy()
            if self._chainids is not None:
                newmol._chainids = self._chainids[indices].copy()
            if self._bfactors is not None:
                newmol._bfactors = self._bfactors[indices].copy()
            if self._occupancies is not None:
                newmol._occupancies = self._occupancies[indices].copy()
            if self._hetero is not None:
                newmol._hetero = self._hetero[indices].copy()
            if self._altlocs is not None:
                newmol._altlocs = self._altlocs[indices].copy()
            if self._elements is not None:
                newmol._elements = self._elements[indices].copy()
            if self._segnames is not None:
                newmol._segnames = self._segnames[indices].copy()
            if self._secondary is not None:
                newmol._secondary = self._secondary[indices].copy()
            if self._anisou is not None:
                newmol._anisou = self._anisou[indices].copy()
            if self._siguij is not None:
                newmol._siguij = self._siguij[indices].copy()
            if self._charges is not None:
                newmol._charges = self._charges[indices].copy()
            if self._masses is not None:
                newmol._masses = self._masses[indices].copy()
            if self._radii is not None:
                newmol._radii = self._radii[indices].copy()
            if self._atomtypes is not None:
                newmol._atomtypes = self._atomtypes[indices].copy()
            if self._icodes is not None:
                newmol._icodes = self._icodes[indices].copy()
        
        return newmol
    
    def getHierView(self):
        """Return a hierarchical view of the atom group."""
        return HierView(self)
