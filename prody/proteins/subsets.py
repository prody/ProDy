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
""":mod:`subsets` module defines classes for manipulating subsets of atomic 
data.

Classes:
    
    * :class:`AtomSubset`
    * :class:`Chain`
    * :class:`Residue`
    * :class:`Selection`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np

import prody
from prody import ProDyLogger as LOGGER
from . import Atom
from .select import ProDyAtomSelect as SELECT

__all__ = ['AtomSubset', 'Selection', 'Chain', 'Residue']

class AtomSubset(object):
    """A class for manipulating subset of atomic data in an :class:`AtomGroup`.
    
    This class stores a reference to an :class:`AtomGroup` instance, a set of 
    atom indices, and active coordinate set index for the atom group.
    
    """
    
    __slots__ = ['_ag', '_indices', '_acsi']
    
    def __init__(self, atomgroup, indices, acsi=None):
        """Instantiate atom group base class. 
        
        :arg atomgroup: an atom group
        :type atomgroup: :class:`AtomGroup`
        
        :arg indices: list of indices of atoms in the subset
        :type indices: list of integers
        
        :arg acsi: active coordinate set index
        :type acsi: integer
        
        """
        
        if not isinstance(atomgroup, prody.AtomGroup):
            raise TypeError('atomgroup must be AtomGroup, not {0:s}'
                            .format(type(atomgroup)))
        self._ag = atomgroup

        if not isinstance(indices, np.ndarray):
            indices = np.array(indices, np.int64)
        elif not indices.dtype == np.int64:
            indices = indices.astype(np.int64)
        else:
            indices = indices
        self._indices = np.unique(indices)

        if acsi is None:
            self._acsi = atomgroup.getActiveCoordsetIndex()
        else:
            self._acsi = int(acsi)

    
    def __iter__(self):
        """Iterate over atoms."""
        acsi = self._acsi
        ag = self._ag 
        for index in self._indices:
            yield Atom(ag, index, acsi)
    
    def __len__(self):
        return len(self._indices)
    
    def __invert__(self):
        
        arange = range(self._ag.getNumOfAtoms())
        indices = list(self._indices)
        while indices:
            arange.pop(indices.pop())
        sel = Selection(self._ag, arange, "not ({0:s}) ".format(
                                                self.getSelectionString()),
                        self._acsi)        
        return sel
    
    def __or__(self, other):
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
        if self is other:
            return self
        acsi = self._acsi
        if acsi != other._acsi:
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0
        if isinstance(other, Atom):
            other_indices = np.array([other._index])
        else:
            other_indices = other._indices
        indices = np.unique(np.concatenate((self._indices, other_indices)))
        return Selection(self._ag, indices, 
                         '({0:s}) or ({1:s})'.format(self.getSelectionString(), 
                                                    other.getSelectionString()),
                          acsi)

    def __and__(self, other):
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
        if self is other:
            return self
        acsi = self._acsi
        if acsi != other._acsi:
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0
        indices = set(self._indices)
        if isinstance(other, Atom):
            other_indices = set([other._index])
        else:
            other_indices = set(other._indices)
        indices = indices.intersection(other_indices)
        indices = np.unique(indices)
        return Selection(self._ag, indices, 
                         '({0:s}) and ({1:s})'.format(self.getSelectionString(), 
                                                    other.getSelectionString()),
                         acsi)    
    
    def getAtomGroup(self):
        """Return associated atom group."""
        return self._ag

    def getIndices(self):
        """Return the indices of atoms."""
        return self._indices.copy()
    
    def getNumOfAtoms(self):
        """Return number of atoms."""
        return self._indices.__len__()

    def getCoordinates(self):
        """Return coordinates from the active coordinate set."""
        if self._ag._coordinates is None:
            return None
        return self._ag._coordinates[self._acsi, self._indices].copy()
    
    def setCoordinates(self, coordinates):
        """Set coordinates in the active coordinate set."""
        self._ag._coordinates[self._acsi, self._indices] = coordinates
        
    def getCoordsets(self, indices):
        """Return coordinate sets at given *indices*.
        
        *indices* may be an integer or a list of integers.
        
        """
        if self._ag._coordinates is None:
            return None
        if indices is None:
            indices = slice(None)
        try: 
            return self._ag._coordinates[indices, self._indices].copy()
        except IndexError:
            raise IndexError('indices may be an integer or a list of integers')

    def getNumOfCoordsets(self):
        """Return number of coordinate sets."""
        return self._ag._n_coordsets
    
    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each coordinate set."""
        for i in range(self._ag._n_coordsets):
            yield self._ag._coordinates[i, self._indices].copy()

    def getActiveCoordsetIndex(self):
        """Return the index of the active coordinate set."""
        return self._acsi
    
    def setActiveCoordsetIndex(self, index):
        """Set the index of the active coordinate set."""
        if self._ag._coordinates is None:
            return None
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if self._ag._n_coordsets <= index or \
           self._ag._n_coordsets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += self._ag._n_coordsets
        self._acsi = index
    
    def getAtomNames(self):
        """Return atom names of the atoms."""
        if self._ag._atomnames is None:
            return None
        return self._ag._atomnames[self._indices]
    
    def setAtomNames(self, atom_names):
        """Set atom names of the atoms."""
        if self._ag._atomnames is None:
            raise AttributeError('atom names are not set')
        self._ag._atomnames[self._indices] = atom_names
    
    def getAlternateLocationIndicators(self):
        """Return alternate location indicators of the atoms."""
        if self._ag._altlocs is None:
            return None
        return self._ag._altlocs[self._indices]
    
    def setAlternateLocationIndicators(self, alternate_location_indicators):
        """Set alternate location indicators of the atoms."""
        if self._ag._altlocs is None:
            raise AttributeError('alternate location indicators are not set')
        self._ag._altlocs[self._indices] = alternate_location_indicators
    
    def getAnisotropicTemperatureFactors(self):
        """Return anisotropic temperature factors of the atoms."""
        if self._ag._anisou is None:
            return None
        return self._ag._anisou[self._indices]
    
    def setAnisotropicTemperatureFactors(self, anisotropic_temperature_factors):
        """Set anisotropic temperature factors of the atoms."""
        if self._ag._anisou is None:
            raise AttributeError('anisotropic temperature factors are not set')
        self._ag._anisou[self._indices] = anisotropic_temperature_factors
    
    def getChainIdentifiers(self):
        """Return chain identifiers of the atoms."""
        if self._ag._chainids is None:
            return None
        return self._ag._chainids[self._indices]
    
    def setChainIdentifiers(self, chain_identifiers):
        """Set chain identifiers of the atoms."""
        if self._ag._chainids is None:
            raise AttributeError('chain identifiers are not set')
        self._ag._chainids[self._indices] = chain_identifiers
    
    def getElementSymbols(self):
        """Return element symbols of the atoms."""
        if self._ag._elements is None:
            return None
        return self._ag._elements[self._indices]
    
    def setElementSymbols(self, element_symbols):
        """Set element symbols of the atoms."""
        if self._ag._elements is None:
            raise AttributeError('element symbols are not set')
        self._ag._elements[self._indices] = element_symbols
    
    def getHeteroFlags(self):
        """Return hetero flags of the atoms."""
        if self._ag._hetero is None:
            return None
        return self._ag._hetero[self._indices]
    
    def setHeteroFlags(self, hetero_flags):
        """Set hetero flags of the atoms."""
        if self._ag._hetero is None:
            raise AttributeError('hetero flags are not set')
        self._ag._hetero[self._indices] = hetero_flags
    
    def getOccupancies(self):
        """Return occupancies of the atoms."""
        if self._ag._occupancies is None:
            return None
        return self._ag._occupancies[self._indices]
    
    def setOccupancies(self, occupancies):
        """Set occupancies of the atoms."""
        if self._ag._occupancies is None:
            raise AttributeError('occupancies are not set')
        self._ag._occupancies[self._indices] = occupancies
    
    def getResidueNames(self):
        """Return residue names of the atoms."""
        if self._ag._resnames is None:
            return None
        return self._ag._resnames[self._indices]
    
    def setResidueNames(self, residue_names):
        """Set residue names of the atoms."""
        if self._ag._resnames is None:
            raise AttributeError('residue names are not set')
        self._ag._resnames[self._indices] = residue_names
    
    def getResidueNumbers(self):
        """Return residue numbers of the atoms."""
        if self._ag._resnums is None:
            return None
        return self._ag._resnums[self._indices]
    
    def setResidueNumbers(self, residue_numbers):
        """Set residue numbers of the atoms."""
        if self._ag._resnums is None:
            raise AttributeError('residue numbers are not set')
        self._ag._resnums[self._indices] = residue_numbers
    
    def getSecondaryStructureAssignments(self):
        """Return secondary structure assignments of the atoms."""
        if self._ag._secondary is None:
            return None
        return self._ag._secondary[self._indices]
    
    def setSecondaryStructureAssignments(self, secondary_structure_assignments):
        """Set secondary structure assignments of the atoms."""
        if self._ag._secondary is None:
            raise AttributeError('secondary structure assignments are not set')
        self._ag._secondary[self._indices] = secondary_structure_assignments
    
    def getSegmentNames(self):
        """Return segment names of the atoms."""
        if self._ag._segnames is None:
            return None
        return self._ag._segnames[self._indices]
    
    def setSsegmentNames(self, segment_names):
        """Set segment names of the atoms."""
        if self._ag._segnames is None:
            raise AttributeError('segment names are not set')
        self._ag._segnames[self._indices] = segment_names
    
    def getAnisotropicStandardDeviations(self):
        """Return standard deviations for the anisotropic temperature
        factors of the atoms."""
        if self._ag._siguij is None:
            return None
        return self._ag._siguij[self._indices]
    
    def setAnisotropicStandardDeviations(self, anisotropic_standard_deviations):
        """Set standard deviations for the anisotropic temperature factors of 
        the atoms."""
        if self._ag._siguij is None:
            raise AttributeError('standard deviations for the anisotropic temperature factors are not set')
        self._ag._siguij[self._indices] = anisotropic_standard_deviations
    
    def getTemperatureFactors(self):
        """Return temperature (B) factors of the atoms."""
        if self._ag._bfactors is None:
            return None
        return self._ag._bfactors[self._indices]
    
    def setTemperatureFactors(self, temperature_factors):
        """Set temperature (B) factors of the atoms."""
        if self._ag._bfactors is None:
            raise AttributeError('temperature (B) factors are not set')
        self._ag._bfactors[self._indices] = temperature_factors
    
    def getRadii(self):
        """Return atomic radii of the atoms."""
        if self._ag._radii is None:
            return None
        return self._ag._radii[self._indices]
    
    def setRadii(self, radii):
        """Set atomic radii of the atoms."""
        if self._ag._radii is None:
            raise AttributeError('atomic radii are not set')
        self._ag._radii[self._indices] = radii
    
    def getMasses(self):
        """Return atomic masses of the atoms."""
        if self._ag._masses is None:
            return None
        return self._ag._masses[self._indices]
    
    def setMasses(self, masses):
        """Set atomic masses of the atoms."""
        if self._ag._masses is None:
            raise AttributeError('atomic masses are not set')
        self._ag._masses[self._indices] = masses
    
    def getCharges(self):
        """Return atomic partial charges of the atoms."""
        if self._ag._charges is None:
            return None
        return self._ag._charges[self._indices]
    
    def setCharges(self, charges):
        """Set atomic partial charges of the atoms."""
        if self._ag._atomtypes is None:
            raise AttributeError('atomic partial charges are not set')
        self._ag._charges[self._indices] = charges
    
    def getAtomTypes(self):
        """Return atom types of the atoms."""
        if self._ag._atomtypes is None:
            return None
        return self._ag._atomtypes[self._indices]
    
    def setAtomTypes(self, atom_types):
        """Set atom types of the atoms."""
        if self._ag._atomtypes is None:
            raise AttributeError('atom types are not set')
        self._ag._atomtypes[self._indices] = atom_types
    
    def getInsertionCodes(self):
        """Return insertion codes of the atoms."""
        if self._ag._icodes is None:
            return None
        return self._ag._icodes[self._indices]
    
    def setInsertionCodes(self, insertion_codes):
        """Set insertion codes of the atoms."""
        if self._ag._icodes is None:
            raise AttributeError('insertion codes are not set')
        self._ag._icodes[self._indices] = insertion_codes

    def select(self, selstr):
        """Return a selection matching the given selection criteria."""
        return SELECT.select(self, selstr)

class Chain(AtomSubset):
    
    __slots__ = AtomSubset.__slots__ + ['_chid', '_seq', '_dict']
    
    def __init__(self, atomgroup, indices, chainid, acsi=None):
        AtomSubset.__init__(self, atomgroup, indices, acsi)
        self._chid = chainid
        self._seq = None
        self._dict = dict()
        
    def __repr__(self):
        return ('<Chain: {0:s} from {1:s} ({2:d} atoms; '
                '{3:d} coordinate sets, active set index: {4:d})>').format(
                self._chid, self._ag.getName(), len(self), 
                self._ag._n_coordsets, self._acsi)

    def __str__(self):
        return ('Chain {0:s} from {1:s}').format(self._chid, self._ag.getName())
        return ('Chain {0:s} from {1:s} ({2:d} atoms; '
                '{3:d} coordinate sets, active set index: {4:d})').format(
                self._chid, self._ag._name, len(self), 
                self._ag._n_coordsets, self._acsi)

    def __getitem__(self, number):
        """Returns the residue with given number, if it exists. Assumes
        the insertion code is an empty string."""
        return self.getResidue(number)
    
    def getResidue(self, number, insertcode=''):
        """Return residue with given number."""
        return self._dict.get((number, insertcode), None)

    def iterResidues(self):
        """Iterate residues in the chain."""
        keys = self._dict.keys()
        keys.sort()
        for key in keys:
            yield self._dict[key]
    
    def getNumOfResidues(self):
        """Return number of residues."""
        return len(self._dict)

    def getIdentifier(self):
        """Return chain identifier."""
        return self._chid
    
    def setIdentifier(self, identifier):
        """Set chain identifier."""
        self.chain_identifiers = identifier
    
    def getSequence(self):
        """Return sequence, if chain is a polypeptide."""
        if self._seq:
            return self._seq
        CAs = self.select('name CA').select('protein')
        if len(CAs) > 0:
            self._seq = prody.proteins.compare._getSequence(CAs.residue_names)
        else:
            self._seq = ''
        return self._seq

    def getSelectionString(self):
        """Return selection string that selects this atom group."""
        return 'chain {0:s}'.format(self._chid)


class Residue(AtomSubset):
    
    __slots__ = AtomSubset.__slots__ + ['_number', '_name', '_icode', '_chain']
    
    def __init__(self, atomgroup, indices, number, name, icode, chain, acsi=None):
        AtomSubset.__init__(self, atomgroup, indices, acsi)
        self._number = number
        self._name = name
        self._icode = icode
        self._chain = chain

    def __repr__(self):
        return '<Residue: {0:s}>'.format(str(self))
        
    def __str__(self):
        return ('{0:s} {1:d}{2:s} from Chain {3:s} from {4:s} '
                '({5:d} atoms; {6:d} coordinate sets, '
                'active set index: {7:d})').format(self._name, self._number, 
                self._icode, self._chain._chid, self._ag._name, len(self), 
                self._ag._n_coordsets, self._acsi)


    def __getitem__(self, name):
        return self.getAtom(name)
    
    def getAtom(self, name):
        """Return atom with given *name*.
        
        Assumes that atom names in a residue are unique. If more than one atoms 
        with the given *name* exists, the one with the smaller index will be 
        returned.
        
        """
        nz = (self.getAtomNames() == name).nonzero()[0]
        if len(nz) > 0:
            return Atom(self._ag, self._indices[nz[0]], self._acsi)
            
    
    def getChain(self):
        """Return the chain that the residue belongs to."""
        return self._chain
    
    def getNumber(self):
        """Return residue number."""
        return self._number
    
    def setNumber(self, number):
        """Set residue number."""
        self.setResidueNumbers(number)
        self._number = number
    
    def getName(self):
        """Return residue name."""
        return self._name
    
    def setName(self, name):
        """Set residue name."""
        self.setResidueNames(name)
        self._name = name

    def getInsertionCode(self):
        """Return residue insertion code."""
        return self._icode
        
    def setInsertionCode(self, icode):
        """Set residue insertion code."""
        self.setInsertionCodes(icode)
        self._icode = icode
    
    def getChainIdentifier(self):
        return self._chain.getIdentifier()
    
    def getSelectionString(self):
        """Return selection string that will select this residue."""
        return 'chain {0:s} and resnum {1:d}{2:s}'.format(
                self._chain._chid, self._number, self._icode)


class Selection(AtomSubset):
    """A class for accessing and manipulating attributes of select of atoms 
    in an :class:`AtomGroup` instance."""
    
    __slots__ = AtomSubset.__slots__ + ['_selstr']
    
    def __init__(self, atomgroup, indices, selstr, acsi=None):
        AtomSubset.__init__(self, atomgroup, indices, acsi)
        self._selstr = str(selstr)
        
    def __repr__(self):
        selstr = self._selstr
        if len(selstr) > 33:
            selstr = selstr[:15] + '...' + selstr[-15:]  
        return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms; '
                '{3:d} coordinate sets, active set index: {4:d})>').format(
                selstr, self._ag.getName(), len(self), 
                         self._ag._n_coordsets, self._acsi)
        return '<Selection: {0:s}>'.format(str(self))
                   
    def __str__(self):
        selstr = self._selstr
        if len(selstr) > 33:
            selstr = selstr[:15] + '...' + selstr[-15:]  
        return 'Selection "{0:s}" from {1:s}'.format(selstr, self._ag.getName())
        
    
    def getSelectionString(self):
        """Return selection string that selects this atom subset."""
        return self._selstr

    def getHierView(self):
        """Return a hierarchical view of the atom subset."""
        LOGGER.warning('HierView will be disabled for selections.')
        return prody.proteins.HierView(self)
