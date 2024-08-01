# -*- coding: utf-8 -*-
"""This module defines base class :class:`Atomic` that all other
:mod:`~prody.atomic` classes are derived from."""

from numpy import all, arange
from os import path
from prody import LOGGER, __path__
from prody.utilities import openData

from . import flags
from .bond import trimBonds
from .fields import READONLY

__all__ = ['Atomic', 'AAMAP']

NOTALLNONE = set(['not', 'all', 'none', 'index', 'sequence', 'x', 'y', 'z'])

MODMAP = {}
with openData('mod_res_map.dat') as f:
    for line in f:
        try:
            mod, aa = line.strip().split(' ')
            MODMAP[mod] = aa
        except:
            continue

AAMAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
    'ASX': 'B', 'GLX': 'Z', 'SEC': 'U', 'PYL': 'O', 'XLE': 'J', '': '-'
}

# add bases
AAMAP.update({'ADE': 'a', 'THY': 't', 'CYT': 'c',
              'GUA': 'g', 'URA': 'u'})

# add reversed mapping
_ = {}
for aaa, a in AAMAP.items():
    _[a] = aaa
AAMAP.update(_)

# add modified AAs
MODAAMAP = {}
for mod, aa in MODMAP.items():
    if aa in AAMAP:
        MODAAMAP[mod] = AAMAP[aa]
AAMAP.update(MODAAMAP)

class Atomic(object):

    """Base class for all atomic classes that can be used for type checking."""

    __slots__ = []

    def __getattribute__(self, name):

        try:
            return object.__getattribute__(self, name)

        except AttributeError:
            if name.startswith('is') and self.isFlagLabel(name[2:]):
                return all(self._getFlags(name[2:]))
            else:
                if name == 'all':
                    try:
                        ag = self.getAtomGroup()
                    except AttributeError:
                        ag = self
                        selstr = name
                        return Selection(ag, arange(self.numAtoms()), 'all',
                                         self._acsi, unique=True)
                    else:
                        try:
                            dummies = self.numDummies()
                        except AttributeError:
                            return Selection(ag, self.getIndices(),
                                             self.getSelstr(),
                                             self._acsi, unique=True)
                        else:
                            return AtomMap(ag, self.getIndices(), self._acsi,
                                           intarrays=True, dummies=dummies,
                                           title=self.getTitle())
                elif name == 'none':
                    return None
                elif self.isFlagLabel(name):
                    try:
                        ag = self.getAtomGroup()
                    except AttributeError:
                        ag = self
                        selstr = name
                    else:
                        selstr = '({0}) and ({1})'.format(name,
                                                          self.getSelstr())
                    try:
                        dummies = self.numDummies()
                    except AttributeError:
                        indices = self._getSubset(name)
                        if len(indices):
                            return Selection(ag, indices, selstr,
                                             self._acsi, unique=True)
                        else:
                            return None
                    else:
                        indices = self._getSubset(name)
                        if len(indices):
                            return AtomMap(ag, indices, self._acsi,
                                           intarrays=True, dummies=dummies,
                                           title='Selection ' + repr(name) +
                                                 ' from ' + str(self))
                        else:
                            return None
                elif name == '_anisous':
                    return None
                else:
                    selstr = name
                    items = name.split('_')
                    word = items[0]
                    if (self.isFlagLabel(word) or self.isDataLabel(word) or
                       word in NOTALLNONE or isSelectionMacro(word)):
                        selstr = ' '.join(items)
                        return SELECT.select(self, selstr)

        raise AttributeError('{0} object has no attribute `{1}` and {2} '
                             'is not a valid selection string'
                             .format(self.__class__.__name__, name,
                                     repr(selstr)))

    def __getstate__(self):

        return dict([(slot, getattr(self, slot))
                     for slot in self.__class__.__slots__])

    def __setstate__(self, state):

        for slot in self.__class__.__slots__:
            try:
                value = state[slot]
            except KeyError:
                pass
            else:
                setattr(self, slot, value)

    def copy(self):
        """Returns a copy of atoms (and atomic data) in an :class:`.AtomGroup`
        instance."""

        dummies = None
        indices = None
        readonly = False
        try:
            ag = self.getAtomGroup()
        except AttributeError:
            ag = self
            readonly = True
            new = AtomGroup(ag.getTitle())
        else:
            indices = self.getIndices()
            new = AtomGroup(ag.getTitle() + ' ' + str(self))
            try:
                dummies = self.numDummies()
            except AttributeError:
                pass
            else:
                if dummies:
                    dummy = self.getFlags('dummy')
                    mapped = self.getFlags('mapped')

        try:
            self.getIndex()
        except AttributeError:
            this = self
        else:
            this = self.all

        if self.numCoordsets():
            new.setCoords(this.getCoordsets(), label=ag.getCSLabels())

        for label in ag.getDataLabels():
            if label in READONLY:
                if readonly:
                    new._data[label] = this.getData(label)
            else:
                new.setData(label, this.getData(label))

        #if readonly:
        #    for label in READONLY:
        #        data = this.getData(label)
        #        if data is not None:
        #            new._data[label] = data

        skip_flags = set()
        for label in ag.getFlagLabels():
            if label in skip_flags:
                continue
            else:
                new._setFlags(label, this.getFlags(label))
                skip_flags.update(flags.ALIASES.get(label, [label]))

        if dummies:
            new._setFlags('dummy', dummy)
            new._setFlags('mapped', mapped)

        bonds = ag._bonds
        bmap = ag._bmap
        if bonds is not None and bmap is not None:
            if indices is None:
                new._bonds = bonds.copy()
                new._bmap = bmap.copy()
                new._data['numbonds'] = ag._data['numbonds'].copy()
            elif dummies:
                if dummies:
                    indices = indices[self._getMapping()]
                if len(set(indices)) == len(indices):
                    new.setBonds(trimBonds(bonds, indices))
                else:
                    LOGGER.warn('Duplicate atoms in mapping, bonds are '
                                'not copied.')
            else:
                bonds = trimBonds(bonds, indices)
                if bonds is not None:
                    new.setBonds(bonds)
        return new

    __copy__ = copy
    toAtomGroup = copy

    def select(self, selstr, **kwargs):
        """Returns atoms matching *selstr* criteria.  See :mod:`~.select` module
        documentation for details and usage examples."""

        return SELECT.select(self, selstr, **kwargs)

    def getTitle(self):
        """Returns title of the instance."""
        try:
            ag = self.getAtomGroup()
        except AttributeError:
            ag = self
        return ag._title
    
    def getSequence(self, **kwargs):
        """Returns one-letter sequence string for amino acids.
        When *allres* keyword argument is **True**, sequence will include all
        residues (e.g. water molecules) in the chain and **X** will be used for
        non-standard residue names."""

        get = AAMAP.get
        if hasattr(self, 'getResnames'):
            seq = ''.join([get(res, 'X') for res in self.getResnames()])
        else:
            res = self.getResname()
            seq = get(res, 'X')
        
        return seq

    def getHierView(self, **kwargs):
        """Returns a hierarchical view of the atom selection."""

        return HierView(self, **kwargs)

    def numSegments(self):
        """Returns number of segments."""

        return self.getHierView().numSegments()

    def numChains(self):
        """Returns number of chains."""

        return self.getHierView().numChains()

    def numResidues(self):
        """Returns number of residues."""

        return self.getHierView().numResidues()

    def toTEMPyAtoms(self):
        """Returns a :class:`TEMPy.protein.prot_rep_biopy.Atom` or list of them as appropriate"""
        try:
            from TEMPy.protein.prot_rep_biopy import Atom as TEMPyAtom
        except ImportError:
            raise ImportError('TEMPy is needed for this functionality')

        if hasattr(self, 'getResnums'):
            return [atom.toTEMPyAtom() for atom in self if atom is not None]
        else:
            return [self.toTEMPyAtom()]

    def toTEMPyStructure(self):
        """Returns a :class:`.protein.prot_rep_biopy.Structure` object""" 
        try:
            from TEMPy.protein.prot_rep_biopy import BioPy_Structure
        except ImportError:
            raise ImportError('TEMPy is needed for this functionality')

        return BioPy_Structure(self.toTEMPyAtoms())

    def toBioPythonStructure(self, header=None, **kwargs):
        """Returns a :class:`Bio.PDB.Structure` object

        :arg atoms: an object with atom and coordinate data
        :type atoms: :class:`.Atomic`

        :arg csets: coordinate set indices, default is all coordinate sets
        """ 
        try:
            from Bio.PDB.Structure import Structure
            from Bio.PDB.StructureBuilder import StructureBuilder
            from Bio.PDB.PDBParser import PDBParser
            from Bio.PDB.PDBExceptions import PDBConstructionException
        except ImportError:
            raise ImportError('Bio StructureBuilder could not be imported. '
                'Reinstall ProDy or install Biopython '
                'to solve the problem.')

        origACSI = self.getACSIndex()

        csets = kwargs.get('csets', None)
        if csets is None:
            csets = range(self.numCoordsets())

        structure_builder = StructureBuilder()
        structure_builder.init_structure(self.getTitle())
        if header is not None:
            structure_builder.set_header(header)

        result = structure_builder.get_structure()
        result.is_pqr = (self.getCharges() is not None 
                         and self.getRadii() is not None)
        
        for i in csets:
            self.setACSIndex(i)
            structure_builder.init_model(i)

            current_segid = None
            current_chain_id = None
            current_residue_id = None

            for global_line_counter, atom in enumerate(self):
                segid = atom.getSegname()
                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)

                chainid = atom.getChid()
                resname = atom.getResname()

                if atom.getFlag('hetatm'):
                    if atom.getFlag('water'):
                        hetero_flag = 'W'
                    else:
                        hetero_flag = 'H'
                else:
                    hetero_flag = ' '

                resseq = atom.getResnum()
                icode = atom.getIcode()
                if len(icode) == 0:
                    icode = ' '
                residue_id = (hetero_flag, resseq, icode)

                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        result._handle_PDB_exception(message, global_line_counter)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        result._handle_PDB_exception(message, global_line_counter)

                name = atom.getName()
                coord = atom.getCoords()
                altloc = atom.getAltloc()
                fullname = atom.getName()
                serial_number = atom.getSerial()
                element = atom.getElement()

                if not result.is_pqr:
                    # init atom with pdb fields
                    try:
                        structure_builder.init_atom(
                            name,
                            coord,
                            atom.getBeta(),
                            atom.getOccupancy(),
                            altloc,
                            fullname,
                            serial_number,
                            element,
                        )
                    except PDBConstructionException as message:
                        result._handle_PDB_exception(message, global_line_counter)
                else:
                    try:
                        structure_builder.init_atom(
                            name,
                            coord,
                            atom.getCharge(),
                            atom.getRadius(),
                            altloc,
                            fullname,
                            serial_number,
                            element,
                            atom.getCharge(),
                            atom.getRadius(),
                            result.is_pqr,
                        )
                    except PDBConstructionException as message:
                        result._handle_PDB_exception(message, global_line_counter)

                if atom.getAnisou() is not None:
                    structure_builder.set_anisou(atom.getAnisou())

                if atom.getAnistd() is not None:
                    structure_builder.set_siguij(atom.getAnistd())

        self.setACSIndex(origACSI)

        return result