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

""":mod:`selector` module defines a class for selecting subsets of atoms based 
on a string.

Classes:

  * :class:`Select`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import time
#import re

import numpy as np
from . import pyparsing as pp
pp.ParserElement.enablePackrat()

import prody
from prody import ProDyLogger as LOGGER
from . import keywords

DEBUG = True

__all__ = ['Select']


""" TODO
- make classmethods from methods not related to instances
- evalonly in within and sameas 
"""

class SelectionError(Exception):    
    pass
    
class Select(object):

    """Select subsets of atoms based on a selection string.
    
    Definitions of single word keywords, such as :term:`protein`, 
    :term:`backbone`, :term:`polar`, etc., are stored in :mod:`definitions`
    module and may be altered using functions in this module. 

    """
    # Numbers and ranges
    PLUSORMINUS = pp.Literal('+') | pp.Literal('-')
    NUMBER = pp.Word(pp.nums) 
    INTEGER = pp.Combine( pp.Optional(PLUSORMINUS) + NUMBER )
    E = pp.CaselessLiteral('E')
    FLOAT = pp.Combine(INTEGER + 
                       pp.Optional(pp.Literal('.') + pp.Optional(NUMBER)) +
                       pp.Optional( E + INTEGER )
                       )
    
    NOT_READY = ('helix', 'alpha_helix', 'helix_3_10', 'pi_helix',
                'sheet', 'extended_beta', 'bridge_beta', 'turn', 'coil', 'purine', 'pyrimidine')
    KEYWORDS_BOOLEAN = ('all', 'none', 'protein', 'nucleic', 'hetero', 'water',
                        'backbone', 'sidechain', 'calpha', 'acidic', 'basic', 'polar', 'charged',
                        'neutral', 'aliphatic', 'hydrophobic', 'aromatic', 'cyclic', 'acyclic', 
                        'noh', 'hydrogen', 'large', 'medium', 'small')
    KEYWORDS_FLOAT = ('x', 'y', 'z', 'beta', 'mass', 'occupancy', 'mass', 'radius', 'charge')
    KEYWORDS_INTEGER = ('serial', 'index', 'resnum', 'resid')
    KEYWORDS_STRING = ('name', 'type', 'resname', 'chain', 'element', 'segment')
    KEYWORDS_NUMERIC = KEYWORDS_FLOAT + KEYWORDS_INTEGER    
    KEYWORDS_VALUE_PAIRED = KEYWORDS_NUMERIC + KEYWORDS_STRING 

    def isFloatKeyword(keyword):
        return keyword in Select.KEYWORDS_FLOAT
    isFloatKeyword = staticmethod(isFloatKeyword)

    def isNumericKeyword(keyword):
        return keyword in Select.KEYWORDS_NUMERIC
    isNumericKeyword = staticmethod(isNumericKeyword)

    def isAlnumKeyword(keyword):
        return keyword in Select.KEYWORDS_STRING
    isAlnumKeyword = staticmethod(isAlnumKeyword)

    def isValuePairedKeyword(keyword):
        return keyword in Select.KEYWORDS_VALUE_PAIRED
    isValuePairedKeyword = staticmethod(isValuePairedKeyword)

    def isBooleanKeyword(keyword):
        return keyword in Select.KEYWORDS_BOOLEAN
    isBooleanKeyword = staticmethod(isBooleanKeyword)
        
    def isKeyword(keyword):
        return (Select.isBooleanKeyword(keyword) or 
                Select.isValuePairedKeyword(keyword) or
                Select.isNumericKeyword(keyword))
    isKeyword = staticmethod(isKeyword)
    
    def __init__(self):
        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._selstr = None
        self._evalonly = None
        
        self._coordinates = None
        self._atomnames = None
        self._resnames = None
        self._chainids = None
        self._elements = None
        self._resnames = None
        self._bfactors = None
        self._occupancies = None
        self._radius = None
        self._charges = None
        self._masses = None
        self._segnames = None
        self._resnums = None
        self._icodes = None
        self._atomtypes = None
        self._kdtree = None
        
        self._tokenizer = pp.operatorPrecedence(
             pp.OneOrMore(pp.Word(pp.alphanums+'.:_') | 
             pp.Group(pp.oneOf('" \'') + pp.Word(pp.alphanums+'~!@#$%^&*_+`=:;,.<>?|') + pp.oneOf('" \' "r \'r'))),
             [(pp.oneOf('+ -'), 1, pp.opAssoc.RIGHT, self._sign),
              (pp.oneOf('sqrt pow sq abs floor ceil sin cos tan atan asin acos sinh cosh tanh exp log log10'), 1, pp.opAssoc.RIGHT, self._func),
              (pp.oneOf('** ^'), 2, pp.opAssoc.LEFT, self._pow),
              (pp.oneOf('* / %'), 2, pp.opAssoc.LEFT, self._mul),
              (pp.oneOf('+ -'), 2, pp.opAssoc.LEFT, self._add),
              (pp.oneOf('< > <= >= == = !='), 2, pp.opAssoc.LEFT, self._comp),
              (pp.Regex('same [a-z]+ as') | pp.Regex('(ex)?within [0-9]+\.?[0-9]* of'), 1, pp.opAssoc.RIGHT, self._special),
              (pp.Keyword('!!!'), 1, pp.opAssoc.RIGHT, self._not),
              (pp.Keyword('&&&'), 2, pp.opAssoc.LEFT, self._and),
              (pp.Keyword('||'), 2, pp.opAssoc.LEFT, self._or),]
            )
        self._tokenizer.setParseAction(self._action)
        
    def _reset(self):
        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._selstr = None
        self._evalonly = None

        self._coordinates = None
        self._atomnames = None
        self._resnames = None
        self._chainids = None
        self._elements = None
        self._resnames = None
        self._bfactors = None
        self._occupancies = None
        self._radius = None
        self._charges = None
        self._masses = None
        self._segnames = None
        self._resnums = None
        self._icodes = None
        self._atomtypes = None
        self._kdtree = None
                    
    def _getAtomNames(self):
        if self._atomnames is None:
            if self._ag._atomnames is None:
                raise SelectionError('atom names are not set')                
            if self._indices is None:
                self._atomnames = self._ag._atomnames
            else:
                self._atomnames = self._atoms.getAtomNames()
        return self._atomnames

    def _getResidueNames(self):
        if self._resnames is None:
            if self._ag._resnames is None:
                raise SelectionError('residue names are not set')                
            if self._indices is None:
                self._resnames = self._ag._resnames
            else:
                self._resnames = self._atoms.getResidueNames()
        return self._resnames
    
    def _getResidueNumbers(self):
        if self._resnums is None:
            if self._ag._resnums is None:
                raise SelectionError('residue numbers are not set')                
            if self._indices is None:
                self._resnums = self._ag._resnums
            else:
                self._resnums = self._atoms.getResidueNumbers()
        return self._resnums
    
    def _getChainIdentifiers(self):
        if self._chainids is None:
            if self._ag._chainids is None:
                raise SelectionError('chain identifiers are not set')                
            if self._indices is None:
                self._chainids = self._ag._chainids
            else:
                self._chainids = self._atoms.getChainIdentifiers()
        return self._chainids

    def _getInsertionCodes(self):
        if self._icodes is None:
            if self._ag._icodes is None:
                raise SelectionError('insertion codes are not set')                
            if self._indices is None:
                self._icodes = self._ag._icodes
            else:
                self._icodes = self._atoms.getInsertionCodes()
        return self._icodes

    
    def _getSegmentNames(self):
        if self._segnames is None:
            if self._ag._segnames is None:
                raise SelectionError('segment names are not set')                
            if self._indices is None:
                self._segnames = self._ag._segnames
            else:
                self._segnames = self._atoms.getSegmentNames()
        return self._segnames
    
    def _getCoordinates(self):
        if self._coordinates is None:
            if self._indices is None:
                self._coordinates = self._ag._coordinates[self._ag._acsi]
            else:
                self._coordinates = self._atoms.getCoordinates()
        return self._coordinates
    
    def _getTemperatureFactors(self):
        if self._bfactors is None:
            if self._ag._bfactors is None:
                raise SelectionError('temperature factors are not set')                
            if self._indices is None:
                self._bfactors = self._ag._bfactors
            else:
                self._bfactors = self._atoms.getTemperatureFactors()
        return self._bfactors

    def _getOccupancies(self):
        if self._occupancies is None:
            if self._ag._occupancies is None:
                raise SelectionError('occupancies are not set')                
            if self._indices is None:
                self._occupancies = self._ag._occupancies
            else:
                self._occupancies = self._atoms.getOccupancies()
        return self._occupancies

    def _getMasses(self):
        if self._masses is None:
            if self._ag._masses is None:
                raise SelectionError('masses are not set')                
            if self._indices is None:
                self._masses = self._ag._masses
            else:
                self._masses = self._atoms.getMasses()
        return self._masses

    def _getRadii(self):
        if self._radii is None:
            if self._ag._radii is None:
                raise SelectionError('radii are not set')                
            if self._indices is None:
                self._radii = self._ag._radii
            else:
                self._radii = self._atoms.getRadii()
        return self._radii

    def _getCharges(self):
        if self._charges is None:
            if self._ag._charges is None:
                raise SelectionError('charges are not set')                
            if self._indices is None:
                self._charges = self._ag._charges
            else:
                self._charges = self._atoms.getCharges()
        return self._charges
    
    def _getElementSymbols(self):
        if self._elements is None:
            if self._ag._elements is None:
                raise SelectionError('segment names are not set')
            if self._indices is None:
                self._elements = self._ag._elements
            else:
                self._elements = self._atoms.getElementSymbols()
        return self._elements
        
    def _getAtomTypes(self):
        if self._atomtypes is None:
            if self._ag._atomtypes is None:
                raise SelectionError('atom names are not set')                
            if self._indices is None:
                self._atomtypes = self._ag._atomtypes
            else:
                self._atomtypes = self._atoms.getAtomTypes()
        return self._atomtypes

   
    def _getKDTree(self):
        if prody.KDTree is None:
            prody.importBioKDTree()
        if self._kdtree is None:
            kdtree = prody.KDTree(3)
            kdtree.set_coords(self._getCoordinates())
            self._kdtree = kdtree
            return kdtree
        return self._kdtree

    def _evalBoolean(self, keyword):
        if DEBUG: print '_evalBoolean', keyword
        
        if self._evalonly is None:
            n_atoms = self._n_atoms
        else:        
            n_atoms = len(self._evalonly)
        
        if keyword == 'calpha':
            return self._and([['name', 'CA', '&&&', 'protein']])
        elif keyword == 'noh':
            return self._not([['!!!', 'name', (['"', keywords.HYDROGEN_REGEX,'"r'])]])
        elif keyword == 'all':
            return np.ones(n_atoms, np.bool)
        elif keyword == 'none':
            return np.zeros(n_atoms, np.bool)
        elif keyword == 'hydrogen':
            return self._evaluate(['name', (['"', keywords.HYDROGEN_REGEX,'"r'])])
            
        
        atom_names = None
        atom_names_not = False
        residue_names = None
        invert = False
        
        if keyword == 'protein':
            residue_names = keywords.PROTEIN_RESIDUE_NAMES
        elif keyword == 'backbone':
            atom_names = keywords.BACKBONE_ATOM_NAMES
            residue_names = keywords.PROTEIN_RESIDUE_NAMES
        elif keyword == 'acidic':
            residue_names = keywords.ACIDIC_RESIDUE_NAMES
        elif keyword == 'basic':
            residue_names = keywords.BASIC_RESIDUE_NAMES 
        elif keyword == 'charged':
            residue_names = keywords.ACIDIC_RESIDUE_NAMES + keywords.BASIC_RESIDUE_NAMES
        elif keyword == 'aliphatic':
            residue_names = keywords.ALIPHATIC_RESIDUE_NAMES
        elif keyword == 'aromatic':
            residue_names = keywords.AROMATIC_RESIDUE_NAMES
        elif keyword == 'small':
            residue_names = keywords.SMALL_RESIDUE_NAMES
        elif keyword == 'medium':
            residue_names = keywords.MEDIUM_RESIDUE_NAMES
        elif keyword == 'cyclic':
            residue_names = keywords.CYCLIC_RESIDUE_NAMES  
        elif keyword == 'large':
            residue_names = tuple(set(keywords.PROTEIN_RESIDUE_NAMES).difference( 
                    set(keywords.SMALL_RESIDUE_NAMES + keywords.MEDIUM_RESIDUE_NAMES)))
        elif keyword == 'neutral':
            residue_names = keywords.ACIDIC_RESIDUE_NAMES + keywords.BASIC_RESIDUE_NAMES
            invert = True
        elif keyword == 'acyclic':
            residue_names = keywords.CYCLIC_RESIDUE_NAMES
            invert = True
        elif keyword in ('water', 'waters'):
            residue_names = keywords.WATER_RESIDUE_NAMES
        elif keyword == 'nucleic':
            residue_names = keywords.NUCLEIC_RESIDUE_NAMES
        elif keyword == 'hetero':
            residue_names = keywords.NUCLEIC_RESIDUE_NAMES + keywords.PROTEIN_RESIDUE_NAMES
            invert = True
        elif keyword == 'sidechain':
            atom_names = keywords.BACKBONE_ATOM_NAMES
            residue_names = keywords.PROTEIN_RESIDUE_NAMES
            atom_names_not = True
        else:
            raise SelectionError('"{0:s}" is not a valid keyword.'.format(keyword))
            
        resnames = self._getResidueNames()
        if self._evalonly is not None:
            resnames = resnames[self._evalonly]
       
        torf = np.zeros(n_atoms, np.bool)

        if atom_names is None:
            for i in xrange(n_atoms):
                torf[i] = (resnames[i] in residue_names)
        else:
            atomnames = self._getAtomNames()
            if self._evalonly is not None:
                atomnames = atomnames[self._evalonly]
            if atom_names_not:
                    torf[i] = (not atomnames[i] in atom_names and
                               resnames[i] in residue_names)                
            else:
                for i in xrange(n_atoms):
                    torf[i] = (atomnames[i] in atom_names and
                               resnames[i] in residue_names)
            
        if invert:
            torf = np.invert(torf, torf)
        
        return torf
    
    def _numrange(self, token):
        tknstr = ' '.join(token)
        while '  ' in tknstr:
            tknstr = tknstr.replace('  ', ' ')
        tknstr = tknstr.replace(' to ', 'to').replace('to ', 'to').replace(' to', 'to')
        tknstr = tknstr.replace(' : ', ':').replace(': ', ':').replace(' :', ':')
        token = []
        for item in tknstr.split():
            if 'to' in item:
                items = item.split('to')
                if len(items) != 2:
                    raise SelectionError('"{0:s}" is not understood.'.format(' to '.join(items)))
                try:
                    token.append( [float(items[0]), float(items[1])] )
                except:
                    raise SelectionError('"{0:s}" is not understood, "to" must be surrounded by numbers.'.format(' to '.join(items)))
            elif ':' in item:
                items = item.split(':')
                if not len(items) in (2, 3):
                    raise SelectionError('"{0:s}" is not understood.'.format(':'.join(items)))
                try:
                    if len(items) == 2:
                        token.append( (int(items[0]), int(items[1])) )
                    else:
                        token.append( (int(items[0]), int(items[1]), int(items[2])) )
                except:
                    raise SelectionError('"{0:s}" is not understood, ":" must be surrounded by integers.'.format(':'.join(items)))
            elif '.' in item:
                try:
                    token.append( float(item) )
                except:
                    raise SelectionError('"{0:s}" is not understood.'.format(item))
            elif item.isdigit():
                try:
                    token.append( int(item) )
                except:
                    raise SelectionError('"{0:s}" is not understood.'.format(item))
            else:
                token.append( item )
        if DEBUG: print '_numrange', token            
        return token
    
    def _resnum(self, token=None):
        if DEBUG: print '_resnum', token
        if token is None:
            return self._getResidueNumbers() 
        icodes = None
        if self._evalonly is None:
            resids = self._getResidueNumbers()
            n_atoms = self._n_atoms
        else:
            evalonly = self._evalonly
            resids = self._getResidueNumbers()[evalonly]
            n_atoms = len(evalonly)
        torf = np.zeros(n_atoms, np.bool)
        
        for item in self._numrange(token):
            if isinstance(item, str):
                if icodes is None:
                    if self._evalonly is None:
                        icodes = self._getInsertionCodes()
                    else:
                        icodes = self._getInsertionCodes()[evalonly]
                icode = str(item[-1])
                if icode == '_':
                    icode = ''
                torf[(resids == int(item[:-1])) * (icodes == icode)] = True
            elif isinstance(item, list):
                fr = item[0] 
                to = item[1]
                for i in xrange(n_atoms):
                    if fr <= resids[i] <= to:
                        torf[i] = True            
            elif isinstance(item, tuple):
                if len(item) == 2:
                    fr = item[0] 
                    to = item[1]
                    for i in xrange(n_atoms):
                        if fr <= resids[i] < to:
                            torf[i] = True
                else:
                    arange = range(item[0], item[1], item[2])
                    for i in xrange(n_atoms):
                            torf[i] = resids[i] in arange
            else:
                torf[resids == item] = True
        return torf
    
    def _index(self, token=None, add=0):
        if token is None:
            if self._indices is not None:
                return self._indices + add
            else:
                return np.arange(add, self._ag._n_atoms + add)
        torf = np.zeros(self._ag._n_atoms, np.bool)
        
        for item in self._numrange(token):
            if isinstance(item, str):
                raise SelectionError('"index/serial {0:s}" is not understood.'.format(item))
            elif isinstance(item, tuple):
                if len(item) == 2:
                    torf[item[0]-add:item[1]-add] = True
                else:
                    torf[item[0]-add:item[1]-add:item[2]-add] = True
            elif isinstance(item, list):
                torf[int(np.ceil(item[0]-add)):int(np.floor(item[1]-add))+1] = True
            else:
                try:
                    torf[int(item)-add] = True
                except IndexError:
                    pass

        if self._indices is not None:
            return torf[self._indices]
        return torf

    def _evalAlnum(self, keyword, values):
        if keyword == 'chain':
            data = self._getChainIdentifiers()
            for i, value in enumerate(values):
                if value == '_':
                    values[i] = ' '
        elif keyword == 'element':
            data = self._getElementSymbols()
        elif keyword == 'name':
            data = self._getAtomNames()
        elif keyword == 'resname':
            data = self._getResidueNames()
        elif keyword == 'segment':
            data = self._getSegmentNames()
        elif keyword == 'type':
            data = self._getAtomTypes()
        else:
            raise SelectionError('"{0:s}" is not a valid keyword.'.format(keyword))
            
        if self._evalonly is not None:
            data = data[self._evalonly]
        n_atoms = len(data)
        torf = np.zeros(n_atoms, np.bool)
        
        for value in values:
            if not isinstance(value, str):
                if len(value[2]) == 1:
                    value = value[1]
                else:
                    if prody.re is None: prody.importRE()
                    value = prody.re.compile(value[1])
                    for i in xrange(n_atoms):
                        torf[i] = (value.match(data[i]) is not None)
                    continue
            torf[ data == value ] = True
        return torf
    
    def _evalFloat(self, keyword, values=None):
        if keyword == 'x':
            data = self._getCoordinates()[:,0]
        elif keyword == 'y':
            data = self._getCoordinates()[:,1]
        elif keyword == 'z':
            data = self._getCoordinates()[:,2]
        elif keyword == 'beta':
            data = self._getTemperatureFactors()   
        elif keyword == 'occupancy':
            data = self._getOccupancies()
        elif keyword == 'mass':
            data = self._getMasses()
        elif keyword == 'radius':
            data = self._getRadii()
        elif keyword == 'charge':
            data = self._getCharges()
        else:
            raise SelectionError('"{0:s}" is not a valid keyword.'.format(keyword))
        
        if values is None:
            return data
    
        if self._evalonly is not None:
            data = data[self._evalonly]
        n_atoms = len(data)
        torf = np.zeros(n_atoms, np.bool)

        for item in self._numrange(values):
            if isinstance(item, str):
                pass
            elif isinstance(item, list):
                fr = item[0] 
                to = item[1]
                for i in xrange(n_atoms):
                    if fr <= data[i] <= to:
                        torf[i] = True            
            elif isinstance(item, tuple):
                if len(item) == 2:
                    fr = item[0] 
                    to = item[1]
                    for i in xrange(n_atoms):
                        if fr <= data[i] < to:
                            torf[i] = True
                else:
                    raise SelectionError('"{0:s}" is not valid for keywords expecting floating values.'.format(':'.join(item)))
            else:
                torf[data == item] = True
        return torf
    
    def select(self, atoms, selstr):
        """Return a Selection (or an AtomMap) of atoms matching *selstr*."""
        
        if not isinstance(atoms, (prody.AtomGroup, prody.AtomSubset, prody.AtomMap)):
            raise TypeError('atoms must be an atom container, not {0:s}'.format(type(atoms)))
        elif not isinstance(selstr, str):
            raise TypeError('selstr must be a string, not a {0:s}'.format(type(selstr)))
        self._reset()
        self._selstr = selstr
        if isinstance(atoms, prody.AtomGroup): 
            self._ag = atoms
            self._atoms = atoms
            self._indices = None
            self._n_atoms = atoms._n_atoms
        else:
            self._ag = atoms.getAtomGroup()
            self._indices = atoms.getIndices()
            if isinstance(atoms, prody.AtomMap):
                self._atoms = prody.proteins.Selection(self._ag, self._indices, '')
                self._atoms._indices = self._indices
            else: 
                self._atoms = atoms
            self._n_atoms = len(self._indices)
        if DEBUG:
            print '_select', selstr
        torf = self._parseSelStr()[0]
        if DEBUG:
            print '_select', torf
        if isinstance(atoms, prody.AtomGroup):
            indices = torf.nonzero()[0]
        else:
            indices = self._indices[torf]
            
        ag = self._ag
        self._reset()
        if isinstance(atoms, prody.AtomMap):
            return prody.AtomMap(ag, indices, np.arange(len(indices)), 
                                 np.array([]),
                                 'Selection "{0:s}" from AtomMap {1:s}'.format(
                                 selstr, atoms.getName()),
                                 atoms.getActiveCoordsetIndex())
        else:
            if isinstance(atoms, prody.AtomSubset):
                selstr = '({0:s}) and ({1:s})'.format(selstr, atoms.getSelectionString())
            return prody.Selection(ag, indices, selstr, 
                                   atoms.getActiveCoordsetIndex())

    def _getStdSelStr(self):
        selstr = self._selstr
        selstr = ' ' + selstr + ' '
        #selstr = selstr.replace('(', ' ( ').replace(')', ' ) ')
        while ' and ' in selstr:
            selstr = selstr.replace(' and ', ' &&& ')
        while ' or ' in selstr:
            selstr = selstr.replace(' or ', ' || ')
        while ' not ' in selstr:
            selstr = selstr.replace(' not ', ' !!! ')
        while '  ' in selstr:
            selstr = selstr.replace('  ', ' ')
        selstr = selstr.strip()
        return selstr

    def _parseSelStr(self):
        selstr = self._getStdSelStr()
        if DEBUG: print '_parseSelStr', selstr
        start = time.time()
        try: 
            tokens = self._tokenizer.parseString(selstr, parseAll=True).asList()
            if DEBUG: print '_parseSelStr', tokens
            return tokens
        except pp.ParseException, err:
            print 'Parse Failure'
            print self._selstr #err.line
            print " "*(err.column-1) + "^"
            raise pp.ParseException(str(err))

    def _special(self, token):
        token = token[0]
        if token[0].startswith('same'):
            return self._sameas(token)
        else:
            return self._within(token, token[0].startswith('exwithin'))
    
    def _within(self, token, exclude):
        terms = token
        if DEBUG: print '_within', terms
        within = float(terms[0].split()[1])
        which = terms[1]
        if not isinstance(which, np.ndarray):
            which = self._evaluate(terms[1:])
        kdtree = self._getKDTree()
        coordinates = self._getCoordinates()
        which = which.nonzero()[0]
        result = []
        for index in which:
            kdtree.search(coordinates[index], within)
            result.append(kdtree.get_indices())
        unique = np.unique(np.concatenate(result))
        
        if self._indices is None:
            torf = np.zeros(self._n_atoms, np.bool)
            torf[unique] = True
        else:
            torf = np.zeros(self._ag._n_atoms, np.bool)
            torf[unique] = True
            torf = torf[self._indices]
        if exclude:
            torf[which] = False
        return torf

    def _sameas(self, token):
        terms = token
        if DEBUG: print '_sameas', terms
        what = token[0].split()[1]
        which = token[1]
        if not isinstance(which, np.ndarray):
            which = self._evaluate(token[1:])
        
        if what == 'residue':
            chainids = self._getChainIdentifiers()
            resids =  self._getResidueNumbers()
            resnum = list(np.unique(resids[which]).astype(np.str))
            torf = np.all(
                [self._evalAlnum('chain', list(np.unique(chainids[which]))),
                 self._resnum(resnum)], 0)
        elif what == 'chain':
            chainids = self._getChainIdentifiers()
            torf = self._evalAlnum('chain', list(np.unique(chainids[which])))        
        elif what == 'segment':
            segnames = self._getSegmentNames()
            torf = self._evalAlnum('segment', list(np.unique(segnames[which]))) 
        return torf

    def _not(self, token):
        if DEBUG: print '_not', token
        torf = self._evaluate(token[0][1:])
        np.invert(torf, torf)
        return torf
    
    def _and(self, tokens):
        if DEBUG: print '_and', tokens
        temp = tokens[0]
        tokenlist = []
        token = []
        while temp:
            tkn = temp.pop(0)
            if tkn == '&&&':
                tokenlist.append(token)
                token = []
            else:
                if Select.isBooleanKeyword(tkn) and not token and temp and temp[0] != '&&&':
                    if DEBUG: print '_and inserting &&&'
                    token.append(tkn)
                    tokenlist.append(token)
                    token = []
                else:
                    token.append(tkn)
        tokenlist.append(token)
        
        if DEBUG: print '_and tokenlist', tokenlist

        for token in tokenlist:
            zero = token[0]
            if isinstance(zero, np.ndarray):                    
                if self._evalonly is None: 
                    self._evalonly = zero.nonzero()[0]
                else:        
                    self._evalonly = self._evalonly[zero[self._evalonly].nonzero()[0]]
            else:
                torf = self._evaluate(token)
                self._evalonly = torf.nonzero()[0]
            if DEBUG: print '_and', self._evalonly
        torf = np.zeros(self._n_atoms, np.bool)
        torf[self._evalonly] = True
        self._evalonly = None
        return torf
    
    def _or(self, tokens):
        if DEBUG: print '_or', tokens
        temp = tokens[0]
        tokenlist = []
        token = []
        while temp:
            tkn = temp.pop(0)
            if tkn == '||':
                tokenlist.append(token)
                token = []
            else:
                token.append(tkn)
        tokenlist.append(token)

        if DEBUG: print '_or tokenlist', tokenlist

        for token in tokenlist:
            zero = token[0]
            if isinstance(zero, np.ndarray):                    
                if self._evalonly is None: 
                    self._evalonly = np.invert(zero).nonzero()[0]
                else:        
                    self._evalonly = self._evalonly[np.invert(zero[self._evalonly]).nonzero()[0]]
            else:
                torf = self._evaluate(token)
                self._evalonly = np.invert(torf).nonzero()[0]
            if DEBUG: print '_or', self._evalonly
        torf = np.ones(self._n_atoms, np.bool)
        torf[self._evalonly] = False
        self._evalonly = None
        return torf
        
    def _evaluate(self, token):
        if DEBUG: print '_evaluate', token
        keyword = token[0]
        if len(token) == 1:
            if Select.isBooleanKeyword(keyword):
                return self._evalBoolean(keyword)
            elif Select.isNumericKeyword(keyword):
                return self._getnum(keyword)
            else:
                raise SelectionError('"{0:s}" is not a valid keyword'.format(keyword))
        elif Select.isBooleanKeyword(keyword):
            return self._and([token])
        elif Select.isAlnumKeyword(keyword):
            return self._evalAlnum(keyword, token[1:])
        elif Select.isFloatKeyword(keyword):
            return self._evalFloat(keyword, token[1:])
        elif keyword in ('resnum', 'resid'):
            return self._resnum(token[1:])
        elif keyword == 'index':
            return self._index(token[1:])
        elif keyword == 'serial':
            return self._index(token[1:], 1)
        for item in token[1:]:
            if Select.isKeyword(item):
                raise SelectionError('"{0:s}" in "{1:s}" is not understood, use "" to escape'.format(item, ' '.join(token)))
        raise SelectionError('{0:s} is not yet implemented'.format(' '.join(token)))
    
    def _action(self, token):
        if DEBUG: print '_action', token
        if isinstance(token[0], np.ndarray):
            return token[0]
        else:
            return self._evaluate(token)

    def _comp(self, token):
        if DEBUG: print '_comp', token
        token = token[0]
        comp = token[1]
        left = self._getnum(token[0])
        if DEBUG: print '_comp', left
        right = self._getnum(token[2])
        if DEBUG: print '_comp', right
        
        if comp == '>':
            return left > right
        elif comp == '<':
            return left < right
        elif comp == '<=':
            return left <= right
        elif comp == '>=':
            return left >= right
        elif comp == '==':
            return left == right
        elif comp == '=':
            return left == right
        elif comp == '!=':
            return left != right
        else:
            raise SelectionError('Unknown error in "{0:s}".'.format(' '.join(token)))

    def _pow(self, token):
        if DEBUG: print '_pow', token
        items = token[0]
        return self._getnum(items[0]) ** self._getnum(items[2])

    def _add(self, token):
        if DEBUG: print '_add', token
        items = token[0]
        left = self._getnum(items[0])
        op = items[1]
        right = self._getnum(items[2])
        if op == '+':
            return left + right
        else:
            return left - right
 
    def _mul(self, token):
        if DEBUG: print '_mul', token
        items = token[0]
        left = self._getnum(items[0])
        op = items[1]
        right = self._getnum(items[2])
        if op == '*':
            return left * right
        elif op == '/':
            if right == 0.0:
                raise ZeroDivisionError(' '.join(items))
            return left / right
        else:
            return left % right

    def _getnum(self, token):
        if DEBUG: print '_getnum', token
        if isinstance(token, np.ndarray):
            return token
        elif Select.isFloatKeyword(token):
            return self._evalFloat(token)
        elif token in ('resnum', 'resid'):
            return self._resnum()
        elif token == 'index':
            return self._index()    
        elif token == 'serial':
            return self._index(None, 1)
        else:
            try:
                num = float(token)
            except ValueError:
                raise SelectionError('"{0:s}" must be a number or a valid keyword'.format(token))
            else:
                return num

    def _sign(self, token):
        token = token[0]
        if DEBUG: print '_sign', token
        num = self._getnum(token[1])
        if token[0] == '-':
            return -num
        return num

    def _func(self, token):
        token = token[0]
        if DEBUG: print '_func', token
        fun = token[0]
        num = token[1] 
        if fun == 'sqrt':
            return np.sqrt(num)
        elif fun == 'sq':
            return np.power(num, 2)
        elif fun == 'abs':
            return np.abs(num)
        elif fun == 'floor':
            return np.floor(num)
        elif fun == 'ceil':
            return np.ceil(num)
        elif fun == 'sin':
            return np.sin(num)
        elif fun == 'cos':
            return np.cos(num)
        elif fun == 'tan':
            return np.tan(num)
        elif fun == 'atan':
            return np.arctan(num)
        elif fun == 'asin':
            return np.arcsin(num)
        elif fun == 'acos':
            return np.arccos(num)
        elif fun == 'sinh':
            return np.sinh(num)
        elif fun == 'cosh':
            return np.cosh(num)
        elif fun == 'tanh':
            return np.tanh(num)
        elif fun == 'exp':
            return np.exp(num)
        elif fun == 'log':
            return np.log(num)
        elif fun == 'log10':
            return np.log10(num)

ProDyAtomSelect = Select()
