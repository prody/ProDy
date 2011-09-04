# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

""" This module defines classes for selecting subsets of atoms and identifying 
contacts, and functions to learn and change definitions of selection keywords.

.. _selections:


Atom selections
===============================================================================


ProDy offers a powerful atom selector. The keywords, selection grammar,
and capabilities of the selector are similar to those found in VMD (|vmd|). 
Small differences between the two should not affect most practical uses of 
atom selections. ProDy selection engine also enables the identification of 
intermolecular contacts. This section describes the keywords and selection 
syntax.

|more| See :ref:`contacts` and :ref:`selection-operations` for more usage
examples.

Getting interactive help
-------------------------------------------------------------------------------

The contents of this web page can be viewed in an interactive session as 
follows:
    
>>> from prody import *
>>> # help(select)

    
Keywords with arguments
-------------------------------------------------------------------------------

Below is the list of keywords that can be used when paired with atomic
attributes as arguments.

================ ================ =============================================
Keyword          Arguments        Description
================ ================ =============================================
name             string           atom name
element          string           element symbol
type [*]         string           atom type
altloc [†‡]      string           one-character alternate location identifier
resname          string           residue name
chain [‡]        string           one-character chain identifier
segment [‡]      string           segment name
secondary [*‡]   string           one-character secondary structure identifier
index            integer, range   internal atom number (starts from 0) 
serial           integer, range   atom serial number (parsed from file)
resnum [§]       integer, range   residue number
resid [§]        integer, range   residue number
x                float, range     x coordinate
y                float, range     y coordinate
z                float, range     z coordinate
beta             float, range     β (temperature) factor
occupancy        float, range     atomic occupancy value
charge [*]       float, range     atomic charge
mass [*]         float, range     atomic mass
radius [*]       float, range     atomic radius
================ ================ =============================================

**[*]** These atomic attributes are not set by the PDB parser when a PDB file 
is parsed. Using them before they are set will raise selection error. 
Secondary structure assignments can be made using 
:func:`~prody.proteins.assignSecondaryStructure` function.

**[†]** Alternate locations are parsed as alternate coordinate sets. This
keyword will work for alternate location specified by "A". This to work for
alternate locations indicated by other letters, they must be parsed 
specifically by passing the identifier to the :func:`~prody.proteins.parsePDB`.

**[‡]** Atoms with unspecified alternate location/chain/segment/secondary 
structure identifiers can be selected using "_". This character is replaced 
with a whitespace.

.. versionchanged:: 0.7.1 
   *segment* keyword is added to the above list.

**[§]** If there are multiple residues with the same number but 
distinguished with insertion codes, the insertion code can be appended
to the residue number. "_" stands for empty insertion code. For example:
    
  * ``"resnum 5"`` selects residue 5 (all insertion codes)
  * ``"resnum 5A"`` selects residue 5 with insertion code A
  * ``"resnum 5_"`` selects residue 5 with no insertion code

**Strings (with special characters)**

Strings can be any combination of the following::

  abcdefghijklmnopqrstuvwxyz
  ABCDEFGHIJKLMNOPQRSTUVWXYZ
  0123456789
  ~@#$.:;_',
  
For example ``"name C' N` O~ C$ C#"`` is a valid selection string. 

.. versionchanged:: 0.7
   Special characters are allowed.

**Integers and floats**

Numbers can be provided as integers or floats, and they will be converted to
appropriate type. For example ``"resnum 10 11.0"`` will select residues
with number 10 and 11, but ``"resnum 10.5"`` will not select anything.

**Number ranges**

Number ranges can be passed as follows:
    
  * ``"resnum 5 10 to 15"`` selects residues 5, 10, 11, 12, 13, 14, and 15
  * ``"resnum 5 10:15"`` selects residues 5, 10, 11, 12, 13, and 14 
    (:, colon, works as it does in Python slicing operations)
  * ``"resnum 1:10:2"`` selects residues 1, 3, 5, 7, and 9
  * ``"x 1 to 10"`` selects atoms whose x coordinates are greater or equal to 1
    or smaller or equal to 10  
  * ``"x 1:10"`` selects atoms whose x coordinates are greater or equal to 1
    or smaller or equal to 10
    

**More special characters (``)**

.. versionadded:: 0.7
   Strings can include the following characters (including whitespace) as well 
   when they are surrounded by grave accent character (``):

::
  
  ~!@#$%^&*()-_=+[{}]\|;:,<>./?()'"

For example ``"name `CA*` `C *`"`` will work.

**Regular expressions ("")**

.. versionadded:: 0.7
   Strings surrounded by double quotes ("") will be treated as regular 
   expressions. The following character set can be used between double 
   quotes:

::
  
  ~!@#$%^&*()-_=+[{}]\|;:,<>./?()'`

For example ``'resname "A.."'`` will select residues whose names start with 
letter A and are three-characters long.

For more information on regular expressions see :mod:`re`. 

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import time
import re as RE

import numpy as np
from . import pyparsing as pp
pp.ParserElement.enablePackrat()
KDTree = None

import prody
from prody import ProDyLogger as LOGGER
from prody.atomic import *
DEBUG = False

__all__ = ['Select', 'Contacts',
           'getProteinResidueNames', 'setProteinResidueNames',
           'getKeywordResidueNames', 'setKeywordResidueNames',
           'getBackboneAtomNames', 'setBackboneAtomNames',
           'getAtomNameRegex', 'setAtomNameRegex',
           'defSelectionMacro', 'delSelectionMacro', 'getSelectionMacro',
           ]

KEYWORDS_STRING = set(('name', 'type', 'resname', 'chain', 'element', 
                       'segment', 'altloc', 'secondary'))
KEYWORDS_INTEGER = set(('serial', 'index', 'resnum', 'resid'))
KEYWORDS_FLOAT = set(('x', 'y', 'z', 'beta', 'mass', 'occupancy', 'mass', 
                      'radius', 'charge'))
KEYWORDS_NUMERIC = KEYWORDS_FLOAT.union(KEYWORDS_INTEGER)    

KEYWORDS_VALUE_PAIRED = KEYWORDS_NUMERIC.union(KEYWORDS_STRING) 


# 21st and 22nd amino acids	    3-Letter	1-Letter
# Selenocysteine	            Sec	        U
# Pyrrolysine	                Pyl	        O

# Ambiguous Amino Acids	                3-Letter	1-Letter
# Asparagine or aspartic acid	        Asx	        B
# Glutamine or glutamic acid	        Glx	        Z
# Leucine or Isoleucine	                Xle	        J
# Unspecified or unknown amino acid     Xaa         X

KEYWORD_RESNAMES = {
    'protein': ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 
                'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
                'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD', 'HSE', 'HSP', 
                'GLX', 'ASX', 'SEC', 'PYL', 'XLE', 'CSO'],
    'nucleic': ['GUA', 'ADE', 'CYT', 'THY', 'URA', 'DA', 'DC', 'DG', 'DT', 
                'A', 'C', 'G', 'T', 'U'],

    'acidic': ['ASP', 'GLU'],
    'aliphatic': ['ALA', 'GLY', 'ILE', 'LEU', 'VAL', 'XLE'],
    'aromatic': ['HIS', 'PHE', 'TRP', 'TYR', 'HSD', 'HSE', 'HSP'],
    'basic': ['LYS', 'ARG', 'HIS', 'HSP', 'HSD'],
    'buried': 'ALA LEU VAL ILE XLE PHE CYS MET TRP'.split(),
    'cyclic': ['HIS', 'PHE', 'PRO', 'TRP', 'TYR', 'HSD', 'HSE', 'HSP'],
    'hydrophobic': ['ALA', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 'VAL', 'XLE'],
    'small': ['ALA', 'GLY', 'SER'],
    'medium': ['VAL', 'THR', 'ASP', 'ASN', 'ASX', 'PRO', 'CYS', 'SEC'],
    
    'water': ['HOH', 'WAT', 'TIP3', 'H2O',
              'HH0', 'OHH', 'OH2', 'SOL', 'TIP', 'TIP2', 'TIP4'],
    'lipid': 'DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE'.split(),
    'heme': 'HEM HEME'.split(),
    'ion': ('AL BA CA CAL CD CES CLA CL CO CS CU CU1 CUA HG IN IOD K MG MN3 '
            'MO3 MO4 MO5 MO6 NA NAW OC7 PB POT PT RB SOD TB TL WO4 YB ZN ZN1 '
            'ZN2').split(),
    'sugar': ['AGLC'],
    
    'at': 'ADE A THY T'.split(),
    'cg': 'CYT C GUA G'.split(),
    'purine': 'ADE A GUA G'.split(),
    'pyrimidine': 'CYT C THY T URA U'.split(),
}

KEYWORD_RESNAMES_READONLY = {
    'acyclic': 'protein and not cyclic',
    'charged': 'acidic or basic',
    'large': 'not (small or medium)',
    'neutral': 'protein and not charged',
    'polar': 'protein and not hydrophobic',
    'surface': 'protein and not buried',
}

def _setReadonlyResidueNames():
    protein = set(KEYWORD_RESNAMES['protein'])
    KEYWORD_RESNAMES['acyclic'] = list(protein.difference(set(KEYWORD_RESNAMES['cyclic'])))
    KEYWORD_RESNAMES['charged'] = list(set(KEYWORD_RESNAMES['acidic'] + KEYWORD_RESNAMES['basic']))
    KEYWORD_RESNAMES['large'] = list(protein.difference(set(KEYWORD_RESNAMES['small'] + KEYWORD_RESNAMES['medium'])))
    KEYWORD_RESNAMES['neutral'] = list(protein.difference(set(KEYWORD_RESNAMES['charged'])))
    KEYWORD_RESNAMES['polar'] = list(protein.difference(set(KEYWORD_RESNAMES['hydrophobic'])))
    KEYWORD_RESNAMES['surface'] = list(protein.difference(set(KEYWORD_RESNAMES['buried'])))
    
_setReadonlyResidueNames()

__doc__ += """

Keywords without arguments
-------------------------------------------------------------------------------

Below is the list of keywords defined based on residue type and/or property.
These definitions can be retrieved or altered using :func:`getKeywordResidueNames` 
and :func:`setKeywordResidueNames`, respectively.

============= =================================================================
Keyword       Description
============= =================================================================
"""
keys = KEYWORD_RESNAMES.keys()
keys.sort()
for key in keys:
    if key in KEYWORD_RESNAMES_READONLY:
        __doc__ += '{0:13s} resname {1:s}\n'.format(key + ' [#]', ' '.join(KEYWORD_RESNAMES[key]))
    else:
        __doc__ += '{0:13s} resname {1:s}\n'.format(key, ' '.join(KEYWORD_RESNAMES[key]))

__doc__ += """============= =================================================================

**[#]** Definitions of these keywords cannot be changed directly, as they 
are defined based on others as follows: 
    
"""
keys = KEYWORD_RESNAMES_READONLY.keys()
keys.sort()
for key in keys:
    __doc__ += '  * "{0:s}" is "{1:s}"\n'.format(key, KEYWORD_RESNAMES_READONLY[key])

__doc__ += """
.. versionchanged:: 0.6.1
   ASX (asparagine or aspartic acid), GLX (lutamine or glutamic acid),
   XLE (leucine or isoleucine), SEC (selenocysteine), and PYL
   (pyrrolysine) has been added to the standard definition of "protein" 
   keyword. Note that this list of residue names can be changed using
   :func:`setProteinResidueNames` function.

The following are additional keywords whose definitions are more restricted:

=============== ===============================================================
Keyword         Description
=============== ===============================================================
all             all atoms
none            nothing (returns ``None``)
hetero          non-protein/nucleic atoms, same as ``"not (protein or nucleic)"``
calpha (ca)     Cα atoms of protein residues, same as ``"name CA and protein"``
backbone (bb)   backbone atoms of protein residues, same as ``"name CA C O N H and protein"``
sidechain (sc)  side-chain atoms of protein residues, same as ``"not name CA C O N H and protein"``
carbon          carbon atoms, same as ``'name "C.*" and not resname ion'``
hydrogen        hydrogen atoms, same as ``'name "[1-9]?H.*"'``
noh             non hydrogen atoms, same as ``'not name "[1-9]?H.*"'``
nitrogen        nitrogen atoms, same as ``'name "N.*"'``
oxygen          oxygen atoms, same as ``'name "O.*"'``
sulfur          sulfur atoms, same as ``'name "S.*"'``
extended        residue in extended conformation, same as ``"secondary E"``
helix           residue in α-helix conformation, same as ``"secondary H"``
helix_3_10      residue in 3_10-helix conformation, same as ``"secondary G"``
helix_pi        residue in π-helix conformation, same as ``"secondary I"``
turn            residue in hydrogen bonded turn conformation, same as ``"secondary T"``
bridge          residue in isolated beta-bridge conformation, same as ``"secondary B"``
bend            residue in bend conformation, same as ``"secondary S"``
coil            residue not in one of above conformations, same as ``"secondary C"``
=============== ===============================================================

.. versionchanged:: 0.6.2
   H is added to the list of backbone atoms.

.. versionadded:: 0.7 
   New keywords are defined: ``lipid, heme, ion, buried, surface, at, 
   cg, purine, pyrimidine, carbon, nitrogen, oxygen, sulfur, extended, helix, 
   helix_pi, helix_3_10, turn, bridge, bend, coil`` 


Among these list of backbone atom names can be changed using 
:func:`setBackboneAtomNames`  and regular expressions for element types
can be changed using :func:`setAtomNameRegex`.

"""

SECONDARY_STRUCTURE_MAP = {
    'extended': 'E',
    'helix': 'H',
    'helix_pi': 'G',
    'helix_3_10': 'I',
    'turn': 'T',
    'bridge': 'B',
    'bend': 'S',
    'coil': 'C',

}
    
KEYWORD_NAME_REGEX = {
    'carbon': 'C.*',
    'hydrogen': '[0-9]?H.*',
    'nitrogen': 'N.*',
    'oxygen': 'O.*',
    'sulfur': 'S.*',
}

BACKBONE_ATOM_NAMES = set(('CA', 'N', 'C', 'O', 'H')) 

KEYWORD_MAP = {}
def _buildKeywordMap():
    global KEYWORD_MAP
    
    protein = KEYWORD_RESNAMES['protein']
    #'keyword' : (residue_names, invert, atom_names, atom_names_not),
    for keyword, resnames in KEYWORD_RESNAMES.iteritems():
        KEYWORD_MAP[keyword] = (resnames, False, None, False)

    KEYWORD_MAP['alpha'] = (protein, False, ['CA'], False)
    KEYWORD_MAP['calpha'] = (protein, False, ['CA'], False)
    KEYWORD_MAP['ca'] = (protein, False, ['CA'], False)
    KEYWORD_MAP['backbone'] = (protein, False, BACKBONE_ATOM_NAMES, False)
    KEYWORD_MAP['bb'] = (protein, False, BACKBONE_ATOM_NAMES, False)
    KEYWORD_MAP['sidechain'] = (protein, False, BACKBONE_ATOM_NAMES, True)
    KEYWORD_MAP['sc'] = (protein, False, BACKBONE_ATOM_NAMES, True)

    KEYWORD_MAP['hetero'] = (protein + KEYWORD_RESNAMES['nucleic'], True, None, False) 

    for name, regex in KEYWORD_NAME_REGEX.iteritems():
        KEYWORD_MAP[name] = (None, False, [['"', regex, '"']], False)
    
    KEYWORD_MAP['carbon'] = (KEYWORD_RESNAMES['ion'], True, [['"', KEYWORD_NAME_REGEX['carbon'], '"']], False)
    KEYWORD_MAP['noh'] = (None, False, [['"', KEYWORD_NAME_REGEX['hydrogen'], '"']], True)
    
_buildKeywordMap()
KEYWORDS_BOOLEAN = set(['all', 'none'] + KEYWORD_MAP.keys() + 
                       SECONDARY_STRUCTURE_MAP.keys())

__doc__ += """

Numerical comparisons
-------------------------------------------------------------------------------

Following keywords can be used in numerical comparisons, as operands of 
arithmetic operations or as arguments to functions:  
 
 * index, serial   
 * resnum, resid
 * x, y, z
 * beta, occupancy
 * charge, mass, radius (these must be set by user before they can be used) 

Numerical attributes of atoms can be used wit the following comparison 

========== =================================
Comparison Description
========== =================================
   <       less than
   >       greater than
   <=      less than or equal
   >=      greater than or equal
   ==      equal
   =       equal
   !=      not equal
========== =================================

Numerical attributes of atoms can be used as operands to the following 
operators:

========= ==================================
Operation Description
========= ==================================
x ** y    x to the power y
x ^ y     x to the power y
x * y     x times y
x / y     x divided by y
x // y    x divided by y (floor division)
x % y     x modulo y
x + y     x plus y 
x - y     x minus y
========= ==================================
   
Numerical attributes of atoms can be used as arguments to the following 
functions:
   
======== ===================================
Function Description
======== ===================================
abs(x)   absolute value of x 
acos(x)  arccos of x
asin(x)  arcsin of x
atan(x)  arctan of x
ceil(x)  smallest integer not less than x
cos(x)   cosine of x"
cosh(x)  hyperbolic cosine of x
floor(x) largest integer not greater than x 
exp(x)   e to the power x
log(x)   natural logarithm of x
log10(x) base 10 logarithm of x
sin(x)   sine of x
sinh(x)  hyperbolic sine of x
sq(x)    square of x
sqrt(x)  square-root of x
tan(x)   tangent of x
tanh(x)  yperbolic tangent of x
======== ===================================

**Examples**
  
  * ``"sqrt(x**2 + y**2 + z**2) < 10"`` selects atoms within 10 Å of the 
    origin
  * ``"resnum <= 100"`` selects atoms with residue numbers less than or equal 
    to 100  


Selection macros
-------------------------------------------------------------------------------

.. versionadded:: 0.7
   Finally, any valid selection string can be used to define selection macros
   using the :func:`defSelectionMacro` function. Macros are saved in ProDy
   configuration and loaded in later sessions automatically. 



:mod:`prody.select`
===============================================================================

Classes
-------------------------------------------------------------------------------

  * :class:`Select`
  * :class:`Contacts`
  
Functions
-------------------------------------------------------------------------------

Below functions can be used to learn and change the definitions of 
some selection keywords:

  * Learn keyword definitions:
    
    * :func:`getAtomNameRegex`
    * :func:`getBackboneAtomNames`
    * :func:`getKeywordResidueNames` 
    * :func:`getProteinResidueNames`
    
  * Change keyword definitions:
    
    * :func:`setAtomNameRegex`
    * :func:`setBackboneAtomNames`
    * :func:`setKeywordResidueNames`
    * :func:`setProteinResidueNames`

Below functions are for manipulating selection macros:
    
  * :func:`defSelectionMacro`
  * :func:`delSelectionMacro`
  * :func:`getSelectionMacro` 
  

"""

FUNCTION_MAP = {
    'sqrt'  : np.sqrt,
    'sq'    : lambda num: np.power(num, 2),
    'abs'   : np.abs,
    'floor' : np.floor,
    'ceil'  : np.ceil,
    'sin'   : np.sin,
    'cos'   : np.cos,
    'tan'   : np.tan,
    'asin'  : np.arcsin,
    'acos'  : np.arccos,
    'atan'  : np.arctan,
    'sinh'  : np.sinh,
    'cosh'  : np.cosh,
    'tahn'  : np.tanh,
    'exp'   : np.exp,
    'log'   : np.log,
    'log10' : np.log10,
}
    
BINARY_OPERATOR_MAP = {
    '+'  : lambda a, b: a + b,
    '-'  : lambda a, b: a - b,
    '*'  : lambda a, b: a * b,
    '/'  : lambda a, b: a / b,
    '%'  : lambda a, b: a % b,
    '>'  : lambda a, b: a > b,
    '<'  : lambda a, b: a < b,
    '>=' : lambda a, b: a >= b,
    '<=' : lambda a, b: a <= b,
    '='  : lambda a, b: a == b,
    '==' : lambda a, b: a == b,
    '!=' : lambda a, b: a != b,
}

COMPARISONS = set(('<', '>', '>=', '<=', '==', '=', '!='))

n_atoms = 10
ATOMGROUP = prody.AtomGroup('Test')
ATOMGROUP.setCoordinates(np.random.random((n_atoms,3)))
ATOMGROUP.setAtomNames(['CA']*n_atoms)
ATOMGROUP.setResidueNames(['GLY']*n_atoms)
ATOMGROUP.setResidueNumbers(np.arange(1,n_atoms+1))
ATOMGROUP.setChainIdentifiers(['A']*n_atoms)
ATOMGROUP.setAltLocIndicators([' ']*n_atoms)
ATOMGROUP.setElementSymbols(['C']*n_atoms)
ATOMGROUP.setHeteroFlags([False]*n_atoms)
ATOMGROUP.setOccupancies([1]*n_atoms)
ATOMGROUP.setSecondaryStrs(['H']*n_atoms)
ATOMGROUP.setSegmentNames(['PDB']*n_atoms)
ATOMGROUP.setAnisoTempFactors(np.random.random((n_atoms,6)))
ATOMGROUP.setAnisoStdDevs(np.random.random((n_atoms,6)))
ATOMGROUP.setInsertionCodes([' ']*n_atoms)
ATOMGROUP.setAtomTypes(['CH2']*n_atoms)
ATOMGROUP.setTempFactors([0]*n_atoms)
ATOMGROUP.setCharges([0]*n_atoms)
ATOMGROUP.setMasses([12]*n_atoms)
ATOMGROUP.setRadii([1.4]*n_atoms)

try:
    MACROS = prody._ProDySettings['selection_macros']
except KeyError:
    MACROS = {}

def defSelectionMacro(name, selstr):
    """Define selection macro *selstr* with name *name*.
    
    .. versionadded:: 0.7
    
    Both *name* and *selstr* must be string. An existing keyword cannot be 
    used as a macro name. If a macro with given *name* exists, it will be 
    overwritten.
    
    >>> defSelectionMacro('cbeta', 'name CB and protein')
    
    """
    
    if not isinstance(name, str) or not isinstance(selstr, str):
        raise TypeError('both name and selstr must be strings')
    elif isKeyword(name):
        raise ValueError('"{0:s}" is an existing keyword, cannot be used as a '
                         'macro name'.format(name))
    elif not (name.isalpha() and name.islower()):
        raise ValueError('macro names must be all lower case letters, {0:s} '
                         'is not a valid macro name'.format(name))
    
    LOGGER.info('Testing validity of selection string:')
    try:
        ATOMGROUP.select(selstr)
    except:
        LOGGER.warning('"{0:s}" is not a valid selection string, macro "{1:s}"'
                    ' is not defined.'.format(selstr, name))
    else:
        LOGGER.info('Macro "{0:s}" is defined as "{1:s}".'
                    .format(name, selstr))
        MACROS[name] = selstr
        prody._ProDySettings['selection_macros'] = MACROS
        prody._saveProDySettings()

def delSelectionMacro(name):
    """Delete the macro *name*.
    
    .. versionadded:: 0.7
    
    >>> delSelectionMacro('cbeta')
    
    """
    
    try:
        MACROS.pop(name)
    except:
        LOGGER.warning('Macro "{0:s}" is not found.'.format(name))
    else:
        LOGGER.info('Macro "{0:s}" is deleted.'.format(name))
        prody._ProDySettings['selection_macros'] = MACROS
        prody._saveProDySettings()

def getSelectionMacro(name=None):
    """Return the definition of the macro *name*. 
        
    .. versionadded:: 0.7
    
    If *name* is not given, returns a copy of the selection macros dictionary.
    
    """
    
    if name is None:        
        return MACROS.copy()
    try:
        return MACROS[name]
    except KeyError:
        LOGGER.info('"{0:s}" is not a user defined macro name.'.format(name))

mapField2Var = {}
for field in ATOMIC_DATA_FIELDS.values():
    mapField2Var[field.name] = field.var

def getKeywordResidueNames(keyword):
    """Return residue names associated with a keyword.
    
    .. versionadded:: 0.7
    
    >>> getKeywordResidueNames('acidic')
    ['ASP', 'GLU']
    
    """
    
    try:
        resnames = KEYWORD_RESNAMES[keyword]
        resnames.sort()
        return resnames  
    except KeyError:
        if keyword in KEYWORD_RESNAMES_READONLY:
            LOGGER.warning('{0:s} is defined as "{1:s}"'.format(keyword, 
                                        KEYWORD_RESNAMES_READONLY[keyword]))
        else:
            LOGGER.warning('{0:s} is not a keyword'.format(keyword))

def setKeywordResidueNames(keyword, resnames):
    """Change the list of residue names associated with a keyword.
    
    .. versionadded:: 0.7
    
    *keyword* must be a string, and *resnames* may be a list, tuple, or set of
    strings. The existing list of residue names will be overwritten with the
    given residue names. Note that changes in keyword definitions are not 
    saved permanently.
    
    >>> setKeywordResidueNames('acidic', ['ASP', 'GLU'])
    
    """
    
    if not isinstance(keyword, str):
        raise TypeError('keyword must be a string')
    if not isinstance(resnames, (list, tuple, set)):
        raise TypeError('resnames must be a list, set, or tuple')
    if keyword in KEYWORD_RESNAMES_READONLY:
        LOGGER.warning('{0:s} is defined as "{1:s}" and cannot be changed '
                           'directly'.format(keyword, 
                                        KEYWORD_RESNAMES_READONLY[keyword]))
        return
    if keyword in KEYWORD_RESNAMES:
        for rn in resnames:
            if not isinstance(rn, str):
                raise TypeError('all items in resnames must be strings')
        KEYWORD_RESNAMES[keyword] = list(set(resnames))
        _setReadonlyResidueNames()
    else:
        raise ValueError('{0:s} is not a valid keyword'.format(keyword))

def getAtomNameRegex(name):
    """Return regular expression used for selecting common elements.
    
    .. versionadded:: 0.7
    
    >>> getAtomNameRegex('nitrogen')
    'N.*'
    
    """
    
    try:
        return KEYWORD_NAME_REGEX[name]   
    except KeyError:
        LOGGER.warning('{0:s} is not a valid element'.format(keyword))

def setAtomNameRegex(name, regex):
    """Set regular expression used for selecting common elements.
    
    Note that changes in keyword definitions are not saved permanently.
    
    .. versionadded:: 0.7
    
    >>> setAtomNameRegex('nitrogen', 'N.*')
    
    """
    
    if not name in KEYWORD_NAME_REGEX:
        raise ValueError('{0:s} is not a valid keyword'.format(name))
    try:
        RE.compile(regex)
    except:
        raise ValueError('{0:s} is not a valid regular expression'
                         .format(regex))
    else:
        KEYWORD_NAME_REGEX[name] = regex

def getBackboneAtomNames():
    """Return protein backbone atom names.
    
    >>> getBackboneAtomNames()
    ['C', 'CA', 'H', 'N', 'O']
    
    """
    
    bban = list(BACKBONE_ATOM_NAMES)
    bban.sort()
    return bban 

def setBackboneAtomNames(backbone_atom_names):
    """Set protein backbone atom names.
    
    Note that changes in keyword definitions are not saved permanently.
    
    """
    
    if not isinstance(backbone_atom_names, (list, tuple, set)):
        raise TypeError('backbone_atom_names must be a list, tuple, or set')
    global BACKBONE_ATOM_NAMES
    BACKBONE_ATOM_NAMES = set(backbone_atom_names)
    _buildKeywordMap()

def getProteinResidueNames():
    """Return list of protein residue names."""
    
    return KEYWORD_RESNAMES['protein']

def setProteinResidueNames(resnames):
    """Set list of protein residue names.
    
    Note that changes in keyword definitions are not saved permanently.
    
    """

    setKeywordResidueNames('protein', resnames)

class SelectionError(Exception):    
    pass


def isFloatKeyword(keyword):
    return keyword in KEYWORDS_FLOAT

def isNumericKeyword(keyword):
    return keyword in KEYWORDS_NUMERIC

def isAlnumKeyword(keyword):
    return keyword in KEYWORDS_STRING

def isValuePairedKeyword(keyword):
    return keyword in KEYWORDS_VALUE_PAIRED

def isBooleanKeyword(keyword):
    return keyword in KEYWORDS_BOOLEAN
    
def isKeyword(keyword):
    return (isBooleanKeyword(keyword) or isValuePairedKeyword(keyword) or
            isNumericKeyword(keyword))

def isComparison(tokens):
    for tkn in tokens:
        if tkn in COMPARISONS:
            return True
    return False

_specialKeywords = set(['secondary', 'chain', 'altloc', 'segment'])

class Select(object):

    """Select subsets of atoms based on a selection string.
    
    Definitions of single word keywords, such as protein, 
    backbone, polar, etc., may be altered using functions in 
    :mod:`~prody.select` module. 
    
    This class makes use of |pyparsing| module.

    """

    def __init__(self):
        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._selstr = None
        self._evalonly = None
        self._acsi = None
        self._timestamp = None
        
        self._coordinates = None
        self._kdtree = None
        self._kwargs  = None
        self._selstr2indices = False
        self._data = dict()
        for var in mapField2Var.values():
            self._data[var] = None        
        
        shortlist = pp.alphanums + '''~@#$.:;_','''
        longlist = pp.alphanums + '''~!@#$%^&*()-_=+[{}]\|;:,<>./?()' '''
        
        self._tokenizer = pp.operatorPrecedence(
             pp.OneOrMore(pp.Word(shortlist) | 
             pp.Group(pp.Literal('"') + pp.Word(longlist + '`') + 
                      pp.Literal('"')) | 
             pp.Group(pp.Literal('`') + pp.Word(longlist + '"') + 
                      pp.Literal('`'))
                      ),
             [(pp.oneOf('sqrt sq abs floor ceil sin cos tan atan '
                        'asin acos sinh cosh tanh exp log log10'), 
                        1, pp.opAssoc.RIGHT, self._func),
              (pp.oneOf('** ^'), 2, pp.opAssoc.LEFT, self._pow),
              (pp.oneOf('+ -'), 1, pp.opAssoc.RIGHT, self._sign),
              (pp.oneOf('* / %'), 2, pp.opAssoc.LEFT, self._mul),
              (pp.oneOf('+ -'), 2, pp.opAssoc.LEFT, self._add),
              (pp.oneOf('< > <= >= == = !='), 2, pp.opAssoc.LEFT, self._comp),
              (pp.Keyword('!!!') | 
               pp.Regex('same [a-z]+ as') | 
               pp.Regex('(ex)?within [0-9]+\.?[0-9]* of'), 
                        1, pp.opAssoc.RIGHT, self._special),
              (pp.Keyword('&&&'), 2, pp.opAssoc.LEFT, self._and),
              (pp.Keyword('||'), 2, pp.opAssoc.LEFT, self._or),]
            )

        self._tokenizer.setParseAction(self._defaultAction)

    def getBoolArray(self, atoms, selstr, **kwargs):
        """Return a boolean array with ``True`` values for *atoms* matching 
        *selstr*.
        
        .. versionadded:: 0.5
        
        .. note:: The length of the boolean :class:`numpy.ndarray` will be
           equal to the number of atoms in *atoms* argument.
            
        """
        
        if not isinstance(atoms, prody.Atomic):
            raise TypeError('atoms must be an Atomic instance, not {0:s}'
                            .format(type(atoms)))
        elif not isinstance(selstr, str):
            raise TypeError('selstr must be a string, not a {0:s}'
                            .format(type(selstr)))
        if self._atoms is atoms:
            if DEBUG: print('atoms is the same')
            if self._acsi != atoms.getActiveCoordsetIndex():
                self._coordinates = None
                self._kdtree = None
            elif self._timestamp != atoms._getTimeStamp(self._acsi):
                self._kdtree = None
        else:
            self._reset()
            if isinstance(atoms, prody.AtomGroup): 
                self._ag = atoms
                self._atoms = atoms
                self._indices = None
                self._n_atoms = atoms.getNumOfAtoms()
            else:
                self._ag = atoms.getAtomGroup()
                self._indices = atoms.getIndices()
                if isinstance(atoms, prody.AtomMap):
                    self._atoms = prody.Selection(self._ag, self._indices, '')
                    self._atoms._indices = self._indices
                else: 
                    self._atoms = atoms
                self._n_atoms = len(self._indices)
        self._selstr = selstr
        self._acsi = atoms.getActiveCoordsetIndex()
        self._timestamp = atoms._getTimeStamp(self._acsi)
            
        self._kwargs = kwargs
        if DEBUG:
            print('getBoolArray', selstr)
        torf = self._evalSelstr()
        if not isinstance(torf, np.ndarray):
            raise SelectionError('{0:s} is not a valid selection string.'
                                 .format(selstr))
        elif torf.dtype != np.bool:
            if DEBUG:
                print('_select torf.dtype', torf.dtype, isinstance(torf.dtype, 
                                                                   np.bool))
            raise SelectionError('{0:s} is not a valid selection string.'
                                 .format(selstr))
        if DEBUG:
            print('_select', torf)
        return torf
    
    def getIndices(self, atoms, selstr, **kwargs):
        """Return indices of atoms matching *selstr*.
        
        .. versionadded:: 0.5
        
        """
        
        torf = self.getBoolArray(atoms, selstr, **kwargs)        
        return torf.nonzero()[0]
        
    def select(self, atoms, selstr, **kwargs):
        """Return a subset of atoms matching *selstr* as a :class:`Selection`.
        
        :arg atoms: atoms to select from which    
        :type atoms: :class:`~prody.atomic.Atomic`
        
        :arg selstr: selection string
        :type selstr: str
        
        :keyword cache: cache atomic data and KDTree, default is ``True``
        :type cache: bool
        
        If type of *atoms* is :class:`~prody.atomic.AtomMap`, an 
        :class:`~prody.atomic.AtomMap` instance is returned. Otherwise,
        :class:`~prody.atomic.Selection` instances are returned.

        .. note:

            * If selection string does not match any atoms, ``None`` is 
              returned.
              
            * :meth:`select` accepts arbitrary keyword arguments which enables 
              identification of intermolecular contacts. See :ref:`contacts` 
              for details.
        
            * :meth:`select` accepts a keyword argument that enables caching
              atomic data and KDTree from previous select operation. It works
              if *atoms* objects in two consecutive selections are the same.
        
            * A special case for making atom selections is passing an
              :class:`~prody.atomic.AtomMap` instance as *atoms* argument. 
              Unmapped atoms will not be included in the returned 
              :class:`~prody.atomic.AtomMap` instance. The order of atoms 
              will be preserved.

        .. warning:: ``cache=True`` should be used if attributes of *atoms* 
           object have not changed since the previous selection.
           
        """
        
        self._selstr2indices = False
        indices = self.getIndices(atoms, selstr, **kwargs)
        if not isinstance(atoms, prody.AtomGroup):
            indices = self._indices[indices]
        ag = self._ag
        if not kwargs.get('cache', True):
            self._reset()
        self._kwargs = None
        if len(indices) == 0:
            return None
        elif isinstance(atoms, prody.AtomMap):
            return prody.AtomMap(ag, indices, np.arange(len(indices)), 
                                 np.array([]),
                                 'Selection "{0:s}" from AtomMap {1:s}'.format(
                                 selstr, atoms.getName()),
                                 atoms.getActiveCoordsetIndex())
        else:
            
            if self._selstr2indices:
                selstr = 'index {0:s}'.format(prody.rangeString(indices))
            elif isinstance(atoms, prody.AtomPointer):
                selstr = '({0:s}) and ({1:s})'.format(selstr, 
                                                    atoms.getSelectionString())
            
            return prody.Selection(ag, indices, selstr, 
                                                atoms.getActiveCoordsetIndex())
        
    def _reset(self):
        if DEBUG: print('_reset')
        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._evalonly = None
        self._acsi = None
        self._timestamp = None
        self._coordinates = None
        self._kdtree = None
        for var in mapField2Var.values():
            self._data[var] = None        

    def _standardizeSelectionString(self):
        selstr = ' ' + self._selstr + ' '

        selstr = selstr.replace(')and(', ')&&&(')
        selstr = selstr.replace(' and(', ' &&&(')
        selstr = selstr.replace(')and ', ')&&& ')
        selstr = selstr.replace(' and ', ' &&& ')
            
        selstr = selstr.replace(')or(', ')||(')
        selstr = selstr.replace(' or(', ' ||(')
        selstr = selstr.replace(')or ', ')|| ')
        selstr = selstr.replace(' or ', ' || ')
        
        selstr = selstr.replace('(not ', '(!!! ')
        selstr = selstr.replace(' not(', ' !!!(')
        selstr = selstr.replace(' not ', ' !!! ')
        
        if MACROS:
            for macro in MACROS.iterkeys():
                selstr = selstr.replace(macro, '(' + MACROS[macro] + ')')
        
        return selstr.strip()

    def _evalSelstr(self):
        selstr = self._selstr.strip() 
        if len(selstr.split()) == 1 and '(' not in selstr and \
           ')' not in selstr and selstr not in MACROS:
            return self._evalBoolean(selstr)
        selstr = self._standardizeSelectionString()
        if DEBUG: print('_evalSelstr', selstr)

        try: 
            tokens = self._tokenizer.parseString(selstr, 
                                                 parseAll=True).asList()
            if DEBUG: print('_evalSelstr', tokens)
            return tokens[0]
        except pp.ParseException, err:
            print('Parse Failure')
            print(self._selstr) #err.line
            print(" "*(err.column-1) + "^")
            raise pp.ParseException(str(err))
    
    def _defaultAction(self, token):
        if DEBUG: print('_defaultAction', token)
        if isinstance(token[0], (np.ndarray, float)):
            return token[0]
        else:
            return self._evaluate(token)        
    
    def _evaluate(self, token):
        if DEBUG: print('_evaluate', token)

        keyword = token[0]
        if len(token) == 1:
            if isBooleanKeyword(keyword):
                return self._evalBoolean(keyword)
            elif isNumericKeyword(keyword):
                return self._getNumArray(keyword)
            elif self._kwargs is not None and keyword in self._kwargs:
                return keyword
            elif self._ag.isAttribute():
                return self._evalAttribute(keyword)
            else:
                try:
                    return float(keyword)
                except ValueError:
                    raise SelectionError('"{0:s}" is not a valid keyword or a '
                                         'number.'.format(keyword))
        elif isAlnumKeyword(keyword):
            return self._evalAlnum(keyword, token[1:])
        elif isFloatKeyword(keyword):
            return self._evalFloat(keyword, token[1:])
        elif keyword in ('resnum', 'resid'):
            return self._resnum(token[1:])
        elif keyword == 'index':
            return self._index(token[1:])
        elif keyword == 'serial':
            return self._serial(token[1:])
        elif keyword == 'within':
            return self._within([' '.join(token[:3])] + token[3:], False)
        elif keyword == 'exwithin':
            return self._within([' '.join(token[:3])] + token[3:], True)
        elif keyword == 'same':
            return self._sameas([' '.join(token[:3])] + token[3:])
        elif keyword == '!!!':
            return self._not(token)
        elif self._ag.isAttribute(keyword):
            return self._evalAttribute(keyword, token[1:])
        elif isBooleanKeyword(keyword):
            raise SelectionError('Single word keywords must be followed with '
                                 'and operator.')            
            return self._and([token])
        raise SelectionError('{0:s} is not understood. Please report this if '
                             'you think there is a bug.'
                             .format(' '.join(token)))

    def _or(self, tokens):
        if DEBUG: print('_or', tokens)
        temp = tokens[0]
        tokenlist = []
        token = []
        while temp:
            tkn = temp.pop(0)
            if isinstance(tkn, str) and isBooleanKeyword(tkn):
                tkn = self._evalBoolean(tkn)
            if tkn == '||':
                tokenlist.append(token)
                token = []
            else:
                token.append(tkn)
        tokenlist.append(token)

        if DEBUG: print('_or tokenlist', tokenlist)

        for token in tokenlist:
            zero = token[0]
            if isinstance(zero, np.ndarray):                    
                if self._evalonly is None: 
                    self._evalonly = np.invert(zero).nonzero()[0]
                else:        
                    self._evalonly = self._evalonly[np.invert(zero[
                                                self._evalonly]).nonzero()[0]]
            else:
                torf = self._evaluate(token)
                if self._evalonly is None:
                    self._evalonly = np.invert(torf).nonzero()[0]
                else:
                    self._evalonly = self._evalonly[np.invert(torf)]
            if DEBUG: print('_or evalonly', self._evalonly)
        torf = np.ones(self._n_atoms, np.bool)
        torf[self._evalonly] = False
        self._evalonly = None
        return torf

    def _and(self, tokens):
        if DEBUG: print('_and', tokens)
        temp = tokens[0]
        tokenlist = []
        token = []
        while temp:
            tkn = temp.pop(0)
            if isinstance(tkn, str) and isBooleanKeyword(tkn):
                tkn = self._evalBoolean(tkn, True)
                if isinstance(tkn, list):
                    tkn.extend(temp)
                    temp = tkn
                    continue
            if tkn == '&&&':
                tokenlist.append(token)
                token = []
            else:
                token.append(tkn)
        tokenlist.append(token)
        if DEBUG: print('_and tokenlist', tokenlist)
        for token in tokenlist:
            zero = token[0]
            if isinstance(zero, np.ndarray):                    
                if self._evalonly is None: 
                    self._evalonly = zero.nonzero()[0]
                else:        
                    self._evalonly = self._evalonly[zero[self._evalonly
                                                                ].nonzero()[0]]
            else:
                torf = self._evaluate(token)
                if self._evalonly is None:
                    self._evalonly = torf.nonzero()[0]
                else:
                    self._evalonly = self._evalonly[torf]
            if DEBUG: print('_and evalonly', self._evalonly)
        torf = np.zeros(self._n_atoms, np.bool)
        torf[self._evalonly] = True
        self._evalonly = None
        return torf
    
    def _special(self, token):
        if DEBUG: print('_special', token)
        token = token[0]
        if token[0] == '!!!':
            return self._not(token)
        elif token[0].startswith('same'):
            return self._sameas(token)
        else:
            return self._within(token, token[0].startswith('exwithin'))

    def _not(self, token):
        if DEBUG: print('_not', token)
        if isinstance(token[1], np.ndarray):
            torf = token[1]
        else:
            torf = self._evaluate(token[1:])
        np.invert(torf, torf)
        return torf
    
    def _within(self, token, exclude):
        terms = token
        if DEBUG: print('_within', terms)
        within = float(terms[0].split()[1])
        which = terms[1]
        if not isinstance(which, np.ndarray):
            which = self._evaluate(terms[1:])
        result = []
        append = result.append
        kdtree = self._getKDTree()
        get_indices = kdtree.get_indices
        search = kdtree.search
        if isinstance(which, np.ndarray):
            which = which.nonzero()[0]
            if len(which) == 0:
                return np.zeros(self._n_atoms, np.bool)
            coordinates = self._getCoordinates()
            for index in which:
                search(coordinates[index], within)
                append(get_indices())
        else:
            if self._kwargs is not None and which in self._kwargs:
                if DEBUG: print('_kwargs', which)
                which = self._kwargs[which]
                exclude=False
                self._selstr2indices = True
            if isinstance(which, np.ndarray):
                if which.ndim == 1 and len(which) == 3:
                    which = [which]
                elif not (which.ndim == 2 and which.shape[1] == 3):
                    raise SelectionError('{0:s} must be a coordinate array, '
                                         'shape (N, 3) or (3,)'.format(kw))
                for xyz in which:
                    if DEBUG: print('xyz', xyz)
                    search(xyz, within)
                    append(get_indices())
            else:
                try:
                    coordinates = which.getCoordinates()
                except:
                    raise SelectionError('{0:s} must have a getCoordinates() '
                                         'method.'.format(kw))
                if not isinstance(coordinates, np.ndarray):
                    raise SelectionError('{0:s}.getCoordinates() method must '
                                         'return a numpy.ndarray instance.'
                                         .format(kw))
                for xyz in coordinates:
                    search(xyz, within)
                    append(get_indices())
                
        unique = np.unique(np.concatenate(result))
        if len(unique) == 0:
            return np.zeros(self._n_atoms, np.bool)
        if self._indices is None:
            torf = np.zeros(self._n_atoms, np.bool)
            torf[unique] = True
        else:
            torf = np.zeros(self._ag._n_atoms, np.bool)
            torf[unique] = True
            torf = torf[self._indices]
        if exclude:
            torf[which] = False
        if self._evalonly is None:
            return torf
        else:
            return torf[self._evalonly]
    
    def _sameas(self, token):
        terms = token
        if DEBUG: print('_sameas', terms)
        what = token[0].split()[1]
        which = token[1]
        if not isinstance(which, np.ndarray):
            which = self._evaluate(token[1:])
        
        if what == 'residue':
            chainids = self._getAtomicData('chain')
            resids =  self._getAtomicData('resnum')
            resnum = np.unique(resids[which]).astype('|S6')
            torf = np.all(
                [self._evalAlnum('chain', list(np.unique(chainids[which]))),
                 self._resnum(resnum)], 0)
        elif what == 'chain':
            chainids = self._getAtomicData('chain')
            torf = self._evalAlnum('chain', list(np.unique(chainids[which])))        
        elif what == 'segment':
            segnames = self._getAtomicData('segment')
            torf = self._evalAlnum('segment', list(np.unique(segnames[which]))) 
        return torf
     
    def _comp(self, token):
        if DEBUG: print('_comp', token)
        token = token[0]
        if len(token) > 3:
            if isBooleanKeyword(token[0]):
                return self._and([[token.pop(0), '&&&', self._comp([token])]])
            elif isBooleanKeyword(token[-1]):
                return self._and([[token.pop(-1), '&&&', self._comp([token])]])
            else:
                raise SelectionError('{0:s} is not valid'
                                     .format(' '.join(token)))
        comp = token[1]
        left = self._getNumArray(token[0])
        if DEBUG: print('_comp left', left)
        right = self._getNumArray(token[2])
        if DEBUG: print('_comp right', right)

        try:
            return BINARY_OPERATOR_MAP[comp](left, right)
        except KeyError:
            raise SelectionError('Unknown error in "{0:s}".'
                                 .format(' '.join(token)))

    def _pow(self, token):
        if DEBUG: print('_pow', token)
        items = token[0]
        return self._getNumArray(items[0]) ** self._getNumArray(items[2])

    def _add(self, token):
        if DEBUG: print('_add', token)
        items = token[0]
        left = self._getNumArray(items.pop(0))
        while items:
            left = BINARY_OPERATOR_MAP[items.pop(0)](
                                        left, self._getNumArray(items.pop(0)))
        if DEBUG: print('_add total', left)
        return left
 
    def _mul(self, token):
        if DEBUG: print('_mul', token)
        items = token[0]
        left = self._getNumArray(items[0])
        i = 1
        while i < len(items):
            op = items[i]
            i += 1
            right = self._getNumArray(items[i])
            i += 1
            if op == '/' and right == 0.0: 
                raise ZeroDivisionError(' '.join(items))
            left = BINARY_OPERATOR_MAP[op](left, right)
        return left
    
    def _getNumArray(self, token):
        if DEBUG: print('_getNumArray', token)
        if isinstance(token, (np.ndarray, float)):
            return token
        elif isFloatKeyword(token):
            return self._evalFloat(token)
        elif token in ('resnum', 'resid'):
            return self._resnum()
        elif token == 'index':
            return self._index()    
        elif token == 'serial':
            return self._serial()
        elif self._ag.isAttribute(token):
            data = self._getAtomicData(token)
            if data.dtype.type in (np.float64, np.int64):
                return data
            else:
                raise SelectionError('attribute "{0:s}" is not a numeric type'
                                     .format(token))
        else:
            try:
                token = float(token)
            except ValueError:
                raise SelectionError('"{0:s}" must be a number or a valid '
                                     'numeric keyword'.format(token))
            else:
                return token

    def _sign(self, tokens):
        if DEBUG: print('_sign', tokens)
        tokens = tokens[0]
        token = self._getNumArray(tokens[1])
        if tokens[0] == '-':
            return -token
        return token

    def _func(self, token):
        if DEBUG: print('_func', token)
        token = token[0]
        return FUNCTION_MAP[token[0]](token[1])

    def _evalAttribute(self, keyword, values=None):
        if values is None:
            torf = self._atoms.getAttribute(keyword)
            if isinstance(torf.dtype, np.bool):
                return torf
            else:
                raise SelectionError('attribute {0:s} is not boolean'
                                     .format(keyword))
        else:
            data = self._getAtomicData(keyword)
            if data.dtype.type in (np.int64, np.float64):
                return self._evalFloat(keyword, values)
            elif data.dtype.type == np.string_:
                return self._evalAlnum(keyword, values)
            else:
                raise SelectionError('type of attribute {0:s} is not valid'
                                     .format(keyword))

    def _evalBoolean(self, keyword, _and=False):
        if DEBUG: print('_evalBoolean', keyword)
        
        if self._evalonly is None:
            n_atoms = self._n_atoms
        else:        
            n_atoms = len(self._evalonly)
        
        if keyword == 'all':
            return np.ones(n_atoms, np.bool)
        elif keyword == 'none':
            return np.zeros(n_atoms, np.bool)
        elif keyword in SECONDARY_STRUCTURE_MAP:
            if _and:
                return ['secondary', SECONDARY_STRUCTURE_MAP[keyword]]
            else:
                return self._evalAlnum('secondary', 
                                       [SECONDARY_STRUCTURE_MAP[keyword]])
        try:
            (residue_names, rn_invert, atom_names, 
                                            an_invert) = KEYWORD_MAP[keyword]
            if DEBUG:
                print('_evalBoolean', residue_names, rn_invert, atom_names, 
                      an_invert)
        except KeyError:
            raise SelectionError('"{0:s}" is not a valid keyword.'
                                 .format(keyword))

        compound = []
        if atom_names is not None:
            if an_invert:
                compound.append('!!!')
            compound.append('name')
            compound.extend(atom_names)
            if residue_names is not None:
                compound.append('&&&')
        if residue_names is not None:
            if rn_invert:
                compound.append('!!!')
            compound.append('resname')
            compound.extend(residue_names)
        if DEBUG:
            print('_evalBoolean compound', compound)
        if _and:
            return compound
        else:
            return self._and([compound])
    
    def _evalAlnum(self, keyword, values):
        if DEBUG: print('_evalAlnum', keyword, values)
        data = self._getAtomicData(keyword)
        if keyword in _specialKeywords:
            for i, value in enumerate(values):
                if value == '_':
                    values[i] = ' '
                    values.append('')
                    break
            
        if self._evalonly is not None:
            data = data[self._evalonly]
        n_atoms = len(data)
        
        regexps = []
        strings = []
        for value in values:
            if isinstance(value, str):
                strings.append(value)
            elif value[0] == value[2] == '"':
                regexps.append(value[1])
            else:
                strings.append(value[1])
                
        if len(strings) == 1:
            torf = data == strings[0]
        elif len(strings) > 4:
            torf = np.zeros(n_atoms, np.bool)
            strings = set(strings)
            for i, datum in enumerate(data):        
                torf[i] = datum in strings
        elif strings: 
            torf = [(data == value).reshape((n_atoms, 1)) for value in strings]
            torf = np.concatenate(torf, 1).sum(1).astype(np.bool) 
        else:
            torf = np.zeros(n_atoms, np.bool)

        for value in regexps:
            value = RE.compile('^' + value + '$')
            for i in xrange(n_atoms):
                torf[i] = (value.match(data[i]) is not None)

        return torf
    
    def _evalFloat(self, keyword, values=None):
        if DEBUG: print('_evalFloat', keyword, values)
        if keyword == 'x':
            data = self._getCoordinates()[:,0]
        elif keyword == 'y':
            data = self._getCoordinates()[:,1]
        elif keyword == 'z':
            data = self._getCoordinates()[:,2]
        else:
            data = self._getAtomicData(keyword)
        
        if values is None:
            return data
    
        if self._evalonly is not None:
            data = data[self._evalonly]
        n_atoms = len(data)
        torf = np.zeros(n_atoms, np.bool)

        for item in self._getNumRange(values):
            if isinstance(item, str):
                pass
            elif isinstance(item, list):
                torf[(item[0] <= data) * (data <= item[1])] = True
            elif isinstance(item, tuple):
                if len(item) == 2:
                    torf[(item[0] <= data) * (data < item[1])] = True
                else:
                    raise SelectionError('"{0:s}" is not a valid range for '
                                         'keywords expecting floating values, '
                                         'such as {1:s}.'.format(':'.join(
                                         [str(i) for i in item]), keyword))
            else:
                torf[data == item] = True
        return torf

    def _resnum(self, token=None):
        if DEBUG: print('_resnum', token)
        if token is None:
            return self._getAtomicData('resnum') 
        icodes = None
        if self._evalonly is None:
            resids = self._getAtomicData('resnum')
            n_atoms = self._n_atoms
        else:
            evalonly = self._evalonly
            resids = self._getAtomicData('resnum')[evalonly]
            n_atoms = len(evalonly)
        torf = np.zeros(n_atoms, np.bool)
        
        for item in self._getNumRange(token):
            if isinstance(item, str):
                if icodes is None:
                    if self._evalonly is None:
                        icodes = self._getAtomicData('icode')
                    else:
                        icodes = self._getAtomicData('icode')[evalonly]
                icode = str(item[-1])
                if icode == '_':
                    icode = ''
                torf[(resids == int(item[:-1])) * (icodes == icode)] = True
            elif isinstance(item, list):
                torf[(item[0] <= resids) * (resids <= item[1])] = True
            elif isinstance(item, tuple):
                if len(item) == 2:
                    torf[(item[0] <= resids) * (resids < item[1])] = True
                else:
                    for i in range(item[0], item[1], item[2]):
                        torf[resids == i] = True
            else:
                torf[resids == item] = True
        return torf

    def _serial(self, token=None):
        if DEBUG: print('_serial', token)
        if token is None:
            return self._getAtomicData('serial') 
        if self._evalonly is None:
            serials = self._getAtomicData('serial')
            n_atoms = self._n_atoms
        else:
            evalonly = self._evalonly
            serials = self._getAtomicData('serial')[evalonly]
            n_atoms = len(evalonly)
        torf = np.zeros(n_atoms, np.bool)
        
        for item in self._getNumRange(token):
            if isinstance(item, list):
                torf[(item[0] <= serials) * (serials <= item[1])] = True
            elif isinstance(item, tuple):
                if len(item) == 2:
                    torf[(item[0] <= serials) * (serials < item[1])] = True
                else:
                    for i in range(item[0], item[1], item[2]):
                        torf[serials == i] = True
            else:
                torf[serials == item] = True
        return torf
    
    def _index(self, token=None):
        if DEBUG: print('_index', token)
        if token is None:
            if self._indices is not None:
                return self._indices
            else:
                return np.arange(self._ag._n_atoms)
        torf = np.zeros(self._ag._n_atoms, np.bool)
        
        for item in self._getNumRange(token):
            if isinstance(item, str):
                raise SelectionError('"index/serial {0:s}" is not understood.'
                                     .format(item))
            elif isinstance(item, tuple):
                if len(item) == 2:
                    torf[item[0]:item[1]] = True
                else:
                    torf[item[0]:item[1]:item[2]] = True
            elif isinstance(item, list):
                torf[int(np.ceil(item[0])):int(
                                    np.floor(item[1]))+1] = True
            else:
                try:
                    torf[item] = True
                except IndexError:
                    pass
        if self._evalonly is not None:
            return torf[self._evalonly]
        if self._indices is not None:
            return torf[self._indices]
        if DEBUG: print('_index return ', len(torf), torf)
        return torf

    def _getNumRange(self, token):
        tknstr = ' '.join(token)
        while '  ' in tknstr:
            tknstr = tknstr.replace('  ', ' ')
        tknstr = tknstr.replace(' to ', 'to').replace(
                                            'to ', 'to').replace(' to', 'to')
        tknstr = tknstr.replace(' : ', ':').replace(
                                            ': ', ':').replace(' :', ':')
        token = []
        for item in tknstr.split():
            if 'to' in item:
                items = item.split('to')
                if len(items) != 2:
                    raise SelectionError('"{0:s}" is not understood.'
                                         .format(' to '.join(items)))
                try:
                    token.append( [float(items[0]), float(items[1])] )
                except:
                    raise SelectionError('"{0:s}" is not understood, "to" '
                                         'must be surrounded by numbers.'
                                         .format(' to '.join(items)))
            elif ':' in item:
                items = item.split(':')
                if not len(items) in (2, 3):
                    raise SelectionError('"{0:s}" is not understood.'
                                         .format(':'.join(items)))
                try:
                    if len(items) == 2:
                        token.append( (int(items[0]), int(items[1])) )
                    else:
                        token.append( (int(items[0]), int(items[1]), 
                                       int(items[2])) )
                except:
                    raise SelectionError('"{0:s}" is not understood, ":" must '
                                         'be surrounded by integers.'
                                         .format(':'.join(items)))
            elif '.' in item:
                try:
                    token.append( float(item) )
                except:
                    raise SelectionError('"{0:s}" is not understood.'
                                         .format(item))
            elif item.isdigit():
                try:
                    token.append( int(item) )
                except:
                    raise SelectionError('"{0:s}" is not understood.'
                                         .format(item))
            else:
                token.append( item )
        if DEBUG: print('_getNumRange', token)            
        return token
    
    def _getAtomicData(self, keyword):
        field = ATOMIC_DATA_FIELDS.get(keyword, None)
        indices = self._indices
        if field is None:
            data = self._atoms.getAttribute(keyword)
            if data is None:
                raise SelectionError('"{0:s}" is not a valid keyword or '
                                     'attribute.'.format(keyword))
            elif data.ndim == 1:
                if indices is None:                
                    return data
                else:
                    return data[indices]
            else:
                raise SelectionError('attribute "{0:s}" is 1-dimensional data '
                                     .format(keyword))                
        else:
            var = field.var
            data = self._data[var]
            if data is None:
                data = self._ag._data[var] 
                if data is None:
                    raise SelectionError('{0:s} are not set.'
                                         .format(field.doc_pl))
                self._data[var] = data
            if indices is None:                
                return data
            else:
                return data[indices]
    
    def _getCoordinates(self):
        if self._coordinates is None:
            if self._indices is None:
                self._coordinates = self._ag._coordinates[self._ag._acsi]
            else:
                self._coordinates = self._atoms._getCoordinates()
        return self._coordinates
   
    def _getKDTree(self):
        if KDTree is None: prody.importBioKDTree()
        if not KDTree:
            raise ImportError('Bio.KDTree is required for distance based '
                              'selections.')
        if self._kdtree is None:
            if DEBUG: print('kdtree')
            kdtree = KDTree(3)
            kdtree.set_coords(self._getCoordinates())
            self._kdtree = kdtree
            return kdtree
        return self._kdtree


class Contacts(object):
    """A class for identification of intermolecular contacts."""
    
    def __init__(self, atoms):
        """*atoms* for which contacts will be identified. *atoms* can be 
        :class:`~prody.atomic.AtomGroup` or :class:`~prody.atomic.AtomSubset`.
        """

        if not isinstance(atoms, (AtomGroup, AtomSubset)):                
            raise TypeError('{0:s} is not a valid type for atoms'
                            .format(type(atoms)))
        self._atoms = atoms
        self._acsi = atoms.getActiveCoordsetIndex()
        self._timestamps = np.zeros(atoms.getNumOfCoordsets()) 
        self._kdtrees = [None] * atoms.getNumOfCoordsets()
        if not isinstance(atoms, AtomGroup):
            self._indices = atoms.getIndices()
            self._ag = atoms.getAtomGroup()
        else:
            self._ag = atoms 
            self._indices = None
        if KDTree is None: prody.importBioKDTree()
        if not KDTree:
            raise ImportError('Bio.KDTree is required for distance based '
                              'selections.')

    def __repr__(self):
        return '<Contacts: {0:s} (active coordset index: {1:d})>'.format(
                                                str(self._atoms), self._acsi)
    

    def _getKDTree(self):

        acsi = self._acsi
        ag = self._ag
        if ag._getTimeStamp(acsi) != self._timestamps[acsi]:    
            kdtree = KDTree(3)
            if self._indices == None:
                kdtree.set_coords(self._ag._getCoordinates())
            else:
                kdtree.set_coords(self._ag._getCoordinates()[self._indices])
            self._kdtrees[acsi] = kdtree
            self._timestamps[acsi] = ag._getTimeStamp(acsi) 
            return kdtree
        else:
            return self._kdtrees[acsi]

    def getActiveCoordsetIndex(self):
        """Return active coordinate set index."""
        
        return self._acsi
    
    def setActiveCoordsetIndex(self, acsi):
        """Set active coordinate set index.
        
        .. note:: Changing active coordinate set index effects only 
           :class:`Contacts` instance. The active coordinate set index
           of the associated :class:`~prody.atomic.Atomic` instance 
           remains the same. 
        """
        
        ts = self._timestamps
        ag = self._ag
        if acsi >= len(ts):
            n_csets = ag.getNumOfCoordsets() 
            diff = n_csets - len(ts)
            self._kdtrees += [None] * diff
            self._timestamps = np.zeros(n_csets)
            self._timestamps[:len(ts)] = ts
        self._acsi = acsi
   

    def select(self, within, what):
        """Select atoms *within* of *what*. *within* is distance in Å and 
        *what* can be point(s) in 3-d space (:class:`~numpy.ndarray` with 
        shape N,3) or a set of atoms (:class:`~prody.atomic.Atomic` instances).
        """
        
        if isinstance(what, np.ndarray):
            if what.ndim == 1 and len(what) == 3:
                what = [what]
            elif not (what.ndim == 2 and what.shape[1] == 3):
                raise SelectionError('*what* must be a coordinate array, '
                                     'shape (N, 3) or (3,).')
        else:
            try:
                what = what._getCoordinates()
            except:
                raise SelectionError('*what* must have a getCoordinates() '
                                     'method.')
            if not isinstance(what, np.ndarray):
                raise SelectionError('what.getCoordinates() method must '
                                     'return a numpy.ndarray instance.')
        kdtree = self._getKDTree()
        search = kdtree.search
        get_indices = kdtree.get_indices
        indices = []
        append = indices.append
        for xyz in what:
            search(xyz, float(within))
            append(get_indices())
        if self._indices is not None:        
            indices = self._indices[indices]
        indices = np.unique(np.concatenate(indices))
        if len(indices) != 0:
            return Selection(self._ag, np.array(indices), 
                'index {0:s}'.format(' '.join(np.array(indices, '|S'))), 
                                     self._acsi)
        return None
