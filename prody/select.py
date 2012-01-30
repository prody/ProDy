# -*- coding: utf-8 -*-
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
chid [‡]         string           same as *chain*
icode [‡]        string           single letter insertion code
segment [‡]      string           segment name
segname [‡]      string           same as *segment*
secondary [\*‡]  string           one-character secondary structure identifier
secstr [\*‡]     string           same as *secondary*
index            integer, range   internal atom number (starts from 0) 
serial           integer, range   atom serial number (parsed from file)
resnum [§]       integer, range   residue number
resid [§]        integer, range   same as *resnum*
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
:func:`~prody.proteins.assignSecstr` function.

**[†]** Alternate locations are parsed as alternate coordinate sets. This
keyword will work for alternate location specified by "A". This to work for
alternate locations indicated by other letters, they must be parsed 
specifically by passing the identifier to the :func:`~prody.proteins.parsePDB`.

**[‡]** Atoms with unspecified alternate location/chain/segment/icode/secondary 
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

Negative numbers must be entered between grave accent symbols, 
e.g. ``"resnum `-3`"``

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
    
Number ranges involving negative numbers must be entered between grave accent 
symbols, e.g. ``"resnum `-3 to 10`"``, ``"resnum `-3:10:2`"``

**More special characters (``)**

.. versionadded:: 0.7
   Strings can include the following characters (including whitespace) as well 
   when they are surrounded by grave accent character (``):

::
  
  ~!@#$%^&*()-_=+[{}]\|;:,<>./?()'"

For example ``"name `CA#` `C #`"`` will work.

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
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import time
import re as RE

import numpy as np
from . import pyparsing as pp
pp.ParserElement.enablePackrat()

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS
from atomic import *
from tools import *
import measure
DEBUG = False

__all__ = ['Select', 'Contacts',
           'getProteinResidueNames', 'setProteinResidueNames',
           'getKeywordResnames', 'setKeywordResnames',
           'getKeywordResidueNames', 'setKeywordResidueNames',
           'getBackboneAtomNames', 'setBackboneAtomNames',
           'getBackboneAtomNames', 'setBackboneAtomNames',
           'getAtomNameRegex', 'setAtomNameRegex',
           'defSelectionMacro', 'delSelectionMacro', 'getSelectionMacro',
           'getReservedWords']

KEYWORDS_STRING = set(['name', 'type', 'resname', 'chain', 'element', 
                       'segment', 'altloc', 'secondary', 'icode',
                       'chid', 'secstr', 'segname'])
KEYWORDS_INTEGER = set(['serial', 'index', 'resnum', 'resid', 
                        'segindex', 'chindex', 'resindex'])
KEYWORDS_FLOAT = set(['x', 'y', 'z', 'beta', 'mass', 'occupancy', 'mass', 
                      'radius', 'charge'])
KEYWORDS_NUMERIC = KEYWORDS_FLOAT.union(KEYWORDS_INTEGER)    

KEYWORDS_VALUE_PAIRED = KEYWORDS_NUMERIC.union(KEYWORDS_STRING)
KEYWORDS_SYNONYMS = {}
for key, field in ATOMIC_DATA_FIELDS.iteritems(): 
    if field.synonym:
        KEYWORDS_SYNONYMS[field.synonym] = key
ATOMIC_ATTRIBUTES = prody.atomic.ATOMIC_ATTRIBUTES
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
    'hydrophobic': ['ALA', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 'VAL', 
                    'XLE'],
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
    KEYWORD_RESNAMES['acyclic'] = list(protein.difference(
        set(KEYWORD_RESNAMES['cyclic'])))
    KEYWORD_RESNAMES['charged'] = list(set(KEYWORD_RESNAMES['acidic'] + 
        KEYWORD_RESNAMES['basic']))
    KEYWORD_RESNAMES['large'] = list(protein.difference(
        set(KEYWORD_RESNAMES['small'] + KEYWORD_RESNAMES['medium'])))
    KEYWORD_RESNAMES['neutral'] = list(protein.difference(
        set(KEYWORD_RESNAMES['charged'])))
    KEYWORD_RESNAMES['polar'] = list(protein.difference(
        set(KEYWORD_RESNAMES['hydrophobic'])))
    KEYWORD_RESNAMES['surface'] = list(protein.difference(
        set(KEYWORD_RESNAMES['buried'])))
    
_setReadonlyResidueNames()

__doc__ += """

Keywords without arguments
-------------------------------------------------------------------------------

Below is the list of keywords defined based on residue type and/or property.
These definitions can be retrieved or altered using 
:func:`getKeywordResnames` and :func:`setKeywordResnames`, 
respectively.

============= =================================================================
Keyword       Description
============= =================================================================
"""
keys = KEYWORD_RESNAMES.keys()
keys.sort()
for key in keys:
    if key in KEYWORD_RESNAMES_READONLY:
        __doc__ += '{0:13s} resname {1:s}\n'.format(
            key + ' [#]', ' '.join(KEYWORD_RESNAMES[key]))
    else:
        __doc__ += '{0:13s} resname {1:s}\n'.format(
            key, ' '.join(KEYWORD_RESNAMES[key]))

__doc__ += """\
============= =================================================================

**[#]** Definitions of these keywords cannot be changed directly, as they 
are defined based on others as follows: 
    
"""
keys = KEYWORD_RESNAMES_READONLY.keys()
keys.sort()
for key in keys:
    __doc__ += '  * ``"{0:s}"`` is ``"{1:s}"``\n'.format(
        key, KEYWORD_RESNAMES_READONLY[key])

__doc__ += """
.. versionchanged:: 0.6.1
   ASX (asparagine or aspartic acid), GLX (lutamine or glutamic acid),
   XLE (leucine or isoleucine), SEC (selenocysteine), and PYL
   (pyrrolysine) has been added to the standard definition of "protein" 
   keyword. Note that this list of residue names can be changed using
   :func:`setKeywordResnames` function.

The following are additional keywords whose definitions are more restricted:

=============== ===============================================================
Keyword         Description
=============== ===============================================================
all             all atoms
none            nothing (returns ``None``)
hetero          non-protein/nucleic atoms, same as \
``"not (protein or nucleic)"``
calpha (ca)     Cα atoms of protein residues, same as ``"name CA and protein"``
backbone (bb)   backbone atoms of protein residues, same as \
``"name CA C O N and protein"``
backbonefull    backbone atoms of protein residues, same as \
``"name CA C O N H H1 H2 H3 OXT and protein"``
bbful           same as ``backbonefull`` 
sidechain (sc)  side-chain atoms of protein residues, same as \
``"not name CA C O N H and protein"``
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
turn            residue in hydrogen bonded turn conformation, same as \
``"secondary T"``
bridge          residue in isolated beta-bridge conformation, same as \
``"secondary B"``
bend            residue in bend conformation, same as ``"secondary S"``
coil            residue not in one of above conformations, same as \
``"secondary C"``
=============== ===============================================================

.. versionchanged:: 0.6.2
   H is added to the list of backbone atoms.

.. versionadded:: 0.7 
   New keywords are defined: ``lipid, heme, ion, buried, surface, at, 
   cg, purine, pyrimidine, carbon, nitrogen, oxygen, sulfur, extended, helix, 
   helix_pi, helix_3_10, turn, bridge, bend, coil`` 

.. versionchanged:: 0.8.2
   H is added to the list of backbone atoms. ``backbone`` now selects maximum 
   of four heavy atoms per residue, to return consistent results for NMR and
   X-ray structures. 
   
.. versionadded:: 0.8.2
   ``backbonefull`` (or ``bbfull``) keyword is added and selects all backbone 
   atoms including hydrogens.     


Among these list of backbone atom names can be changed using 
:func:`setBackboneAtomNames`  and regular expressions for element types
can be changed using :func:`setAtomNameRegex`.

"""

SECSTR_MAP = {
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
    'carbon': RE.compile('C.*'),
    'hydrogen': RE.compile('[0-9]?H.*'),
    'nitrogen': RE.compile('N.*'),
    'oxygen': RE.compile('O.*'),
    'sulfur': RE.compile('S.*'),
}

BACKBONE_ATOM_NAMES = set(('CA', 'N', 'C', 'O'))
BACKBONE_FULL_ATOM_NAMES = set(('CA', 'N', 'C', 'O', 
                                'H', 'H1', 'H2', 'H3', 'OXT'))

KEYWORD_MAP = {}
def _buildKeywordMap():
    global KEYWORD_MAP
    
    protein = KEYWORD_RESNAMES['protein']
    #'keyword' : (residue_names, invert, atom_names, atom_names_not),
    for keyword, resnames in KEYWORD_RESNAMES.iteritems():
        KEYWORD_MAP[keyword] = (resnames, False, None, False)

    KEYWORD_MAP['alpha'] = (protein, False, ['CA'], False)
    KEYWORD_MAP['calpha'] = (protein, False, ['CA'], False)
    KEYWORD_MAP['ca'] = KEYWORD_MAP['calpha']
    KEYWORD_MAP['backbone'] = (protein, False, BACKBONE_ATOM_NAMES, False)
    KEYWORD_MAP['bb'] = KEYWORD_MAP['backbone']
    KEYWORD_MAP['backbonefull'] = (protein, False, 
                                   BACKBONE_FULL_ATOM_NAMES, False)
    KEYWORD_MAP['bbfull'] = KEYWORD_MAP['backbonefull']
    KEYWORD_MAP['sidechain'] = (protein, False, BACKBONE_FULL_ATOM_NAMES, True)
    KEYWORD_MAP['sc'] = KEYWORD_MAP['sidechain']

    KEYWORD_MAP['hetero'] = (protein + KEYWORD_RESNAMES['nucleic'], True, 
                             None, False) 

    for name, regex in KEYWORD_NAME_REGEX.iteritems():
        KEYWORD_MAP[name] = (None, False, [regex], False)
    
    KEYWORD_MAP['carbon'] = (KEYWORD_RESNAMES['ion'], True, 
                             [KEYWORD_NAME_REGEX['carbon']], False)
    KEYWORD_MAP['noh'] = (None, False, [KEYWORD_NAME_REGEX['hydrogen']], True)
    
_buildKeywordMap()
KEYWORDS_BOOLEAN = set(['all', 'none'] + KEYWORD_MAP.keys() + 
                       SECSTR_MAP.keys())

__doc__ += """

Numerical comparisons
-------------------------------------------------------------------------------

Following keywords can be used in numerical comparisons, as operands of 
arithmetic operations or as arguments to functions:  
 
 * index, serial   
 * resnum, resid
 * x, y, z
 * beta, occupancy
 * charge, mass, radius (these must be set by the user before they can be used) 

Numerical attributes of atoms can be used with the following comparison 

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

*Examples:* ``"x < 0"``, ``"occupancy != 1"``

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
   
These operations must be used with a numerical comparison, e.g. 
``"x ** 2 < 10"``, ``"x ** 2 ** 2 < 10"``, ``"occupancy != 1"``
   
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
cos(x)   cosine of x
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
tanh(x)  hyperbolic tangent of x
======== ===================================

**Examples**
  
  * ``"sqrt(x**2 + y**2 + z**2) < 10"`` selects atoms within 10 Å of the 
    origin
  * ``"resnum <= 100"`` selects atoms with residue numbers less than or equal 
    to 100  


Distance based selections
-------------------------------------------------------------------------------

Atoms within a user specified distance (Å) from a set of user specified atoms
can be selected using ``within . of ..`` keyword, e.g. ``within 5 of water``
selects atoms that are within 5 Å of water molecules. This setting will
results selecting water atoms as well.

User can avoid selecting specified atoms using ``exwithin . of ..`` setting,
e.g. ``exwithin 5 of water`` will not select water molecules and is equivalent
to ``within 5 of water and not water``


Expanding selections
-------------------------------------------------------------------------------

A selection can be expanded to include the atoms in the same *residue*, 
*chain*, or *segment* using ``same .. as ..`` setting, e.g.
``same residue as exwithin 4 of water`` will select residues that have
at least an atom within 4 Å of any water molecule.    

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
    * :func:`getKeywordResnames` 
    
  * Change keyword definitions:
    
    * :func:`setAtomNameRegex`
    * :func:`setBackboneAtomNames`
    * :func:`setKeywordResnames`

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
    
BINOP_MAP = {
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
ATOMGROUP.setCoords(np.random.random((n_atoms,3)))
ATOMGROUP.setNames(['CA']*n_atoms)
ATOMGROUP.setResnames(['GLY']*n_atoms)
ATOMGROUP.setResnums(np.arange(1,n_atoms+1))
ATOMGROUP.setChids(['A']*n_atoms)
ATOMGROUP.setAltlocs([' ']*n_atoms)
ATOMGROUP.setElements(['C']*n_atoms)
ATOMGROUP.setHeteros([False]*n_atoms)
ATOMGROUP.setOccupancies([1]*n_atoms)
ATOMGROUP.setSecstrs(['H']*n_atoms)
ATOMGROUP.setSegnames(['PDB']*n_atoms)
ATOMGROUP.setAnisous(np.random.random((n_atoms,6)))
ATOMGROUP.setAnistds(np.random.random((n_atoms,6)))
ATOMGROUP.setIcodes([' ']*n_atoms)
ATOMGROUP.setTypes(['CH2']*n_atoms)
ATOMGROUP.setBetas([0]*n_atoms)
ATOMGROUP.setCharges([0]*n_atoms)
ATOMGROUP.setMasses([12]*n_atoms)
ATOMGROUP.setRadii([1.4]*n_atoms)

MACROS = SETTINGS.get('selection_macros', {})

def areAllStrings(container):
    """Return ``True`` if all items in *container* are instances of 
    :func:`str`."""
    
    for item in container: 
        if not isinstance(item, str):
            return False
    return True

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
        SETTINGS['selection_macros'] = MACROS
        SETTINGS.save()

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
        SETTINGS['selection_macros'] = MACROS
        SETTINGS.save()

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
    """Deprecated, use :func:`getKeywordResnames`."""
    
    prody.deprecate('getKeywordResidueNames', 'getKeywordResnames')
    return getKeywordResnames(keyword) 
    
def getKeywordResnames(keyword):
    """Return residue names associated with a keyword.
    
    .. versionadded:: 0.7
    
    >>> getKeywordResnames('acidic')
    ['ASP', 'GLU']
    
    """
    
    assert isinstance(keyword, str), 'keyword must be a string instance'
    try:
        resnames = KEYWORD_RESNAMES[keyword]
        resnames.sort()
        return resnames  
    except KeyError:
        if keyword in KEYWORD_RESNAMES_READONLY:
            LOGGER.warning('"{0:s}" is defined as "{1:s}"'.format(keyword, 
                                        KEYWORD_RESNAMES_READONLY[keyword]))
        else:
            LOGGER.warning('"{0:s}" is not a keyword'.format(keyword))

def setKeywordResidueNames(keyword, resnames):
    """Deprecated, use :func:`setKeywordResnames`."""
    
    prody.deprecate('setKeywordResidueNames', 'setKeywordResnames')
    return setKeywordResnames(keyword, resnames)
    
def setKeywordResnames(keyword, resnames):
    """Change the list of residue names associated with a keyword.  *keyword* 
    must be a string, and *resnames* may be a list, tuple, or set of strings. 
    The existing list of residue names will be overwritten with the given 
    residue names.  Note that changes in keyword definitions are not saved 
    permanently.
    
    >>> setKeywordResnames('acidic', ['ASP', 'GLU'])
    """
    
    if not isinstance(keyword, str):
        raise TypeError('keyword must be a string')
    if not isinstance(resnames, (list, tuple, set)):
        raise TypeError('resnames must be a list, set, or tuple')
    if not areAllStrings(resnames):
        raise TypeError('all items in resnames must be string instances')
    
    if keyword in KEYWORD_RESNAMES_READONLY:
        LOGGER.warning('"{0:s}" is defined as "{1:s}" and cannot be changed '
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
        raise ValueError('"{0:s}" is not a valid keyword'.format(keyword))

def getAtomNameRegex(name):
    """Return regular expression used for selecting common elements.
    
    .. versionadded:: 0.7
    
    >>> getAtomNameRegex('nitrogen')
    'N.*'
    
    """
    
    assert isinstance(name, str), 'name must be a string instance'
    try:
        return KEYWORD_NAME_REGEX[name].pattern   
    except KeyError:
        LOGGER.warning('{0:s} is not a valid element'.format(keyword))

def setAtomNameRegex(name, regex):
    """Set regular expression used for selecting common elements.
    
    Note that changes in keyword definitions are not saved permanently.
    
    .. versionadded:: 0.7
    
    >>> setAtomNameRegex('nitrogen', 'N.*')
    
    """
    
    assert isinstance(name, str), 'name must be a string instance'
    if not name in KEYWORD_NAME_REGEX:
        raise ValueError('"{0:s}" is not a valid keyword'.format(name))
    if not isinstance(regex, str):
        raise TypeError('regex must be a string instance')
    try:
        regex = RE.compile(regex)
    except:
        raise ValueError('"{0:s}" is not a valid regular expression'
                         .format(regex))
    else:
        KEYWORD_NAME_REGEX[name] = regex

def getBackboneAtomNames(full=False):
    """Return protein backbone [full] atom names.
    
    >>> getBackboneAtomNames()
    ['C', 'CA', 'N', 'O']
    
    .. versionchanged:: 0.8.2
       User can get the atom names for *backbonefull* keyword by passing
       ``full=True`` argument.

    """
    
    assert isinstance(full, bool), 'full must be a boolean instance'
    if full:
        bban = list(BACKBONE_FULL_ATOM_NAMES)
    else:
        bban = list(BACKBONE_ATOM_NAMES)
    bban.sort()
    return bban 

def setBackboneAtomNames(backbone_atom_names, full=False):
    """Set protein backbone [full] atom names.
    
    Note that changes in keyword definitions are not saved permanently.
    
    .. versionchanged:: 0.8.2
       User can set the atom names for *backbonefull* keyword by passing
       ``full=True`` argument.
    
    """
    
    if not isinstance(backbone_atom_names, (list, tuple, set)):
        raise TypeError('backbone_atom_names must be a list, tuple, or set')
    if not areAllStrings(backbone_atom_names):
        raise TypeError('all items in backbone_atom_names must be string '
                        'instances')
    assert isinstance(full, bool), 'full must be a boolean instance'
    if full:    
        global BACKBONE_FULL_ATOM_NAMES
        BACKBONE_FULL_ATOM_NAMES = set(backbone_atom_names)
    else:
        global BACKBONE_ATOM_NAMES
        BACKBONE_ATOM_NAMES = set(backbone_atom_names)
    _buildKeywordMap()

def getProteinResidueNames():
    """Deprecated, use :func:`getKeywordResnames`."""
    
    prody.deprecate('getProteinResidueNames', 'getKeywordResnames')
    return getKeywordResnames('protein')
    
def setProteinResidueNames(resnames):
    """Deprecated, use :func:`setKeywordResnames`."""
    
    prody.deprecate('setProteinResidueNames', 'setKeywordResnames')
    return setKeywordResidueNames('protein', resnames)
    
class SelectionError(Exception):    
    
    def __init__(self, selstr, *args):
        selstr = selstr.replace(AND, 'and').replace(OR, 'or').replace(NOT, 
                                                                     'not')
        Exception.__init__(self, '"{0:s}" is not a valid selection string. '
                                 .format(selstr) + ' '.join(args) )


def isFloatKeyword(keyword):
    return keyword in KEYWORDS_FLOAT

def isIntKeyword(keyword):
    return keyword in KEYWORDS_INTEGER

def isNumericKeyword(keyword):
    return keyword in KEYWORDS_NUMERIC

def isAlnumKeyword(keyword):
    return keyword in KEYWORDS_STRING

def isValuePairedKeyword(keyword):
    return keyword in KEYWORDS_VALUE_PAIRED

def isBooleanKeyword(keyword):
    return keyword in KEYWORDS_BOOLEAN
    
def isKeyword(keyword):
    return isBooleanKeyword(keyword) or isValuePairedKeyword(keyword)

AND = '&&&'
NOT = '!!!'
OR  = '||'



RESERVED = set(ATOMIC_DATA_FIELDS.keys() + ATOMIC_ATTRIBUTES.keys() +
               ['and', 'or', 'not', 'within', 'of', 'exwithin', 'same', 'as'] +
               KEYWORDS_SYNONYMS.keys() + 
               ['n_atoms', 'n_csets', 'cslabels', 'title', 'coordinates',
                'bonds', 'bmap', 'numbonds'])

def isReserved(word):
    return word in RESERVED or isKeyword(word) or word in FUNCTION_MAP
        
        
def getReservedWords():
    """Return a list of words reserved for atom selections and internal 
    variables. These words are: """

    words = list(set(list(RESERVED) + FUNCTION_MAP.keys() + 
                     list(KEYWORDS_BOOLEAN) + list(KEYWORDS_VALUE_PAIRED)))
    
    words.sort()
    return words

getReservedWords.__doc__ += "*{0:s}*.".format('*, *'.join(getReservedWords()))

_specialKeywords = set(['secondary', 'chain', 'altloc', 'segment', 'icode'])

def tkn2str(token):
    
    if isinstance(token, str):
        return token
    else:
        return ' '.join(token)


def expandBoolean(keyword):
    
    if keyword in KEYWORD_MAP:
        (residue_names, rn_invert, atom_names, 
                                        an_invert) = KEYWORD_MAP[keyword]
        tokens = []
        if atom_names is not None:
            if an_invert:
                tokens.append(NOT)
            tokens.append('name')
            tokens.extend(atom_names)
            if residue_names is not None:
                tokens.append(AND)
        if residue_names is not None:
            
            if rn_invert:
                tokens.append(NOT)
            tokens.append('resname')
            tokens.extend(residue_names)
        return tokens
    elif keyword in SECSTR_MAP:
        return ['secondary', SECSTR_MAP[keyword]]
    else:
        return keyword

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
        self._acsi = None
        self._timestamp = None
        
        self._coords = None
        self._kdtree = None
        self._kwargs  = None
        self._selstr2indices = False
        self._data = dict()
        for var in mapField2Var.values():
            self._data[var] = None        
        
        shortlist = pp.alphanums + '''~@#$.:;_','''
        longlist = pp.alphanums + '''~!@#$%^&*()-_=+[{}]\|;:,<>./?()' '''
        specialchars = pp.Group(pp.Literal('`') + 
                                pp.Optional(pp.Word(longlist + '"')) + 
                                pp.Literal('`'))
        def specialCharsParseAction(token):
            if len(token[0]) == 2:
                return '_'
            else:
                return token[0][1]
        specialchars.setParseAction(specialCharsParseAction)
        regularexp = pp.Group(pp.Literal('"') + 
                              pp.Optional(pp.Word(longlist + '`')) + 
                              pp.Literal('"'))
        def regularExpParseAction(token): 
            token = token[0]
            if len(token[0]) == 2:
                return RE.compile('^()$')
            else:
                try:
                    regexp = RE.compile('^(' + token[1] + ')$')
                except:
                    raise SelectionError('failed to compile regular expression'
                                         ' "{0:s}"'.format(token[1]))
                else:
                    return regexp  
        regularexp.setParseAction(regularExpParseAction)
        oneormore = pp.OneOrMore(pp.Word(shortlist) | regularexp | 
                                 specialchars)
        funcnames = FUNCTION_MAP.keys()
        functions = pp.Keyword(funcnames[0])
        for func in funcnames[1:]:
            functions = functions | pp.Keyword(func)
        self._tokenizer = pp.operatorPrecedence(
             oneormore,
             [(functions, 1, pp.opAssoc.RIGHT, self._func),
              (pp.oneOf('+ -'), 1, pp.opAssoc.RIGHT, self._sign),
              (pp.oneOf('** ^'), 2, pp.opAssoc.LEFT, self._pow),
              (pp.oneOf('* / %'), 2, pp.opAssoc.LEFT, self._mul),
              (pp.oneOf('+ -'), 2, pp.opAssoc.LEFT, self._add),
              (pp.oneOf('< > <= >= == = !='), 2, pp.opAssoc.LEFT, self._comp),
              (pp.Keyword(NOT) | 
               pp.Regex('same [a-z]+ as') | 
               pp.Regex('(ex)?within [0-9]+\.?[0-9]* of'), 
                        1, pp.opAssoc.RIGHT, self._unary),
              (pp.Keyword(AND), 2, pp.opAssoc.LEFT, self._and),
              (pp.Keyword(OR), 2, pp.opAssoc.LEFT, self._or),]
            )

        self._tokenizer.setParseAction(self._defaultAction)
        self._tokenizer.leaveWhitespace()
        
        
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
            if self._acsi != atoms.getACSIndex():
                self._coords = None
                self._kdtree = None
            elif self._timestamp != atoms._getTimeStamp(self._acsi):
                self._kdtree = None
        else:
            self._reset()
            if isinstance(atoms, prody.AtomGroup): 
                self._ag = atoms
                self._atoms = atoms
                self._indices = None
                self._n_atoms = atoms.numAtoms()
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
        self._acsi = atoms.getACSIndex()
        self._timestamp = atoms._getTimeStamp(self._acsi)
            
        self._kwargs = kwargs
        if DEBUG:
            print('getBoolArray', selstr)
        torf = self._evalSelstr()
        if not isinstance(torf, np.ndarray):
            raise SelectionError(selstr)
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
                                 selstr, atoms.getTitle()),
                                 atoms.getACSIndex())
        else:
            
            if self._selstr2indices:
                selstr = 'index {0:s}'.format(rangeString(indices))
            elif isinstance(atoms, prody.AtomPointer):
                selstr = '({0:s}) and ({1:s})'.format(selstr, 
                                                      atoms.getSelstr())
            
            return prody.Selection(ag, indices, selstr, atoms.getACSIndex(),
                                   unique=True)
        
    def _reset(self):
        if DEBUG: print('_reset')
        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._acsi = None
        self._timestamp = None
        self._coords = None
        self._kdtree = None
        for var in mapField2Var.values():
            self._data[var] = None        

    def _prepareSelstr(self):
        if DEBUG: print('_prepareSelstr', self._selstr) 
        selstr = ' ' + self._selstr + ' '
        selstr = selstr.replace(')and(', ')&&&(')
        selstr = selstr.replace(' and(', ' &&&(')
        selstr = selstr.replace(')and ', ')&&& ')
        while ' and ' in selstr:
            selstr = selstr.replace(' and ', ' &&& ')
            
        selstr = selstr.replace(')or(', ')||(')
        selstr = selstr.replace(' or(', ' ||(')
        selstr = selstr.replace(')or ', ')|| ')
        while ' or ' in selstr:
            selstr = selstr.replace(' or ', ' || ')
        
        #if selstr.startswith('not '):
        #    selstr = selstr.replace('not ', '!!! ')
        selstr = selstr.replace('(not ', '(!!! ')
        selstr = selstr.replace(' not(', ' !!!(')
        while ' not ' in selstr:
            selstr = selstr.replace(' not ', ' !!! ')
        
        if MACROS:
            for macro in MACROS.iterkeys():
                selstr = selstr.replace(' ' + macro + ' ', 
                                        ' (' + MACROS[macro] + ') ')
                selstr = selstr.replace('(' + macro + ' ', 
                                        '((' + MACROS[macro] + ') ')
                selstr = selstr.replace(' ' + macro + ')', 
                                        ' (' + MACROS[macro] + '))')
        
        if DEBUG: print('_prepareSelstr', selstr) 
        return selstr.strip()

    def _evalSelstr(self):
        selstr = self._selstr.strip() 
        if DEBUG: print('_evalSelstr', selstr)
        if len(selstr.split()) == 1 and '(' not in selstr and \
           ')' not in selstr and selstr not in MACROS:
            if isBooleanKeyword(selstr):
                return self._evalBoolean(selstr)
            elif self._ag.isData(selstr):
                return self._evalUserdata(selstr)
            elif isValuePairedKeyword(selstr):
                raise SelectionError(selstr, '"{0:s}" must be followed by '
                                     'values.'.format(selstr))
            else:
                raise SelectionError(selstr, '"{0:s}" is not a user set atom '
                                     'group attribute either.'.format(selstr))
        
        selstr = self._prepareSelstr()
        try:
            if DEBUG: print('_evalSelstr using Pyparsing')
            tokens = self._tokenizer.parseString(selstr, 
                                             parseAll=True).asList()
            if DEBUG: print('_evalSelstr', tokens)
            return tokens[0]
        except pp.ParseException as err:
            raise SelectionError(selstr, '\n' + ' ' * (err.column + 16) + 
                                         '^ parsing the rest failed.')
    
    def _isValid(self, token):
        """Check the validity of part of a selection string. Expects a Python
        :func:`str` or a :func:`list`."""
        
        if DEBUG: print('_isValid', token)
        
        if isinstance(token, str):
            return isBooleanKeyword(token) or \
                   self._atoms.getDataType(token) == bool   
        elif isinstance(token, list):
            tkn = token[0]
            return isValuePairedKeyword(tkn) or tkn in self._kwargs or \
                    self._atoms.isData(tkn) or tkn == NOT
        return False
        
    def _defaultAction(self, tokens):
        if DEBUG: print('_evaluate', tokens)

        if isinstance(tokens[0], np.ndarray):
            return tokens[0]
        torf = self._evaluate(tokens)
        if torf is None:
            try:
                selstr = ' '.join(tokens)
            except:
                selstr = tokens[0]
            raise SelectionError(selstr)
        return torf
    
    def _evaluate(self, tokens, evalonly=None):
        if DEBUG: print('_evaluate', tokens)
        
        if isinstance(tokens, str):
            if isBooleanKeyword(tokens):
                return self._evalBoolean(tokens, evalonly=evalonly)
            elif self._ag.isData(tokens):
                return self._evalUserdata(tokens, evalonly=evalonly)
            else:
                return None
        elif isinstance(tokens, (np.ndarray, float)):
            return tokens
    
        keyword = tokens[0]
        if len(tokens) == 1:
            if isBooleanKeyword(keyword):
                return self._evalBoolean(keyword, evalonly=evalonly)
            elif isNumericKeyword(keyword):
                return self._evalNumeric(keyword)
            elif self._ag.isData(keyword):
                return self._evalUserdata(keyword, evalonly=evalonly)
            elif isinstance(keyword, np.ndarray):
                return keyword
            else:
                try:
                    return float(keyword)
                except ValueError:
                    pass
        elif isAlnumKeyword(keyword):
            return self._evalAlnum(keyword, tokens[1:], evalonly=evalonly)
        elif keyword in ('resnum', 'resid'):
            return self._resnum(tokens[1:], evalonly=evalonly)
        elif keyword == 'index':
            return self._index(tokens[1:], evalonly=evalonly)
        elif keyword == 'serial':
            return self._serial(tokens[1:], evalonly=evalonly)
        elif isFloatKeyword(keyword) or isIntKeyword(keyword):
            return self._evalFloat(keyword, tokens[1:], evalonly=evalonly)
        elif keyword == NOT:
            return self._not(tokens, evalonly=evalonly)
        elif self._ag.isData(keyword):
            return self._evalUserdata(keyword, tokens[1:], evalonly=evalonly)
        return None

    def _or(self, selstr, location, tokens):
        if DEBUG: print('_or\n_or tokens '+str(tokens))
        previous = None
        evalonly = None
        selection = None
        for current in tokens[0]:
            if previous is None:
                previous = current
                continue
            if current == OR:
                if isinstance(previous, np.ndarray):
                    if evalonly is None:
                        torf = previous
                    else:
                        torf = previous[evalonly]
                else:
                    if not self._isValid(previous):
                        raise SelectionError(selstr)
                    torf = self._evaluate(previous, evalonly=evalonly)
                    if torf is None:
                        raise SelectionError(selstr)
                if evalonly is None:
                    selection = torf
                    evalonly = np.invert(selection).nonzero()[0]
                else:
                    selection[evalonly[torf]] = True
                    evalonly = evalonly[torf.nonzero()[0]]
                previous = None
            else:
                if isinstance(previous, str):
                    previous = [previous, current]
                else:
                    previous.append(current)
        if DEBUG: print('_or last item', previous)
        if isinstance(previous, np.ndarray):
            if evalonly is None:
                torf = previous
            else:
                torf = previous[evalonly]
        else:
            if not self._isValid(previous):
                raise SelectionError(selstr)
            torf = self._evaluate(previous, evalonly=evalonly)
            if torf is None:
                raise SelectionError(selstr)
        selection[evalonly[torf]] = True
        return selection

    def _and(self, selstr, location, tokens):
        if DEBUG: print('_and\n_and tokens '+str(tokens))
        evalonly = None
        if DEBUG and evalonly is not None: print('_and evalonly ', len(evalonly))
        previous = None
        selection = None
        for current in tokens[0]:
            if previous is None:
                previous = current
                continue
            if current == AND:
                if isinstance(previous, np.ndarray):
                    if evalonly is None:
                        torf = previous
                    else:
                        torf = previous[evalonly]
                else:
                    if not self._isValid(previous):
                        raise SelectionError(selstr)
                    
                    torf = self._evaluate(previous, evalonly=evalonly)
                    if torf is None:
                        raise SelectionError(selstr)
                if selection is None:
                    selection = torf
                if evalonly is None:
                    evalonly = selection.nonzero()[0]
                else:
                    selection[evalonly] = torf
                    evalonly = evalonly[torf]
                previous = None
            else:
                if isinstance(previous, str):
                    previous = [previous, current]
                else:
                    previous.append(current)
        if DEBUG: print('_and last item', previous)
        if isinstance(previous, np.ndarray):
            if evalonly is None:
                torf = previous
            else:
                torf = previous[evalonly]
        else:
            if not self._isValid(previous):
                raise SelectionError(selstr)
            torf = self._evaluate(previous, evalonly=evalonly)
            if torf is None:
                raise SelectionError(selstr)
        if selection is None:
            selection = torf
        else:
            selection[evalonly] = torf
        return selection
    
    def _unary(self, selstr, location, tokens):
        """Perform the unary operation."""
        
        if DEBUG: print('_unary', tokens)
        tokens = tokens[0]
        if tokens[0] == NOT:
            result = self._not(tokens)
        elif tokens[0].startswith('same'):
            result = self._sameas(tokens)
        else:
            result = self._within(tokens, tokens[0].startswith('exwithin'))
        if result is None:
            raise SelectionError(selstr)
        return result

    def _not(self, tokens, evalonly=None):
        """Negate selection."""
        
        if DEBUG: print('_not', tokens)
        if isinstance(tokens[1], np.ndarray):
            torf = tokens[1]
        else:
            torf = self._evaluate(tokens[1:], evalonly=evalonly)
            if torf is None:
                return None
        np.invert(torf, torf)
        return torf
    
    def _within(self, tokens, exclude):
        """Perform distance based selection."""

        if DEBUG: print('_within', tokens)
        try:
            within = float(tokens[0].split()[1])
        except:
            return None
        if DEBUG: print('_within', within)
        which = tokens[1]
        if not isinstance(which, np.ndarray):
            which = self._evaluate(tokens[1:])

        if DEBUG: print('_within', which)
        if which is None and self._kwargs is not None:
            try:
                coordinates = self._kwargs[tokens[1]]
            except:
                return None
            if DEBUG: print('_kwargs', tokens[1])
            if isinstance(coordinates, np.ndarray):
                if coordinates.ndim == 1 and len(coordinates) == 3:
                    coordinates = np.array([coordinates])
                elif not (coordinates.ndim == 2 and coordinates.shape[1] == 3):
                    return None
            else:
                try:
                    coordinates = coordinates.getCoords()
                except:
                    return None
                if not isinstance(coordinates, np.ndarray):
                    return None
            exclude=False
            self._selstr2indices = True
            which = np.arange(len(coordinates))
        elif isinstance(which, np.ndarray) and which.dtype == np.bool: 
            which = which.nonzero()[0]
            coordinates = self._getCoords()
        else:
            return None
        selection = []
        append = selection.append
        kdtree = self._getKDTree()
        get_indices = kdtree.get_indices
        search = kdtree.search
        for index in which:
            search(coordinates[index], within)
            append(get_indices())
        unique = np.unique(np.concatenate(selection))
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
        return torf
    
    def _sameas(self, token):
        """Expand selection."""
        
        terms = token
        if DEBUG: print('_sameas', terms)
        what = token[0].split()[1]
        which = token[1]
        if not isinstance(which, np.ndarray):
            which = self._evaluate(token[1:])
            if not isinstance(which, np.ndarray):
                return None
        self._ag.getHierView()
        if what == 'residue':
            resindex = self._getAtomicData('resindex')[which]
            torf = self._evalFloat('resindex', np.unique(resindex))
            #chainids = self._getAtomicData('chain')[which]
            #resnum =  self._getAtomicData('resnum')[which]
            #icodes = self._getAtomicData('icode')[which]
            #if icodes is None:
            #    resnum = [str(resnum[i]) + code 
            #                        if code else resids[i] 
            #                        for i, code in enumerate(icodes)]
            #torf = np.zeros(self._n_atoms, np.bool)
            #for chain in np.unique(chainids):
            #    torf[np.all([self._evalAlnum('chain', [chain]),
            #                 self._resnum(np.unique(resnum[chainids == chain]), 
            #                              numRange=False)], 0)] = True
        elif what == 'chain':
            chindex = self._getAtomicData('chindex')[which]
            torf = self._evalFloat('chindex', np.unique(chindex))        
            #chainids = self._getAtomicData('chain')
            #torf = self._evalAlnum('chain', list(np.unique(chainids[which])))        
        elif what == 'segment':
            segindex = self._getAtomicData('segindex')[which]
            torf = self._evalFloat('segindex', np.unique(segindex))
            #segnames = self._getAtomicData('segment')
            #torf = self._evalAlnum('segment', list(np.unique(segnames[which])))
        else: 
            raise SelectionError('"{0:s}" is not valid, selections can be '
                                 'expanded to same "chain", "residue", or ' 
                                 '"segment"'.format(token[0]))
        return torf
     
    def _comp(self, selstr, location, tokens):
        """Perform numeric comparisons. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_comp', tokens)
        tokens = tokens[0]
        if len(tokens) >= 3 and len(tokens) % 2 != 1:
            raise SelectionError(selstr)
        i = 1
        left = self._evalNumeric(tokens[0])
        if DEBUG: print('_comp left', left)
        if left is None:
            raise SelectionError(selstr)
        result = None
        while i < len(tokens): 
            comp = tokens[i]
            right = self._evalNumeric(tokens[i + 1])
            if DEBUG: print('_comp right', right)
            if right is None:
                raise SelectionError(selstr)
            if result is None:
                result = BINOP_MAP[comp](left, right)
            else:
                result *= BINOP_MAP[comp](left, right)
            left = right
            i += 2
        return result

    def _pow(self, selstr, location, tokens):
        """Perform power operation. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_pow', tokens)
        tokens = tokens[0]
        base = self._evalNumeric(tokens.pop(0))
        if base is None:
            raise SelectionError(selstr)
        power = self._evalNumeric(tokens.pop())
        if power is None:
            raise SelectionError(selstr)
        tokens.pop()
        while tokens:
            number = self._evalNumeric(tokens.pop()) 
            if number is None:
                raise SelectionError(selstr)
            power = number * power
            tokens.pop()
        return base ** power

    def _add(self, selstr, location, tokens):
        """Perform addition operations. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_add', tokens)
        tokens = tokens[0]
        left = self._evalNumeric(tokens.pop(0))
        if left is None:
            raise SelectionError(selstr)
        while tokens:
            op = tokens.pop(0)
            right = self._evalNumeric(tokens.pop(0))
            if right is None:
                raise SelectionError(selstr)
            left = BINOP_MAP[op](left, right)
        if DEBUG: print('_add total', left)
        return left
 
    def _mul(self, selstr, location, tokens):
        """Perform multiplication operations. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_mul', tokens)
        tokens = tokens[0]
        left = self._evalNumeric(tokens[0])
        if left is None:
            raise SelectionError(selstr)
        i = 1
        while i < len(tokens):
            op = tokens[i]
            i += 1
            right = self._evalNumeric(tokens[i])
            if right is None:
                raise SelectionError(selstr)
            i += 1
            if op == '/' and right == 0.0: 
                raise SelectionError(selstr, 
                                     'This leads to zero division error.')
            left = BINOP_MAP[op](left, right)
        return left
    
    def _evalNumeric(self, token):
        """Evaluate a number operand or a numeric keyword."""
        
        if DEBUG: print('_evalNumeric', token)
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
        elif isNumericKeyword(token): 
            return self._getAtomicData(token)
        elif self._ag.isData(token):
            data = self._getAtomicData(token)
            if data.dtype in (np.float, np.int):
                return data
            else:
                return None
        else:
            try:
                token = float(token)
            except ValueError:
                return None
            else:
                return token
    
    def _sign(self, selstr, location, tokens):
        tokens = tokens[0]
        if DEBUG: print('_sign', tokens)
        token = self._evalNumeric(tokens[1])
        if token is None:
            raise SelectionError(selstr)
        if tokens[0] == '-':
            return -token
        return token

    def _func(self, selstr, location, tokens):
        if DEBUG: print('_func', tokens)
        tokens = tokens[0]
        if len(tokens) != 2:
            raise SelectionError(selstr)
        arg = tokens[1]
        if not isinstance(arg, (np.ndarray, float)):
            arg = self._evaluate(arg)

        if isinstance(arg, float) or \
            isinstance(arg, np.ndarray) and \
            arg.dtype in (np.float, np.int):
            return FUNCTION_MAP[tokens[0]](arg)
        else:
            raise SelectionError(selstr)

    def _evalUserdata(self, keyword, values=None, evalonly=None):
        if DEBUG: print('_evalAttribute', keyword, values)
        data = self._atoms.getData(keyword)
        if values is None:
            if data.dtype == bool:
                if evalonly is None:
                    return data
                else:
                    return data[evalonly]
            else:
                return None
        else:
            if data.dtype in (np.int, np.float):
                return self._evalFloat(keyword, values, evalonly=evalonly)
            elif data.dtype.type == np.string_:
                return self._evalAlnum(keyword, values, evalonly=evalonly)
            else:
                return None


    def _evalBoolean(self, keyword, evalonly=None):
        """Evaluate a boolean keyword."""
    
        if DEBUG: print('_evalBoolean', keyword)
        if evalonly is None:
            n_atoms = self._n_atoms
        else:        
            n_atoms = len(evalonly)
        
        if keyword == 'all':
            return np.ones(n_atoms, np.bool)
        elif keyword == 'none':
            return np.zeros(n_atoms, np.bool)
        else:
            torf = self._and(keyword, 0, [expandBoolean(keyword)])
            if evalonly is None:
                return torf
            else:
                return torf[evalonly] 

    def _evalAlnum(self, keyword, values, evalonly=None):
        """Evaluate keywords associated with alphanumeric data, e.g. residue 
        names, atom names, etc."""
        
        if DEBUG: print('_evalAlnum', keyword, values)
        keyword = KEYWORDS_SYNONYMS.get(keyword, keyword)
        data = self._getAtomicData(keyword)
        if keyword in _specialKeywords:
            for i, value in enumerate(values):
                if value == '_':
                    values[i] = ' '
                    values.append('')
                    break
            
        if evalonly is not None:
            data = data[evalonly]
        n_atoms = len(data)
        
        regexps = []
        strings = []
        for value in values:
            if isinstance(value, str):
                strings.append(value)
            else:
                regexps.append(value)
                
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
            for i in xrange(n_atoms):
                torf[i] = (value.match(data[i]) is not None)

        return torf
    
    def _evalFloat(self, keyword, values=None, evalonly=None):
        """Evaluate a keyword associated with atom attributes of type float. 
        If *values* is not passed, return the attribute array."""
        
        if DEBUG: print('_evalFloat', keyword, values)
        if keyword == 'x':
            data = self._getCoords()[:,0]
        elif keyword == 'y':
            data = self._getCoords()[:,1]
        elif keyword == 'z':
            data = self._getCoords()[:,2]
        else:
            data = self._getAtomicData(keyword)
        
        if values is None:
            return data
    
        if evalonly is not None:
            data = data[evalonly]
        n_atoms = len(data)
        torf = np.zeros(n_atoms, np.bool)

        numbers = self._getNumRange(values)
        if numbers is None:
            return None
        for item in numbers:
            if isinstance(item, str):
                pass
            elif isinstance(item, list):
                torf[(item[0] <= data) & (data <= item[1])] = True
            elif isinstance(item, tuple):
                if len(item) == 2:
                    torf[(item[0] <= data) & (data < item[1])] = True
                else:
                    return None
            else:
                torf[data == item] = True
        return torf

    def _resnum(self, token=None, numRange=True, evalonly=None):
        """Evaluate "resnum" keyword."""
        
        if DEBUG: print('_resnum', token)
        if token is None:
            return self._getAtomicData('resnum') 
        icodes = None
        if evalonly is None:
            resids = self._getAtomicData('resnum')
            n_atoms = self._n_atoms
        else:
            resids = self._getAtomicData('resnum')[evalonly]
            n_atoms = len(evalonly)
        torf = np.zeros(n_atoms, np.bool)
        
        if numRange:
            token = self._getNumRange(token, False)
            if token is None:
                return None
        
        for item in token:
            if isinstance(item, str):
                if icodes is None:
                    if evalonly is None:
                        icodes = self._getAtomicData('icode')
                    else:
                        icodes = self._getAtomicData('icode')[evalonly]
                icode = str(item[-1])
                if icode == '_':
                    icode = ''
                try:
                    number = int(item[:-1])
                except ValueError:
                    return None
                torf[(resids == number) * (icodes == icode)] = True
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

    def _serial(self, token=None, evalonly=None):
        """Evaluate "serial" keyword."""
        
        if DEBUG: print('_serial', token)
        if token is None:
            return self._getAtomicData('serial') 
        if evalonly is None:
            serials = self._getAtomicData('serial')
            n_atoms = self._n_atoms
        else:
            serials = self._getAtomicData('serial')[evalonly]
            n_atoms = len(evalonly)
        torf = np.zeros(n_atoms, np.bool)
        
        numbers = self._getNumRange(token)
        if numbers is None:
            return None
        for item in numbers:
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
        if DEBUG: print('_serial n_selected', torf.sum())
        return torf
    
    def _index(self, token=None, evalonly=None):
        """Evaluate "index" keyword."""
        
        if DEBUG: print('_index', token)
        if token is None:
            return self._indices or np.arange(self._ag._n_atoms)
        torf = np.zeros(self._ag._n_atoms, np.bool)
        
        numbers = self._getNumRange(token)
        if numbers is None:
            return None
        for item in numbers:
            try:
                if isinstance(item, tuple):
                    if len(item) == 2:
                        torf[item[0]:item[1]] = True
                    else:
                        torf[item[0]:item[1]:item[2]] = True
                elif isinstance(item, list):
                    torf[int(np.ceil(item[0])):int(
                                        np.floor(item[1]))+1] = True
                else:
                    torf[item] = True
            except IndexError:
                pass
        if DEBUG: print('_index n_selected', torf.sum())
        if self._indices is None:
            if evalonly is None:
                return torf
            else:
                return torf[evalonly]
        else:
            if evalonly is None:
                return torf[self._indices]
            else:
                return torf[self._indices][evalonly]

    def _getNumRange(self, token, intfloat=True):
        """Evaluate numeric values. Identify ranges, integers, and floats,
        put them in a list and return."""
        
        if DEBUG: print('_getNumRange', type(token), token)
        if isinstance(token, np.ndarray):
            return token
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
                # to means upper bound is included in the range
                # boundaries are placed in a LIST
                items = item.split('to')
                if len(items) != 2:
                    raise SelectionError('"{0:s}" is not understood.'
                                         .format(' to '.join(items)))
                try:
                    token.append([float(items[0]), float(items[1])])
                except:
                    raise SelectionError('"{0:s}" is not understood, "to" '
                                         'must be surrounded by numbers.'
                                         .format(' to '.join(items)))
            elif ':' in item:
                # : means upper bound is NOT included in the range
                # boundaries are placed in a TUPLE
                items = item.split(':')
                if not len(items) in (2, 3):
                    raise SelectionError('"{0:s}" is not understood.'
                                         .format(':'.join(items)))
                try:
                    if len(items) == 2:
                        token.append( (int(items[0]), int(items[1])) )
                    else:
                        token.append((int(items[0]), int(items[1]),
                                      int(items[2])))
                except:
                    raise SelectionError('"{0:s}" is not understood, ":" must '
                                         'be surrounded by integers.'
                                         .format(':'.join(items)))
            else:
                try: 
                    item = int(item)
                except ValueError:
                    try:
                        item = float(item)
                    except ValueError:
                        if intfloat:
                            return None
                token.append(item)
        if DEBUG: print('_getNumRange', token)            
        return token
    
    def _getAtomicData(self, keyword):
        """Return atomic data."""
        
        field = ATOMIC_DATA_FIELDS.get(keyword, None)
        indices = self._indices
        if field is None:
            data = self._atoms.getData(keyword)
            if data is None:
                raise SelectionError('"{0:s}" is not a valid keyword or '
                                     'attribute.'.format(keyword))
            elif isinstance(data, np.ndarray) and data.ndim == 1:
                if indices is None:                
                    return data
                else:
                    return data[indices]
            else:
                raise SelectionError('attribute "{0:s}" must be a 1-dimensional'
                                     ' numpy array'.format(keyword))
        else:
            var = field.var
            data = self._data[var]
            if data is None:
                data = getattr(self._ag, '_get' + field.meth_pl)() 
                if data is None:
                    raise SelectionError('{0:s} are not set.'
                                         .format(field.doc_pl))
                self._data[var] = data
            if indices is None:                
                return data
            else:
                return data[indices]
    
    def _getCoords(self):
        """Return atomic coordinates."""
        
        if self._coords is None:
            if self._indices is None:
                self._coords = self._ag._getCoords()
            else:
                self._coords = self._atoms._getCoords()
        if self._coords is None:
            raise AttributeError('coordinates are not set')
        return self._coords
   
    def _getKDTree(self):
        """Return KDTree."""
        
        if self._indices is None:
            return self._ag._getKDTree(self._acsi)
        else:
            if self._kdtree is None:
                kdtree = measure.getKDTree(self._atoms._getCoords())
                self._kdtree = kdtree
            return kdtree


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
        self._acsi = atoms.getACSIndex()
        self._timestamps = np.zeros(atoms.numCoordsets()) 
        self._kdtrees = [None] * atoms.numCoordsets()
        if not isinstance(atoms, AtomGroup):
            self._indices = atoms.getIndices()
            self._ag = atoms.getAtomGroup()
        else:
            self._ag = atoms 
            self._indices = None

    def __repr__(self):
        return '<Contacts: {0:s} (active coordset index: {1:d})>'.format(
                                                str(self._atoms), self._acsi)

    def _getKDTree(self):

        acsi = self._acsi
        ag = self._ag
        indices = self._indices
        if indices == None:
            kdtree = self._ag._getKDTree(acsi)
        else:
            if ag._getTimeStamp(acsi) != self._timestamps[acsi]:    
                kdtree = measure.getKDTree(self._ag._getCoords()[indices])
            self._kdtrees[acsi] = kdtree
            self._timestamps[acsi] = ag._getTimeStamp(acsi) 
        return kdtree

    def getACSIndex(self):
        """Return active coordinate set index."""
        
        return self._acsi
    
    def setACSIndex(self, acsi):
        """Set active coordinate set index.
        
        .. note:: Changing active coordinate set index effects only 
           :class:`Contacts` instance. The active coordinate set index
           of the associated :class:`~prody.atomic.Atomic` instance 
           remains the same. 
        """
        
        ts = self._timestamps
        ag = self._ag
        if acsi >= len(ts):
            n_csets = ag.numCoordsets() 
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
                what = what._getCoords()
            except:
                raise SelectionError('*what* must have a getCoordinates() '
                                     'method.')
            if not isinstance(what, np.ndarray):
                raise SelectionError('what.getCoords() method must '
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
                                     self._acsi, unique=True)
        return None
