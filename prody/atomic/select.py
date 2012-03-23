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

""" This module defines a class for selecting subsets of atoms and functions 
to learn and change definitions of selection keywords.

.. _selections:


Atom selections
===============================================================================


ProDy offers a powerful atom selector for :class:`~.AtomGroup` and other 
:mod:`~prody.atomic` classes.  The keywords, selection grammar, and features 
of the selector are similar to those found in VMD (|vmd|).  Small differences 
between the two should not affect most practical uses of atom selections. 
ProDy selection engine also enables the identification of intermolecular 
contacts.  This section describes the keywords and selection syntax.

|more| See :ref:`contacts` and :ref:`selection-operations` for more usage
examples.

The contents of this web page can be viewed in an interactive session as 
follows:
    
>>> from prody import *
>>> # help(select)

    
Atom attributes 
===============================================================================

Below is the list of atomic attributes that can be used in atom selections:

===============  ==============  ============================================
Keyword          Arguments         Description
===============  ==============  ============================================
name             string          atom name
element          string          element symbol
type [*]         string          atom type
altloc [†‡]      string          one-character alternate location identifier
resname          string          residue name
chain [‡]        string          one-character chain identifier
chid [‡]         string          same as *chain*
icode [‡]        string          single letter insertion code
segment [‡]      string          segment name
segname [‡]      string          same as *segment*
secondary [\*‡]  string          one-character secondary structure identifier
secstr [\*‡]     string          same as *secondary*
sequence         string          one-letter amino acid sequence
index            integer, range  internal atom number (starts from 0) 
serial           integer, range  atom serial number (parsed from file)
resnum [§]       integer, range  residue number
resid [§]        integer, range  same as *resnum*
resindex [¶]     integer, range  unique index number for distinct residues  
chindex [¶]      integer, range  unique index number for distinct chains
segindex [¶]     integer, range  unique index number for distinct segments
x                float, range    x coordinate
y                float, range    y coordinate
z                float, range    z coordinate
beta             float, range    β (temperature) factor
occupancy        float, range    atomic occupancy value
charge [*]       float, range    atomic charge
mass [*]         float, range    atomic mass
radius [*]       float, range    atomic radius
===============  ==============  ============================================

**[*]** These atomic attributes are not set by the PDB parser when a PDB file 
is parsed. Using them before they are set will raise selection error. 
Secondary structure assignments can be made using :func:`~.assignSecstr` 
function.

**[†]** Alternate locations are parsed as alternate coordinate sets. This
keyword will work for alternate location specified by "A". This to work for
alternate locations indicated by other letters, they must be parsed 
specifically by passing the identifier to the :func:`~.parsePDB`.

**[‡]** Atoms with unspecified alternate location/chain/segment/icode/secondary 
structure identifiers can be selected using "_". This character is replaced 
with a whitespace.

**[§]** If there are multiple residues with the same number but 
distinguished with insertion codes, the insertion code can be appended
to the residue number. "_" stands for empty insertion code. For example:
    
  * ``'resnum 5'`` selects residue 5 (all insertion codes)
  * ``'resnum 5A'`` selects residue 5 with insertion code A
  * ``'resnum 5_'`` selects residue 5 with no insertion code

**[¶]** Distinct residues, chains, and segments can be selected using 
*resindex*, *chindex*, and *segindex* keywords, respectively.  Unique
numbers to these entitites are assigned by :class:`~.HierView` class
upon building of a hierarchical view for an :class:`~.AtomGroup`.
Note that hierarchical views are build automatically when needed.

**Strings (with special characters)**

Strings can be any combination of the following::

  abcdefghijklmnopqrstuvwxyz
  ABCDEFGHIJKLMNOPQRSTUVWXYZ
  0123456789
  ~@#$.:;_',
  
For example ``"name C' N` O~ C$ C#"`` is a valid selection string. 


**Integers and floats**

Numbers can be provided as integers or floats, and they will be converted to
appropriate type. For example ``'resnum 10 11.0'`` will select residues
with number 10 and 11, but ``'resnum 10.5'`` will not select anything.

Negative numbers must be entered between grave accent symbols, 
e.g. ``'resnum `-3`'``

**Number ranges**

Number ranges can be passed as follows:
    
  * ``'resnum 5 10 to 15'`` selects residues 5, 10, 11, 12, 13, 14, and 15
  * ``'resnum 5 10:15'`` selects residues 5, 10, 11, 12, 13, and 14 
    (:, colon, works as it does in Python slicing operations)
  * ``'resnum 1:10:2'`` selects residues 1, 3, 5, 7, and 9
  * ``'x 1 to 10'`` selects atoms whose x coordinates are greater or equal to 1
    or smaller or equal to 10  
  * ``'x 1:10'`` selects atoms whose x coordinates are greater or equal to 1
    or smaller or equal to 10
    
Number ranges involving negative numbers must be entered between grave accent 
symbols, e.g. ``'resnum `-3 to 10`'``, ``'resnum `-3:10:2`'``

**More special characters (``)**

Strings can include the following characters (including whitespace) as well 
when they are surrounded by grave accent character (``):
  
  ~!@#$%^&*()-_=+[{}]\|;:,<>./?()'"

For example ``'name `CA#` `C #`'`` will work.

**Regular expressions ("")**

Strings surrounded by double quotes ("") will be treated as regular 
expressions. The following character set can be used between double 
quotes:
  
  ~!@#$%^&*()-_=+[{}]\|;:,<>./?()'`

For example ``'resname "A.."'`` will select residues whose names start with 
letter A and are three-characters long.

For more information on regular expressions see :mod:`re`. 

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import re as RE
from types import NoneType

import numpy as np
from numpy import ndarray, ones, zeros, invert, unique, concatenate

import pyparsing as pp
pp.ParserElement.enablePackrat()

from prody import LOGGER, SETTINGS

from atomic import Atomic
from fields import ATOMIC_ATTRIBUTES, ATOMIC_FIELDS

from atomgroup import AtomGroup 
from chain import Chain, getSequence, AAMAP
from pointer import AtomPointer
from selection import Selection
from segment import Segment
from atommap import AtomMap

from prody.tools import rangeString
from prody.KDTree import getKDTree

DEBUG = 0

__all__ = ['Select', 'SelectionError', 'TypoWarning',
           'getKeywordResnames', 'setKeywordResnames',
           'getBackboneAtomNames', 'setBackboneAtomNames',
           'getBackboneAtomNames', 'setBackboneAtomNames',
           'getAtomNameRegex', 'setAtomNameRegex',
           'defSelectionMacro', 'delSelectionMacro', 'getSelectionMacro',
           'getReservedWords']

KEYWORDS_STRING = set(['name', 'type', 'resname', 'chain', 'element', 
                       'segment', 'altloc', 'secondary', 'icode',
                       'chid', 'secstr', 'segname', 'sequence'])

ALNUM_VALLEN = {}
for key, field in ATOMIC_FIELDS.iteritems():
    if isinstance(field.dtype, str) and field.dtype.startswith('|S'):
        itemsize = np.dtype(field.dtype).itemsize
        ALNUM_VALLEN[key] = itemsize
        if field.synonym:
            ALNUM_VALLEN[field.synonym] = itemsize

KEYWORDS_INTEGER = set(['serial', 'index', 'resnum', 'resid', 
                        'segindex', 'chindex', 'resindex'])
KEYWORDS_FLOAT = set(['x', 'y', 'z', 'beta', 'mass', 'occupancy', 'mass', 
                      'radius', 'charge'])
KEYWORDS_NUMERIC = KEYWORDS_FLOAT.union(KEYWORDS_INTEGER)    

KEYWORDS_VALUE_PAIRED = KEYWORDS_NUMERIC.union(KEYWORDS_STRING)
KEYWORDS_SYNONYMS = {}
for key, field in ATOMIC_FIELDS.iteritems(): 
    if field.synonym:
        KEYWORDS_SYNONYMS[field.synonym] = key
ATOMIC_ATTRIBUTES = ATOMIC_ATTRIBUTES
# 21st and 22nd amino acids	    3-Letter	1-Letter
# Selenocysteine	            Sec	        U
# Pyrrolysine	                Pyl	        O

# Ambiguous Amino Acids	                3-Letter	1-Letter
# Asparagine or aspartic acid	        Asx	        B
# Glutamine or glutamic acid	        Glx	        Z
# Leucine or Isoleucine	                Xle	        J
# Unspecified or unknown amino acid     Xaa         X

# Phosphorylated amino acids           3-Letter	    1-Letter
# Phosphothreonine                     TPO          T
# O-Phosphotyrosine                    PTR          Y
# Phosphoserine                        SEP          S

# Other amino acids           3-Letter	    1-Letter
# S-Hydroxycysteine           CSO           C

KEYWORD_RESNAMES = {
    'protein': ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 
                'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
                'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD', 'HSE', 'HSP', 
                'GLX', 'ASX', 'SEC', 'PYL', 'XLE', 'CSO', 'TPO', 'PTR', 'SEP'],
    'nucleic': ['GUA', 'ADE', 'CYT', 'THY', 'URA', 'DA', 'DC', 'DG', 'DT', 
                'A', 'C', 'G', 'T', 'U'],

    'acidic': ['ASP', 'GLU', 'TPO', 'PTR', 'SEP'],
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

Keywords arguments
===============================================================================

Below is the list of keywords defined based on residue type and/or property.
These definitions can be retrieved or altered using :func:`getKeywordResnames` 
and :func:`setKeywordResnames`, respectively.

=============  ================================================================
Keyword        Description
=============  ================================================================
"""
keys = KEYWORD_RESNAMES.keys()
keys.sort()
from textwrap import wrap
for key in keys:
    lines = wrap('resname ' + ' '.join(KEYWORD_RESNAMES[key]))
    if key in KEYWORD_RESNAMES_READONLY:
        __doc__ += '{0:13s}  {1:s}\n'.format(key + ' [#]', lines[0])
    else:
        __doc__ += '{0:13s}  {1:s}\n'.format(key, lines[0])
    for line in lines[1:]:
        __doc__ += '{0:13s}  {1:s}\n'.format('', line)

__doc__ += """\
=============  ================================================================

**[#]** Definitions of these keywords are based on others and cannot be changed
directly: 
    
"""
keys = KEYWORD_RESNAMES_READONLY.keys()
keys.sort()
for key in keys:
    __doc__ += '  * ``{0:s}`` is ``{1:s}``\n'.format(
        repr(key), repr(KEYWORD_RESNAMES_READONLY[key]))

__doc__ += """

Following are additional keywords whose definitions are more restricted:

===============  ==============================================================
Keyword          Description
===============  ==============================================================
all              all atoms
none             nothing (returns **None**)
hetero           non-protein/nucleic atoms, same as
                 ``'not (protein or nucleic)'``
calpha (ca)      Cα atoms of protein residues, same as 
                 ``'name CA and protein'``
backbone (bb)    backbone atoms of protein residues, same as
                 ``'name CA C O N and protein'``
backbonefull     backbone atoms of protein residues, same as
                 ``'name CA C O N H H1 H2 H3 OXT and protein'``
bbful            same as ``'backbonefull'`` 
sidechain (sc)   side-chain atoms of protein residues, same as
                 ``'not name CA C O N H and protein'``
carbon           carbon atoms, same as ``'name "C.*" and not resname ion'``
hydrogen         hydrogen atoms, same as ``'name "[1-9]?H.*"'``
noh              non hydrogen atoms, same as ``'not name "[1-9]?H.*"'``
nitrogen         nitrogen atoms, same as ``'name "N.*"'``
oxygen           oxygen atoms, same as ``'name "O.*"'``
sulfur           sulfur atoms, same as ``'name "S.*"'``
extended         residue in extended conformation, same as ``'secondary E'``
helix            residue in α-helix conformation, same as ``'secondary H'``
helix_3_10       residue in 3_10-helix conformation, same as ``'secondary G'``
helix_pi         residue in π-helix conformation, same as ``'secondary I'``
turn             residue in hydrogen bonded turn conformation, same as
                 ``'secondary T'``
bridge           residue in isolated beta-bridge conformation, same as
                 ``'secondary B'``
bend             residue in bend conformation, same as ``'secondary S'``
coil             residue not in one of above conformations, same as
                 ``'secondary C'``
===============  ==============================================================

Among these list of backbone atom names can be changed using 
:func:`setBackboneAtomNames`  and regular expressions for element types
can be changed using :func:`setAtomNameRegex`.

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

SRE_Pattern = type(KEYWORD_NAME_REGEX['carbon'])

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
===============================================================================

Following keywords can be used in numerical comparisons, as operands of 
arithmetic operations or as arguments to functions:  
 
 * index, serial   
 * resnum, resid
 * x, y, z
 * beta, occupancy
 * charge, mass, radius (these must be set by the user before they can be used) 

Numerical attributes of atoms can be used with the following comparison 

==========  =================================
Comparison  Description
==========  =================================
   <        less than
   >        greater than
   <=       less than or equal
   >=       greater than or equal
   ==       equal
   =        equal
   !=       not equal
==========  =================================

*Examples:* ``'x < 0'``, ``'occupancy != 1'``

Numerical attributes of atoms can be used as operands to the following 
operators:

=========  ==================================
Operation  Description
=========  ==================================
x ** y     x to the power y
x ^ y      x to the power y
x * y      x times y
x / y      x divided by y
x // y     x divided by y (floor division)
x % y      x modulo y
x + y      x plus y 
x - y      x minus y
=========  ==================================
   
These operations must be used with a numerical comparison, e.g. 
``'x ** 2 < 10'``, ``'x ** 2 ** 2 < 10'``, ``'occupancy != 1'``
   
Numerical attributes of atoms can be used as arguments to the following 
functions:
   
========  ===================================
Function  Description
========  ===================================
abs(x)    absolute value of x 
acos(x)   arccos of x
asin(x)   arcsin of x
atan(x)   arctan of x
ceil(x)   smallest integer not less than x
cos(x)    cosine of x
cosh(x)   hyperbolic cosine of x
floor(x)  largest integer not greater than x 
exp(x)    e to the power x
log(x)    natural logarithm of x
log10(x)  base 10 logarithm of x
sin(x)    sine of x
sinh(x)   hyperbolic sine of x
sq(x)     square of x
sqrt(x)   square-root of x
tan(x)    tangent of x
tanh(x)   hyperbolic tangent of x
========  ===================================

**Examples**
  
  * ``'sqrt(x**2 + y**2 + z**2) < 10'`` selects atoms within 10 Å of the 
    origin
  * ``'resnum <= 100'`` selects atoms with residue numbers less than or equal 
    to 100  


Distance based selections
===============================================================================

Atoms within a user specified distance (Å) from a set of user specified atoms
can be selected using ``within . of ..`` keyword, e.g. ``within 5 of water``
selects atoms that are within 5 Å of water molecules. This setting will
results selecting water atoms as well.

User can avoid selecting specified atoms using ``exwithin . of ..`` setting,
e.g. ``exwithin 5 of water`` will not select water molecules and is equivalent
to ``within 5 of water and not water``

Sequence selections
===============================================================================

One-letter amino acid sequences can be used to make atom selections. 
``'sequence SAR'`` will select **SER-ALA-ARG** residues in a chain.  Note
that the selection does not consider connectivity within a chain.  Regular 
expressions can also be used to make selections: ``'sequence S..R'`` will
select **SER-XXX-XXX-ARG** pattern, if  present. 
    

Expanding selections
===============================================================================

A selection can be expanded to include the atoms in the same *residue*, 
*chain*, or *segment* using ``same .. as ..`` setting, e.g.
``same residue as exwithin 4 of water`` will select residues that have
at least an atom within 4 Å of any water molecule.

Additionally, a selection may be expanded to the immediately bonded atoms using
``bonded to ...`` method, e.f. ``bonded to calpha`` will select atoms bonded 
to Cα.  For this to work, bonds must be set by the user using :meth:`.AtomGroup
.setBonds` method.  It is also possible to select bonded atoms by excluding the
atoms from which the bonds will originate, i.e. ``exbonded to ...``.  

Selection macros
===============================================================================

Any valid selection string can be used to define selection macros using the 
:func:`defSelectionMacro` function.  Macros are saved in ProDy configuration 
and loaded in later sessions automatically.  Below functions are for 
manipulating selection macros:
    
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

ATOMGROUP = None

MACROS = SETTINGS.get('selection_macros', {})

def isMacro(word):
    
    return word in MACROS

def areAllStrings(container):
    """Return ``True`` if all items in *container* are instances of 
    :func:`str`."""
    
    for item in container: 
        if not isinstance(item, str):
            return False
    return True

def defSelectionMacro(name, selstr):
    """Define selection macro *selstr* with name *name*.  Both *name* and 
    *selstr* must be string.  An existing keyword cannot be used as a macro 
    name. If a macro with given *name* exists, it will be overwritten.
    
    >>> defSelectionMacro('cbeta', 'name CB and protein')"""
    
    if not isinstance(name, str) or not isinstance(selstr, str):
        raise TypeError('both name and selstr must be strings')
    elif isKeyword(name):
        raise ValueError("{0:s} is an existing keyword, cannot be used as a "
                         "macro name".format(repr(name)))
    elif not (name.isalpha() and name.islower()):
        raise ValueError('macro names must be all lower case letters, {0:s} '
                         'is not a valid macro name'.format(repr(name)))
    
    LOGGER.info('Testing validity of selection string:')
    try:
        ATOMGROUP.select(selstr)
    except SelectionError:
        LOGGER.warn('{0:s} is not a valid selection string, macro {1:s} is not'
                    'defined.'.format(repr(selstr), repr(name)))
    else:
        LOGGER.info("Macro {0:s} is defined as {1:s}."
                    .format(repr(name), repr(selstr)))
        MACROS[name] = selstr
        SETTINGS['selection_macros'] = MACROS
        SETTINGS.save()

def delSelectionMacro(name):
    """Delete the macro *name*.
    
    >>> delSelectionMacro('cbeta')"""
    
    try:
        MACROS.pop(name)
    except:
        LOGGER.warn("Macro {0:s} is not found.".format(repr(name)))
    else:
        LOGGER.info("Macro {0:s} is deleted.".format(repr(name)))
        SETTINGS['selection_macros'] = MACROS
        SETTINGS.save()

def getSelectionMacro(name=None):
    """Return the definition of the macro *name*.  If *name* is not given, 
    returns a copy of the selection macros dictionary."""
    
    if name is None:        
        return MACROS.copy()
    try:
        return MACROS[name]
    except KeyError:
        LOGGER.info("{0:s} is not a user defined macro name."
                    .format(repr(name)))

mapField2Var = {}
for field in ATOMIC_FIELDS.values():
    mapField2Var[field.name] = field.var

def getKeywordResnames(keyword):
    """Return residue names associated with a keyword.
    
    >>> getKeywordResnames('acidic')
    ['ASP', 'GLU', 'PTR', 'SEP', 'TPO']"""
    
    assert isinstance(keyword, str), 'keyword must be a string instance'
    try:
        resnames = KEYWORD_RESNAMES[keyword]
        resnames.sort()
        return resnames  
    except KeyError:
        if keyword in KEYWORD_RESNAMES_READONLY:
            LOGGER.warn('{0:s} is defined as {1:s}'.format(repr(keyword), 
                                    repr(KEYWORD_RESNAMES_READONLY[keyword])))
        else:
            LOGGER.warn("{0:s} is not a keyword".format(repr(keyword)))

def setKeywordResnames(keyword, resnames):
    """Change the list of residue names associated with a keyword.  *keyword* 
    must be a string, and *resnames* may be a list, tuple, or set of strings. 
    The existing list of residue names will be overwritten with the given 
    residue names.  Note that changes in keyword definitions are not saved 
    permanently.
    
    >>> setKeywordResnames('acidic', ['ASP', 'GLU', 'PTR', 'SEP', 'TPO'])"""
    
    if not isinstance(keyword, str):
        raise TypeError('keyword must be a string')
    if not isinstance(resnames, (list, tuple, set)):
        raise TypeError('resnames must be a list, set, or tuple')
    if not areAllStrings(resnames):
        raise TypeError('all items in resnames must be string instances')
    
    if keyword in KEYWORD_RESNAMES_READONLY:
        LOGGER.warn("{0:s} is defined as {1:s} and cannot be changed directly"
            .format(repr(keyword), repr(KEYWORD_RESNAMES_READONLY[keyword])))
        return
    if keyword in KEYWORD_RESNAMES:
        for rn in resnames:
            if not isinstance(rn, str):
                raise TypeError('all items in resnames must be strings')
        KEYWORD_RESNAMES[keyword] = list(set(resnames))
        _setReadonlyResidueNames()
    else:
        raise ValueError("{0:s} is not a valid keyword".format(repr(keyword)))

def getAtomNameRegex(name):
    """Return regular expression used for selecting common elements.
    
    >>> getAtomNameRegex('nitrogen')
    'N.*'"""
    
    assert isinstance(name, str), 'name must be a string instance'
    try:
        return KEYWORD_NAME_REGEX[name].pattern   
    except KeyError:
        LOGGER.warn('{0:s} is not a valid element'.format(name))

def setAtomNameRegex(name, regex):
    """Set regular expression used for selecting common elements.  Note that 
    changes in keyword definitions are not saved permanently.
    
    >>> setAtomNameRegex('nitrogen', 'N.*')"""
    
    assert isinstance(name, str), 'name must be a string instance'
    if not name in KEYWORD_NAME_REGEX:
        raise ValueError("{0:s} is not a valid keyword".format(repr(name)))
    if not isinstance(regex, str):
        raise TypeError("regex must be a string instance")
    try:
        regex = RE.compile(regex)
    except:
        raise ValueError("{0:s} is not a valid regular expression"
                         .format(repr(regex)))
    else:
        KEYWORD_NAME_REGEX[name] = regex

def getBackboneAtomNames(full=False):
    """Return protein backbone atom names.  ``full=True`` argument returns 
    atom names for *backbonefull* keyword.
    
    >>> getBackboneAtomNames()
    ['C', 'CA', 'N', 'O']"""
    
    assert isinstance(full, bool), 'full must be a boolean instance'
    if full:
        bban = list(BACKBONE_FULL_ATOM_NAMES)
    else:
        bban = list(BACKBONE_ATOM_NAMES)
    bban.sort()
    return bban 

def setBackboneAtomNames(backbone_atom_names, full=False):
    """Set protein backbone atom names.  Atom names for *backbonefull* keyword 
    can be set by passing ``full=True`` argument.  Note that changes in keyword
    definitions are not saved permanently."""
    
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


class SelectionError(Exception):    
    
    """Exception raised when there are errors in the selection string."""
    
    def __init__(self, sel, loc=0, msg=''):
        
        msg = ("An invalid selection string is encountered:\n{0:s}\n"
               .format(repr(sel)) + 
               ' ' * (loc + 1) + '^ ' + msg)
        Exception.__init__(self, msg)

class TypoWarning(object):
    
    """A class used for issuing warning messages when potential typos are
    detected in a selection string.  Warnings are issued to ``sys.stderr`` 
    via ProDy package logger. :func:`~.confProDy` function can be used to 
    turn typo warnings *on* or *off*, e.g. ``confProDy(typo_warnings=False)``.
    """
    
    def __init__(self, sel='', loc=0, msg='', typo=None):

        if SETTINGS['typo_warnings']:        
            shift = sel.find(typo, loc)
            msg = ("{0:s} might contain typo(s)\n".format(repr(sel)) +
                   ' ' * (shift + 12) + '^ ' + msg)
            LOGGER.warn(msg)

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



RESERVED = set(ATOMIC_FIELDS.keys() + ATOMIC_ATTRIBUTES.keys() +
               ['and', 'or', 'not', 'within', 'of', 'exwithin', 'same', 'as',
                'bonded', 'exbonded', 'to'] +
               KEYWORDS_SYNONYMS.keys() + 
               ['n_atoms', 'n_csets', 'cslabels', 'title', 'coordinates',
                'bonds', 'bmap', 'numbonds'])

def isReserved(word):
    return (word in RESERVED or isKeyword(word) or word in FUNCTION_MAP)
        
        
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

SAMEAS_MAP = {'residue': 'resindex', 'chain': 'chindex', 'segment': 'segindex'}

def expandBoolean(keyword):
    
    if keyword in KEYWORD_MAP:
        (resnames, rn_invert, names, an_invert) = KEYWORD_MAP[keyword]
        tokens = []
        if names is not None:
            if an_invert:
                tokens.append(NOT)
            tokens.append('name')
            tokens.extend(names)
            if resnames is not None:
                tokens.append(AND)
        if resnames is not None:
            if rn_invert:
                tokens.append(NOT)
            tokens.append('resname')
            tokens.extend(resnames)
        return tokens
    elif keyword in SECSTR_MAP:
        return ['secondary', SECSTR_MAP[keyword]]
    else:
        return keyword

def splitList(alist, sep):
    """Return list of lists obtained by splitting *alist* at the position of 
    *sep*."""
    
    result = [[]]
    for item in alist:
        if item == sep:
            result.append([])
        else:
            result[-1].append(item)
    return result

class Select(object):

    """Select subsets of atoms based on a selection string.  
    
    See :mod:`~.select` module documentation for detailed documentation.  
    Definitions of single word keywords, such as *protein*, *backbone*, 
    or *polar*, etc., may be altered using functions in :mod:`~.select` 
    module. 
    
    This class makes use of |pyparsing| module."""

    def __init__(self):
        
        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._selstr = None
        
        self._coords = None
        self._kwargs  = None
        # set True when selection string alone cannot reproduce the selection
        self._ss2idx = False  
        self._data = dict()
        self._replace = False
        
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
        def regularExpParseAction(sel, loc, token): 
            token = token[0]
            if len(token[0]) == 2:
                return RE.compile('^()$')
            else:
                try:
                    regexp = RE.compile(token[1])
                except:
                    raise SelectionError(sel, loc, 'failed to compile regular '
                                    'expression {0:s}'.format(repr(token[1])))
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
               pp.Regex('(ex)?bonded to') |
               pp.Regex('(ex)?bonded [0-9]+ to') |
               pp.Regex('same [a-z]+ as') | 
               pp.Regex('(ex)?within [0-9]+\.?[0-9]* of'), 
                        1, pp.opAssoc.RIGHT, self._unary),
              (pp.Keyword(AND), 2, pp.opAssoc.LEFT, self._and),
              (pp.Keyword(OR), 2, pp.opAssoc.LEFT, self._or),]
            )

        self._tokenizer.setParseAction(self._defaultAction)
        self._tokenizer.leaveWhitespace()
        
        
    def _reset(self):

        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._coords = None
        self._data.clear()
        
    def getBoolArray(self, atoms, selstr, **kwargs):
        """Return a boolean array with ``True`` values for *atoms* matching 
        *selstr*.
        
        .. note:: The length of the boolean :class:`numpy.ndarray` will be
           equal to the number of atoms in *atoms* argument."""
        
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance, not {0:s}'
                            .format(type(atoms)))
        elif not isinstance(selstr, str):
            raise TypeError('selstr must be a string, not a {0:s}'
                            .format(type(selstr)))
        self._reset()

        if any([isReserved(key) for key in kwargs.iterkeys()]):   
            for key in kwargs.iterkeys():
                if isReserved(key):
                    raise SelectionError(selstr, selstr.find(key), "{0:s} "
                        "is a reserved word and cannot be used as a keyword "
                        "argument".format(repr(key)))

        if isinstance(atoms, AtomGroup): 
            self._ag = atoms
            self._atoms = atoms
            self._indices = None
        else:
            self._ag = atoms.getAtomGroup()
            self._indices = atoms._getIndices()
            if isinstance(atoms, AtomMap):
                self._atoms = Selection(self._ag, self._indices, '', )
                self._atoms._indices = self._indices
            else: 
                self._atoms = atoms
        self._n_atoms = atoms.numAtoms()
        self._selstr = selstr
        self._kwargs = kwargs
        
        if DEBUG:
            print('getBoolArray', selstr)
        torf = self._evalSelstr()
        if not isinstance(torf, ndarray):
            if DEBUG: print(torf)
            raise SelectionError(selstr)
        elif torf.dtype != bool:
            if DEBUG:
                print('_select torf.dtype', torf.dtype, isinstance(torf.dtype, 
                                                                   bool))
            raise SelectionError(selstr)
        if DEBUG:
            print('_select', torf)
        return torf
    
    def getIndices(self, atoms, selstr, **kwargs):
        """Return indices of atoms matching *selstr*.
        
        .. note:: The indices correspond to indices in *atoms* argument.  When
           *atoms* is not an :class:`~.AtomGroup` instance, indexing the 
           :class:`~.AtomGroup` may return a different set of atoms."""
        
        torf = self.getBoolArray(atoms, selstr, **kwargs)        
        return torf.nonzero()[0]
        
    def select(self, atoms, selstr, **kwargs):
        """Return a subset of atoms matching *selstr* as a :class:`Selection`.
        
        :arg atoms: atoms to be evaluated    
        :type atoms: :class:`~.Atomic`
        
        :arg selstr: selection string
        :type selstr: str
        
        If type of *atoms* is :class:`~.AtomMap`, an :class:`~.AtomMap` 
        instance is returned. Otherwise, :class:`~.Selection` instances 
        are returned.

        .. note:

            * If selection string does not match any atoms, ``None`` is 
              returned.
              
            * :meth:`select` accepts arbitrary keyword arguments which enables 
              identification of intermolecular contacts. See :ref:`contacts` 
              for details.
        
            * A special case for making atom selections is passing an
              :class:`~.AtomMap` instance as *atoms* argument.  Dummy 
              atoms will not be included in the result but, the order 
              of atoms will be preserved."""
        
        self._ss2idx = False
        self._replace = False
        
        indices = self.getIndices(atoms, selstr, **kwargs)
        
        if not isinstance(atoms, AtomGroup):
            indices = self._indices[indices]
            
        ag = self._ag

        self._kwargs = None

        if len(indices) == 0:
            return None
            
        elif isinstance(atoms, AtomMap):
            return AtomMap(ag, indices, np.arange(len(indices)), 
                     np.array([]), 'Selection {0:s} from AtomMap {1:s}'
                    .format(repr(selstr), atoms.getTitle()), 
                            atoms.getACSIndex())
        else:
            if self._ss2idx:
                selstr = 'index {0:s}'.format(rangeString(indices))
            else:
                if self._replace:
                    for key, value in kwargs.iteritems():
                        if (isinstance(value, AtomPointer) and
                            not isinstance(value, AtomMap) and
                            key in selstr):
                            selstr = selstr.replace(key, 
                                                '(' + value.getSelstr() + ')')
                if isinstance(atoms, AtomPointer):
                    selstr = '({0:s}) and ({1:s})'.format(selstr, 
                                                      atoms.getSelstr())
            
            return Selection(ag, indices, selstr, atoms.getACSIndex(),
                             unique=True)
        
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
        if len(selstr.split()) == 1 and selstr not in MACROS:
            if isBooleanKeyword(selstr):
                return self._evalBoolean(self._selstr, 0, selstr)
            elif self._ag.isData(selstr):
                return self._evalUserdata(self._selstr, 0, selstr)
            elif isValuePairedKeyword(selstr):
                raise SelectionError(selstr, 0, 'must be followed by values')
            else:
                raise SelectionError(selstr, 0, 'is not a valid selection or '
                                     'user data label')
        
        selstr = self._prepareSelstr()
        try:
            if DEBUG: print('_evalSelstr using Pyparsing')
            tokens = self._tokenizer.parseString(selstr, 
                                             parseAll=True).asList()
            if DEBUG: print('_evalSelstr', tokens)
            return tokens[0]
        except pp.ParseException as err:
            raise SelectionError(self._selstr, err.column, 
                                 'parsing failed here')
    
    def _isValid(self, token):
        """Check the validity of part of a selection string. Expects a Python
        :func:`str` or a :func:`list`."""
        
        if DEBUG: print('_isValid', token)
        
        if isinstance(token, str):
            return (isBooleanKeyword(token) or
                   self._atoms.getDataType(token) == bool)   
        elif isinstance(token, list):
            tkn = token[0]
            return (isBooleanKeyword(tkn) or
                    isValuePairedKeyword(tkn) or tkn in self._kwargs or
                    self._atoms.isData(tkn) or tkn == NOT)
        return False
        
    def _defaultAction(self, sel, loc, tkns):
        """Called at the end of parsing action, when evaluating a simple 
        selection string, e.g. ``'index 5'``, or part of a selection string 
        in parentheses, e.g. for ``'index 5'`` in ``'within 5 of (index 5)'``.
        """
        
        if DEBUG: print('_defaulAction', tkns)
        
        if isinstance(tkns[0], ndarray):
            return tkns[0]
        torf = self._evaluate(sel, loc, tkns)
        if isinstance(torf, SelectionError):
            raise torf
        return torf
    
    def _evaluate(self, sel, loc, tkns, evalonly=None):
        """Evaluates statements in a selection string, e.g. ``'calpha'``,
        ``'index 5'``."""
        
        if DEBUG: print('_evaluate', tkns)
        
        if isinstance(tkns, str):
            if isBooleanKeyword(tkns):
                return self._evalBoolean(sel, loc, tkns, evalonly=evalonly, 
                                         rtrn=True)
            elif self._ag.isData(tkns):
                return self._evalUserdata(sel, loc, tkns, evalonly=evalonly)
            else:
                return SelectionError(sel, loc, "{0:s} is not understood"
                                      .format(repr(tkns)))
        elif isinstance(tkns, (ndarray, float)):
            return tkns
    
        keyword = tkns[0]
        if len(tkns) == 1:
            if isBooleanKeyword(keyword):
                return self._evalBoolean(sel, loc, keyword, evalonly=evalonly, 
                                         rtrn=True)
            elif isNumericKeyword(keyword):
                return self._evalNumeric(sel, loc, keyword)
            elif self._ag.isData(keyword):
                return self._evalUserdata(sel, loc, keyword, evalonly=evalonly)
            elif isinstance(keyword, ndarray):
                return keyword
            elif keyword in self._kwargs:
                arg = self._kwargs[keyword]
                if (isinstance(arg, AtomPointer) and
                    arg.getAtomGroup() is self._ag and
                    not isinstance(arg, AtomMap)):
                    torf = zeros(self._ag.numAtoms(), bool)
                    torf[arg._getIndices()] = True
                    if self._indices is not None:
                        torf = torf[self._indices]
                    self._replace = True
                    if DEBUG: print('_evaluate', keyword, torf)
                    return torf
                return keyword
                    
            else:
                try:
                    return float(keyword)
                except ValueError:
                    pass
        elif isAlnumKeyword(keyword):
            return self._evalAlnum(sel, loc, keyword, tkns[1:], 
                                   evalonly=evalonly)
        elif keyword in ('resnum', 'resid'):
            return self._resnum(sel, loc, tkns[1:], evalonly=evalonly)
        elif keyword == 'index':
            return self._index(sel, loc, tkns[1:], evalonly=evalonly)
        elif keyword == 'serial':
            return self._serial(sel, loc, tkns[1:], evalonly=evalonly)
        elif isNumericKeyword(keyword):
            return self._evalFloat(sel, loc, keyword, tkns[1:], 
                                   evalonly=evalonly)
        elif keyword == NOT:
            return self._not(sel, loc, tkns, evalonly=evalonly)
        elif self._ag.isData(keyword):
            return self._evalUserdata(sel, loc, keyword, tkns[1:], 
                                      evalonly=evalonly)
        
        return SelectionError(sel, loc, "{0:s} is not understood"
                              .format(repr(' '.join(tkns))))

    def _or(self, sel, loc, tokens):
        """Evaluate statements containing ``'or'`` operator."""
        
        if DEBUG: print('_or\n_or tokens '+str(tokens))
        previous = None
        evalonly = None
        selection = None
        
        tokens = splitList(tokens[0], OR)
        n_tokens = len(tokens) - 1
        
        arrays = []
        for i in xrange(n_tokens, -1, -1):
            torf = tokens[i][0]
            if isinstance(torf, ndarray):
                tokens.pop(i)
                arrays.append(torf)
        if arrays:
            if len(arrays) == 1:
                selection = arrays[0]
            else:
                selection = np.any(arrays, 0, arrays[0])
            if tokens:
                evalonly = invert(selection).nonzero()[0]
        
        n_tokens = len(tokens) - 1
        for i, tokens in enumerate(tokens):
            if not self._isValid(tokens):
                raise SelectionError(sel, loc, "{0:s} is not a valid usage"
                         .format(repr(' '.join([str(tkn) for tkn in tokens]))))
            torf = self._evaluate(sel, loc, tokens, evalonly=evalonly)
            if isinstance(torf, SelectionError):
                raise torf
           
            if evalonly is None:
                selection = torf
                evalonly = invert(selection).nonzero()[0]
            else:
                selection[evalonly[torf]] = True
                if i < n_tokens:
                    evalonly = evalonly[invert(torf, torf).nonzero()]
        return selection

    def _and(self, sel, loc, tokens, rtrn=False):
        """Evaluate statements containing ``'and'`` operator."""
        
        if DEBUG: print('_and\n_and tokens '+str(tokens))
        evalonly = None
        previous = None
        selection = None
        
        tokens = splitList(tokens[0], AND)
        n_tokens = len(tokens) - 1
        
        arrays = []
        for i in xrange(n_tokens, -1, -1):
            torf = tokens[i][0]
            if isinstance(torf, ndarray):
                tokens.pop(i)
                arrays.append(torf)
        if arrays:
            if len(arrays) == 1:
                selection = arrays[0]
            else:
                selection = np.all(arrays, 0, arrays[0])
            if tokens:
                evalonly = selection.nonzero()[0]

        n_tokens = len(tokens) - 1
        for i, token in enumerate(tokens):
            
            if not self._isValid(token):
                raise SelectionError(sel, loc, "{0:s} is not a valid usage"
                        .format(repr(' '.join([str(tkn) for tkn in token]))))
            torf = self._evaluate(sel, loc, token, evalonly=evalonly)
            if isinstance(torf, SelectionError):
                if rtrn:
                    return torf
                else:
                    raise torf
            if evalonly is None:
                selection = torf
                evalonly = selection.nonzero()[0]
            else:
                selection[evalonly] = torf
                if i < n_tokens:
                    evalonly = evalonly[torf]
        return selection
    
    def _unary(self, sel, loc, tokens):
        """Perform the unary operation."""

        if DEBUG: print('_unary', tokens)
        tokens = tokens[0]
        what = tokens[0]
        if what == NOT:
            torf = self._not(sel, loc, tokens)
        elif what.startswith('same'):
            torf = self._sameas(sel, loc, tokens)
        elif what.endswith('to'):
            torf = self._bondedto(sel, loc, tokens, what.startswith('ex'))
        else:
            torf = self._within(sel, loc, tokens, what.startswith('ex'))
        if isinstance(torf, SelectionError):
            raise torf
        if torf is None:
            raise SelectionError(sel, loc)
        return torf

    def _not(self, sel, loc, tokens, evalonly=None):
        """Negate selection."""
        
        if DEBUG: print('_not', tokens)
        if isinstance(tokens[1], ndarray):
            torf = tokens[1]
        else:
            torf = self._evaluate(sel, loc, tokens[1:], evalonly=evalonly)
            if isinstance(torf, (SelectionError, NoneType)):
                return torf
        invert(torf, torf)
        return torf
    
    def _within(self, sel, loc, tokens, exclude):
        """Perform distance based selection."""

        if DEBUG: print('_within', tokens)
        within = tokens[0].split()[1]
        try:
            within = float(within)
        except:
            return SelectionError(sel, loc, '{0:s} must be a number'
                                  .format(repr(within)))
        if DEBUG: print('_within', within)
        which = tokens[1]
        if not isinstance(which, ndarray):
            which = self._evaluate(sel, loc, tokens[1:])

        if DEBUG: print('_within', which)
        other = False
        if isinstance(which, str) and which in self._kwargs:
            coords = self._kwargs[which]
            if DEBUG: print('_kwargs', which)
            if isinstance(coords, ndarray):
                if coords.ndim == 1 and len(coords) == 3:
                    coords = np.array([coords])
                elif not (coords.ndim == 2 and coords.shape[1] == 3):
                    return SelectionError(sel, loc + sel.find(which), 
                        "{0:s} is not a coordinate array".format(repr(which)))
            else:
                try:
                    coords = coords.getCoords()
                except:
                    return SelectionError(sel, loc + sel.find(which), 
                        "{0:s} does not have `getCoords` method"
                        .format(repr(which)))
                if not isinstance(coords, ndarray):
                    return SelectionError(sel, loc +  + sel.find(which), 
                        "coordinates of {0:s} are not set".format(repr(which)))
            exclude=False
            self._ss2idx = True
            which = np.arange(len(coords))
            other = True
        elif isinstance(which, ndarray) and which.dtype == np.bool: 
            which = which.nonzero()[0]
            coords = self._getCoords(sel, loc)
        else:
            return SelectionError(sel, loc)
        if isinstance(coords, SelectionError): 
            return coords

        if other or len(which) < 20:
            kdtree = self._atoms._getKDTree()
            get_indices = kdtree.get_indices
            search = kdtree.search
            torf = zeros(self._ag.numAtoms(), bool)
            for index in which:
                search(coords[index], within)
                torf[get_indices()] = True
            if self._indices is not None:
                torf = torf[self._indices]
            if exclude:
                torf[which] = False
        else:
            torf = ones(self._n_atoms, bool)
            torf[which] = False
            check = torf.nonzero()[0]
            torf = zeros(self._n_atoms, bool)
            
            cxyz = coords[check]
            kdtree = getKDTree(coords[which])
            get_indices = kdtree.get_indices
            search = kdtree.search
            select = []
            append = select.append
            for i, xyz in enumerate(cxyz):
                search(xyz, within)
                if len(get_indices()):
                    append(i)

            torf[check[select]] = True
            if not exclude:
                torf[which] = True
        return torf
    
    def _sameas(self, sel, loc, token):
        """Evaluate ``'same entity as ...'`` expression."""
        
        if DEBUG: print('_sameas', token)
        what = token[0].split()[1]
        
        which = token[1]
        if not isinstance(which, ndarray):
            which = self._evaluate(sel, loc, token[1:])
            if isinstance(which, SelectionError):
                return which
        
        self._ag.getHierView()
        xindex = SAMEAS_MAP.get(what)
        if xindex is None:
            return SelectionError(sel, loc, 'entity in "same entity as" must '
                                  'be one of "chain", "residue", or "segment"')
        
        indices = self._getData(sel, loc, xindex)
        if isinstance(indices, SelectionError):
            return indices
        return self._evalFloat(sel, loc, xindex, unique(indices[which]))
     
    def _bondedto(self, sel, loc, token, exclude):
        """Expand selection to immediately bonded atoms."""
        
        if DEBUG: print('_bondedto', token)
        
        keyword = token[0].split()
        if len(keyword) == 2:
            repeat = 1
        else: 
            repeat = int(keyword[1])
            if repeat == 0:
                return SelectionError(sel, loc + len('bonded '), 
                                      'number must be a greater than zero')

        which = token[1]
        if not isinstance(which, ndarray):
            which = self._evaluate(sel, loc, token[1:])
            if isinstance(which, SelectionError):
                return which
        bmap = self._ag._bmap
        if bmap is None:
            return SelectionError(sel, loc, 'bonds are not set')
        
        indices = self._indices
        if indices is not None:
            bmap = bmap[indices]
        n_atoms = self._ag.numAtoms()
        for i in xrange(repeat):
            torf = zeros(n_atoms, bool)
            bonded = unique(bmap[which])
            if bonded[0] == -1:
                torf[bonded[1:]] = True
            else: 
                torf[bonded] = True
            if indices is not None:
                torf = torf[indices]
            if exclude:
                torf[which] = False
            else:
                torf[which] = True
            if i + 1 < repeat: 
                which = torf.nonzero()[0]
            if DEBUG: print('_bondedto repeat', i+1, 'selected', len(which))
        return torf
    
    def _comp(self, sel, loc, tokens):
        """Perform numeric comparisons. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_comp', tokens)
        tokens = tokens[0]
        if len(tokens) >= 3 and len(tokens) % 2 != 1:
            raise SelectionError(sel, loc, 'each comparison must have exactly '
                                 'two operands')
        i = 1
        left = self._evalNumeric(sel, loc, tokens[0])
        if DEBUG: print('_comp left', left)
        if isinstance(left, SelectionError):
            raise left
        result = None
        while i < len(tokens): 
            comp = tokens[i]
            right = self._evalNumeric(sel, loc, tokens[i + 1])
            if DEBUG: print('_comp right', right)
            if isinstance(right, SelectionError):
                raise right
            if result is None:
                result = BINOP_MAP[comp](left, right)
            else:
                result *= BINOP_MAP[comp](left, right)
            left = right
            i += 2
        return result

    def _pow(self, sel, loc, tokens):
        """Perform power operation. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_pow', tokens)
        tokens = tokens[0]
        base = self._evalNumeric(sel, loc, tokens.pop(0))
        if isinstance(base, SelectionError):
            raise base
        power = self._evalNumeric(sel, loc, tokens.pop())
        if isinstance(power, SelectionError):
            raise power
        tokens.pop()
        while tokens:
            number = self._evalNumeric(sel, loc, tokens.pop()) 
            if isinstance(number, SelectionError):
                raise number
            power = number * power
            tokens.pop()
        return base ** power

    def _add(self, sel, loc, tokens):
        """Perform addition operations. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_add', tokens)
        tokens = tokens[0]
        left = self._evalNumeric(sel, loc, tokens.pop(0))
        if isinstance(left, SelectionError):
            raise left
        while tokens:
            op = tokens.pop(0)
            right = self._evalNumeric(sel, loc, tokens.pop(0))
            if isinstance(right, SelectionError):
                raise right
            left = BINOP_MAP[op](left, right)
        if DEBUG: print('_add total', left)
        return left
 
    def _mul(self, sel, loc, tokens):
        """Perform multiplication operations. Expected operands are numbers 
        and numeric atom attributes."""
        
        if DEBUG: print('_mul', tokens)
        tokens = tokens[0]
        left = self._evalNumeric(sel, loc, tokens[0])
        if isinstance(left, SelectionError):
            raise left
        i = 1
        while i < len(tokens):
            op = tokens[i]
            i += 1
            right = self._evalNumeric(sel, loc, tokens[i])
            if isinstance(right, SelectionError):
                raise right
            i += 1
            if op == '/' and right == 0.0: 
                raise SelectionError(sel, loc, 'zero division error')
            left = BINOP_MAP[op](left, right)
        return left
    
    def _evalNumeric(self, sel, loc, token):
        """Evaluate a number operand or a numeric keyword."""
        
        if DEBUG: print('_evalNumeric', token)
        if isinstance(token, (ndarray, float)):
            return token
        elif isFloatKeyword(token):
            return self._evalFloat(sel, loc, token)
        elif token in ('resnum', 'resid'):
            return self._resnum(sel, loc)
        elif token == 'index':
            return self._index(sel, loc)
        elif token == 'serial':
            return self._serial(sel, loc)
        elif isNumericKeyword(token): 
            return self._getData(sel, loc, token)
        elif self._ag.isData(token):
            data = self._getData(sel, loc, token)
            if data.dtype in (float, int):
                return data
            else:
                return SelectionError(sel, loc, "data type of {0:s} must be "
                                      "int or float".format(repr(token)))
                
        elif isinstance(token, SRE_Pattern):
            return SelectionError(sel, loc, 'regular expressions cannot be '
                                  'used as numeric values')
        else:
            try:
                token = float(token)
            except ValueError:
                pass
            else:
                return token
        return SelectionError(sel, loc, "{0:s} does not have a numeric value"
                              .format(repr(token)))
    
    def _sign(self, sel, loc, tokens):
        """Change the sign of a selection argument."""
        
        tokens = tokens[0]
        if DEBUG: print('_sign', tokens)
        
        if len(tokens) != 2:
            raise SelectionError(sel, loc, "sign operators (+/-) must be "
                             "followed by single keyword, e.g. '-x', '-beta'")
        arg = self._evalNumeric(sel, loc, tokens[1])
        if arg is None:
            raise SelectionError(sel, loc, "sign operators (+/-) must be "
                         "followed by a numeric keyword, e.g. '-x', '-beta'")
        if tokens[0] == '-':
            return -arg
        return arg

    def _func(self, sel, loc, tokens):
        """Evaluate functions used in selection strings."""
        
        if DEBUG: print('_func', tokens)
        tokens = list(tokens[0])
        
        if len(tokens) != 2:
            raise SelectionError(sel, loc, "functions (sin/abs/etc.) must have"
                    " a single numeric keyword as argument, e.g. 'sin(x)'")
        arg = tokens[1]
        if not isinstance(arg, (ndarray, float)):
            arg = self._evaluate(sel, loc, arg)
         
        if (isinstance(arg, float) or isinstance(arg, ndarray) and 
            arg.dtype in (float, int)):
            return FUNCTION_MAP[tokens[0]](arg)
        else:
            raise SelectionError(sel, loc, "functions (sin/abs/etc.) must have"
                                        " numeric arguments, e.g. 'sin(x)'")

    def _evalUserdata(self, sel, loc, keyword, values=None, 
                      evalonly=None):
        
        if DEBUG: print('_evalUserdata', keyword, values)
        data = self._atoms.getData(keyword)
        if values is None:
            if data.dtype == bool:
                if evalonly is None:
                    return data
                else:
                    return data[evalonly]
            else:
                return SelectionError(sel, loc, "data type of {0:s} must be "
                                      "bool".format(repr(keyword)))
        else:
            if data.dtype in (int, float):
                return self._evalFloat(sel, loc, keyword, values, 
                                       evalonly=evalonly)
            elif data.dtype.type == np.string_:
                return self._evalAlnum(sel, loc, keyword, values, 
                                       evalonly=evalonly)
            else:
                return SelectionError(sel, loc, "data type of {0:s} must be "
                                  "int, float, or str".format(repr(keyword)))

    def _evalBoolean(self, sel, loc, keyword, evalonly=None, rtrn=False):
        """Evaluate a boolean keyword."""
    
        if DEBUG: print('_evalBoolean', keyword)
        if evalonly is None:
            n_atoms = self._n_atoms
        else:        
            n_atoms = len(evalonly)
        
        if keyword == 'all':
            return ones(n_atoms, bool)
        elif keyword == 'none':
            return zeros(n_atoms, bool)
        else:
            keyword = expandBoolean(keyword)
            torf = self._and(sel, loc, [keyword], rtrn=rtrn)
            if evalonly is None or isinstance(torf, SelectionError):
                return torf
            else:
                return torf[evalonly] 

    def _evalAlnum(self, sel, loc, keyword, values, evalonly=None):
        """Evaluate keywords associated with alphanumeric data, e.g. residue 
        names, atom names, etc."""
        
        if DEBUG: print('_evalAlnum', keyword, values)
        
        if keyword == 'sequence':
            return self._sequence(sel, loc, keyword, values, evalonly)
    
        keyword = KEYWORDS_SYNONYMS.get(keyword, keyword)
        data = self._getData(sel, loc, keyword)
        if isinstance(data, SelectionError):
            return data
            
        if keyword in _specialKeywords:
            for i, value in enumerate(values):
                if value == '_':
                    values[i] = ' '
                    values.append('')
                    break
        if evalonly is not None:
            data = data[evalonly]
        if DEBUG: print('_evalAlnum set(data)', set(data))
        n_atoms = len(data)
        vallen = ALNUM_VALLEN[keyword]

        regexps = []
        strings = []
        for value in values:
            if isinstance(value, str):
                strings.append(value)
                if len(value) > vallen:
                    return SelectionError(sel, sel.find(value, loc), 'invalid '
                        'value, maximum {0:d} char(s) allowed for {1:s}'
                        .format(vallen, repr(keyword)))
            else:
                regexps.append(value)
        
        if len(strings) == 1:
            torf = data == strings[0]
        elif len(strings) > 4:
            torf = zeros(n_atoms, np.bool)
            strings = set(strings)
            for i, datum in enumerate(data):        
                torf[i] = datum in strings
        elif strings: 
            torf = [(data == value).reshape((n_atoms, 1)) for value in strings]
            torf = np.concatenate(torf, 1).sum(1).astype(np.bool) 
        else:
            torf = zeros(n_atoms, np.bool)

        for value in regexps:
            for i in xrange(n_atoms):
                torf[i] = (value.match(data[i]) is not None)

        return torf
    
    def _sequence(self, sel, loc, keyword, values, evalonly):
        
        if DEBUG: print('_sequence', values)
        atoms = self._atoms
        if isinstance(atoms, AtomGroup):
            citer = atoms.iterChains()
        elif isinstance(atoms, Chain):
            citer = iter([atoms])
        elif isinstance(atoms, Segment):
            citer = iter(atoms)
        else:
            citer = iter(HierView(self._atoms))
        
        matches = []
        for chain in citer:
            sequence = chain.getSequence()
            residues = [res for res in chain if res.getResname() in AAMAP]
            if len(sequence) != len(residues):
                return SelectionError(sel, loc, 'amino acid sequence length '
                    'does not match number of amino acid residues')
            for value in values: 
                if isinstance(value, str):
                    if not value.isalpha() or not value.isupper():
                        return SelectionError(sel, sel.find(value, loc),
                            'sequence string must be upper case letters or '
                            'a regular expression')
                    value = RE.compile(value)
                for match in value.finditer(sequence):  
                    matches.extend(residues[match.start():match.end()])
        if matches:
            torf = zeros(len(self._ag), bool)
            indices = concatenate([res._getIndices() for res in matches])
            torf[indices] = True
            if self._indices is not None:
                torf = torf[self._indices]
            if evalonly is not None:
                torf = torf[evalonly]
        elif evalonly is not None:
            torf = zeros(len(evalonly), bool)
        else:
            torf = zeros(self._n_atoms, bool)
        return torf
    
    def _evalFloat(self, sel, loc, keyword, values=None, evalonly=None):
        """Evaluate a keyword associated with atom attributes of type float. 
        If *values* is not passed, return the attribute array."""
        
        if DEBUG: print('_evalFloat', keyword, values)
        if keyword == 'x':
            data = self._getCoords(sel, loc)[:,0]
        elif keyword == 'y':
            data = self._getCoords(sel, loc)[:,1]
        elif keyword == 'z':
            data = self._getCoords(sel, loc)[:,2]
        else:
            data = self._getData(sel, loc, keyword)
    
        if values is None or isinstance(data, SelectionError): 
            return data
    
        if evalonly is not None:
            data = data[evalonly]
        n_atoms = len(data)
        torf = zeros(n_atoms, np.bool)

        numbers = self._getNumRange(sel, loc, values)
        if isinstance(numbers, SelectionError):
            return numbers
    
        for item in numbers:
            if isinstance(item, str):
                pass
            elif isinstance(item, list):
                torf[(item[0] <= data) & (data <= item[1])] = True
            elif isinstance(item, tuple):
                if len(item) == 2:
                    torf[(item[0] <= data) & (data < item[1])] = True
                else:
                    item = [str(i) for i in item]
                    return SelectionError(sel, loc, repr(':'.join(item)) + 
                          ' cannot be used with ' + repr(keyword) + ', but ' + 
                          repr(':'.join(item[:2])) + ' is accepted')
            else:
                torf[data == item] = True
        return torf

    def _resnum(self, sel, loc, token=None, numRange=True, evalonly=None):
        """Evaluate "resnum" keyword."""
        
        if DEBUG: print('_resnum', token)
        resids = self._getData(sel, loc, 'resnum')
        if token is None or isinstance(resids, SelectionError):
            return resids
    
        icodes = None
        n_atoms = self._n_atoms
        if evalonly is not None:
            resids = resids[evalonly]
            n_atoms = len(evalonly)
        torf = zeros(n_atoms, np.bool)
        
        if numRange:
            token = self._getNumRange(sel, loc, token, False)
            if isinstance(token, SelectionError):
                return token
        
        for item in token:
            if isinstance(item, str):
                if icodes is None:
                    if evalonly is None:
                        icodes = self._getData(sel, loc, 'icode')
                    else:
                        icodes = self._getData(sel, loc, 'icode')[evalonly]
                icode = str(item[-1])
                if icode == '_':
                    icode = ''
                try:
                    number = int(item[:-1])
                except ValueError:
                    return SelectionError(sel, loc, 'all values must start '
                                          'with a number')
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

    def _serial(self, sel, loc, token=None, evalonly=None):
        """Evaluate *serial* keyword."""
        
        if DEBUG: print('_serial', token)
        
        if token is None:
            return self._getData(sel, loc, 'serial')

        try:
            sn2i = self._ag._getSN2I()
        except (AttributeError, ValueError) as err:
            return SelectionError(sel, loc, str(err))
        
        numbers = self._getNumRange(sel, loc, token)
        if isinstance(numbers, SelectionError):
            return numbers
        indices = []
        for item in numbers:
            if isinstance(item, list):
                indices.append(sn2i[int(np.ceil(item[0])):
                                    int(np.floor(item[1]))+1])
            elif isinstance(item, tuple):
                if len(item) == 2:
                    indices.append(sn2i[item[0]:item[1]])
                else:
                    indices.append(sn2i[item[0]:item[1]:item[2]])
            else:
                indices.append([sn2i[item]])

        indices = unique(concatenate(indices))
        torf = zeros(self._ag.numAtoms(), bool)
        if indices[0] == -1:
            torf[indices[1:]] = True
        else: 
            torf[indices] = True
        if self._indices is not None:
            torf = torf[self._indices]
        if evalonly is None:
            return torf
        else:
            return torf[evalonly]
    
    def _index(self, sel, loc, token=None, evalonly=None):
        """Evaluate "index" keyword."""
        
        if DEBUG: print('_index', token)
        
        if token is None:
            return self._indices or np.arange(self._ag._n_atoms)
        torf = zeros(self._ag._n_atoms, np.bool)
        
        numbers = self._getNumRange(sel, loc, token)
        if isinstance(numbers, SelectionError):
            return numbers
    
        for item in numbers:
            try:
                if isinstance(item, tuple):
                    if len(item) == 2:
                        torf[item[0]:item[1]] = True
                    else:
                        torf[item[0]:item[1]:item[2]] = True
                elif isinstance(item, list):
                    torf[int(np.ceil(item[0])):int(np.floor(item[1]))+1] = True
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

    def _getNumRange(self, sel, loc, token, intfloat=True):
        """Evaluate numeric values. Identify ranges, integers, and floats,
        put them in a list and return."""
        
        if DEBUG: print('_getNumRange', type(token), token)
        if isinstance(token, ndarray):
            return token
        if any([isinstance(tkn, SRE_Pattern) for tkn in token]):
            return SelectionError(sel, loc, 'values must be numbers or ranges')
            
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
                    return SelectionError(sel, loc, "{0:s} is not understood"
                                     .format(repr(' to '.join(items))))
                try:
                    token.append([float(items[0]), float(items[1])])
                except:
                    return SelectionError(sel, loc, "'to' must be surrounded "
                                          "by numbers")
            elif ':' in item:
                # : means upper bound is NOT included in the range
                # boundaries are placed in a TUPLE
                items = item.split(':')
                if not len(items) in (2, 3):
                    return SelectionError(sel, loc, "{0:s} is not understood"
                                         .format(repr(':'.join(items))))
                try:
                    if len(items) == 2:
                        token.append( (int(items[0]), int(items[1])) )
                    else:
                        token.append((int(items[0]), int(items[1]),
                                      int(items[2])))
                except:
                    return SelectionError(sel, loc, "':' must be surrounded "
                                          "by integers")
            else:
                try: 
                    item = int(item)
                except ValueError:
                    try:
                        item = float(item)
                    except ValueError:
                        if intfloat:
                            return SelectionError(sel, loc, "all values must "
                                                  "be numbers")
                token.append(item)
        if DEBUG: print('_getNumRange', token)            
        return token
    
    def _getData(self, sel, loc, keyword):
        """Return atomic data."""
        
        data = self._data.get(keyword)
        if data is None:        
            field = ATOMIC_FIELDS.get(keyword)
            if field is None:
                data = self._ag._getData(keyword)
                if data is None:
                    return SelectionError(sel, loc, "{0:s} is not a valid "
                          "keyword or user data label".format(repr(keyword)))
                elif not isinstance(data, ndarray) and data.ndim == 1:
                    return SelectionError(sel, loc, "{0:s} must be a 1d "
                                          "array".format(repr(keyword)))
            else:
                data = getattr(self._ag, '_get' + field.meth_pl)() 
                if data is None:
                    return SelectionError(sel, loc, "{0:s} is not set by "
                                          "user".format(repr(keyword)))
            self._data[keyword] = data
        indices = self._indices
        if indices is None:               
            return data
        else:
            return data[indices]
    
    def _getCoords(self, sel, loc):
        """Return coordinates of selected atoms.  This method is reduces array 
        copy operations when :class:`~.AtomPointer` subclasses are used for
        making atom selections."""
        
        if self._coords is None:
            self._coords = self._atoms._getCoords()
            if self._coords is None:
                return SelectionError(sel, loc, 'coordinates are not set')
        return self._coords
