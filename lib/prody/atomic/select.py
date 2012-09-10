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

""" This module defines a class for selecting subsets of atoms.  Read this
page in interactive sessions using ``help(select)``.

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


Atom flags 
===============================================================================

All :ref:`flags` can be used as keywords in atom selections:
   
>>> from prody import * 
>>> p = parsePDB('1ubi')
>>> p.select('protein')
<Selection: 'protein' from 1ubi (602 atoms)>
>>> p.select('water')
<Selection: 'water' from 1ubi (81 atoms)>
    
Atom data fields 
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
fragindex [\\\\]   integer, range  unique index number for distinct fragments
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
is parsed. Using them before they are set will raise a selection error. 
Secondary structure assignments can be made using :func:`.assignSecstr` 
function.

**[†]** Alternate locations are parsed as alternate coordinate sets. This
keyword will work for alternate location specified by "A". This to work for
alternate locations indicated by other letters, they must be parsed 
specifically by passing the identifier to the :func:`.parsePDB`.

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
numbers to these entitites are assigned by :class:`.HierView` class
upon building of a hierarchical view for an :class:`.AtomGroup`.
Note that hierarchical views are build automatically when needed.

**[\\\\]** Distinct fragments are connected subsets of atoms.  Fragments are 
determined automatically when needed and bond information is available.

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
  * ``'x 1 to 10'`` selects atoms whose x coordinates are greater than or 
    equal to 1 and less than or equal to 10  
  * ``'x 1:10'`` selects atoms whose x coordinates are greater than or 
    equal to 1 and less than or equal to 10
    
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
``bonded to ...`` setting, e.f. ``bonded to calpha`` will select atoms bonded 
to Cα atoms.  For this setting to work, bonds must be set by the user using the
:meth:`.AtomGroup.setBonds` method.  It is also possible to select bonded atoms
by excluding the originating atoms using ``exbonded to ...`` setting.  

Selection macros
===============================================================================

Any valid selection string can be used to define selection macros using the 
:func:`defSelectionMacro` function.  Macros are saved in ProDy configuration 
and loaded in later sessions automatically.  Below functions are for 
manipulating selection macros:
    
  * :func:`defSelectionMacro`
  * :func:`delSelectionMacro`
  * :func:`getSelectionMacro`"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from re import compile as re_compile
from code import interact
from types import NoneType

import numpy as np
from numpy import array, ndarray, ones, zeros, arange
from numpy import invert, unique, concatenate, all, any
from numpy import logical_and, logical_or 

import pyparsing as pp

from prody import LOGGER, SETTINGS

from .atomic import Atomic
from .fields import ATOMIC_FIELDS
from .flags import PLANTERS as FLAG_PLANTERS

from .atomgroup import AtomGroup 
from .chain import Chain, getSequence, AAMAP
from .pointer import AtomPointer
from .selection import Selection
from .segment import Segment
from .atommap import AtomMap

from prody.utilities import rangeString
from prody.kdtree import KDTree

DEBUG = 0
NUMB = 0 # Select instance will not really evaluate string for the atoms 

def debug(sel, loc, *args):
    
    if DEBUG:
        print ''
        if args:
            print(args[0], args[1:])
        print(repr(sel))
        print(' ' * (loc + 1) + '^')

__all__ = ['Select', 'SelectionError', 'TypoWarning',
           'defSelectionMacro', 'delSelectionMacro', 'getSelectionMacro',
           'isSelectionMacro']
           
SRE_Pattern = type(re_compile('C.*'))

DTYPES_NUMERIC = set([np.int, np.int32, np.int64, 
                      np.float, np.float32, np.float64])
try:
    DTYPES_NUMERIC.add(np.float16)
except:
    pass
try:
    DTYPES_NUMERIC.add(np.int16)
except:
    pass
try:
    DTYPES_NUMERIC.add(np.int8)
except:
    pass

EXCMSG_NUMERIC = ('{0:s} must be a number or a label for numeric data that is '
                  'present or that can be calculated')

KEYWORDS_STRING = set(['name', 'type', 'resname', 'chain', 'element', 
                       'segment', 'altloc', 'secondary', 'icode',
                       'chid', 'secstr', 'segname', 'sequence'])
FIELDS_ALNUM = dict()
FIELDS_NUMERIC = set()
FIELDS_SYNONYMS = dict()
for name, field in ATOMIC_FIELDS.iteritems():
    if field.ndim == 1:
        if field.dtype in DTYPES_NUMERIC:
            FIELDS_NUMERIC.add(name)
            if field.synonym:
                FIELDS_NUMERIC.add(field.synonym)
                FIELDS_SYNONYMS[field.synonym] = name
        else:
            itemsize = np.dtype(field.dtype).itemsize
            FIELDS_ALNUM[name] = itemsize
            if field.synonym:
                FIELDS_ALNUM[field.synonym] = itemsize
                FIELDS_SYNONYMS[field.synonym] = name 


KEYWORDS_INTEGER = set(['serial', 'index', 'resnum', 'resid', 
                        'segindex', 'chindex', 'resindex', 
                        'fragindex', 'fragment', 'numbonds'])
KEYWORDS_FLOAT = set(['x', 'y', 'z', 'beta', 'mass', 'occupancy', 'mass', 
                      'radius', 'charge'])
KEYWORDS_NUMERIC = KEYWORDS_FLOAT.union(KEYWORDS_INTEGER)    

KEYWORDS_VALUE_PAIRED = KEYWORDS_NUMERIC.union(KEYWORDS_STRING)

XYZMAP = {'x': 0, 'y': 1, 'z': 2}

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
    '+'  : np.add,
    '-'  : np.subtract,
    '*'  : np.multiply,
    '/'  : np.divide,
    '%'  : np.remainder,
    '>'  : np.greater,
    '<'  : np.less,
    '>=' : np.greater_equal,
    '<=' : np.less_equal,
    '='  : np.equal,
    '==' : np.equal,
    '!=' : np.not_equal,
}

_TOSPACE = set(['secondary', 'chain', 'altloc', 'segment', 'icode'])
    
ATOMGROUP = None

MACROS = SETTINGS.get('selection_macros', {})
MACROS_REGEX = None

def isSelectionMacro(word):
    """Return **True** if *word* is a user defined selection macro."""
    
    try:
        return word in MACROS
    except:
        return False


def defSelectionMacro(name, selstr):
    """Define selection macro *selstr* with name *name*.  Both *name* and 
    *selstr* must be string.  An existing keyword cannot be used as a macro 
    name. If a macro with given *name* exists, it will be overwritten.
    
    >>> defSelectionMacro('cbeta', 'name CB and protein')"""
    
    if not isinstance(name, str) or not isinstance(selstr, str):
        raise TypeError('both name and selstr must be strings')
    elif isReserved(name):
        raise ValueError('{0:s} is a reserved word and cannot be used as a '
                         'macro name'.format(repr(name)))
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

def evalMacros(selstr):
    
    if MACROS:
        global MACROS_REGEX
        if MACROS_REGEX is None: MACROS_REGEX = {}
        selstr = ' ' + selstr + ' '
        for key, macro in MACROS.items():
            
            selstr = (' (' + macro + ') ').join(MACROS_REGEX.setdefault(key, 
                        re_compile('[( )]' + key + '[( )]')).split(selstr))
        selstr = selstr[1:-1]
    return selstr   


def areAllStrings(container):
    """Return ``True`` if all items in *container* are instances of 
    :func:`str`."""
    
    for item in container: 
        if not isinstance(item, str):
            return False
    return True


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

KEYWORDS_BOOLEAN = set(FLAG_PLANTERS.keys() + ['all', 'none'])

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
    return isBooleanKeyword(keyword) or isValuePairedKeyword(keyword)

def tkn2str(token):
    
    if isinstance(token, str):
        return token
    else:
        return ' '.join(token)


SAMEAS_MAP = {'residue': 'resindex', 'chain': 'chindex', 
              'segment': 'segindex', 'fragment': 'fragindex'}


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


XYZDIST = set(['x', 'y', 'z', 'within', 'exwithin'])

def checkSelstr(selstr, what, error=ValueError):
    """Check *selstr* if it satisfies a selected condition.  For now, only
    whether coordinate/distance based selection are checked.  If *error* is 
    a subclass of :class:`Exception`, an exception will be raised, otherwise 
    return **True** or **False** will be returned."""
    
    selstr = selstr.replace('(', ' ( ')
    selstr = selstr.replace(')', ' ) ')
    
    if what in set(['dist']):
        for item in selstr.split():
            if item in XYZDIST:
                if issubclass(error, Exception):
                    raise error('invalid selection {0:s}, coordinate '
                                'based selections are not accepted'
                                .format(repr(selstr)))
                else:
                    return False

def evalNumeric(sel, loc, token, allnum=True):
    """Return a list of ranges, integers, and floats extracted from the token.
    If *allnum* is **True**, a :class:`SelectionError` will be returned.
    """
    
    if isinstance(token, ndarray):
        return token
    for tkn in token:
        try:
            _ = tkn.pattern, tkn.finditer
        except: 
            pass
        else:
            return SelectionError(sel, loc, 'values must be numbers or ranges')
        
    tknstr = ' '.join(' '.join(token).split())
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
                    if allnum:
                        return SelectionError(sel, loc, 
                                              'all values must be numbers')
            token.append(item)
    return token


AND = pp.Keyword('and')
OR = pp.Keyword('or')
NOT = pp.Keyword('not')

SHORTLIST = pp.alphanums + '''~@#$.:;_','''

LONGLIST = pp.alphanums + '''~!@#$%^&*()-_=+[{}]\|;:,<>./?()' '''
SPECIALCHARS = pp.Group(pp.Literal('`') + 
                        pp.Optional(pp.Word(LONGLIST + '"')) + 
                        pp.Literal('`'))
def specialCharsParseAction(token):
    if len(token[0]) == 2: # meaning `` was used
        return '_'
    else:
        return token[0][1]
SPECIALCHARS.setParseAction(specialCharsParseAction)

REGULAREXP = pp.Group(pp.Literal('"') + 
                      pp.Optional(pp.Word(LONGLIST + '`')) + 
                      pp.Literal('"'))
def regularExpParseAction(sel, loc, token):
    
    token = token[0]
    if len(token) == 2:
        return re_compile('^()$')
    else:
        try:
            regexp = re_compile(token[1])
        except:
            raise SelectionError(sel, loc, 'failed to compile regular '
                            'expression {0:s}'.format(repr(token[1])))
        else:
            return regexp
REGULAREXP.setParseAction(regularExpParseAction)

# number ranges

FLOAT = pp.Regex(r'\d+(\.\d*)?([eE]\d+)?')
RANGEPY = pp.Group(FLOAT + pp.Literal(':') + FLOAT +
                   pp.Optional(pp.Group(pp.Literal(':') + FLOAT)))
RANGETO = pp.Group(FLOAT + pp.Literal('to') + FLOAT)

def rangeParseAction(sel, loc, tokens):
    
    tokens = tokens[0]
    debug(sel, loc, '_range', tokens)
    arange = [None, None, None]
    start = float(tokens[0])
    stop = float(tokens[2])

    if start > stop:
        raise SelectionError(sel, loc, 'range start value ({0:s}) is greater '
            'than and stop value ({1:s})'.format(tokens[0], tokens[2]))
    elif start == stop:
        return tokens[0]

    if tokens[1] == 'to':
        comp = '<='
    elif len(tokens) == 5:
        comp = float(tokens[4])
    else:
        comp = '<'
    return 'range', start, stop, comp
    
RANGEPY.setParseAction(rangeParseAction)
RANGETO.setParseAction(rangeParseAction)

WORDS = ~AND + ~OR + ~NOT

BONDED = pp.Regex('(ex)?bonded to')
BONDED2 = pp.Regex('(ex)?bonded [0-9]+ to')
SAMEAS = pp.Regex('same [a-z]+ as')
WITHIN = pp.Regex('(ex)?within [0-9]+\.?[0-9]* of')

WORDS += ~BONDED + ~BONDED2 + ~SAMEAS + ~WITHIN 

FUNCNAMES = FUNCTION_MAP.keys()
FUNCTIONS = pp.Keyword(FUNCNAMES[0])
for func in FUNCNAMES[1:]:
    kwfunc = pp.Keyword(func)
    FUNCTIONS = FUNCTIONS | kwfunc 
    WORDS += ~kwfunc
WORDS += pp.Word(SHORTLIST)
ONEORMORE = pp.OneOrMore(RANGEPY | RANGETO | REGULAREXP | SPECIALCHARS | WORDS)


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
        
        self._parser = pp.operatorPrecedence(
             ONEORMORE,
             [(FUNCTIONS, 1, pp.opAssoc.RIGHT, self._func),
              (pp.oneOf('+ -'), 1, pp.opAssoc.RIGHT, self._sign),
              (pp.oneOf('** ^'), 2, pp.opAssoc.LEFT, self._pow),
              (pp.oneOf('* / %'), 2, pp.opAssoc.LEFT, self._binop),
              (pp.oneOf('+ -'), 2, pp.opAssoc.LEFT, self._binop),
              (pp.oneOf('< > <= >= == = !='), 2, pp.opAssoc.LEFT, self._comp),
              (NOT | BONDED | BONDED2 | SAMEAS | WITHIN, 
               1, pp.opAssoc.RIGHT, self._unary),
              (pp.Optional(AND), 2, pp.opAssoc.LEFT, self._and),
              (OR, 2, pp.opAssoc.LEFT, self._or),]
            )

        self._parser.setParseAction(self._default)
        self._parser.leaveWhitespace()
        self._parser.enablePackrat()
        
        self._map = {'resnum': self._resnum, 'resid': self._resnum, 
            'serial': self._serial, 'index': self._index,
            'x': self._generic, 'y': self._generic, 'z': self._generic,
            'sequence': self._sequence}
        
    def _reset(self):

        self._ag = None
        self._atoms = None
        self._indices = None
        self._n_atoms = None
        self._coords = None
        self._data.clear()
        
    def _evalAtoms(self, atoms):
        
        self._atoms = atoms
        try:
            self._ag = atoms.getAtomGroup()
        except AttributeError:
            self._ag = atoms
            self._indices = None
        else:
            self._indices = atoms._getIndices()
            if len(self._indices) == 1:
                try:                
                    index = atoms.getIndex()
                except AttributeError:
                    pass
                else:
                    self._atoms = Selection(self._ag, array([index]), 
                                    'index ' + str(index), atoms.getACSIndex())
        
    def select(self, atoms, selstr, **kwargs):
        """Return a subset of atoms matching *selstr* as a :class:`Selection`.
        If selection string does not match any atoms, **None** is returned.
        
        :arg atoms: atoms to be evaluated    
        :type atoms: :class:`.Atomic`
        
        :arg selstr: selection string
        :type selstr: str
        
        If type of *atoms* is :class:`.AtomMap`, an :class:`.AtomMap` instance
        is returned.  Otherwise, :class:`.Selection` instances are returned.

        .. note:
              
            * :meth:`select` accepts arbitrary keyword arguments which enables 
              identification of intermolecular contacts. See :ref:`contacts` 
              for details.
        
            * A special case for making atom selections is passing an
              :class:`.AtomMap` instance as *atoms* argument.  Dummy 
              atoms will not be included in the result, but the order 
              of atoms will be preserved."""
        
        self._ss2idx = False
        self._replace = False
        
        self._selstr = selstr
        indices = self.getIndices(atoms, selstr, **kwargs)
        
        self._kwargs = None

        if len(indices) == 0:
            return None

        if not isinstance(atoms, AtomGroup):
            indices = self._indices[indices]
            
        ag = self._ag

        try:
            dummies = atoms.numDummies()
        except AttributeError:
            if self._ss2idx:
                selstr = 'index {0:s}'.format(rangeString(indices))
            else:
                if self._replace:
                    for key, value in kwargs.iteritems():
                        if value in self._ag and key in selstr:
                            selstr = selstr.replace(key, 
                                                '(' + value.getSelstr() + ')')
                if isinstance(atoms, AtomPointer):
                    selstr = '({0:s}) and ({1:s})'.format(selstr, 
                                                      atoms.getSelstr())
            
            return Selection(ag, indices, selstr, atoms.getACSIndex(),
                             unique=True)
        else:
            return AtomMap(ag, indices, atoms.getACSIndex(), dummies=dummies,
               title='Selection {0:s} from '.format(repr(selstr)) + str(atoms))
               
    def getIndices(self, atoms, selstr, **kwargs):
        """Return indices of atoms matching *selstr*.  Indices correspond to 
        the order in *atoms* argument.  If *atoms* is a subset of atoms, they 
        should not be used for indexing the corresponding :class:`.AtomGroup`
        instance."""
        
        selstr = selstr.strip() 
        if (len(selstr.split()) == 1 and selstr.isalnum() and 
            selstr not in MACROS):
            self._evalAtoms(atoms)
            if selstr == 'none':
                return array([])
            elif selstr == 'all':
                return arange(atoms.numAtoms())
            elif atoms.isFlagLabel(selstr):
                return atoms._getFlags(selstr).nonzero()[0]
            elif isValuePairedKeyword(selstr):
                raise SelectionError(selstr, 0, 'must be followed by values')
            else:
                raise SelectionError(selstr, 0, 'is not a valid selection or '
                                     'user data label')
        else:
            torf = self.getBoolArray(atoms, selstr, **kwargs)        
            return torf.nonzero()[0]
        
    def getBoolArray(self, atoms, selstr, **kwargs):
        """Return a boolean array with **True** values for *atoms* matching 
        *selstr*.  The length of the boolean :class:`numpy.ndarray` will be
        equal to the number of atoms in *atoms* argument."""
        
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance, not {0:s}'
                            .format(type(atoms)))

        self._reset()

        for key in kwargs.iterkeys():
            if not key.isalnum():
                raise TypeError('{0:s} is not a valid keyword argument, '
                                  'keywords must be all alpha numeric '
                                  'characters'.format(repr(key)))
            if isReserved(key):
                loc = selstr.find(key)
                if loc > -1:
                    raise SelectionError(selstr, loc, '{0:s} is a reserved '
                                'word and cannot be used as a keyword argument'
                                .format(repr(key)))

        self._n_atoms = atoms.numAtoms()
        self._selstr = selstr
        self._kwargs = kwargs
        
        if DEBUG: print('getBoolArray', selstr)
            
        self._evalAtoms(atoms)
            
        selstr = selstr.strip() 
        if (len(selstr.split()) == 1 and selstr.isalnum() and 
            selstr not in MACROS):
            if selstr == 'none':
                return zeros(atoms.numAtoms(), bool)
            elif selstr == 'all':
                return ones(atoms.numAtoms(), bool)
            elif atoms.isFlagLabel(selstr):
                return atoms.getFlags(selstr)
            elif isValuePairedKeyword(selstr):
                raise SelectionError(selstr, 0, 'must be followed by values')
            else:
                raise SelectionError(selstr, 0, 'is not a valid selection or '
                                     'user data label')
            

        try:
            tokens = self._parser.parseString(selstr, parseAll=True).asList()
        except pp.ParseException as err:
            pass
            which = selstr.rfind(' ', 0, err.column)
            if which > -1:
                if selstr[which + 1] == '(':
                    msg = ('an arithmetic, comparison, or logical operator '
                           'must precede the opening parenthesis')
                elif selstr[which - 1] == ')':
                    msg = ('an arithmetic, comparison, or logical operator '
                           'must follow the closing parenthesis')
                else:
                    msg = 'parsing failed here'
            else:
                msg = 'parsing failed here'
                
            raise SelectionError(selstr, err.column, msg + '\n' + str(err))
        else:    
            if DEBUG: print('_evalSelstr', tokens)
            torf = tokens[0]

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

        
    def _default(self, sel, loc, tokens, subset=None):
        
        debug(sel, loc, '_default', tokens)
        if NUMB: return
        
        if subset != None:
            # is neeeded
            raise SelectionError(sel, loc, 'susbet is not None')

        if len(tokens) == 1:
            token = tokens[0]
            try:
                dtype = token.dtype
            except AttributeError:
                pass
            else:
                return token
            
            if token == 'none':
                return zeros(self._n_atoms if subset is None else len(subset), 
                              bool)
            elif token == 'all':
                return ones(self._n_atoms if subset is None else len(subset), 
                            bool)
            elif self._atoms.isFlagLabel(token):
                return self._getFlags(token, subset)
            try: 
                arg = self._kwargs[first]
            except KeyError:
                pass
            else:
                if arg in self._ag:
                    
                    try:
                        dummies = arg.numDummies()
                    except AttributeError:
                        indices = arg._getIndices()
                    else:
                        if dummies:
                            indices = arg.getIndices()[arg.getFlags('mapped')]
                        else:
                            indices = arg.getIndices()
                            
                    torf = zeros(self._ag.numAtoms(), bool)
                    torf[indices] = True
                    if self._indices is not None:
                        torf = torf[self._indices]
                    self._replace = True
                    if DEBUG: print('_evaluate', first, torf)
                    return torf
                return first
            try:
                return float(first)
            except ValueError:
                pass
        else:
            # is neeeded
            raise SelectionError(sel, loc, '_default len(tokens) > 1')
    
    
    def _eval(self, sel, loc, tokens, subset=None):
        """"""
        
        debug(sel, loc, '_eval', tokens)
        if NUMB: return
        #if isinstance(tokens, ndarray):
        #    return tokens

        if len(tokens) == 1:
            token = tokens[0]
            if token == 'none':
                return zeros(self._n_atoms if subset is None else len(subset), 
                              bool), False
            elif token == 'all':
                return ones(self._n_atoms if subset is None else len(subset), 
                            bool), False
            elif self._atoms.isFlagLabel(token):
                return self._getFlags(token, subset), False
            try: 
                arg = self._kwargs[first]
            except KeyError:
                pass
            else:
                if arg in self._ag:
                    
                    try:
                        dummies = arg.numDummies()
                    except AttributeError:
                        indices = arg._getIndices()
                    else:
                        if dummies:
                            indices = arg.getIndices()[arg.getFlags('mapped')]
                        else:
                            indices = arg.getIndices()
                            
                    torf = zeros(self._ag.numAtoms(), bool)
                    torf[indices] = True
                    if self._indices is not None:
                        torf = torf[self._indices]
                    self._replace = True
                    if DEBUG: print('_evaluate', first, torf)
                    return torf, False
                return first, False
            try:
                return float(first), False
            except ValueError:
                pass
        else:
            return self._map.get(tokens[0], self._generic)(sel, loc, tokens, 
                                                           subset)
            
            #if self._atoms.isFlagLabel(tokens[0]):
            #    torf = self._getFlags(tokens.pop(0), subset)
            #    return logical_and(torf, self._eval(sel, loc, tokens, subset), 
            #                       torf), False
                
    def _getFlags(self, label, subset):
        
        if subset is None:
            # get a copy to avoid alterations
            return self._atoms.getFlags(label)
        else:
            return self._atoms._getFlags(label)[subset]
        
    
    def _dontevaluate(self, sel, loc, tkns, evalonly=None):
        """Evaluates statements in a selection string, e.g. ``'calpha'``,
        ``'index 5'``."""
        
        debug(sel, loc, '_evaluate', tkns)
        
        if isinstance(tkns, str):
            # NOT ENCOUNTERED
            #with open('/home/abakan/evaluate.txt', 'a') as out:
            #    out.write('STRING ' + sel + '\n')
            if DEBUG: print('_evaluate', 'STRING STRING STRING STRING', tkns)
            if tkns == 'none':
                return zeros(self._n_atoms, bool)
            if self._atoms.isFlagLabel(tkns):
                if evalonly is None:
                    return self._atoms.getFlags(tkns) # get a copy of flags
                else:
                    return self._atoms._getFlags(tkns)[evalonly]
            elif self._ag.isDataLabel(tkns):
                return self._evalUserdata(sel, loc, tkns, evalonly=evalonly)
            else:
                return SelectionError(sel, loc, '{0:s} is not understood'
                                      .format(repr(tkns)))
        elif isinstance(tkns, (ndarray, float)):
            # NOT ENCOUNTERED
            #with open('/home/abakan/evaluate.txt', 'a') as out:
            #    out.write('NDARRAY ' + sel + '\n')
            return tkns
    
        keyword = tkns[0]
        if len(tkns) == 1:
            #with open('/home/abakan/evaluate.txt', 'a') as out:
            #    out.write('LIST ' + sel + ' - ' + str(keyword) + '\n')
            if keyword == 'none':
                return zeros(self._n_atoms, bool)
            if self._atoms.isFlagLabel(keyword):
                if evalonly is None:
                    return self._atoms.getFlags(keyword) # get a copy of flags
                else:
                    return self._atoms._getFlags(keyword)[evalonly]
            elif isNumericKeyword(keyword):
                return self._evalNumeric(sel, loc, keyword)
            elif self._ag.isDataLabel(keyword):
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
                #with open('/home/abakan/evaluate.txt', 'a') as out:
                #    out.write('FLOAT ' + sel + ' - ' + str(keyword) + '\n')
                try:
                    return float(keyword)
                except ValueError:
                    pass
        elif self._atoms.isFlagLabel(keyword):
            if evalonly is None:
                return all((self._atoms._getFlags(keyword), 
                            self._evaluate(sel, loc, tkns[1:])), 0)
            else:
                return all((self._atoms._getFlags(keyword)[evalonly], 
                            self._evaluate(sel, loc, tkns[1:], evalonly)), 0)
        #try: 
        #    return self._map[keyword](sel, loc, keyword, tkns[1:], evalonly)
        #except KeyError:
        #    pass
        
        elif isAlnumKeyword(keyword):
            return self._evalAlnum(sel, loc, keyword, tkns[1:], evalonly=evalonly)
        elif keyword in ('resnum', 'resid'):
            return self._resnum(sel, loc, tkns[1:], evalonly=evalonly)
        elif keyword == 'index':
            return self._index(sel, loc, tkns[1:], evalonly=evalonly)
        elif keyword == 'serial':
            return self._serial(sel, loc, tkns[1:], evalonly=evalonly)
        elif isNumericKeyword(keyword):
            return self._evalFloat(sel, loc, keyword, tkns[1:], evalonly=evalonly)
        elif keyword == NOT:
            return self._not(sel, loc, tkns, evalonly=evalonly)
        elif self._ag.isDataLabel(keyword):
            return self._evalUserdata(sel, loc, keyword, tkns[1:], 
                                      evalonly=evalonly)
        
        return SelectionError(sel, loc, "{0:s} is not understood"
                              .format(repr(' '.join(tkns))))

    def _or(self, sel, loc, tokens):
        """Evaluate statements containing ``'or'`` operator."""
        
        debug(sel, loc, '_or', tokens)
        if NUMB: return
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
                selection = any(arrays, 0, arrays[0])
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

    def _and(self, sel, loc, tokens):
        """Evaluate statements containing ``'and'`` operator."""
        
        tokens = tokens[0]
        debug(sel, loc, '_and', tokens)
        if NUMB: return
        
        flags = []
        torfs = []
        evals = []
        atoms = self._atoms
        isFlagLabel = atoms.isFlagLabel
        isDataLabel = atoms.isDataLabel

        for token in tokens:
            # check whether token is an array to avoid array == str comparison
            try:
                dtype = token.dtype
            except AttributeError:
                if token == 'and':
                    continue
                elif isFlagLabel(token):
                    flags.append(token)
                elif isDataLabel(token):
                    evals.append([])
                    evals[-1].append(token)
                else:
                    try:
                        evals[-1].append(token)
                    except IndexError:
                        raise SelectionError(sel, loc) 
            else:
                if dtype == bool:
                    torfs.append(token)
                else:
                    try:
                        evals[-1].append(token)
                    except IndexError:
                        raise SelectionError(sel, loc)

        torf = None
        if torfs:
            torf = torfs.pop(0)
            while torfs:
                subset = torf.nonzero()[0]
                if not len(subset): return torf
                torf[subset] = torfs.pop(0)[subset]

        if flags:
            if torf is None:
                torf = atoms.getFlags(flags.pop(0))
            while flags:
                subset = torf.nonzero()[0]
                if not len(subset): return torf
                torf[subset] = atoms._getFlags(flags.pop(0))[subset]

        if evals:
            if torf is None:
                torf, err = self._eval(sel, loc, evals.pop(0))
                if err: raise err
            while evals:
                subset = torf.nonzero()[0]
                if not len(subset): return torf
                arr, err = self._eval(sel, loc, evals.pop(0), subset=subset)
                if err: raise err
                torf[subset] = arr
        # ?? check torf.shape/ndim
        return torf
    
    def _unary(self, sel, loc, tokens):
        """Perform the unary operation."""

        debug(sel, loc, '_unary', tokens)
        if NUMB: return
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
        which = int(tokens[0] == NOT)
        if isinstance(tokens[which], ndarray):
            torf = tokens[which]
        else:
            torf = self._evaluate(sel, loc, tokens[which:], evalonly=evalonly)
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
            which = arange(len(coords))
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
            get_indices = kdtree.getIndices
            search = kdtree.search
            get_count = kdtree.getCount
            torf = zeros(self._ag.numAtoms(), bool)
            for index in which:
                search(within, coords[index])
                if get_count():
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
            kdtree = KDTree(coords[which])
            search = kdtree.search
            get_count = kdtree.getCount
            select = []
            append = select.append
            for i, xyz in enumerate(cxyz):
                search(within, xyz)
                if get_count():
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

        torf = token[1]
        if not isinstance(torf, ndarray):
            torf = self._evaluate(sel, loc, token[1:])
            if isinstance(torf, SelectionError):
                return torf
        which = torf.nonzero()[0]
        if not len(which):
            return torf
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
    
    def _getNumeric(self, sel, loc, arg, copy=False):
        """Return numeric data or a number."""

        debug(sel, loc, '_getNumeric', arg)

        # arg may be an array, a string, or a regular expression 
        try:
            dtype, ndim = arg.dtype, arg.ndim
        except AttributeError:
            pass
        else:
            # i don't expect that a string array may show up
            if dtype == bool:
                return None, SelectionError(sel, loc, 'operands and function '
                    'arguments must be numbers or numeric data labels')
            else:
                return arg, False
        
        # no regular expressions
        try:
            pattern = arg.pattern
        except AttributeError:
            pass
        else:
            return None, SelectionError(sel, sel.index(pattern, loc),
               'operands and function arguments cannot be regular expressions')
        
        if arg in XYZMAP:
            coords = self._getCoordsNew() # how about atoms._getCoords() ?
            if coords is None:
                return None, SelectionError(sel, loc, 
                    'coordinates are not set')
            else:
                if copy:
                    return coords[:, XYZMAP[arg]].copy(), False
                else:
                    return coords[:, XYZMAP[arg]], False
        
        try:
            if copy:
                data = self._atoms.getData(arg)
            else:
                data = self._atoms._getData(arg)
        except Exception as err:
            return None, SelectionError(sel, loc, 'following exception '
                        'occurred when evaluating {0:s}: {1:s}'
                        .format(repr(arg), str(err)))

        if data is not None:
            if str(data.dtype).startswith('|S'):
                return None, SelectionError(sel, loc, '{0:s} is not a numeric '
                        'data label'.format(repr(arg)))
            else:
                return data, False

        if arg == 'index':
            try:
                if copy:
                    return self._atoms.getIndices(), False
                else:
                    return self._atoms._getIndices(), False
            except AttributeError:
                return arange(self._atoms.numAtoms()), False
            
        try:
            return float(arg), False
        except Exception as err:
            return None, SelectionError(sel, loc, '{0:s} is not a number or a '
                    'numeric data label'
                    .format(repr(arg), ))
            
    def _comp(self, sel, loc, tokens):
        """Perform comparison."""
        
        tokens = tokens[0]
        debug(sel, loc, '_comp', tokens)
        if NUMB: return
        
        if len(tokens) >= 3 and len(tokens) % 2 != 1:
            raise SelectionError(sel, loc, 
                                 'invalid number of operators and operands')
        token = tokens.pop(0)
        left, err = self._getNumeric(sel, loc, token)
        if err: raise err
        
        torf = None
        while tokens:
            try:
                binop = BINOP_MAP[tokens.pop(0)]
            except KeyError:
                raise SelectionError(sel, loc, 'invalid operator encountered')
            
            right, err = self._getNumeric(sel, loc, tokens.pop(0))
            if err: raise err

            if torf is None:
                torf = binop(left, right)
            else:
                logical_and(binop(left, right), torf, torf)
        
        # check whether atomic data was contained in comparison
        # i.e. len(atoms) == len(torf)
        try:
            ndim, shape = torf.ndim, torf.shape
        except AttributeError:
            raise SelectionError(sel, loc, 
                                'comparison must contain atomic data')
        else:
            if ndim != 1 or shape[0] != self._atoms.numAtoms():
                raise SelectionError(sel, loc, 
                                    'comparison must contain atomic data')
            else:
                return torf
     
    def _binop(self, sel, loc, tokens):
        """Perform binary operation."""
        
        tokens = tokens[0]
        debug(sel, loc, '_binop', tokens)
        if NUMB: return
        
        if len(tokens) >= 3 and len(tokens) % 2 != 1:
            raise SelectionError(sel, loc, 'invalid number of items')
        left, err = self._getNumeric(sel, loc, tokens.pop(0), copy=True)
        if err: raise err

        while tokens:
            binop = tokens.pop(0)
            if binop not in BINOP_MAP:
                raise SelectionError(sel, loc, 'invalid operator encountered')
            
            right, err = self._getNumeric(sel, loc, tokens.pop(0))
            if err: raise err
            
            if DEBUG: print(binop, left, right)
            if binop == '/' and any(right == 0.0):
                raise SelectionError(sel, loc, 'zero division error')
            binop = BINOP_MAP[binop]

            try:
                ndim = left.ndim
            except:
                left = binop(left, right)
            else:
                # ndim must not be zero for in place operation
                if ndim:
                    binop(left, right, left)
                else:
                    left = binop(left, right)
        return left

    def _pow(self, sel, loc, tokens):
        """Perform power operation. Expected operands are numbers 
        and numeric atom attributes."""
        
        tokens = tokens[0]
        debug(sel, loc, '_pow', tokens)
        if NUMB: return

        base, err = self._getNumeric(sel, loc, tokens.pop(0))
        if err: raise err
        power, err = self._getNumeric(sel, loc, tokens.pop())
        if err: raise err
                
        if tokens.pop() not in ('^', '**'):
            raise SelectionError(sel, loc, 'invalid power operator')
        while tokens:
            number, err = self._getNumeric(sel, loc,   tokens.pop())
            if err: raise err
            power = number ** power
            if tokens.pop() not in ('^', '**'):
                raise SelectionError(sel, loc, 'invalid power operator')
        
        return base ** power

    def _sign(self, sel, loc, tokens):
        """Change the sign of a selection argument."""
        
        tokens = tokens[0]
        debug(sel, loc, '_sign', tokens)
        if NUMB: return
        
        if len(tokens) != 2:
            raise SelectionError(sel, loc, 
                                 'sign operators (+/-) must be followed '
                                 'by single keyword, e.g. "-x", "-beta"')
        arg, err = self._getNumeric(sel, loc, tokens[1])
        if err: raise err
        
        if tokens[0] == '-':
            arg =  -arg
        return arg

    def _func(self, sel, loc, tokens):
        """Evaluate functions used in selection strings."""
        
        tokens = list(tokens[0])
        debug(sel, loc, '_func', tokens)
        if NUMB: return
        
        if len(tokens) != 2:
            raise SelectionError(sel, loc, '{0:s} accepts a single numeric '
                             'argument, e.g. {0:s}(x)'.format(repr(tokens[0])))
        arg, err = self._getNumeric(sel, loc, tokens[1])
        if err: raise err
        debug(sel, loc, tokens[0], arg)
        return FUNCTION_MAP[tokens[0]](arg)

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
        elif self._ag.isDataLabel(token):
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

    def _generic(self, sel, loc, tokens, subset=None):
        
        debug(sel, loc, '_generic', tokens)

        label = tokens.pop(0)
        data = self._getData(sel, loc, label)
        if subset is not None:
            data = data[subset]
            subset = None
        dtype = data.dtype
        type_ = dtype.type
        torf = None
        regexp = []
        values = []
        ranges = []
        
        for token in tokens:

            # check for regular expressions which are only compatible with 
            # string data type
            try:
                token.pattern
            except AttributeError:
                pass
            else:
                if not dtype.str.startswith('|S'):
                    ptrn = '"{0:s}"'.format(token.pattern)
                    if sel.count(ptrn, loc) == 1: loc = sel.index(ptrn, loc) 
                    return None, SelectionError(sel, loc, '{0:s} is a regular '
                        'expression and is not compatible with type of {1:s} '
                        '({2:s})'.format(ptrn, repr(label), str(dtype)))
                else:
                    regexp.append(token)
                    continue
            
            
            # check for ranges
            if token[0] == 'range': 
                if dtype.str.startswith('|S'):
                    return None, SelectionError(sel, loc, 'number ranges '
                        'are not compatible with type of {0:s} ({1:s})'
                        .format(repr(label), str(dtype)))
                else:
                    if token[-1] in BINOP_MAP:
                        ranges.append(token)
                    else:                
                        values.extend(arange(*token[1:]))
                    continue
            if token == '_':
                token = ' '
            try:
                value = type_(token)
            except Exception as err:
                if sel.count(label, loc) == 1: loc = sel.index(label, loc)
                if sel.count(token, loc) == 1: loc = sel.index(token, loc)
                return None, SelectionError(sel, loc, '{0:s} could not be '
                    'converted to type of {1:s} ({2:s})'
                    .format(repr(token), repr(label), str(err)))
            values.append(value)

        if values:
            if len(values) > 4:
                torf = array([val in values for val in data], bool)
            else:
                torf = data == values.pop(0)
                for val in values:
                    subset = (torf == False).nonzero()[0]
                    if len(subset) == 0: return torf, False
                    torf[subset] = data[subset] == val
        if ranges:
            while ranges:
                _, start, stop, comp = ranges.pop(0)
                if torf is None:
                    torf = start <= data
                else:
                    subset = (torf == False).nonzero()[0]
                    if len(subset) == 0: return torf, False
                    subdata = data[subset]
                    torf[subset] = start <= subdata
                if subset is None:
                    torf = logical_and(torf, BINOP_MAP[comp](data, stop), torf)
                else: 
                    torf[subset] = logical_and(torf[subset], 
                                       BINOP_MAP[comp](subdata, stop))


        for re in regexp:
            if torf is None:
                torf = array([re.match(val) is not None for val in data], bool)
            else:
                subset = (torf == False).nonzero()[0]
                if len(subset) == 0: return torf, False
                torf[subset] = [re.match(val) is not None 
                                for val in data[subset]]
        return torf, False

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
                    value = re_compile(value)
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
            token = evalNumeric(sel, loc, token, False)
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
        
        numbers = evalNumeric(sel, loc, token)
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
            return (self._indices if self._indices is not None else 
                    arange(self._ag._n_atoms))
        torf = zeros(self._ag._n_atoms, np.bool)
        
        numbers = evalNumeric(sel, loc, token)
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
    
    def _getData(self, sel, loc, keyword):
        """Return atomic data."""
        
        try:
            idx = XYZMAP[keyword]
        except KeyError:
            pass
        else:
            return self._getCoords(sel, loc)[:,idx]
        
        data = self._data.get(keyword)
        if data is None:
            field = ATOMIC_FIELDS.get(FIELDS_SYNONYMS.get(keyword, keyword))
            if field is None:
                data = self._atoms._getData(keyword)
                if data is None:
                    return SelectionError(sel, loc, "{0:s} is not a valid "
                          "keyword or user data label".format(repr(keyword)))
                elif not isinstance(data, ndarray) and data.ndim == 1:
                    return SelectionError(sel, loc, "{0:s} must be a 1d "
                                          "array".format(repr(keyword)))
            else:
                try:
                    data = getattr(self._atoms, '_get' + field.meth_pl)()
                except Exception as err: 
                    return SelectionError(sel, loc, str(err))
                if data is None:
                    return SelectionError(sel, loc, "{0:s} is not set by "
                                          "user".format(repr(keyword)))
            self._data[keyword] = data
        return data
        #indices = self._indices
        #if indices is None or self._dummies is None:               
        #    return data
        #elif self._dummies is None:
        #    return data[indices]
    
    def _getCoords(self, sel, loc):
        """Return coordinates of selected atoms.  This method is reduces array 
        copy operations when :class:`.AtomPointer` subclasses are used for
        making atom selections."""
        
        if self._coords is None:
            self._coords = self._atoms._getCoords()
            if self._coords is None:
                return SelectionError(sel, loc, 'coordinates are not set')
        return self._coords

    def _getCoordsNew(self):
        """Return coordinates of atoms."""
        
        if self._coords is None:
            self._coords = self._atoms._getCoords()
        return self._coords

