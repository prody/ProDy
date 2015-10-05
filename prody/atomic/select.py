# -*- coding: utf-8 -*-
"""This module defines a class for selecting subsets of atoms.  You can read
this page in interactive sessions using ``help(select)``.

ProDy offers a fast and powerful atom selection class, :class:`.Select`.
Selection features, grammar, and keywords are similar to those of VMD.
Small differences, that is described below, should not affect most practical
uses of atom selections.  With added flexibility of Python, ProDy selection
engine can also be used to identify intermolecular contacts.  You may see
this and other usage examples in :ref:`contacts` and
:ref:`selection-operations`.

First, we import everything from ProDy and parse a protein-DNA-ligand
complex structure:

.. ipython:: python

   from prody import *
   p = parsePDB('3mht')

:func:`.parsePDB` returns :class:`.AtomGroup` instances, ``p`` in this case,
that stores all atomic data in the file.  We can count different types of
atoms using :ref:`flags` and :meth:`~.AtomGroup.numAtoms` method as follows:

.. ipython:: python

   p.numAtoms('protein')
   p.numAtoms('nucleic')
   p.numAtoms('hetero')
   p.numAtoms('water')

Last two counts suggest that ligand has 26 atoms, i.e. number of :term:`hetero`
atoms less the number of :term:`water` atoms.


Atom flags
-------------------------------------------------------------------------------

We select subset of atoms by using :meth:`.AtomGroup.select` method.
All :ref:`flags` can be input arguments to this methods as follows:

.. ipython:: python

   p.select('protein')
   p.select('water')

This operation returns :class:`.Selection` instances, which can be an input
to functions that accepts an *atoms* argument.


Logical operators
-------------------------------------------------------------------------------

Flags can be combined using ``'and'`` and ``'or'`` operators:

.. ipython:: python

   p.select('protein and water')

``'protein and water'`` did not result in selection of :term:`protein` and
:term:`water` atoms.  This is because, no atom is flagged as a protein and a
water atom at the same time.

.. note::

   **Interpreting selection strings**

   You may think as if a selection string, such as ``'protein and water'``, is
   evaluated on a per atom basis and an atom is selected if it satisfies the
   given criterion.  To select both water and protein atoms, ``'or'`` logical
   operator should be used instead.  A protein or a water atom would satisfy
   ``'protein or water'`` criterion.

.. ipython:: python

   p.select('protein or water')

We can also use ``'not'`` operator to negate an atom flag.  For example,
the following selection will only select ligand atoms:

.. ipython:: python

   p.select('not water and hetero')

If you omit the ``'and'`` operator, you will get the same result:

.. ipython:: python

   p.select('not water hetero')

.. note::

   **Default operator** between two flags, or other selection tokens that will
   be discussed later, is ``'and'``.  For example, ``'not water hetero'``
   is equivalent to ``'not water and hetero'``.

We can select Cα atoms of acidic residues by omitting the default logical
operator as follows:

.. ipython:: python

   sel = p.select('acidic calpha')
   sel
   set(sel.getResnames())


Quick selections
-------------------------------------------------------------------------------

For simple selections, such as shown above, following may be preferable over
the :meth:`~.AtomGroup.select` method:

.. ipython:: python

   p.acidic_calpha

The result is the same as using ``p.select('acidic calpha')``.  Underscore,
``_``, is considered as a whitespace.  The limitation of this approach is that
special characters cannot be used.


Atom data fields
-------------------------------------------------------------------------------

In addition to :ref:`flags`, :ref:`fields` can be used in atom selections
when combined with some values.  For example, we can select Cα and Cβ atoms
of alanine residues as follows:

.. ipython:: python

   p.select('resname ALA name CA CB')

Note that we omitted the default ``'and'`` operator.

.. note::

   **Whitespace** or **empty string** can be specified using an ``'_'``.
   Atoms with string data fields empty, such as those with no a chain
   identifiers or alternate location identifiers, can be selected using
   an underscore.


.. ipython:: python

   p.select('chain _')  # chain identifiers of all atoms are specified in 3mht
   p.select('altloc _')  # altloc identifiers for all atoms are empty

Numeric data fields can also be used to make selections:

.. ipython:: python

   p.select('ca resnum 1 2 3 4')


A special case for residues is having insertion codes.  Residue numbers and
insertion codes can be specified together as follows:

  * ``'resnum 5'`` selects residue 5 (all insertion codes)
  * ``'resnum 5A'`` selects residue 5 with insertion code A
  * ``'resnum 5_'`` selects residue 5 with no insertion code


Number ranges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A range of numbers using ``'to'`` or Python style slicing with ``':'``:

.. ipython:: python

   p.select('ca resnum 1to4')
   p.select('ca resnum 1:4')
   p.select('ca resnum 1:4:2')

.. note::

   **Number ranges** specify continuous intervals:

     * ``'to'`` is all inclusive, e.g. ``'resnum 1 to 4'`` means
       ``'1 <= resnum <= 4'``

     * ``':'`` is left inclusive, e.g. ``'resnum 1:4'`` means
       ``'1 <= resnum < 4'``

   Consecutive use of ``':'``, however, specifies a discrete range of numbers,
   e.g. ``'resnum 1:4:2'`` means ``'resnum 1 3'``


Special characters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following characters can be specified when using :ref:`fields` for atom
selections::

  abcdefghijklmnopqrstuvwxyz
  ABCDEFGHIJKLMNOPQRSTUVWXYZ
  0123456789
  ~@#$.:;_',

For example, ``"name C' N` O~ C$ C#"`` is a valid selection string.

.. note::

   **Special characters** (``~!@#$%^&*()-_=+[{}]\|;:,<>./?()'"``) must be
   escaped using grave accent characters (``````).


Negative numbers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Negative numbers and number ranges must also be escaped using grave accent
characters, since negative sign ``'-'`` is considered a special character
unless it indicates subtraction operation (see below).

.. ipython:: python

   p.select('x `-25 to 25`')
   p.select('x `-22.542`')


Omitting the grave accent character will cause a :exc:`.SelectionError`.


Regular expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, you can specify regular expressions to select atoms based on
data fields with type string.  Following will select residues whose names
start with capital letter A

.. ipython:: python

   sel = p.select('resname "A.*"')
   set(sel.getResnames())


.. note::

   **Regular expressions** can be specified using double quotes, ``"..."``.
   For more information on regular expressions see :mod:`re`.


Numerical comparisons
-------------------------------------------------------------------------------

:ref:`fields` with numeric types can be used as operands in numerical
comparisons:

.. ipython:: python

   p.select('x < 0')
   p.select('occupancy = 1')

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

It is also possible to chain comparison statements as follows:

.. ipython:: python

   p.select('-10 <= x < 0')

This would be the same as the following selection:

.. ipython:: python

   p.select('-10 <= x and x < 0') == p.select('-10 <= x < 0')

Furthermore, numerical comparisons may involve the following operations:

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

.. ipython:: python

   p.select('x ** 2 < 10')
   p.select('x ** 2 ** 2 < 10')

Finally, following functions can be used in numerical comparisons:

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

.. ipython:: python

   p.select('sqrt(sq(x) + sq(y) + sq(z)) < 100')  # within 100 Å of origin


Distance based selections
-------------------------------------------------------------------------------

Atoms within a user specified distance (Å) from a set of user specified atoms
can be selected using ``'within . of .'`` keyword, e.g. ``'within 5 of water'``
selects atoms that are within 5 Å of water molecules. This setting will
results selecting water atoms as well.

User can avoid selecting specified atoms using ``exwithin . of ..`` setting,
e.g. ``'exwithin 5 of water'`` will not select water molecules and is
equivalent to ``'within 5 of water and not water'``

.. ipython:: python

   p.select('exwithin 5 of water') == p.select('not water within 5 of water')


Sequence selections
-------------------------------------------------------------------------------

One-letter amino acid sequences can be used to make atom selections.
``'sequence SAR'`` will select **SER-ALA-ARG** residues in a chain.  Note
that the selection does not consider connectivity within a chain.  Regular
expressions can also be used to make selections: ``'sequence "MI.*KQ"'`` will
select **MET-ILE-(XXX)n-ASP-LYS-GLN** pattern, if present.

.. ipython:: python

   sel = p.select('ca sequence "MI.*DKQ"')
   sel
   sel.getResnames()


Expanding selections
-------------------------------------------------------------------------------

A selection can be expanded to include the atoms in the same *residue*,
*chain*, or *segment* using ``same .. as ..`` setting, e.g.
``'same residue as exwithin 4 of water'`` will select residues that have
at least an atom within 4 Å of any water molecule.

.. ipython:: python

   p.select('same residue as exwithin 4 of water')

Additionally, a selection may be expanded to the immediately bonded atoms using
``bonded [n] to ...`` setting, e.f. ``bonded 1 to calpha`` will select atoms
bonded to Cα atoms.  For this setting to work, bonds must be set by the user
using the :meth:`.AtomGroup.setBonds` method.  It is also possible to select
bonded atoms by excluding the originating atoms using ``exbonded [n] to ...``
setting.  Number ``'[n]'`` indicates number of bonds to consider from the
originating selection and defaults to 1.


Selection macros
-------------------------------------------------------------------------------

ProDy allows you to define a macro for any valid selection string.  Below
functions are for manipulating selection macros:

  * :func:`defSelectionMacro`
  * :func:`delSelectionMacro`
  * :func:`getSelectionMacro`
  * :func:`isSelectionMacro`

.. ipython:: python

   defSelectionMacro('alanine', 'resname ALA')
   p.select('alanine') == p.select('resname ALA')

You can also use this macro as follows:

.. ipython:: python

   p.alanine

Macros are stored in ProDy configuration file permanently.  You can delete
them if you wish as follows:

.. ipython:: python

   delSelectionMacro('alanine')


Keyword arguments
-------------------------------------------------------------------------------

:meth:`~.Select.select` method also accepts keyword arguments that can simplify
some selections.  Consider the following case where you want to select some
protein atoms that are close to its center:

.. ipython:: python

   protein = p.protein
   calcCenter(protein).round(2)
   sel1 = protein.select('sqrt(sq(x--21.17) + sq(y-35.86) + sq(z-79.97)) < 5')
   sel1

Instead, you could pass a keyword argument and use the keyword in the
selection string:

.. ipython:: python

   sel2 = protein.select('within 5 of center', center=calcCenter(protein))
   sel2
   sel1 == sel2

Note that selection string for *sel2* lists indices of atoms.  This
substitution is performed automatically to ensure reproducibility of the
selection without the keyword *center*.

Keywords cannot be reserved words (see :func:`.listReservedWords`) and must be
all alphanumeric characters."""

import sys
from re import compile as re_compile

import numpy as np
from numpy import array, ndarray, ones, zeros, arange
from numpy import invert, unique, concatenate, all, any
from numpy import logical_and, logical_or, floor, ceil, where

try:
    from . import pyparsing as pp
    from .pyparsing import ParseException
except ImportError:
    import pyparsing as pp
    from pyparsing import ParseException


from prody import LOGGER, SETTINGS, PY2K

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

if PY2K:
    range = xrange

DEBUG = 0
NUMB = 0 # Select instance will not really evaluate string for the atoms
TIMER = 0

def debug(sel, loc, *args):

    if DEBUG:
        print('')
        if args:
            print(args[0], args[1:])
        print(repr(sel))
        print(' ' * (loc + 1) + '^')

__all__ = ['Select', 'SelectionError', 'SelectionWarning',
           'defSelectionMacro', 'delSelectionMacro', 'getSelectionMacro',
           'isSelectionMacro']

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

    .. ipython:: python

       defSelectionMacro('cbeta', 'name CB and protein')"""

    if not isinstance(name, str) or not isinstance(selstr, str):
        raise TypeError('both name and selstr must be strings')
    elif isReserved(name):
        raise ValueError('{0} is a reserved word and cannot be used as a '
                         'macro name'.format(repr(name)))
    elif not (name.isalpha() and name.islower()):
        raise ValueError('macro names must be all lower case letters, {0} '
                         'is not a valid macro name'.format(repr(name)))

    LOGGER.info('Testing validity of selection string.')
    try:
        ATOMGROUP.select(selstr)
    except SelectionError:
        LOGGER.warn('{0} is not a valid selection string, macro {1} is not'
                    'defined.'.format(repr(selstr), repr(name)))
    else:
        LOGGER.info("Macro {0} is defined as {1}."
                    .format(repr(name), repr(selstr)))
        MACROS[name] = selstr
        SETTINGS['selection_macros'] = MACROS
        SETTINGS.save()


def delSelectionMacro(name):
    """Delete the macro *name*.

    .. ipython:: python

       delSelectionMacro('cbeta')"""

    try:
        MACROS.pop(name)
    except:
        LOGGER.warn("Macro {0} is not found.".format(repr(name)))
    else:
        if MACROS_REGEX is not None: MACROS_REGEX.pop(name, None)
        LOGGER.info("Macro {0} is deleted.".format(repr(name)))
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
        LOGGER.info("{0} is not a user defined macro name."
                    .format(repr(name)))


def replaceMacros(selstr):

    global MACROS_REGEX
    if MACROS_REGEX is None: MACROS_REGEX = {}
    selstr = ' ' + selstr + ' '
    if MACROS:
        for name, macro in MACROS.items(): # PY3K: OK
            re = MACROS_REGEX.setdefault(name,
                re_compile('[( )]' + name + '[( )]'))

        for match in re.finditer(selstr):
            start, end = match.start(), match.end()
            match = selstr[start:end]
            selstr = (selstr[:start] + match[0] + '(' + macro+ ')' +
                      match[-1] + selstr[end:])
    return selstr[1:-1]


def checkSelstr(selstr, what, error=ValueError):
    """Check *selstr* if it satisfies a selected condition.  For now, only
    whether coordinate/distance based selection are checked.  If *error* is
    a subclass of :exc:`Exception`, an exception will be raised, otherwise
    return **True** or **False** will be returned."""

    selstr = selstr.replace('(', ' ( ')
    selstr = selstr.replace(')', ' ) ')

    if what in set(['dist']):
        for item in selstr.split():
            if item in XYZDIST:
                if issubclass(error, Exception):
                    raise error('invalid selection {0}, coordinate '
                                'based selections are not accepted'
                                .format(repr(selstr)))
                else:
                    return False


class SelectionError(Exception):

    """Exception raised when there are errors in the selection string."""

    def __init__(self, sel, loc=0, msg='', tkns=None):

        if tkns:
            for tkn in tkns:
                tkn = str(tkn)
                if sel.count(tkn, loc) == 1: loc = sel.index(tkn, loc)

        msg = ("An invalid selection string is encountered:\n{0}\n"
               .format(repr(sel)) +
               ' ' * (loc + 1) + '^ ' + msg)
        Exception.__init__(self, msg)


class SelectionWarning(Warning):

    """A class used for issuing warning messages when potential typos are
    detected in a selection string.  Warnings are issued to ``sys.stderr``
    via ProDy package logger.  Use :func:`.confProDy` to selection warnings
    *on* or *off*, e.g. ``confProDy(selection_warning=False)``."""

    def __init__(self, sel='', loc=0, msg='', tkns=None):

        if SETTINGS.get('selection_warning', True):

            if tkns:
                for tkn in tkns:
                    tkn = str(tkn)
                    if sel.count(tkn, loc) == 1: loc = sel.index(tkn, loc)

            msg = ('Selection string contains typo(s):\n'
                   '{0}\n '.format(repr(sel)) +
                   ' ' * loc + '^ ' + msg)
            LOGGER.warn(msg)

FIELDS_SYNONYMS = {'chid': 'chain',
                   'fragment': 'fragindex',
                   'resid': 'resnum',
                   'secstr': 'secondary',
                   'segname': 'segment'}


XYZ2INDEX = {'x': 0, 'y': 1, 'z': 2}

FUNCTIONS = {
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

OPERATORS = {
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

SAMEAS_MAP = {'residue': 'resindex', 'chain': 'chindex',
              'segment': 'segindex', 'fragment': 'fragindex'}

# maybe a value to a field/data label or field label
XYZ = set(['x', 'y', 'z'])
XYZDIST = set(['x', 'y', 'z', 'within', 'exwithin'])

OR = pp.Keyword('or')
AND = pp.Keyword('and')
WORD = pp.Word(pp.alphanums + '''~@#$.:;_'`,''')

_ = list(FUNCTIONS)
FUNCNAMES = set(_)

kwfunc = pp.Keyword(_[0])
FUNCNAMES_OPLIST = kwfunc
FUNCNAMES_EXPR = ~kwfunc
for func in _[1:]:
    kwfunc = pp.Keyword(func)
    FUNCNAMES_OPLIST = FUNCNAMES_OPLIST | kwfunc
    FUNCNAMES_EXPR += ~kwfunc

RE_SCHARS = re_compile('`[\w\W]*?`')
PP_SCHARS = pp.Regex(RE_SCHARS)


def specialCharsParseAction(sel, loc, token):

    token = token[0][1:-1]
    if not token:
        raise SelectionError(sel, loc, '`` is invalid, no special characters')
    if ':' in token or 'to' in token:
        try:
            token = PP_NRANGE.parseString(token)[0]
        except ParseException:
            pass
    return token

PP_SCHARS.setParseAction(specialCharsParseAction)


RE_REGEXP = re_compile('"[\w\W]*"')
PP_REGEXP = pp.Regex(RE_REGEXP.pattern)


def regularExpParseAction(sel, loc, token):

    token = token[0]
    if token == '""':
        raise SelectionError(sel, loc, '"" is invalid, no regular expression')
    try:
        regexp = re_compile(token[1:-1])
    except:
        raise SelectionError(sel, loc, 'failed to compile regular '
                             'expression {0}'.format(repr(token)))
    else:
        return regexp

PP_REGEXP.setParseAction(regularExpParseAction)

_ = '[-+]?\d+(\.\d*)?([eE]\d+)?'
RE_NRANGE = re_compile(_ + '\ *(to|:)\ *' + _)

PP_NRANGE = pp.Group(pp.Regex(RE_NRANGE.pattern) +
                     pp.Optional(pp.Regex('(\ *:\ *' + _ + ')')))


def rangeParseAction(sel, loc, tokens):

    tokens = tokens[0]
    debug(sel, loc, '_nrange', tokens)

    token = tokens[0]
    sep = ':' if ':' in token else 'to'
    first, last = token.split(sep)
    try:
        start = int(first)
    except ValueError:
        start = float(first)
    try:
        stop = int(last)
    except ValueError:
        stop = float(last)

    if start > stop or (start == stop and sep == ':'):
        raise SelectionError(sel, loc, 'range start value ({0}) is greater '
                             'than stop value ({1})'
                             .format(repr(start),
                                     repr(stop if sep != ':' else stop - 1)))
    elif start == stop:
        return first

    if sep == 'to':
        comp = '<='
    elif len(tokens) == 1:
        comp = '<'
    else:
        try:
            comp = int(tokens[1][1:])
        except ValueError:
            comp = float(tokens[1][1:])
    return 'range', start, stop, comp

PP_NRANGE.setParseAction(rangeParseAction)


UNARY = set(['not', 'bonded', 'exbonded', 'within', 'exwithin', 'same'])


class Select(object):

    """Select subsets of atoms based on a selection string.
    See :mod:`~.atomic.select` module documentation for selection grammar
    and examples.  This class makes use of pyparsing_ module."""

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

        self._parsers = {}


        self._evalmap = {'resnum': self._resnum, 'resid': self._resnum,
            'serial': self._serial, 'index': self._index,
            'x': self._generic, 'y': self._generic, 'z': self._generic,
            'chid': self._generic, 'secstr': self._generic,
            'fragment': self._generic, 'fragindex': self._generic,
            'segment': self._generic, 'sequence': self._sequence, }


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
        """Return a :class:`.Selection` of atoms matching *selstr*, or
        **None**, if selection string does not match any atoms.

        :arg atoms: atoms to be evaluated
        :type atoms: :class:`.Atomic`

        :arg selstr: selection string
        :type selstr: str

        Note that, if *atoms* is an :class:`.AtomMap` instance, an
        :class:`.AtomMap` is returned, instead of a a :class:`.Selection`.

        .. note:

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
                selstr = 'index {0}'.format(rangeString(indices))
            else:
                if self._replace:
                    for key, value in kwargs.items(): # PY3K: OK
                        if value in ag and key in selstr:
                            if value == ag:
                                ss = 'all'
                            else:
                                ss = value.getSelstr()
                            selstr = selstr.replace(key, '(' + ss + ')')
                if isinstance(atoms, AtomPointer):
                    selstr = '({0}) and ({1})'.format(selstr,
                                                      atoms.getSelstr())

            return Selection(ag, indices, selstr, atoms.getACSIndex(),
                             unique=True)
        else:
            return AtomMap(ag, indices, atoms.getACSIndex(), dummies=dummies,
               title='Selection {0} from '.format(repr(selstr)) + str(atoms))

    def getIndices(self, atoms, selstr, **kwargs):
        """Return indices of atoms matching *selstr*.  Indices correspond to
        the order in *atoms* argument.  If *atoms* is a subset of atoms, they
        should not be used for indexing the corresponding :class:`.AtomGroup`
        instance."""

        ss = selstr.strip()
        if (len(ss.split()) == 1 and ss.isalnum() and ss not in MACROS):
            self._evalAtoms(atoms)
            if ss == 'none':
                return array([])
            elif ss == 'all':
                return arange(atoms.numAtoms())
            elif atoms.isFlagLabel(ss):
                return atoms._getFlags(ss).nonzero()[0]
            elif atoms.isDataLabel(ss) or ss in self._evalmap:
                raise SelectionError(selstr, 0, 'must be followed by values',
                                     [ss])
            else:
                raise SelectionError(selstr, 0, 'is not a valid selection '
                                     'string', [ss])
        else:
            torf = self.getBoolArray(atoms, selstr, **kwargs)
            return torf.nonzero()[0]

    def getBoolArray(self, atoms, selstr, **kwargs):
        """Return a boolean array with **True** values for *atoms* matching
        *selstr*.  The length of the boolean :class:`numpy.ndarray` will be
        equal to the length of *atoms* argument."""

        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance, not {0}'
                            .format(type(atoms)))

        self._reset()

        for key in kwargs.keys():
            if not key.isalnum():
                raise TypeError('{0} is not a valid keyword argument, '
                                  'keywords must be all alpha numeric '
                                  'characters'.format(repr(key)))
            if isReserved(key):
                loc = selstr.find(key)
                if loc > -1:
                    raise SelectionError(selstr, loc, '{0} is a reserved '
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
            elif atoms.isDataLabel(selstr):
                raise SelectionError(selstr, 0, 'must be followed by values')
            else:
                raise SelectionError(selstr, 0, 'is not a valid selection or '
                                     'user data label')

        selstr = replaceMacros(selstr)
        try:
            parser = self._getParser(selstr)
            tokens = parser(selstr, parseAll=True)
        except pp.ParseException as err:
            self._parsers.pop(self._parser, None)
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

    def _getParser(self, selstr):
        """Return an efficient parser that can handle *selstr*."""

        alnum = selstr
        alpha = selstr
        for ch in selstr:
            if not ch.isalnum(): alnum = alnum.replace(ch, ' ')
            if not ch.isalpha(): alpha = alpha.replace(ch, ' ')
        items = set(alnum.split())
        chars = set(selstr)


        funcs = 4 if items.intersection(FUNCNAMES) else 0
        opers = 2 if chars.intersection(OPERATORS) else 0
        logic = 1 if 'or' in items or '(' in chars else 0

        schars = 8 if '`' in chars and RE_SCHARS.search(selstr) else 0
        regexp = 16 if '"' in chars and RE_REGEXP.search(selstr) else 0
        nrange = 32 if ((':' in chars or ' to ' in alpha) and
                        RE_NRANGE.search(selstr)) else 0

        self._parser = key = (logic + opers + funcs,
                              logic + funcs + schars + regexp + nrange)

        if key == (0, 0):
            return self._noParser

        try:
            return self._parsers[key][0].parseString
        except KeyError:
            pass


        word = ~AND + ~OR

        oplist = []
        if funcs:
            oplist.append((FUNCNAMES_OPLIST, 1, pp.opAssoc.RIGHT, self._func))
            # following causes 20% slow down
            #word += FUNCNAMES_EXPR

        if funcs or opers:
            oplist.extend([
                (pp.oneOf('+ -'), 1, pp.opAssoc.RIGHT, self._sign),
                (pp.oneOf('** ^'), 2, pp.opAssoc.LEFT, self._pow),
                (pp.oneOf('* / %'), 2, pp.opAssoc.LEFT, self._binop),
                (pp.oneOf('+ -'), 2, pp.opAssoc.LEFT, self._binop),
                (pp.oneOf('< > <= >= == = !='), 2, pp.opAssoc.LEFT,
                 self._comp)])

        oplist.extend([
          (pp.Optional(AND), 2, pp.opAssoc.LEFT, self._and),
          (OR, 2, pp.opAssoc.LEFT, self._or)])

        word += WORD

        expr = word
        if schars: expr = PP_SCHARS | expr
        if regexp: expr = PP_REGEXP | expr
        if nrange: expr = PP_NRANGE | expr

        parser = pp.operatorPrecedence(expr, oplist)
        parser.setParseAction(self._default)
        parser.leaveWhitespace()
        parser.enablePackrat()
        self._parsers[key] = parser, expr, oplist
        return parser.parseString

    def _noParser(self, selstr, parseAll=True):

        debug(selstr, 0, ['_noParser'])
        return [self._default(selstr, 0, selstr.split())]

    def _getZeros(self, subset=None):
        """Return a bool array with zero elements."""

        if subset is None:
            return zeros(self._atoms.numAtoms(), bool)
        else:
            return zeros(len(subset), bool)

    def _default(self, sel, loc, tokens):

        debug(sel, loc, '_default', tokens)
        if NUMB: return

        if len(tokens) == 1:
            torf, err = self._eval(sel, loc, tokens)
        else:
            torf, err = self._and2(sel, loc, tokens)
        if err: raise err
        return torf


    def _eval(self, sel, loc, tokens, subset=None):

        debug(sel, loc, '_eval', tokens)
        if NUMB: return
        #if isinstance(tokens, ndarray):
        #    return tokens

        if len(tokens) == 1:

            token = tokens[0]
            try:
                dtype = token.dtype
            except AttributeError:
                pass
            else:
                return token, False

            if token == 'none':
                return zeros(self._n_atoms if subset is None else len(subset),
                              bool), False
            elif token == 'all':
                return ones(self._n_atoms if subset is None else len(subset),
                            bool), False

            elif self._atoms.isFlagLabel(token):
                return self._getFlags(token, subset), False

            elif self._atoms.isDataLabel(token):
                data, err = self._getData(sel, loc, token)
                if subset is None:
                    return data, False
                else:
                    return None, SelectionError(sel, loc, 'subset??')
                    return data[subset], False

            try:
                arg = self._kwargs[token]
            except KeyError:
                try:
                    return float(token), False
                except ValueError:

                    data, err = self._getData(sel, loc, token)
                    if data is not None:
                        return data, False

                    if token in ATOMIC_FIELDS:
                        return None, SelectionError(sel, loc, '{0} data is '
                            'not found'.format(repr(token)), [token])
                    else:
                        return None, SelectionError(sel, loc, '{0} could '
                            'not be evaluated'.format(repr(token)), [token])
            else:
                if arg in self._ag:
                    try:
                        dummies = arg.numDummies()
                    except AttributeError:
                        try:
                            indices = arg._getIndices()
                        except AttributeError:
                            indices = arange(self._ag.numAtoms())
                    else:
                        if dummies:
                            indices = arg._getIndices()[arg.getFlags('mapped')]
                        else:
                            indices = arg._getIndices()

                    torf = zeros(self._ag.numAtoms(), bool)
                    torf[indices] = True
                    if self._indices is not None:
                        torf = torf[self._indices]
                    self._replace = True
                    return torf, False
                return token, False
        else:
            return self._evalmap.get(tokens[0], self._generic)(
                                            sel, loc, tokens, subset)

    def _getFlags(self, label, subset):

        if subset is None:
            # get a copy to avoid alterations
            return self._atoms.getFlags(label)
        else:
            return self._atoms._getFlags(label)[subset]

    def _or(self, sel, loc, tokens):

        debug(sel, loc, '_or', tokens)
        tokens = tokens[0]
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
                if token == 'or':
                    continue
                elif isFlagLabel(token):
                    flags.append(token)
                elif isDataLabel(token) or token in self._evalmap:
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
                ss = where(torf == 0)[0]
                if len(ss) == 0: return torf
                torf[ss] = torfs.pop(0)[ss]

        if flags:
            if torf is None:
                torf = atoms.getFlags(flags.pop(0))
            while flags:
                ss = where(torf == 0)[0]
                if len(ss) == 0: return torf
                torf[ss] = atoms._getFlags(flags.pop(0))[ss]

        if evals:
            if torf is None:
                tokens = evals.pop(0)
                first = str(tokens[0])
                torf, err = self._eval(sel, loc, tokens)
                if err: raise err
                try:
                    dtype = torf.dtype
                except AttributeError:
                    raise SelectionError(sel, loc, 'a problem '
                        'occurred when evaluating token {0}'
                        .format(repr(first)), [first])

                else:
                    if dtype != bool:
                        raise SelectionError(sel, loc, 'a problem '
                            'occurred when evaluating token {0}'
                            .format(repr(first)), [first])
            while evals:
                ss = where(torf == 0)[0]
                if len(ss) == 0: return torf
                tokens = evals.pop(0)
                first = str(tokens[0])
                arr, err = self._eval(sel, loc, tokens, subset=ss)
                if err: raise err

                try:
                    dtype = arr.dtype
                except AttributeError:
                    raise SelectionError(sel, loc, 'a problem '
                        'occurred when evaluating token {0}'
                        .format(repr(first)), [first])

                else:
                    if dtype != bool:
                        raise SelectionError(sel, loc, 'a problem '
                            'occurred when evaluating token {0}'
                            .format(repr(first)), [first])

                torf[ss] = arr

        return torf

    def _and(self, sel, loc, tokens):

        debug(sel, loc, '_and', tokens)
        tokens = tokens[0]
        torf, err = self._and2(sel, loc, tokens)
        if err: raise err
        return torf

    def _and2(self, sel, loc, tokens, subset=None):

        debug(sel, loc, '_and2', tokens)
        if NUMB: return

        if tokens[0] == 'and' or tokens[-1] == 'and':
            return None, SelectionError(sel, loc, '{0} operator must be '
                'surrounded with arguments'.format(repr('and')), [tokens[0]])

        flags = []
        torfs = []
        evals = []
        unary = []
        atoms = self._atoms
        isFlagLabel = atoms.isFlagLabel
        isDataLabel = atoms.isDataLabel
        append = None
        wasand = False
        while tokens:
            # check whether token is an array to avoid array == str comparison
            token = tokens.pop(0)
            try:
                dtype = token.dtype
            except AttributeError:
                if token == 'and':
                    if wasand:
                        return None, SelectionError(sel, loc, 'incorrect use '
                            'of `and` operator, expected {0}'
                            .format(repr('and ... and')), ['and', 'and'])
                    append = None
                    wasand = True
                    continue

                elif isFlagLabel(token):
                    flags.append(token)
                    append = None

                elif (isDataLabel(token) or token in self._evalmap or
                      token in ATOMIC_FIELDS):
                    # not evals, must start a new list
                    # last evals list must have more than one values
                    #
                    if not evals:
                        evals.append([])
                        append = evals[-1].append
                    elif len(evals[-1]) == 1:
                        if token in XYZ:
                            pass
                        else:
                            return None, SelectionError(sel, loc, '{0} '
                                'must be followed by values'
                                .format(evals[-1][0]), [evals[-1][0]])
                    elif token in XYZ:
                        if (not tokens or tokens[0] in self._evalmap or
                              tokens[0] in ATOMIC_FIELDS or
                              isDataLabel(tokens[0]) or tokens[0] == 'and'):
                            pass
                        else:
                            try:
                                float(tokens[0])
                            except ValueError as err:
                                pass
                            else:
                                evals.append([])
                                append = evals[-1].append
                    elif not tokens:
                        return None, SelectionError(sel, loc, '{0} '
                            'must be followed by values'
                            .format(token), [token])
                    else:
                        evals.append([])
                        append = evals[-1].append
                    append(token)

                elif token in UNARY:
                    unary.append([])
                    append = unary[-1].append

                    if token == 'not':
                        append((token,))

                    elif token == 'same':
                        if len(tokens) < 3 or tokens[1] != 'as':
                            return None, SelectionError(sel, loc, 'incorrect '
                                'use of `same as` statement, expected {0}'
                                .format('same entity as ...'), [token])
                        append((token, tokens.pop(0), tokens.pop(0)))

                    elif token.endswith('within'):
                        if len(tokens) < 3 or tokens[1] != 'of':
                            return None, SelectionError(sel, loc, 'incorrect '
                                'use of `within` statement, expected {0}'
                                .format('[ex]within x.y of ...'), [token])
                        append((token, tokens.pop(0), tokens.pop(0)))

                    elif token.endswith('bonded'):
                        token2 = tokens.pop(0)
                        if len(tokens) < (1 + int(token2 == 'to')):
                            return None, SelectionError(sel, loc, 'incorrect '
                                'use of `bonded` statement, expected {0}'
                                .format('[ex]bonded [n] to ...'), [token2])
                        if token2 == 'to':
                            append((token, 'to'))
                        else:
                            append((token, token2, tokens.pop(0)))

                    anyargs = False
                    while tokens:
                        next = tokens[0]
                        try:
                            dtype = next.dtype
                        except AttributeError:
                            append(tokens.pop(0))
                            if next == 'and' or isFlagLabel(next):
                                break
                            if isDataLabel(next) or next in self._evalmap:
                                if anyargs:
                                    break
                                else:
                                    anyargs = True
                        else:
                            append(tokens.pop(0))
                            break
                else:
                    try:
                        append(token)
                    except TypeError:
                        return None, SelectionError(sel, loc, 'a problem '
                                    'occurred when evaluation token {0}'
                                    .format(repr(token)), [token])
            else:
                if dtype == bool:
                    torfs.append(token)
                else:
                    return None, SelectionError(sel, loc, 'a problem '
                                'occurred when evaluation token {0}'
                                .format(repr(token)), [token])
            wasand = False
        torf = None

        if torfs:
            torf = torfs.pop(0)

            while torfs:
                ss = torf.nonzero()[0]
                if len(ss) == 0: return torf, False
                torf[ss] = torfs.pop(0)[ss]

        if flags:
            if torf is None:
                torf = atoms.getFlags(flags.pop(0))

            while flags:
                ss = torf.nonzero()[0]
                if len(ss) == 0: return torf, False
                torf[ss] = atoms._getFlags(flags.pop(0))[ss]

        if unary:
            if torf is None:
                torf, err = self._unary(sel, loc, unary.pop(0))
                if err: return None, err

            while unary:
                ss = torf.nonzero()[0]
                if len(ss) == 0: return torf, False
                arr, err = self._unary(sel, loc, unary.pop(0))
                if err: return None, err
                torf[ss] = arr[ss]

        if evals:
            if torf is None:
                tokens = evals.pop(0)
                first = str(tokens[0])
                torf, err = self._eval(sel, loc, tokens, subset=subset)
                if err: return None, err
                try:
                    dtype = torf.dtype
                except AttributeError:
                    return None, SelectionError(sel, loc, 'a problem '
                        'occurred when evaluating token {0}'
                        .format(repr(first)), [first])

                else:
                    if dtype != bool:
                        return None, SelectionError(sel, loc, 'a problem '
                            'occurred when evaluating token {0}'
                            .format(repr(first)), [first])
            while evals:
                ss = torf.nonzero()[0]
                if len(ss) == 0: return torf, False
                tokens = evals.pop(0)
                first = str(tokens[0])
                arr, err = self._eval(sel, loc, tokens, subset=ss)
                if err: return None, err

                try:
                    dtype = arr.dtype
                except AttributeError:
                    return None, SelectionError(sel, loc, 'a problem '
                        'occurred when evaluating token {0}'
                        .format(repr(first)), [first])

                else:
                    if dtype != bool:
                        return None, SelectionError(sel, loc, 'a problem '
                            'occurred when evaluating token {0}'
                            .format(repr(first)), [first])

                torf[ss] = arr

        # ?? check torf.shape/ndim
        if subset is None:
            return torf, False
        else:
            return torf[subset], False

    def _unary(self, sel, loc, tokens):

        debug(sel, loc, '_unary', tokens)
        if NUMB: return

        what = tokens[0]
        which = tokens[1:]

        if not which:
            return None, SelectionError(sel, loc, '{0} must be followed by '
                .format(repr(' '.join(what))), what)
        if len(which) == 1:
            which, err = self._eval(sel, loc, which)
        else:
            which, err = self._and2(sel, loc, which)
        if err: raise err

        tokens = [what, which]
        if what[0] == 'not':
            return self._not(sel, loc, tokens)
        elif what[0] == 'same':
            return self._sameas(sel, loc, tokens)
        elif what[-1] == 'to':
            return self._bondedto(sel, loc, tokens)
        else:
            return self._within(sel, loc, tokens)

    def _not(self, sel, loc, tokens):
        """Negate selection."""

        debug(sel, loc, '_not', tokens)
        label, torf = tokens
        return invert(torf, torf), False

    def _within(self, sel, loc, tokens):
        """Perform distance based selection."""

        if DEBUG: print('_within', tokens)
        label, which = tokens
        within = label[1]
        label = ' '.join(label)
        try:
            within = float(within)
        except Exception as err:
            return None, SelectionError('could not convert {0} in {1} to '
                'float ({2})'.format(within, repr(label), str(err)),
                [label, within])
        exclude = label.startswith('ex')
        other = False
        try:
            dtype = which.dtype
        except AttributeError:

            if which in self._kwargs:
                coords = self._kwargs[which]
                try:
                    ndim, shape = coords.ndim, coords.shape
                except AttributeError:
                    try:
                        coords = coords._getCoords()
                    except AttributeError:
                        try:
                            coords = coords.getCoords()
                        except AttributeError:
                            return None, SelectionError(sel, loc,
                                '{0} must be a coordinate array or have '
                                '`getCoords` method'.format(repr(which)),
                                [label, which])
                    if coords is None:
                        return None, SelectionError(sel, loc,
                            'coordinates are not set for {0} ({1})'
                            .format(repr(which), repr(self._kwargs[which])),
                            [label, which])
                    else:
                        ndim, shape = coords.ndim, coords.shape
                if ndim == 1 and shape[0] == 3:
                    coords = array([coords])
                elif not (ndim == 2 and shape[1] == 3):
                    return None, SelectionError(sel, loc,
                        '{0} must be a coordinate array or have '
                        '`getCoords` method'.format(repr(which)),
                        [label, which])
                exclude=False
                self._ss2idx = True
                which = arange(len(coords))
                other = True
        else:
            if dtype == bool:
                which = which.nonzero()[0]
                coords = self._getCoords()
                if coords is None:
                    return None, SelectionError(sel, loc, 'coordinates are '
                                                'not set')
            else:
                return None, SelectionError(sel, loc, 'not understood')

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
            n_atoms = self._ag.numAtoms()
            torf = ones(n_atoms, bool)
            torf[which] = False
            check = torf.nonzero()[0]
            torf = zeros(n_atoms, bool)

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

        return torf, False

    def _sameas(self, sel, loc, tokens):
        """Evaluate ``'same entity as ...'`` expression."""

        debug(self, loc, '_sameas', tokens)
        label, which = tokens
        what = label[1]
        label = ' '.join(label)
        index = SAMEAS_MAP.get(what)
        if index is None:
            return None, SelectionError(sel, loc, 'entity in "same ... as" '
                'must be one of "chain", "residue", "segment", or "fragment",'
                ' not {0}'.format(repr(what)), [label])

        indices, err = self._getData(sel, loc, index)
        iset = set(indices[which])
        torf = array([i in iset for i in indices], bool)

        return torf, False

    def _bondedto(self, sel, loc, tokens):
        """Expand selection to immediately bonded atoms."""

        debug(sel, loc, '_bondedto', tokens)
        label, torf = tokens
        token = label[1]
        label = ' '.join(label)
        if token == 'to':
            repeat = 1
        else:
            try:
                repeat = int(token)
            except TypeError:
                return None, SelectionError(sel, loc, '{0} in {0} could not '
                    'be converted to an integer'.format(token, repr(label)),
                    [label])
            else:
                if float(token) != repeat:
                    SelectionWarning(sel, loc, 'number in {0} should be an '
                        'integer'.format(repr(label)), [label])

            if repeat <= 0:
                SelectionWarning(sel, loc, 'number in {0} should be a '
                    'positive integer'.format(repr(label)), [label])
                return zeros(self._atoms.numAtoms(), bool), False

        bmap = self._ag._bmap
        if bmap is None:
            return None, SelectionError(sel, loc, 'bonds are not set',
                                        [label])
        which = torf.nonzero()[0]
        if not len(which):
            return torf, False

        indices = self._indices
        if indices is not None:
            bmap = bmap[indices]
        n_atoms = self._ag.numAtoms()
        for i in range(repeat):
            torf = zeros(n_atoms, bool)
            bonded = unique(bmap[which])
            if bonded[0] == -1:
                torf[bonded[1:]] = True
            else:
                torf[bonded] = True
            if indices is not None:
                torf = torf[indices]
            if label.startswith('ex'):
                torf[which] = False
            else:
                torf[which] = True
            if i + 1 < repeat:
                which = torf.nonzero()[0]

        return torf, False

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

        if arg in XYZ2INDEX:
            coords = self._getCoords() # how about atoms._getCoords() ?
            if coords is None:
                return None, SelectionError(sel, loc,
                    'coordinates are not set')
            else:
                if copy:
                    return coords[:, XYZ2INDEX[arg]].copy(), False
                else:
                    return coords[:, XYZ2INDEX[arg]], False
        arg = FIELDS_SYNONYMS.get(arg, arg)
        try:
            if copy:
                data = self._atoms.getData(arg)
            else:
                data = self._atoms._getData(arg)
        except Exception as err:
            return None, SelectionError(sel, loc, 'following exception '
                        'occurred when evaluating {0}: {1}'
                        .format(repr(arg), str(err)))

        if data is not None:
            if data.dtype.char in 'US':
                return None, SelectionError(sel, loc, '{0} is not a numeric '
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
            return None, SelectionError(sel, loc, '{0} is not a number or a '
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
        left, err = self._getNumeric(sel, loc, tokens.pop(0))
        if err: raise err

        torf = None
        while tokens:
            try:
                binop = OPERATORS[tokens.pop(0)]
            except KeyError:
                raise SelectionError(sel, loc, 'invalid operator encountered')

            right, err = self._getNumeric(sel, loc, tokens.pop(0))
            if err: raise err

            if torf is None:
                torf = binop(left, right)
            else:
                logical_and(binop(left, right), torf, torf)
            left = right
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
            if binop not in OPERATORS:
                raise SelectionError(sel, loc, 'invalid operator encountered')

            right, err = self._getNumeric(sel, loc, tokens.pop(0))
            if err: raise err

            if DEBUG: print(binop, left, right)
            if binop == '/' and any(right == 0.0):
                raise SelectionError(sel, loc, 'zero division error')
            binop = OPERATORS[binop]

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
            raise SelectionError(sel, loc, '{0} accepts a single numeric '
                             'argument, e.g. {0}(x)'.format(repr(tokens[0])))
        arg, err = self._getNumeric(sel, loc, tokens[1])
        if err: raise err
        debug(sel, loc, tokens[0], arg)
        return FUNCTIONS[tokens[0]](arg)

    def _generic(self, sel, loc, tokens, subset=None):

        debug(sel, loc, '_generic', tokens)

        label = tokens.pop(0)
        data, err = self._getData(sel, loc, label)
        if err: return None, err

        if subset is not None:
            data = data[subset]
            subset = None
        dtype = data.dtype
        type_ = dtype.type
        isstr = dtype.char == 'S' or dtype.char == 'U'
        if isstr:
            maxlen = int(dtype.str[2:])
        torf = None
        regexp = []
        values = []
        ranges = []

        valset = True
        for token in tokens:

            # check for regular expressions which are only compatible with
            # string data type
            try:
                token.pattern
            except AttributeError:
                pass
            else:
                if not isstr:
                    ptrn = '"{0}"'.format(token.pattern)
                    SelectionWarning(sel, loc, '{0} is a regular '
                        'expression and is not evaluated for {1}'
                        .format(ptrn, repr(label)), [label, ptrn])
                else:
                    regexp.append(token)
                continue


            # check for ranges
            if token[0] == 'range':
                if isstr:
                    SelectionWarning(sel, loc, 'number ranges '
                        'are not evaluated with data type of {0}'
                        .format(repr(label)), [label, token[1]])
                else:
                    if token[-1] in OPERATORS:
                        ranges.append(token)
                    else:
                        nrange = arange(*token[1:])
                        # if dtypes are not the same, don't use set method
                        if nrange.dtype != dtype: valset = False
                        values.extend(nrange)
                continue

            if isstr:
                if token == '_':
                    values.append('')
                    values.append(' ')
                else:
                    if len(token) > maxlen:
                        SelectionWarning(sel, loc, '{0} is longer than the '
                            'maximum characters for data field {1}'
                            .format(repr(token), repr(label)), [label, token])
                    values.append(token)
            else:
                try:
                    value = type_(token)
                except Exception as err:
                    SelectionWarning(sel, loc, '{0} could not be '
                        'converted to type of {1} ({2})'
                        .format(repr(token), repr(label), str(err)),
                        [label, token])
                    continue
                try:
                    val2 = float(token)
                except:
                    pass
                else:
                    if val2 != value:
                        SelectionWarning(sel, loc, '{0} has a different '
                            'values when converted to a float and to type of '
                            '{1}'.format(repr(token), repr(label)),
                            [label, token])
                values.append(value)

        if values:
            # use first option only if values and data array has the same dtype
            if valset and len(values) > 10:
                valset = set(values)
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
                    torf = logical_and(torf, OPERATORS[comp](data, stop), torf)
                else:
                    torf[subset] = logical_and(torf[subset],
                                       OPERATORS[comp](subdata, stop))

        if regexp:
            for re in regexp:
                if torf is None:
                    torf = array([re.match(val) is not None
                                  for val in data], bool)
                else:
                    subset = (torf == False).nonzero()[0]
                    if len(subset) == 0: return torf, False
                    torf[subset] = [re.match(val) is not None
                                    for val in data[subset]]
        if torf is None:
            torf = self._getZeros(subset)
        return torf, False

    def _index(self, sel, loc, tokens, subset=None):

        debug(sel, loc, '_index', tokens)
        label = tokens.pop(0)
        torf = zeros(self._ag.numAtoms(), bool)

        for token in tokens:
            try:
                remainder = token % 1.
            except TypeError:
                pass
            else:
                return None, SelectionError(sel, loc, 'it is a number, index')
                if remainder == 0:
                    try:
                        torf[token] = True
                    except IndexError:
                        pass
                else:
                    SelectionWarning(sel, loc, '{0} must be followed by '
                        'integers and/or number ranges'.format(repr(label)),
                        [label, token])
                continue

            try:
                token.pattern
            except AttributeError:
                pass
            else:
                ptrn = '"{0}"'.format(token.pattern)
                SelectionWarning(sel, loc, '{0} is a regular '
                    'expression and is not evaluated for {1}'
                    .format(ptrn, repr(label)), [label, ptrn])
                continue

            if token[0] == 'range':
                _, start, stop, step = token
                if step in OPERATORS:
                    if step == '<=':
                        stop += 1
                    step = 1
                if start % 1.0 != 0 or stop % 1.0 != 0 or step % 1.0 != 0:
                    SelectionWarning(sel, loc, '{0} number ranges should be '
                        'specified by integers'.format(repr(label)),
                        [label, str(start)])
                torf[start:stop:step] = True
                continue

            try:
                val = float(token)
            except TypeError:
                return None, SelectionError(sel, loc, '{0} must be '
                    'followed by integers and/or number ranges',
                    [label, token])
            else:
                if val % 1.0 == 0:
                    try:
                        torf[int(val)] = True
                    except IndexError:
                        pass
                else:
                    SelectionWarning(sel, loc, '{0} must be followed by '
                        'integers and/or number ranges'.format(repr(label)),
                        [label, token])

        try:
            indices = self._atoms._getIndices()
        except AttributeError:
            pass
        else:
            torf = torf[indices]

        if subset is not None:
            torf = torf[subset]
        return torf, False

    def _serial(self, sel, loc, tokens, subset=None):

        debug(sel, loc, '_serial', tokens)
        label = tokens.pop(0)
        sn2i = self._ag._getSN2I()
        if sn2i is None:
            return None, SelectionError(sel, loc, 'serial numbers are not set'
                                        ['serial'])
        torf = zeros(len(sn2i), bool)
        for token in tokens:
            try:
                remainder = token % 1.
            except TypeError:
                pass
            else:
                return None, SelectionError(sel, loc, '??? it is a number, serial')
                if remainder == 0:
                    try:
                        torf[token] = True
                    except IndexError:
                        pass
                else:
                    SelectionWarning(sel, loc, '{0} must be followed by '
                        'integers and/or number ranges'.format(repr(label)),
                        [label, token])
                continue

            try:
                pattern = token.pattern
            except AttributeError:
                pass
            else:
                ptrn = '"{0}"'.format(pattern)
                SelectionWarning(sel, loc, '{0} is a regular '
                    'expression and is not evaluated for {1}'
                    .format(ptrn, repr(label)), [label, ptrn])
                continue

            if token[0] == 'range':
                _, start, stop, step = token
                if step in OPERATORS:
                    if step == '<=':
                        stop += 1
                    step = 1
                if start % 1.0 != 0 or stop % 1.0 != 0 or step % 1.0 != 0:
                    SelectionWarning(sel, loc, '{0} number ranges should be '
                        'specified by integers'.format(repr(label)),
                        [label, str(start)])
                torf[start:stop:step] = True
                continue

            try:
                val = float(token)
            except TypeError:
                return None, SelectionError(sel, loc, '{0} must be '
                    'followed by integers and/or number ranges',
                    [label, token])
            else:
                if val % 1.0 == 0:
                    try:
                        torf[int(val)] = True
                    except IndexError:
                        pass
                else:
                    SelectionWarning(sel, loc, '{0} must be followed by '
                        'integers and/or number ranges'.format(repr(label)),
                        [label, token])

        indices = sn2i[torf]
        indices = indices[where(indices != -1)[0]]
        if len(indices) == 0:
            return self._getZeros(subset), False
        torf = zeros(self._ag.numAtoms(), bool)
        torf[indices] = True
        try:
            indices = self._atoms._getIndices()
        except AttributeError:
            pass
        else:
            torf = torf[indices]

        if subset is not None:
            torf = torf[subset]
        return torf, False

    def _resnum(self, sel, loc, tokens, subset=None):

        debug(sel, loc, '_resnum', tokens)
        label = tokens.pop(0)

        resnums, err = self._getData(sel, loc, 'resnum')
        if err: return None, err
        wicode = set([])
        values = ['resnum']
        for token in tokens:

            try:
                pattern = token.pattern
            except AttributeError:
                pass
            else:
                ptrn = '"{0}"'.format(pattern)
                SelectionWarning(sel, loc, '{0} is a regular expression and '
                                 'its use with {1} is not recommended'
                                 .format(ptrn, repr(label)), [label, ptrn])
                continue

            if token[0] == 'range':
                values.append(token)
                continue

            try:
                value = float(token)
            except (TypeError, ValueError):
                icode = token[-1]
                value = token[:1]
                try:
                    value = int(value)
                except:
                    SelectionWarning(sel, loc, '{0} must be followed by '
                        'integers, number ranges, or integer and insertion '
                        'code combinations, e.g. {1}'
                        .format(repr(label), repr('10A')),
                            [label, token])
                else:
                    wicode.add((value, '' if icode == '_' else icode))
            else:
                if value != int(token):
                    SelectionWarning(sel, loc, '{0} must be followed by '
                        'integers and/or number ranges'.format(repr(label)),
                        [label, token])
                values.append(token)
        torf = None

        if len(values) > 1:
            torf, _ = self._generic(sel, loc, values, subset)

        if wicode:
            icode, err = self._getData(sel, loc, 'icode')
            if err: return None, err
            if subset is None:
                rnic = zip(resnums, icode) # PY3K: OK
            else:
                rnic = zip(resnums, icode[subset]) # PY3K: OK

            if torf is None:
                torf = array([val in wicode for val in rnic], bool)
            else:
                torf = logical_or(torf, array([val in wicode for val in rnic],
                                               bool), torf)
        return torf, False

    def _sequence(self, sel, loc, tokens, subset=None):

        debug(sel, loc, '_sequence', tokens)
        label = tokens.pop(0)

        regexp = []
        for token in tokens:

            try:
                token.pattern
            except AttributeError:
                if not token.isalpha() or not token.isupper():
                    SelectionWarning(sel, loc, '{0} does not look like a '
                        'valid sequence'.format(repr(token)),
                        [label, token])
                try:
                    token = re_compile(token)
                except Exception as err:
                    SelectionWarning(sel, loc, '{0} could not be compiled '
                        'as a regular expression for sequence evaluation'
                        .format(repr(token)), [label, token])
                else:
                    regexp.append(token)
            else:
                regexp.append(token)

        if not regexp:
            return self._getZeros(subset), False

        calpha = self._atoms.calpha
        if calpha is None:
            return self._getZeros(subset), False

        matches = []
        for chain in iter(HierView(calpha)):
            sequence = chain.getSequence()
            indices = chain._getIndices()
            for re in regexp:
                for match in re.finditer(sequence):
                    matches.extend(indices[match.start():match.end()])

        if matches:
            torf = zeros(self._ag.numAtoms(), bool)
            torf[matches] = True
            if self._indices is not None:
                torf = torf[self._indices]
            torf = self._sameas(sel, loc, [('same', 'residue', 'as'), torf])[0]
            if subset is None:
                return torf, False
            else:
                return torf[subset], False
        else:
            return self._getZeros(subset), False

    def _getData(self, sel, loc, keyword):
        """Return atomic data."""

        data = self._data.get(keyword)
        if data is not None:
            return data, False

        try:
            idx = XYZ2INDEX[keyword]
        except KeyError:
            pass
        else:
            data = self._getCoords()
            if data is not None:
                data = data[:,idx]
                self._data[keyword] = data
            return data, False

        if keyword == 'index':
            try:
                data = self._atoms._getIndices()
            except AttributeError:
                data = arange(self._atoms.numAtoms())
            self._data['index'] = data
            return data, False

        field = ATOMIC_FIELDS.get(FIELDS_SYNONYMS.get(keyword, keyword))
        if field is None:
            data = self._atoms._getData(keyword)
            if data is None:
                return None, SelectionError(sel, loc, '{0} is not a valid '
                                       'data label'.format(repr(keyword)))
            elif not isinstance(data, ndarray) and data.ndim == 1:
                return None, SelectionError(sel, loc, '{0} is not a 1d '
                                      'array'.format(repr(keyword)))
        else:
            try:
                data = getattr(self._atoms, '_get' + field.meth_pl)()
            except Exception as err:
                return None, SelectionError(sel, loc, str(err))
            if data is None:
                return None, SelectionError(sel, loc, '{0} is not set'
                                                       .format(repr(keyword)))
        self._data[keyword] = data
        return data, False

    def _getCoords(self):
        """Return coordinates of atoms."""

        if self._coords is None:
            self._coords = self._atoms._getCoords()
        return self._coords
