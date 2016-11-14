"""This module defines miscellaneous utility functions."""

from numpy import unique, linalg, diag, sqrt, dot

__all__ = ['Everything', 'rangeString', 'alnum', 'importLA', 'dictElement',
           'intorfloat', 'startswith', 'showFigure', 'countBytes', 'sqrtm']


class Everything(object):

    """A place for everything."""

    def __contains__(self, what):

        return True


def rangeString(lint, sep=' ', rng=' to ', exc=False, pos=True):
    """Returns a structured string for a given list of integers.

    :arg lint: integer list or array
    :arg sep: range or number separator
    :arg rng: range symbol
    :arg exc: set **True** if range symbol is exclusive
    :arg pos: only consider zero and positive integers

    .. ipython:: python

       from prody.utilities import rangeString
       lint = [1, 2, 3, 4, 10, 15, 16, 17]
       rangeString(lint)
       rangeString(lint, sep=',', rng='-')
       rangeString(lint, ',', ':', exc=True)"""

    ints = unique(lint)
    if pos and ints[0] < 0:
        ints = ints[ints > -1]

    prev = ints[0]
    lint = [[prev]]
    for i in ints[1:]:
        if i - prev > 1:
            lint.append([i])
        else:
            lint[-1].append(i)
        prev = i
    exc = int(exc)
    return sep.join([str(l[0]) if len(l) == 1 else
                     str(l[0]) + rng + str(l[-1] + exc) for l in lint])


def alnum(string, alt='_', trim=False, single=False):
    """Replace non alpha numeric characters with *alt*.  If *trim* is **True**
    remove preceding and trailing *arg* characters.  If *single* is **True**,
    contain only a single joining *alt* character. """

    result = ''
    multi = not bool(single)
    prev = None
    for char in string:
        if char.isalnum():
            result += char
            prev = char
        else:
            if multi or prev != alt:
                result += alt
            prev = alt
    trim = int(bool(trim))
    result = result[trim * (result[0] == alt):
                    len(result) - trim * (result[-1] == alt)]
    return result


def importLA():
    """Returns one of :mod:`scipy.linalg` or :mod:`numpy.linalg`."""

    try:
        import scipy.linalg as linalg
    except ImportError:
        try:
            import numpy.linalg as linalg
        except:
            raise ImportError('scipy.linalg or numpy.linalg is required for '
                              'NMA and structure alignment calculations')
    return linalg


def dictElement(element, prefix=None):
    """Returnss a dictionary built from the children of *element*, which must be
    a :class:`xml.etree.ElementTree.Element` instance.  Keys of the dictionary
    are *tag* of children without the *prefix*, or namespace.  Values depend on
    the content of the child.  If a child does not have any children, its text
    attribute is the value.  If a child has children, then the child is the
    value."""

    dict_ = {}
    length = False
    if isinstance(prefix, str):
        length = len(prefix)
    for child in element:
        tag = child.tag
        if length and tag.startswith(prefix):
            tag = tag[length:]
        if len(child) == 0:
            dict_[tag] = child.text
        else:
            dict_[tag] = child
    return dict_


def intorfloat(x):
    """Returns ``int(x)``, or ``float(x)`` upon :exc:`ValueError`."""

    try:
        return int(x)
    except ValueError:
        return float(x)


def startswith(this, that):
    """Returns **True** if *this* or *that* starts with the other."""

    if len(this) < len(that):
        return that.startswith(this)
    else:
        return this.startswith(that)


def showFigure():
    """Call :func:`~matplotlib.pyplot.show` function with ``block=False``
    argument to avoid blocking behavior in non-interactive sessions.  If
    *block* keyword argument is not recognized, try again without it."""

    from matplotlib.pyplot import show
    try:
        show(block=False)
    except TypeError:
        show()


def countBytes(arrays, base=False):
    """Returns total number of bytes consumed by elements of arrays.  If
    *base* is **True**, use number of bytes from the base array."""

    if base:
        getbase = lambda arr: arr if arr.base is None else getbase(arr.base)
        nbytes = lambda arr: getbase(arr).nbytes
    else:
        nbytes = lambda arr: arr.nbytes

    return sum(nbytes(arr) for arr in arrays)

def sqrtm(matrix):
    """Returns the square root of a matrix."""
    (U,S,VT)=linalg.svd(matrix)
    D = diag(sqrt(S))
    return dot(dot(U,D),VT)

def getMasses(elements):
    """get the mass atom. """
    
    import numpy as np
    mass_dict = {'C':12,'N':14,'S':32,'O':16,'H':1}

    if isinstance(elements, str):
        return mass_dict[elements]
    else:
        masses = np.zeros(len(elements))
        for i,element in enumerate(elements):
            if mass_dict.has_key(element):
                masses[i] = mass_dict[element]
            else:
                masses[i] = 0.
        return masses
