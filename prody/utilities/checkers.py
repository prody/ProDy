"""This module defines functions for type, value, and/or attribute checking."""

from numpy import any, float32, tile

__all__ = ['checkCoords', 'checkWeights', 'checkTypes']

COORDS_NDIM = set([2])
CSETS_NDIMS = set([2, 3])

def checkCoords(coords, csets=False, natoms=None, dtype=(float, float32),
                name='coords'):
    """Returns **True** if shape, dimensionality, and data type of *coords*
    array are as expected.

    :arg coords: coordinate array

    :arg csets: whether multiple coordinate sets (i.e. ``.ndim in (2, 3)``) are
        allowed, default is **False**

    :arg natoms: number of atoms, if **None** number of atoms is not checked

    :arg dtype: allowed data type(s), default is ``(float, numpy.float32)``,
        if **None** data type is not checked

    :arg name: name of the coordinate argument

    :raises: :exc:`TypeError` when *coords* is not an instance of
        :class:`numpy.ndarray`

    :raises: :exc:`ValueError` when wrong shape, dimensionality, or data type
        is encountered"""

    try:
        ndim, shape = coords.ndim, coords.shape
    except AttributeError:
        raise TypeError('coords must be a numpy.ndarray instance')

    ndims = CSETS_NDIMS if csets else COORDS_NDIM
    if ndim not in ndims:
        raise ValueError(str(name) + '.ndim must be ' +
                         ' or '.join([str(d) for d in ndims]))

    elif shape[-1] != 3:
        raise ValueError(str(name) + '.shape[-1] must be 3')

    elif natoms and shape[-2] != natoms:
        raise ValueError(str(name) + '.shape[-2] must match number of atoms')

    elif dtype:
        if isinstance(dtype, type) and coords.dtype != dtype:
            raise ValueError(str(name) + '.dtype must be ' + dtype.__name__)
        elif coords.dtype not in dtype:
            if len(dtype) > 1:
                msg = ', '.join([repr(dt.__name__) for dt in dtype[:-1]]
                                ) + ', or ' + repr(dtype[-1].__name__)
            else:
                msg = dtype[0].__name__
            raise ValueError(str(name) + '.dtype must be ' + msg)

    return True

NDIM12 = set([1, 2])
NDIM123 = set([1, 2, 3])

def checkWeights(weights, natoms, ncsets=None, dtype=float):
    """Returns *weights* if it has correct shape ([ncsets, ]natoms, 1).
    after its shape and data type is corrected. otherwise raise an exception.
    All items of *weights* must be greater than zero."""

    try:
        ndim, shape, wtype = weights.ndim, weights.shape, weights.dtype
    except AttributeError:
        raise TypeError('weights must be a numpy.ndarray instance')

    if csets:
        if ndim not in NDIM123:
            raise ValueError('weights.dim must be 1, 2, or 3')
        if csets > 1:
            if ndim == 3 and shape[0] != csets:
                raise ValueError('weights.shape must be '
                                   '(ncsets, natoms[, 1])')
            weights = tile(weights.reshape((1, natoms, 1)), (ncsets, 1, 1))
        elif ndim < 3:
            weights = weights.reshape((1, natoms, 1))
    else:
        if ndim not in NDIM12:
            raise ValueError('weights.dim must be 1 or 2')
        if shape[0] != natoms:
            raise ValueError('weights.shape must be (natoms[, 1])')
        if ndim == 1:
            weights = weights.reshape((natoms, 1))

    if dtype:
        if wtype != dtype:
            try:
                weights = weights.astype(dtype)
            except ValueError:
                raise ValueError('weights.astype({0}) failed, {0} type '
                                   'could not be assigned'.format(str(dtype)))
    if any(weights < 0):
        raise ValueError('all weights must be greater or equal to 0')

    return weights


def checkTypes(args, **types):
    """Returns **True** if types of all *args* match those given in *types*.

    :raises: :exc:`TypeError` when type of an argument is not one of allowed
        types

    ::

        def incr(n, i):
            '''Return sum of *n* and *i*.'''

            checkTypes(locals(), n=(float, int), i=(float, int))
            return n + i"""

    for arg, allowed in types.items():
        if arg in args and not isinstance(args[arg], types[arg]):

            val = args[arg]
            if isinstance(allowed, (list, tuple)):
                if len(allowed) > 1:
                    tstr = ', '.join([repr(tp.__name__) for tp in allowed[:-1]]
                                     ) + ', or ' + repr(allowed[-1].__name__)

                else:
                    tstr = repr(allowed[0].__name__)

            else:
                tstr = repr(allowed.__name__)

            raise TypeError('{0} must be an instance of {1}, not {2}'
                            .format(repr(arg), tstr, repr(type(val).__name__)))

    return True
