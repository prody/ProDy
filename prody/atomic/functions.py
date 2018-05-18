# -*- coding: utf-8 -*-
"""This module defines some functions for handling atomic classes and data."""

from textwrap import wrap

from numpy import load, savez, ones, zeros, array, argmin, where
from numpy import ndarray, asarray, isscalar, concatenate, arange, ix_

from prody.utilities import openFile, rangeString, getDistance
from prody import LOGGER

from . import flags
from . import select

from .atomic import Atomic
from .atomgroup import AtomGroup
from .atommap import AtomMap
from .bond import trimBonds, evalBonds
from .fields import ATOMIC_FIELDS
from .selection import Selection
from .hierview import HierView

__all__ = ['iterFragments', 'findFragments', 'loadAtoms', 'saveAtoms',
           'isReserved', 'listReservedWords', 'sortAtoms', 'sliceAtoms', 
           'extendAtoms', 'sliceAtomicData', 'extendAtomicData']


SAVE_SKIP_ATOMGROUP = set(['numbonds', 'fragindex'])
SAVE_SKIP_POINTER = set(['numbonds', 'fragindex', 'segindex', 'chindex',
                         'resindex'])


def saveAtoms(atoms, filename=None, **kwargs):
    """Save *atoms* in ProDy internal format.  All :class:`.Atomic` classes are
    accepted as *atoms* argument.  This function saves user set atomic data as
    well.  Note that title of the :class:`.AtomGroup` instance is used as the
    filename when *atoms* is not an :class:`.AtomGroup`.  To avoid overwriting
    an existing file with the same name, specify a *filename*."""

    try:
        atoms.getACSIndex()
    except AttributeError:
        raise TypeError('atoms must be Atomic instance, not {0}'
                        .format(type(atoms)))

    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        ag = atoms
        title = ag.getTitle()
        SKIP = SAVE_SKIP_ATOMGROUP
    else:
        SKIP = SAVE_SKIP_POINTER
        title = str(atoms)

    if filename is None:
        filename = ag.getTitle().replace(' ', '_')
    if '.ag.npz' not in filename:
        filename += '.ag.npz'

    attr_dict = {'title': title}
    attr_dict['n_atoms'] = atoms.numAtoms()
    attr_dict['n_csets'] = atoms.numCoordsets()
    attr_dict['cslabels'] = atoms.getCSLabels()
    attr_dict['flagsts'] = ag._flagsts
    coords = atoms._getCoordsets()
    if coords is not None:
        attr_dict['coordinates'] = coords
    bonds = ag._bonds
    bmap = ag._bmap
    if bonds is not None and bmap is not None:
        if atoms == ag:
            attr_dict['bonds'] = bonds
            attr_dict['bmap'] = bmap
            attr_dict['numbonds'] = ag._data['numbonds']
            frags = ag._data.get('fragindex')
            if frags is not None:
                attr_dict['fragindex'] = frags
        else:
            bonds = trimBonds(bonds, atoms._getIndices())
            attr_dict['bonds'] = bonds
            attr_dict['bmap'], attr_dict['numbonds'] = \
                evalBonds(bonds, len(atoms))

    for label in atoms.getDataLabels():
        if label in SKIP:
            continue
        attr_dict[label] = atoms._getData(label)
    for label in atoms.getFlagLabels():
        if label in SKIP:
            continue
        attr_dict[label] = atoms._getFlags(label)

    ostream = openFile(filename, 'wb', **kwargs)
    savez(ostream, **attr_dict)
    ostream.close()
    return filename


SKIPLOAD = set(['title', 'n_atoms', 'n_csets', 'bonds', 'bmap',
                'coordinates', 'cslabels', 'numbonds', 'flagsts',
                'segindex', 'chindex', 'resindex'])


def loadAtoms(filename):
    """Returns :class:`.AtomGroup` instance loaded from *filename* using
    :func:`numpy.load` function.  See also :func:`saveAtoms`."""

    LOGGER.timeit('_prody_loadatoms')
    attr_dict = load(filename)
    files = set(attr_dict.files)

    if not 'n_atoms' in files:
        raise ValueError('{0} is not a valid atomic data file'
                         .format(repr(filename)))
    title = str(attr_dict['title'])

    if 'coordinates' in files:
        coords = attr_dict['coordinates']
        ag = AtomGroup(title)
        ag._n_csets = int(attr_dict['n_csets'])
        ag._coords = coords
    ag._n_atoms = int(attr_dict['n_atoms'])
    ag._setTimeStamp()
    if 'flagsts' in files:
        ag._flagsts = int(attr_dict['flagsts'])

    if 'bonds' in files and 'bmap' in files and 'numbonds' in files:
        ag._bonds = attr_dict['bonds']
        ag._bmap = attr_dict['bmap']
        ag._data['numbonds'] = attr_dict['numbonds']

    skip_flags = set()

    for label, data in attr_dict.items():
        if label in SKIPLOAD:
            continue
        if data.ndim == 1 and data.dtype == bool:
            if label in skip_flags:
                continue
            else:
                ag._setFlags(label, data)
                skip_flags.update(flags.ALIASES.get(label, [label]))
        else:
            ag.setData(label, data)

    for label in ['segindex', 'chindex', 'resindex']:
        if label in attr_dict:
            ag._data[label] = attr_dict[label]

    if ag.numCoordsets() > 0:
        ag._acsi = 0

    if 'cslabels' in files:
        ag.setCSLabels(list(attr_dict['cslabels']))

    LOGGER.report('Atom group was loaded in %.2fs.', '_prody_loadatoms')
    return ag


def iterFragments(atoms):
    """Yield fragments, connected subsets in *atoms*, as :class:`.Selection`
    instances."""

    try:
        return atoms.iterFragments()
    except AttributeError:
        pass

    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        raise TypeError('atoms must be an Atomic instance')

    bonds = atoms._iterBonds()
    return _iterFragments(atoms, ag, bonds)


def _iterFragments(atoms, ag, bonds):

    fids = zeros((len(ag)), int)
    fdict = {}
    c = 0
    for a, b in bonds:
        af = fids[a]
        bf = fids[b]
        if af and bf:
            if af != bf:
                frag = fdict[af]
                temp = fdict[bf]
                fids[temp] = af
                frag.extend(temp)
                fdict.pop(bf)
        elif af:
            fdict[af].append(b)
            fids[b] = af
        elif bf:
            fdict[bf].append(a)
            fids[a] = bf
        else:
            c += 1
            fdict[c] = [a, b]
            fids[a] = fids[b] = c
    fragments = []
    append = fragments.append
    fidset = set()
    indices = atoms._getIndices()
    for i, fid in zip(indices, fids[indices]):
        if fid in fidset:
            continue
        elif fid:
            fidset.add(fid)
            indices = fdict[fid]
            indices.sort()
            append(indices)
        else:
            # these are non-bonded atoms, e.g. ions
            append([i])

    acsi = atoms.getACSIndex()
    for indices in fragments:
        yield Selection(ag, indices, 'index ' + rangeString(indices), acsi,
                        unique=True)


def findFragments(atoms):
    """Returns list of fragments, connected subsets in *atoms*.  See also
    :func:`iterFragments`."""

    return list(iterFragments(atoms))


RESERVED = set(ATOMIC_FIELDS)
RESERVED.update(['and', 'or', 'not', 'within', 'of', 'exwithin', 'same', 'as',
                 'bonded', 'exbonded', 'to', 'all', 'none',
                 'index', 'sequence', 'x', 'y', 'z'])
RESERVED.update(flags.PLANTERS)
RESERVED.update(select.FUNCTIONS)
RESERVED.update(select.FIELDS_SYNONYMS)
RESERVED.update(['n_atoms', 'n_csets', 'cslabels', 'title', 'coordinates',
                 'bonds', 'bmap'])


def isReserved(word):
    """Returns **True** if *word* is reserved for internal data labeling or atom
    selections.  See :func:`listReservedWords` for a list of reserved words."""

    return word in RESERVED


def listReservedWords():
    """Returns list of words that are reserved for atom selections and internal
    variables. These words are: """

    words = list(RESERVED)
    words.sort()
    return words

_ = listReservedWords.__doc__ + '*' + '*, *'.join(listReservedWords()) + '*.'
listReservedWords.__doc__ = '\n'.join(wrap(_, 79))


def sortAtoms(atoms, label, reverse=False):
    """Returns an :class:`.AtomMap` pointing to *atoms* sorted in ascending
    data *label* order, or optionally in *reverse* order."""

    try:
        data, acsi = atoms.getData(label), atoms.getACSIndex()
    except AttributeError:
        raise TypeError('atoms must be an Atomic instance')
    else:
        if data is None:
            raise ValueError('{0} data is not set for {1}'
                             .format(repr(label), atoms))
    sort = data.argsort()
    if reverse:
        sort = sort[::-1]
    try:
        indices = atoms.getIndices()
    except AttributeError:
        ag = atoms
    else:
        ag = atoms.getAtomGroup()
        sort = indices[sort]
    return AtomMap(ag, sort, acsi)


def sliceAtoms(atoms, select):
    """Slice *atoms* using the selection defined by *select*.

    :arg atoms: atoms to be selected from
    :type atoms: :class:`Atomic`

    :arg select: a :class:`Selection` instance or selection string
    :type select: :class:`Selection`, str

    """

    if atoms == select:
        raise ValueError('atoms and select arguments are the same')
    if select in atoms:
        indices = select._getIndices()
    elif isinstance(select, str):
        select = atoms.select(select)
        indices = select._getIndices()
    else:
        raise TypeError('select must be a string or a Selection instance')

    if isinstance(atoms, AtomGroup):
        which = indices
    else:
        idxset = set(indices)
        which = array([i for i, idx in enumerate(atoms._getIndices())
                          if idx in idxset])

    return which, select

def extendAtoms(nodes, atoms, is3d=False):
    """Returns extended mapping indices and an :class:`.AtomMap`."""

    try:
        i_nodes = nodes.iterAtoms()
    except AttributeError:
        raise ValueError('nodes must be an Atomic instance')

    if not nodes in atoms:
        raise ValueError('nodes must be a subset of atoms')

    atom_indices = []
    real_indices = []   # indices of atoms that are used as nodes (with real mode data)
    indices = []
    get = HierView(atoms).getResidue
    residues = []

    for i, node in enumerate(i_nodes):
        res = get(node.getChid() or None, node.getResnum(),
                  node.getIcode() or None, node.getSegname() or None)
        if res is None:
            raise ValueError('atoms must contain a residue for all atoms')

        res_atom_indices = res._getIndices()
        if res not in residues:
            atom_indices.append(res_atom_indices)

            res_real_indices = ones(len(res_atom_indices)) * -1
            real_indices.append(res_real_indices)

            if is3d:
                nma_indices = list(range(i*3, (i+1)*3)) * len(res)
            else:
                nma_indices = [i] * len(res)

            indices.append(nma_indices)
            residues.append(res)
        else:
            k = where(res_atom_indices==node.getIndex())[0][0]
            if is3d:
                nma_indices[k*3:(k+1)*3] = list(range(i*3, (i+1)*3))
            else:
                nma_indices[k] = i

            res_real_indices = real_indices[residues.index(res)]
        
        # register the real node
        node_index = node.getIndex()
        res_real_indices[res_atom_indices == node_index] = i
    
    def getClosest(a, B):
        D = []
        for b in B:
            d = getDistance(a._getCoords(), b._getCoords())
            D.append(d)
        
        i = argmin(D)
        return B[i]

    for i, res_real_indices in enumerate(real_indices):
        arr = array(res_real_indices)
        # this residue is represented by one node, so no correction is needed
        if sum(arr >= 0) == 1: 
            continue
        # otherwise replace the data of extended atoms by that of 
        # the nearest real node in the residue
        else:
            # get all the atoms in this residue
            res_atoms = array(residues[i])
            # get the real and extended atoms
            real_atoms = res_atoms[arr >= 0]
            for j in range(len(res_real_indices)):
                if res_real_indices[j] >= 0:
                    continue
                else:
                    atom = res_atoms[j]
                    closest_real_atom = getClosest(atom, real_atoms)
                    k = where(real_atoms == closest_real_atom)[0][0]
                    
                    nma_indices = indices[i]
                    if is3d:
                        nma_indices[j*3:(j+1)*3] = nma_indices[k*3:(k+1)*3]
                    else:
                        nma_indices[j] = nma_indices[k]

    atom_indices = concatenate(atom_indices)
    indices = concatenate(indices)

    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        ag = atoms
    atommap = AtomMap(ag, atom_indices, atoms.getACSIndex(),
                      title=str(atoms), intarrays=True)
    return indices, atommap

def sliceAtomicData(data, atoms, select, axis=None):
    """Slice a matrix using indices extracted using :func:`sliceAtoms`.

    :arg data: any data array
    :type data: `~numpy.ndarray`

    :arg atoms: atoms to be selected from
    :type atoms: :class:`Atomic`

    :arg select: a :class:`Selection` instance or selection string
    :type select: :class:`Selection`, str

    :arg axis: the axis along which the data is sliced. See :mod:`~numpy` 
               for details of this parameter. 
               Default is **None** (all axes)
    :type axis: int, list

    """

    if isscalar(data):
        raise TypeError('The data must be array-like.')

    if not isinstance(data, ndarray):
        data = asarray(data)

    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance')

    natoms = atoms.numAtoms()

    is3d = False
    if len(data) != natoms:
        if data.shape[0] == natoms * 3:
            is3d = True
        else:
            raise ValueError('data and atoms must have the same size')

    indices, _ = sliceAtoms(atoms, select)
    if is3d:
        indices = array([[i*3, i*3+1, i*3+2] 
                        for i in indices]
                        ).reshape(3*len(indices))

    if axis is not None:
        I = [arange(s) for s in data.shape] 
        axes = [axis] if isscalar(axis) else axis
        for ax in axes:
            I[ax] = indices
    else:
        I = [indices] * data.ndim
    
    profiles = data[ix_(*I)]
        
    return profiles

def extendAtomicData(data, nodes, atoms, axis=None):
    """Extend a coarse grained data obtained for *nodes* to *atoms*.

    :arg data: any data array
    :type data: `~numpy.ndarray`

    :arg nodes: a set of atoms that has been used
        as nodes in data generation
    :type nodes: :class:`

    :arg atoms: atoms to be selected from
    :type atoms: :class:`Atomic`

    :arg axis: the axis/direction you want to use to slice data from the matrix.
        The options are **0** or **1** or **None** like in :mod:`~numpy`. 
        Default is **None** (all axes)
    :type axis: int

    """
    
    try:
        data = asarray(data)
    except:
        raise TypeError('The data must be array-like.')

    nnodes = nodes.numAtoms()

    is3d = False
    if data.shape[0] != nnodes:
        if data.shape[0] == nnodes * 3:
            is3d = True
        else:
            raise ValueError('data and nodes must have the same size')

    indices, atommap = extendAtoms(nodes, atoms, is3d)
    
    if axis is not None:
        I = [arange(s) for s in data.shape] 
        axes = [axis] if isscalar(axis) else axis
        for ax in axes:
            I[ax] = indices
    else:
        I = [indices] * data.ndim

    data_ext = data[ix_(*I)]
        
    return data_ext, atommap
