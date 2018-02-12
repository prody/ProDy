"""This module defines a functions for handling conformational ensembles."""

import os.path

import numpy as np

from prody.proteins import fetchPDB, parsePDB, writePDB, mapOntoChain
from prody.utilities import openFile, showFigure
from prody import LOGGER, SETTINGS
from prody.atomic import AtomMap, Chain, AtomGroup, Selection, Segment, Select

from .ensemble import *
from .pdbensemble import *
from .conformation import *

__all__ = ['saveEnsemble', 'loadEnsemble', 'trimPDBEnsemble',
           'calcOccupancies', 'showOccupancies', 'alignPDBEnsemble',
           'calcTree', 'buildPDBEnsemble', 'addPDBEnsemble']


def saveEnsemble(ensemble, filename=None, **kwargs):
    """Save *ensemble* model data as :file:`filename.ens.npz`.  If *filename*
    is ``None``, title of the *ensemble* will be used as the filename, after
    white spaces in the title are replaced with underscores.  Extension is
    :file:`.ens.npz`. Upon successful completion of saving, filename is
    returned. This function makes use of :func:`numpy.savez` function."""

    if not isinstance(ensemble, Ensemble):
        raise TypeError('invalid type for ensemble, {0}'
                        .format(type(ensemble)))
    if len(ensemble) == 0:
        raise ValueError('ensemble instance does not contain data')

    dict_ = ensemble.__dict__
    attr_list = ['_title', '_confs', '_weights', '_coords']
    if isinstance(ensemble, PDBEnsemble):
        attr_list.append('_labels')
        attr_list.append('_trans')
    if filename is None:
        filename = ensemble.getTitle().replace(' ', '_')
    attr_dict = {}
    for attr in attr_list:
        value = dict_[attr]
        if value is not None:
            attr_dict[attr] = value
            
    attr_dict['_atoms'] = np.array([dict_['_atoms'], 0])
    filename += '.ens.npz'
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename


def loadEnsemble(filename):
    """Returns ensemble instance loaded from *filename*.  This function makes
    use of :func:`numpy.load` function.  See also :func:`saveEnsemble`"""

    attr_dict = np.load(filename)
    if '_weights' in attr_dict:
        weights = attr_dict['_weights']
    else:
        weights = None   
    isPDBEnsemble = False
    try:
        title = str(attr_dict['_title'])
    except KeyError:
        title = str(attr_dict['_name'])
    if weights is not None and weights.ndim == 3:
        isPDBEnsemble = True
        ensemble = PDBEnsemble(title)
    else:
        ensemble = Ensemble(title)
    ensemble.setCoords(attr_dict['_coords'])
    if '_atoms' in attr_dict:
        atoms = attr_dict['_atoms'][0]
    else:
        atoms = None
    ensemble.setAtoms(atoms)    
    if isPDBEnsemble:
        confs = attr_dict['_confs']
        ensemble.addCoordset(attr_dict['_confs'], weights)
        if '_identifiers' in attr_dict.files:
            ensemble._labels = list(attr_dict['_identifiers'])
        if '_labels' in attr_dict.files:
            ensemble._labels = list(attr_dict['_labels'])
        if '_trans' in attr_dict.files:
            ensemble._trans = attr_dict['_trans']
    else:
        ensemble.addCoordset(attr_dict['_confs'])
        if weights is not None:
            ensemble.setWeights(weights)
    return ensemble


def trimPDBEnsemble(pdb_ensemble, **kwargs):
    """Returns a new PDB ensemble obtained by trimming given *pdb_ensemble*.
    This function helps selecting atoms in a pdb ensemble based on one of the
    following criteria, and returns them in a new :class:`.PDBEnsemble`
    instance.

    **Occupancy**

    Resulting PDB ensemble will contain atoms whose occupancies are greater
    or equal to *occupancy* keyword argument.  Occupancies for atoms will be
    calculated using ``calcOccupancies(pdb_ensemble, normed=True)``.

    :arg occupancy: occupancy for selecting atoms, must satisfy
        ``0 < occupancy <= 1``
    :type occupancy: float

    :arg selstr: The function will trim residues that are NOT specified by 
        the selection string.
    :type selstr: str

    """

    if not isinstance(pdb_ensemble, PDBEnsemble):
        raise TypeError('pdb_ensemble argument must be a PDBEnsemble')
    if pdb_ensemble.numConfs() == 0 or pdb_ensemble.numAtoms() == 0:
        raise ValueError('pdb_ensemble must have conformations')

    if 'occupancy' in kwargs:
        occupancy = float(kwargs['occupancy'])
        assert 0 < occupancy <= 1, ('occupancy is not > 0 and <= 1: '
                                    '{0}'.format(repr(occupancy)))
        n_confs = pdb_ensemble.numConfs()
        assert n_confs > 0, 'pdb_ensemble does not contain any conformations'
        occupancies = calcOccupancies(pdb_ensemble, normed=True)
        #assert weights is not None, 'weights must be set for pdb_ensemble'
        #weights = weights.flatten()
        #mean_weights = weights / n_confs
        torf = occupancies >= occupancy
    elif 'selstr' in kwargs:
        selstr = kwargs['selstr']
        atoms = pdb_ensemble.getAtoms()
        assert atoms is not None, 'atoms are empty'
        selector = Select()
        torf = selector.getBoolArray(atoms, selstr)
    else:
        n_atoms = pdb_ensemble.getCoords().shape[0]
        torf = np.ones(n_atoms, dtype=bool)

    trimmed = PDBEnsemble(pdb_ensemble.getTitle())

    atoms = pdb_ensemble.getAtoms()
    if atoms is not None:
        trim_atoms_idx = [n for n,t in enumerate(torf) if t]
        if type(atoms) is Chain:
            trim_atoms=Chain(atoms.copy(), trim_atoms_idx, atoms._hv)
        elif type(atoms) is AtomGroup:
            trim_atoms= AtomMap(atoms,trim_atoms_idx)
        else:
            trim_atoms= AtomMap(atoms.copy(),trim_atoms_idx)
        trimmed.setAtoms(trim_atoms)

    coords = pdb_ensemble.getCoords()
    if coords is not None:
        trimmed.setCoords(coords[torf])
    confs = pdb_ensemble.getCoordsets()
    if confs is not None:
        weights = pdb_ensemble.getWeights()
        labels = pdb_ensemble.getLabels()
        trimmed.addCoordset(confs[:, torf], weights[:, torf], labels)
    return trimmed


def calcOccupancies(pdb_ensemble, normed=False):
    """Returns occupancy calculated from weights of a :class:`.PDBEnsemble`.
    Any non-zero weight will be considered equal to one.  Occupancies are
    calculated by binary weights for each atom over the conformations in
    the ensemble. When *normed* is ``True``, total weights will be divided
    by the number of atoms.  This function can be used to see how many times
    a residue is resolved when analyzing an ensemble of X-ray structures."""

    if not isinstance(pdb_ensemble, PDBEnsemble):
        raise TypeError('pdb_ensemble must be a PDBEnsemble instance')
    if len(pdb_ensemble) == 0:
        raise ValueError('pdb_ensemble does not contain any conformations')
    assert isinstance(normed, bool), 'normed must be a boolean'
    weights = pdb_ensemble.getWeights()
    if weights is None:
        raise ValueError('pdb_ensemble weights are not set')

    occupancies = weights.astype(bool).sum(0).astype(float).flatten()
    if normed:
        return occupancies / len(pdb_ensemble)
    else:
        return occupancies


def showOccupancies(pdbensemble, *args, **kwargs):
    """Show occupancies for the PDB ensemble using :func:`~matplotlib.pyplot.
    plot`.  Occupancies are calculated using :meth:`calcOccupancies`."""

    import matplotlib.pyplot as plt

    if not isinstance(pdbensemble, PDBEnsemble):
        raise TypeError('pdbensemble must be a PDBEnsemble instance')
    weights = calcOccupancies(pdbensemble)
    if weights is None:
        return None
    show = plt.plot(weights, *args, **kwargs)
    axis = list(plt.axis())
    axis[2] = 0
    axis[3] += 1
    plt.axis(axis)
    plt.xlabel('Atom index')
    plt.ylabel('Sum of weights')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def checkWeights(weights, n_atoms, n_csets=None):
    """Returns weights if checks pass, otherwise raise an exception."""

    assert isinstance(n_atoms, int) or isinstance(n_atoms, long) and n_atoms > 0, \
        'n_atoms must be a positive integer'
    assert n_csets is None or isinstance(n_csets, int) or isinstance(n_atoms, long) and n_csets > 0, \
        'n_csets must be a positive integer'

    if not isinstance(weights, np.ndarray):
        raise TypeError('weights must be a Numpy array')
    elif n_csets is None and weights.ndim not in (1, 2):
        raise ValueError('weights.dim must be 1 or 2')
    elif n_csets is not None:
        if weights.ndim not in (1, 2, 3):
            raise ValueError('weights.dim must be 1, 2, or 3')
        elif weights.ndim == 3 and weights.shape[0] != n_csets:
            raise ValueError('weights.shape must be (n_csets, n_atoms, 1)')
    elif weights.ndim in (2, 3) and weights.shape[-1] != 1:
        raise ValueError('shape of weights must be ([n_csets,] n_atoms, 1)')
    elif weights.dtype != float:
        try:
            weights = weights.astype(float)
        except ValueError:
            raise ValueError('weights.astype(float) failed, float type could '
                             'not be assigned')
    if np.any(weights < 0):
        raise ValueError('all weights must greater or equal to 0')
    if weights.ndim == 1:
        weights = weights.reshape((n_atoms, 1))
    if n_csets is not None and weights.ndim == 2:
        weights = np.tile(weights.reshape((1, n_atoms, 1)), (n_csets, 1, 1))
    return weights


def alignPDBEnsemble(ensemble, suffix='_aligned', outdir='.', gzip=False):
    """Align PDB files using transformations from *ensemble*, which may be
    a :class:`.PDBEnsemble` or a :class:`.PDBConformation` instance. Label of
    the conformation (see :meth:`~.PDBConformation.getLabel`) will be used to
    determine the PDB structure and model number.  First four characters of
    the label is expected to be the PDB identifier and ending numbers to be the
    model number.  For example, the :class:`.Transformation` from conformation
    with label *2k39_ca_selection_'resnum_<_71'_m116* will be applied to 116th
    model of structure **2k39**.  After applicable transformations are made,
    structure will be written into *outputdir* as :file:`2k39_aligned.pdb`.
    If *gzip* is **True**, output files will be compressed.  Return value is
    the output filename or list of filenames, in the order files are processed.
    Note that if multiple models from a file are aligned, that filename will
    appear in the list multiple times."""

    if not isinstance(ensemble, (PDBEnsemble, PDBConformation)):
        raise TypeError('ensemble must be a PDBEnsemble or PDBConformation')
    if isinstance(ensemble, PDBConformation):
        ensemble = [ensemble]
    if gzip:
        gzip = '.gz'
    else:
        gzip = ''
    output = []
    pdbdict = {}
    for conf in ensemble:
        trans = conf.getTransformation()
        if trans is None:
            raise ValueError('transformations are not calculated, call '
                             '`superpose` or `iterpose`')
        label = conf.getLabel()

        pdb = label[:4]
        filename = pdbdict.get(pdb, fetchPDB(pdb))
        if filename is None:
            LOGGER.warning('PDB file for conformation {0} is not found.'
                           .format(label))
            output.append(None)
            continue
        LOGGER.info('Parsing PDB file {0} for conformation {1}.'
                    .format(pdb, label))

        acsi = None
        model = label.rfind('m')
        if model > 3:
            model = label[model+1:]
            if model.isdigit():
                acsi = int(model) - 1
            LOGGER.info('Applying transformation to model {0}.'
                        .format(model))

        if isinstance(filename, str):
            ag = parsePDB(filename)
        else:
            ag = filename

        if acsi is not None:
            if acsi >= ag.numCoordsets():
                LOGGER.warn('Model number {0} for {1} is out of range.'
                            .format(model, pdb))
                output.append(None)
                continue
            ag.setACSIndex(acsi)
        trans.apply(ag)
        outfn = os.path.join(outdir, pdb + suffix + '.pdb' + gzip)
        if ag.numCoordsets() > 1:
            pdbdict[pdb] = ag
        else:
            writePDB(outfn, ag)
        output.append(os.path.normpath(outfn))

    for pdb, ag in pdbdict.items():  # PY3K: OK
        writePDB(os.path.join(outdir, pdb + suffix + '.pdb' + gzip), ag)
    if len(output) == 1:
        return output[0]
    else:
        return output

def calcTree(ensemble, distance_matrix, method='nj'):
    """ Given a distance matrix for an ensemble, it creates an returns a tree structure.
    :arg ensemble: an ensemble with labels. 
    :type ensemble: prody.ensemble.Ensemble or prody.ensemble.PDBEnsemble
    :arg distance_matrix: a square matrix with length of ensemble. If numbers does not mismatch
    it will raise an error. 
    :type distance_matrix: numpy.ndarray 
    """
    try: 
        from Bio import Phylo
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')
    
    names = ensemble.getLabels()
    if len(names) != distance_matrix.shape[0] or len(names) != distance_matrix.shape[1]:
        raise ValueError("The size of matrix and ensemble has a mismatch.")
        
    matrix = []
    k = 1
    for row in distance_matrix:
        matrix.append(list(row[:k]))
        k = k + 1
    from Bio.Phylo.TreeConstruction import _DistanceMatrix
    dm = _DistanceMatrix(names, matrix)
    constructor = Phylo.TreeConstruction.DistanceTreeConstructor()

    method = method.strip().lower()
    if method == 'nj':
        tree = constructor.nj(dm)
    elif method == 'upgma':
        tree = constructor.upgma(dm)
    else:
        raise ValueError('Method can be only either "nj" or "upgma".')

    for node in tree.get_nonterminals():
        node.name = None
    return tree

def buildPDBEnsemble(refpdb, PDBs, title='Unknown', labels=None, seqid=94, coverage=85, mapping_func=mapOntoChain, occupancy=None, unmapped=None):
    """Builds a PDB ensemble from a given reference structure and a list of PDB structures. 
    Note that the reference structure should be included in the list as well.

    :arg refpdb: Reference structure
    :type refpdb: :class:`.Chain`, :class:`.Selection`, or :class:`.AtomGroup`

    :arg PDBs: A list of PDB structures
    :type PDBs: iterable

    :arg title: The title of the ensemble
    :type title: str

    :arg labels: labels of the conformations
    :type labels: list

    :arg seqid: Minimal sequence identity (percent)
    :type seqid: int

    :arg coverage: Minimal sequence overlap (percent)
    :type coverage: int

    :arg occupancy: Minimal occupancy of columns (range from 0 to 1). Columns whose occupancy
    is below this value will be trimmed.
    :type occupancy: float

    :arg unmapped: A list of PDB IDs that cannot be included in the ensemble. This is an 
    output argument. 
    :type unmapped: list
    """

    if not isinstance(refpdb, (Chain, Segment, Selection, AtomGroup)):
        raise TypeError('Refpdb must be a Chain, Segment, Selection, or AtomGroup.')
    
    if labels is not None:
        if len(labels) != len(PDBs):
            raise TypeError('Labels and PDBs must have the same lengths.')

    # obtain the hierarhical view of the reference PDB
    refhv = refpdb.getHierView()
    refchains = list(refhv)

    # obtain the atommap of all the chains combined.
    atoms = refchains[0]
    for i in range(1, len(refchains)):
        atoms += refchains[i]
    
    # initialize a PDBEnsemble with referrence atoms and coordinates
    ensemble = PDBEnsemble(title)
    ensemble.setAtoms(atoms)
    ensemble.setCoords(atoms.getCoords())
    
    # build the ensemble
    if unmapped is None: unmapped = []
    for i,pdb in enumerate(PDBs):
        if not isinstance(pdb, (Chain, Selection, AtomGroup)):
            raise TypeError('PDBs must be a list of Chain, Selection, or AtomGroup.')
        
        if labels is None:
            lbl = pdb.getTitle()
        else:
            lbl = labels[i]

        atommaps = []
        # find the mapping of the pdb to each reference chain
        for chain in refchains:
            mappings = mapping_func(pdb, chain,
                                    seqid=seqid,
                                    coverage=coverage)
            if len(mappings) > 0:
                atommaps.append(mappings[0][0])
            else:
                break

        if len(atommaps) != len(refchains):
            unmapped.append(lbl)
            continue
        
        # combine the mappings of pdb to reference chains
        atommap = atommaps[0]
        for i in range(1, len(atommaps)):
            atommap += atommaps[i]
        
        # add the mappings to the ensemble
        ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'), label = lbl)
    
    if occupancy is not None:
        ensemble = trimPDBEnsemble(ensemble, occupancy=occupancy)
    ensemble.iterpose()

    return ensemble

def addPDBEnsemble(ensemble, PDBs, refpdb=None, labels=None, seqid=94, coverage=85, mapping_func=mapOntoChain, occupancy=None, unmapped=None):  
    """Adds extra structures to a given PDB ensemble. 

    :arg ensemble: The ensemble to which the PDBs are added.
    :type ensemble: :class:`.PDBEnsemble`
    :arg refpdb: Reference structure. If set to `None`, it will be set to `ensemble.getAtoms()` automatically.
    :type refpdb: :class:`.Chain`, :class:`.Selection`, or :class:`.AtomGroup`
    :arg PDBs: A list of PDB structures
    :type PDBs: iterable
    :arg title: The title of the ensemble
    :type title: str
    :arg labels: labels of the conformations
    :type labels: list
    :arg seqid: Minimal sequence identity (percent)
    :type seqid: int
    :arg coverage: Minimal sequence overlap (percent)
    :type coverage: int
    :arg occupancy: Minimal occupancy of columns (range from 0 to 1). Columns whose occupancy
    is below this value will be trimmed.
    :type occupancy: float
    :arg unmapped: A list of PDB IDs that cannot be included in the ensemble. This is an 
    output argument. 
    :type unmapped: list
    """

    if labels is not None:
        if len(labels) != len(PDBs):
            raise TypeError('Labels and PDBs must have the same lengths.')

    # obtain the hierarhical view of the referrence PDB
    if refpdb is None:
        refpdb = ensemble.getAtoms()
    refhv = refpdb.getHierView()
    refchains = list(refhv)

    # obtain the atommap of all the chains combined.
    atoms = refchains[0]
    for i in range(1, len(refchains)):
        atoms += refchains[i]
    
    # add the PDBs to the ensemble
    if unmapped is None: unmapped = []
    for i,pdb in enumerate(PDBs):
        if not isinstance(pdb, (Chain, Selection, AtomGroup)):
            raise TypeError('PDBs must be a list of Chain, Selection, or AtomGroup.')
        
        if labels is None:
            lbl = pdb.getTitle()
        else:
            lbl = labels[i]

        atommaps = []
        # find the mapping of the pdb to each reference chain
        for chain in refchains:
            mappings = mapping_func(pdb, chain,
                                    seqid=seqid,
                                    coverage=coverage)
            if len(mappings) > 0:
                atommaps.append(mappings[0][0])
            else:
                break

        if len(atommaps) != len(refchains):
            unmapped.append(lbl)
            continue
        
        # combine the mappings of pdb to reference chains
        atommap = atommaps[0]
        for i in range(1, len(atommaps)):
            atommap += atommaps[i]
        
        # add the mappings to the ensemble
        ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'), label = lbl)
    
    if occupancy is not None:
        ensemble = trimPDBEnsemble(ensemble, occupancy=occupancy)
    ensemble.iterpose()

    return ensemble

