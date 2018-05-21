"""This module defines a functions for handling conformational ensembles."""

import os.path
import time
from numbers import Integral

import numpy as np

from prody.proteins import fetchPDB, parsePDB, writePDB, mapOntoChain
from prody.utilities import openFile, showFigure, copy, isListLike
from prody import LOGGER, SETTINGS
from prody.atomic import AtomMap, Chain, AtomGroup, Selection, Segment, Select, AtomSubset

from .ensemble import *
from .pdbensemble import *
from .conformation import *

__all__ = ['saveEnsemble', 'loadEnsemble', 'trimPDBEnsemble',
           'calcOccupancies', 'showOccupancies', 'alignPDBEnsemble',
           'buildPDBEnsemble', 'addPDBEnsemble', 'refineEnsemble']


def saveEnsemble(ensemble, filename=None, **kwargs):
    """Save *ensemble* model data as :file:`filename.ens.npz`.  If *filename*
    is **None**, title of the *ensemble* will be used as the filename, after
    white spaces in the title are replaced with underscores.  Extension is
    :file:`.ens.npz`. Upon successful completion of saving, filename is
    returned. This function makes use of :func:`~numpy.savez` function."""

    if not isinstance(ensemble, Ensemble):
        raise TypeError('invalid type for ensemble, {0}'
                        .format(type(ensemble)))
    if len(ensemble) == 0:
        raise ValueError('ensemble instance does not contain data')

    dict_ = ensemble.__dict__
    attr_list = ['_title', '_confs', '_weights', '_coords', '_indices']
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

    atoms = dict_['_atoms']
    if atoms:
        atoms = atoms.copy()
    attr_dict['_atoms'] = np.array([atoms, None])

    if isinstance(ensemble, PDBEnsemble):
        msa = dict_['_msa']
        attr_dict['_msa'] = np.array([msa, None])

    if filename.endswith('.ens'):
        filename += '.npz'
    if not filename.endswith('.npz'):
        filename += '.ens.npz'
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename


def loadEnsemble(filename, **kwargs):
    """Returns ensemble instance loaded from *filename*.  This function makes
    use of :func:`~numpy.load` function.  See also :func:`saveEnsemble`"""

    if not 'encoding' in kwargs:
        kwargs['encoding'] = 'latin1'
    attr_dict = np.load(filename, **kwargs)
    if '_weights' in attr_dict:
        weights = attr_dict['_weights']
    else:
        weights = None   
    isPDBEnsemble = False
    try:
        title = attr_dict['_title']
    except KeyError:
        title = attr_dict['_name']
    if isinstance(title, np.ndarray):
        title = np.asarray(title, dtype=str)
    title = str(title)
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
    if '_indices' in attr_dict:
        indices = attr_dict['_indices']
    else:
        indices = None
    ensemble._indices = indices
    if isPDBEnsemble:
        confs = attr_dict['_confs']
        ensemble.addCoordset(confs, weights)
        if '_identifiers' in attr_dict.files:
            ensemble._labels = list(attr_dict['_identifiers'])
        if '_labels' in attr_dict.files:
            ensemble._labels = list(attr_dict['_labels'])
        if ensemble._labels:
            for i, label in enumerate(ensemble._labels):
                if not isinstance(label, str):
                    try:
                        ensemble._labels[i] = label.decode()
                    except AttributeError:
                        ensemble._labels[i] = str(label)
        if '_trans' in attr_dict.files:
            ensemble._trans = attr_dict['_trans']
        if '_msa' in attr_dict.files:
            ensemble._msa = attr_dict['_msa'][0]
    else:
        ensemble.addCoordset(attr_dict['_confs'])
        if weights is not None:
            ensemble.setWeights(weights)
    return ensemble


def trimPDBEnsemble(pdb_ensemble, occupancy=None, **kwargs):
    """Returns a new PDB ensemble obtained by trimming given *pdb_ensemble*.
    This function helps selecting atoms in a pdb ensemble based on one of the
    following criteria, and returns them in a new :class:`.PDBEnsemble`
    instance.

    Resulting PDB ensemble will contain atoms whose occupancies are greater
    or equal to *occupancy* keyword argument. Occupancies for atoms will be
    calculated using ``calcOccupancies(pdb_ensemble, normed=True)``.

    :arg occupancy: occupancy for selecting atoms, must satisfy
        ``0 < occupancy <= 1``.
        If set to *None* then *hard* trimming will be performed.
    :type occupancy: float

    :arg hard: Whether to perform hard trimming.
        Default is **False**
        If set to **True**, atoms will be completely removed from *pdb_ensemble*.
        If set to **False**, a soft trimming of *pdb_ensemble* will be done
        where atoms will be removed from the active selection. This is useful, 
        for example, when one uses :func:`calcEnsembleENMs` 
        together with :func:`sliceModel` or :func:`reduceModel`
        to calculate the modes from the remaining part while still taking the 
        removed part into consideration (e.g. as the environment).
    :type hard: bool
    """

    hard = kwargs.pop('hard', False) or pdb_ensemble._atoms is None \
           or occupancy is None

    atoms = pdb_ensemble.getAtoms(selected=hard)

    if not isinstance(pdb_ensemble, PDBEnsemble):
        raise TypeError('pdb_ensemble argument must be a PDBEnsemble')
    if pdb_ensemble.numConfs() == 0 or pdb_ensemble.numAtoms() == 0:
        raise ValueError('pdb_ensemble must have conformations')

    if occupancy is not None:
        occupancy = float(occupancy)
        assert 0 < occupancy <= 1, ('occupancy is not > 0 and <= 1: '
                                    '{0}'.format(repr(occupancy)))
        n_confs = pdb_ensemble.numConfs()
        assert n_confs > 0, 'pdb_ensemble does not contain any conformations'
        occupancies = calcOccupancies(pdb_ensemble, normed=True)
        #assert weights is not None, 'weights must be set for pdb_ensemble'
        #weights = weights.flatten()
        #mean_weights = weights / n_confs
        torf = occupancies >= occupancy
    else:
        n_atoms = pdb_ensemble.getCoords().shape[0]
        torf = np.ones(n_atoms, dtype=bool)

    trimmed = PDBEnsemble(pdb_ensemble.getTitle())
    if hard:
        if atoms is not None:
            trim_atoms_idx = [n for n,t in enumerate(torf) if t]
            trim_atoms = atoms[trim_atoms_idx]
            trimmed.setAtoms(trim_atoms)

        coords = pdb_ensemble.getCoords()
        if coords is not None:
            trimmed.setCoords(coords[torf])
        confs = pdb_ensemble.getCoordsets()
        if confs is not None:
            weights = pdb_ensemble.getWeights()
            labels = pdb_ensemble.getLabels()
            msa = pdb_ensemble.getMSA()
            if msa:
                msa = msa[:, torf]
            trimmed.addCoordset(confs[:, torf], weights[:, torf], labels, sequence=msa)
    else:
        indices = np.where(torf)[0]
        selids = pdb_ensemble._indices

        if selids is not None:
            indices = selids[indices]

        select = atoms[indices]
        trimmed.setAtoms(atoms)
        trimmed.setAtoms(select)

        coords = copy(pdb_ensemble._coords)
        if coords is not None:
            trimmed.setCoords(coords)
        confs = copy(pdb_ensemble._confs)
        if confs is not None:
            weights = copy(pdb_ensemble._weights)
            labels = pdb_ensemble.getLabels()
            msa = pdb_ensemble._msa
            trimmed.addCoordset(confs, weights, labels, sequence=msa)

        trimmed.setAtoms(select)

    return trimmed

def calcOccupancies(pdb_ensemble, normed=False):
    """Returns occupancy calculated from weights of a :class:`.PDBEnsemble`.
    Any non-zero weight will be considered equal to one.  Occupancies are
    calculated by binary weights for each atom over the conformations in
    the ensemble. When *normed* is **True**, total weights will be divided
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
    If ``gzip=True``, output files will be compressed.  Return value is
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


def buildPDBEnsemble(PDBs, ref=None, title='Unknown', labels=None, seqid=94, coverage=85, 
                     mapping_func=mapOntoChain, unmapped=None, **kwargs):
    """Builds a PDB ensemble from a given reference structure and a list of PDB structures. 
    Note that the reference structure should be included in the list as well.

    :arg PDBs: A list of PDB structures
    :type PDBs: iterable

    :arg ref: Reference structure or the index to the reference in ``PDBs``. If **None**,
        then the first item in ``PDBs`` will be considered as the reference. 
        Default is **None**
    :type ref: int, :class:`.Chain`, :class:`.Selection`, or :class:`.AtomGroup`

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

    occupancy = kwargs.pop('occupancy', None)
    degeneracy = kwargs.pop('degeneracy', True)
    subset = str(kwargs.get('subset', 'calpha')).lower()

    if len(PDBs) == 1:
        raise ValueError('PDBs should have at least two items')

    if labels is not None:
        if len(labels) != len(PDBs):
            raise ValueError('labels and PDBs must be the same length')

    if ref is None:
        refpdb = PDBs[0]
    elif isinstance(ref, Integral):
        refpdb = PDBs[ref]
    else:
        refpdb = ref
        if refpdb not in PDBs:
            raise ValueError('refpdb should be also in the PDBs')

    # obtain refchains from the hierarhical view of the reference PDB
    if subset != 'all':
        refpdb = refpdb.select(subset)
        
    try:
        refchains = list(refpdb.getHierView())
    except AttributeError:
        raise TypeError('refpdb must have getHierView')

    start = time.time()
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

    LOGGER.progress('Building the ensemble...', len(PDBs), '_prody_buildPDBEnsemble')
    for i, pdb in enumerate(PDBs):
        LOGGER.update(i, 'Mapping %s to the reference...'%pdb.getTitle(), 
                      label='_prody_buildPDBEnsemble')
        try:
            pdb.getHierView()
        except AttributeError:
            raise TypeError('PDBs must be a list of instances having the access to getHierView')
            
        if labels is None:
            lbl = pdb.getTitle()
        else:
            lbl = labels[i]

        atommaps = []
        # find the mapping of the pdb to each reference chain
        for chain in refchains:
            mappings = mapping_func(pdb, chain,
                                    seqid=seqid,
                                    coverage=coverage,
                                    index=i,
                                    **kwargs)
            if len(mappings) > 0:
                atommaps.append(mappings[0][0])
            else:
                break

        if len(atommaps) != len(refchains):
            unmapped.append(lbl)
            continue
        
        # combine the mappings of pdb to reference chains
        atommap = atommaps[0]
        for j in range(1, len(atommaps)):
            atommap += atommaps[j]
        
        # add the mappings to the ensemble
        ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'), 
                             label = lbl, degeneracy=degeneracy)

    LOGGER.finish()

    if occupancy is not None:
        ensemble = trimPDBEnsemble(ensemble, occupancy=occupancy)
    ensemble.iterpose()
    
    LOGGER.info('Ensemble ({0} conformations) were built in {1:.2f}s.'
                     .format(ensemble.numConfs(), time.time()-start))

    return ensemble

def addPDBEnsemble(ensemble, PDBs, refpdb=None, labels=None, seqid=94, coverage=85, 
                   mapping_func=mapOntoChain, occupancy=None, unmapped=None, **kwargs):  
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

    degeneracy = kwargs.pop('degeneracy', True)

    if labels is not None:
        if len(labels) != len(PDBs):
            raise TypeError('Labels and PDBs must have the same lengths.')
    else:
        labels = []
        
        for pdb in PDBs:
            if pdb is None:
                labels.append(None)
            else:
                labels.append(pdb.getTitle())

    # obtain refchains from the hierarhical view of the reference PDB
    if refpdb is None:
        refpdb = ensemble.getAtoms()
    refchains = list(refpdb.getHierView())

    start = time.time()

    # obtain the atommap of all the chains combined.
    atoms = refchains[0]
    for i in range(1, len(refchains)):
        atoms += refchains[i]
    
    # add the PDBs to the ensemble
    if unmapped is None: unmapped = []

    LOGGER.progress('Appending the ensemble...', len(PDBs), '_prody_addPDBEnsemble')
    for i, pdb in enumerate(PDBs):
        lbl = labels[i]
        if pdb is None:
            unmapped.append(labels[i])
            continue

        LOGGER.update(i, 'Mapping %s to the reference...'%pdb.getTitle(), 
                      label='_prody_addPDBEnsemble')
        if not isinstance(pdb, (Chain, Selection, AtomGroup)):
            raise TypeError('PDBs must be a list of Chain, Selection, or AtomGroup.')

        atommaps = []
        # find the mapping of the pdb to each reference chain
        for chain in refchains:
            mappings = mapping_func(pdb, chain,
                                    seqid=seqid,
                                    coverage=coverage,
                                    index=i,
                                    **kwargs)
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
        ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'), 
                             label=lbl, degeneracy=degeneracy)
    LOGGER.finish()

    if occupancy is not None:
        ensemble = trimPDBEnsemble(ensemble, occupancy=occupancy)
    ensemble.iterpose()

    LOGGER.info('{0} PDBs were added to the ensemble in {1:.2f}s.'
                     .format(len(PDBs) - len(unmapped), time.time()-start))

    return ensemble

def refineEnsemble(ens, lower=.5, upper=10.):
    """Refine a PDB ensemble based on RMSD criterions.""" 

    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    ### calculate RMSDs ###
    RMSD = ens.getRMSDs(pairwise=True)
    rmsd = ens.getRMSDs()

    ### imposeing upper bound ###
    I = np.where(rmsd < upper)[0]
    reens = ens[I]
    I = I.reshape(-1, 1)
    reRMSD = RMSD[I, I.T]

    ### hierarchical clustering ###
    v = squareform(reRMSD)
    Z = linkage(v)

    labels = fcluster(Z, lower, criterion='distance')
    uniq_labels = unique(labels)

    clusters = []
    for label in uniq_labels:
        indices = np.where(labels==label)[0]
        clusters.append(indices)

    J = ones(len(clusters), dtype=int) * -1
    for i, cluster in enumerate(clusters):
        if len(cluster) > 0:
            weights = [ens[j].getWeights().sum() for j in cluster]
            j = np.argmax(weights)
            J[i] = cluster[j]
        else:
            J[i] = cluster[0]

    ### refine ensemble ###
    reens = reens[J]

    return reens
