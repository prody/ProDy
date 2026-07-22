# -*- coding: utf-8 -*-
"""This module defines functions for performing perturbation response scanning
from PCA and normal modes."""


import numpy as np
from numpy import isin

from prody import LOGGER
from prody.atomic import AtomGroup, Selection, Atomic, sliceAtomicData
from prody.utilities import div0

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector

__all__ = ['calcPerturbResponse', 'comparePerturbResponses',
           'calcDynamicFlexibilityIndex', 'calcDynamicCouplingIndex']

def calcPerturbResponse(model, **kwargs):

    """This function implements the perturbation response scanning (PRS) method
    described in [CA09]_ and [IG14]_. It returns a PRS matrix, and effectiveness 
    and sensitivity profiles.
    
    Rows of the matrix are the average magnitude of the responses obtained by 
    perturbing the atom/node position at that row index, i.e. ``prs_matrix[i,j]`` 
    will give the response of residue/node *j* to perturbations in residue/node *i*. 
    
    PRS is performed using the covariance matrix from a *model*, e.g. 
    a :class:`.ANM` instance. To use an external matrix, please provide it to 
    a :class:`.PCA` instance using the :meth:`.PCA.setCovariance`.

    When an *atoms* instance is given, the PRS matrix will be added as data, 
    which can be retrieved with ``atoms.getData('prs_matrix')``.  

    *model* and *atoms* must have the same number of atoms. *atoms* must be an
    :class:`.AtomGroup` instance. 

    .. [CA09] Atilgan C, Atilgan AR, Perturbation-Response Scanning
       Reveals Ligand Entry-Exit Mechanisms of Ferric Binding Protein.
       *PLoS Comput Biol* **2009** 5(10):e1000544.

    .. [IG14] General IJ, Liu Y, Blackburn ME, Mao W, Gierasch LM, Bahar I.
        ATPase subdomain IA is a mediator of interdomain allostery in Hsp70
        molecular chaperones. *PLoS Comput. Biol.* **2014** 10:e1003624.

    If *turbo* is **True** (default), then PRS is approximated by the limit of 
    large numbers of forces and no perturbation forces are explicitly applied. 
    If set to **False**, then each residue/node is perturbed *repeats* times (default 100) 
    with a random unit force vector as in ProDy v1.8 and earlier.

    If *return_vectors* is **False** (default), then the average NxN PRS matrix
    is calculated and used to calculate effectiveness and sensitivity.
    If **True**, raw responses are returned as a vector ensemble (using :class:`.ModeEnsemble`)
    with shape n_atoms (perturbed) members x *repeats* vectors x n_atoms (moving).
    This has to run through the explicit forces branch, so overrides *turbo*.
    Each response vector is the displacement ``cov . f`` and is stored with its
    squared magnitude as the eigenvalue, so it can be overlapped with a target
    conformational change (e.g. from :func:`.calcDeformVector`) using
    :func:`.comparePerturbResponses`.

    The directional scan can be steered with two options:

    * *force*: a specific perturbing force to apply instead of the random unit
      forces. Either a single 3-vector or an ``(m, 3)`` array of *m* directions.
      When given, *repeats* is ignored (each node gets exactly *m* responses).
    * *perturb_node*: an index or list of indices restricting the scan to the
      selected node(s) rather than every node.

    Either option routes through the explicit-force branch (ignoring *turbo*).
    They still fill the PRS matrix, but only for the perturbed rows (and, with a
    custom *force*, along that direction rather than the isotropic average), so
    when *return_vectors* is **False** a *partial* matrix and its effectiveness
    and sensitivity are returned. A 3d *model* (ANM or PCA) is required whenever
    *return_vectors* is **True**; GNM responses have no direction.

    By default *return_vectors* returns only the response ensemble. Set
    *return_matrix* to **True** together with *return_vectors* to return both,
    as a 4-tuple ``(prs_matrix, effectiveness, sensitivity, response_ensemble)``.
    """

    if not isinstance(model, (NMA, ModeSet, Mode)):
        raise TypeError('model must be an NMA, ModeSet, or Mode instance')

    if isinstance(model, NMA) and len(model) == 0:
        raise ValueError('model must have normal modes calculated')

    atoms = kwargs.get('atoms', None)
    suppress_diag = kwargs.get('suppress_diag', False)
    no_diag = kwargs.get('no_diag', suppress_diag)

    if atoms is not None:
        if isinstance(atoms, Selection):
            atoms = atoms.copy()
        if not isinstance(atoms, AtomGroup):
            raise TypeError('atoms must be an AtomGroup instance')
        elif atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')

    n_atoms = model.numAtoms()
    # LOGGER.timeit('_prody_prs_all')
    # LOGGER.info('Calculating covariance matrix')
    # LOGGER.timeit('_prody_cov')

    cov = model.getCovariance()

    turbo = kwargs.get('turbo', True)
    return_vectors = kwargs.get('return_vectors', False)

    force = kwargs.pop('force', None)
    perturb_node = kwargs.pop('perturb_node', None)
    directional = force is not None or perturb_node is not None

    if turbo and not return_vectors and not directional:
        if not model.is3d():
            prs_matrix = cov**2

        else:
            cov_squared = cov**2
            n_by_3n_cov_squared = np.zeros((n_atoms, 3 * n_atoms))
            prs_matrix = np.zeros((n_atoms, n_atoms))
            i3 = -3
            i3p3 = 0
            for i in range(n_atoms):
                i3 += 3
                i3p3 += 3
                n_by_3n_cov_squared[i,:] = (cov_squared[i3:i3p3,:]).sum(0)

            j3 = -3
            j3p3 = 0
            for j in range(n_atoms):
                j3 += 3
                j3p3 += 3                
                prs_matrix[:,j] = (n_by_3n_cov_squared[:,j3:j3p3]).sum(1)
    else:
        from .signature import ModeEnsemble
        response_ensemble = ModeEnsemble("response ensemble")

        if return_vectors and not model.is3d():
            raise ValueError('return_vectors=True requires a 3d model (ANM or '
                             'PCA); GNM responses have no direction')

        # resolve which nodes to perturb (default: every node)
        if perturb_node is None:
            nodes = list(range(n_atoms))
        elif np.isscalar(perturb_node):
            nodes = [int(perturb_node)]
        else:
            nodes = [int(i) for i in perturb_node]
        if any(i < 0 or i >= n_atoms for i in nodes):
            raise ValueError('perturb_node indices must be in [0, {0})'.format(n_atoms))

        # resolve forces: a fixed direction (or set of directions), otherwise
        # *repeats* random unit forces drawn per node
        fixed_forces = None
        if force is not None:
            fixed_forces = np.atleast_2d(np.asarray(force, dtype=float))
            if fixed_forces.shape[1] != 3:
                raise ValueError('force must be a 3-vector or an (m, 3) array')

        repeats = kwargs.pop('repeats', 100)
        LOGGER.info('Calculating perturbation response with {0} repeats'.format(repeats))
        LOGGER.timeit('_prody_prs_mat')

        # rows for nodes that are not perturbed stay zero (a partial PRS matrix)
        prs_matrix = np.zeros((n_atoms, n_atoms))
        LOGGER.progress('Calculating perturbation response', len(nodes), '_prody_prs')
        for k, i in enumerate(nodes):
            i3 = i * 3
            if fixed_forces is not None:
                forces = fixed_forces
            else:
                forces = np.random.rand(repeats * 3).reshape((repeats, 3))
                forces /= ((forces**2).sum(1)**0.5).reshape((repeats, 1))

            # response displacement vectors cov . f, shape (3*n_atoms, n_forces)
            responses = np.dot(cov[:, i3:i3+3], forces.T)
            vals = (responses ** 2).sum(0)  # squared response magnitude per force

            # average squared response magnitude over the perturbing forces
            prs_matrix[i] = (responses ** 2).reshape(
                (n_atoms, 3, -1)).sum(axis=1).mean(axis=1)

            responses_nma = NMA("response to perturbing {0}th node".format(i))
            responses_nma.setEigens(responses, vals)
            response_ensemble.addModeSet(responses_nma)
            LOGGER.update(k, label='_prody_prs')

        LOGGER.clear()
        LOGGER.report('Perturbation response matrix calculated in %.1fs.',
                    '_prody_prs_mat')

    norm_prs_matrix = np.zeros((n_atoms, n_atoms))
    self_dp = np.diag(prs_matrix)  
    self_dp = self_dp.reshape(n_atoms, 1)
    re_self_dp = np.repeat(self_dp, n_atoms, axis=1)
    norm_prs_matrix = div0(prs_matrix, re_self_dp)

    if no_diag:
       # suppress the diagonal (self displacement) to facilitate
       # visualizing the response profile
       norm_prs_matrix = norm_prs_matrix - np.diag(np.diag(norm_prs_matrix))
    
    W = 1 - np.eye(n_atoms)
    effectiveness = np.average(norm_prs_matrix, weights=W, axis=1)
    sensitivity = np.average(norm_prs_matrix, weights=W, axis=0)

    # LOGGER.report('Perturbation response scanning completed in %.1fs.',
    #               '_prody_prs_all')

    if atoms is not None:
        try:
            ag = atoms.getAtomGroup()
            defdata = np.zeros(ag.numAtoms(), dtype=float)
            ag.setData('effectiveness', defdata.copy())
            ag.setData('sensitivity', defdata.copy())
        except AttributeError:
            pass
        atoms.setData('effectiveness', effectiveness)
        atoms.setData('sensitivity', sensitivity)

        #atoms.setData('prs_matrix', norm_prs_matrix)

    if return_vectors:
        if atoms is not None:
            response_ensemble.setAtoms(atoms)
        if kwargs.get('return_matrix', False):
            return norm_prs_matrix, effectiveness, sensitivity, response_ensemble
        return response_ensemble

    return norm_prs_matrix, effectiveness, sensitivity


def comparePerturbResponses(response_ensemble, target, stat='max'):
    """Overlap analysis of a directional PRS *response_ensemble* against a *target*.

    The *response_ensemble* is a :class:`.ModeEnsemble` as returned by
    :func:`.calcPerturbResponse` with ``return_vectors=True`` (one modeset per
    perturbed node, holding that node's response vectors). Each response is
    overlapped with *target* and the per-force overlaps are reduced to a single
    value per perturbed node, giving a profile of how well perturbing each node
    drives motion resembling *target*.

    :arg target: what to overlap the responses against. Either

        * an :class:`.NMA`, :class:`.ModeSet`, :class:`.Mode` or :class:`.Vector`
          (e.g. an :class:`.ANM` model, one of its modes, or a deformation vector
          from :func:`.calcDeformVector`); for a multi-mode target the overlap is
          taken against the best-matching target mode per force; or
        * another response :class:`.ModeEnsemble`, compared node-by-node (the two
          ensembles must have the same number of modesets).
    :type target: :class:`.NMA`, :class:`.ModeSet`, :class:`.Mode`,
        :class:`.Vector` or :class:`.ModeEnsemble`

    :arg stat: how to reduce each node's absolute overlaps over the force
        directions. One of ``'max'`` (default, the best-matching force), ``'mean'``
        (the isotropic average over force directions) or ``'sum'``.
    :type stat: str

    Returns a 1D :class:`~numpy.ndarray` with one value per perturbed node.
    """
    from .signature import ModeEnsemble
    from .compare import calcOverlap

    if not isinstance(response_ensemble, ModeEnsemble):
        raise TypeError('response_ensemble must be a ModeEnsemble')

    reducers = {'max': np.max, 'mean': np.mean, 'sum': np.sum}
    if stat not in reducers:
        raise ValueError("stat must be one of 'max', 'mean' or 'sum'")
    reducer = reducers[stat]

    if isinstance(target, ModeEnsemble):
        if target.numModeSets() != response_ensemble.numModeSets():
            raise ValueError('the two ensembles must have the same number of '
                             'modesets (perturbed nodes)')
        profile = [reducer(np.abs(calcOverlap(ms_i, ms_j)))
                   for ms_i, ms_j in zip(response_ensemble, target)]

    elif isinstance(target, (NMA, ModeSet, Mode, Vector)):
        profile = []
        for ms in response_ensemble:
            # (n_forces, n_target_modes); keep the best-matching target mode
            ov = np.abs(calcOverlap(ms, target)).reshape(ms.numModes(), -1)
            profile.append(reducer(ov.max(axis=1)))

    else:
        raise TypeError('target must be an NMA, ModeSet, Mode, Vector or '
                        'ModeEnsemble')

    return np.array(profile)


def calcDynamicFlexibilityIndex(matrix, atoms, select, **kwargs):
    """
    Calculate the dynamic flexibility index for the selected residue(s).
    This function implements the dynamic flexibility index (Dfi) method
    described in [ZNG13]_.

    :arg matrix: PRS or covariance matrix, or a model from which to calculate them
    :type matrix: :class:`~numpy.ndarray`, :class:`.ANM`, :class:`.PCA`, :class:`.GNM`

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg select: a selection string or selection for residues of interest
    :type select: str, :class:`.Selection`

    :arg norm: whether to normalise the covariance to a PRS matrix, default False
        This option is only valid when providing a model.
    :type norm: bool

    .. [ZNG13] Gerek ZN, Kumar S, Ozkan SB, Structural dynamics flexibility 
       informs function and evolution at a proteome scale.
       *Evol Appl.* **2013** 6(3):423-33.

    """
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms should be an Atomic object')

    if not isinstance(select, (str, Selection)):
        raise TypeError('select should be a Selection or selection string')

    if isinstance(matrix, NMA):
        model = matrix
        
        norm = kwargs.get('norm', False)
        if norm:
            matrix, _, _ = calcPerturbResponse(model, atoms=atoms, **kwargs)
        else:
            matrix = model.getCovariance()

    elif not isinstance(matrix, np.ndarray):
        raise TypeError('matrix must be an array, ANM, GNM or PCA, not {0}'
                        .format(type(model)))

    profiles = sliceAtomicData(matrix, atoms, select, axis=0)
    return np.sum(profiles, axis=1)/np.sum(matrix)


def calcDynamicCouplingIndex(matrix, atoms, select, func_sel, **kwargs):
    """
    Calculate the dynamic coupling index for the selected residue(s).
    This function implements the dynamic coupling index (DCI) 
    or functional DFI method described in [AK15]_.

    :arg matrix: PRS or covariance matrix, or a model from which to calculate them
    :type matrix: :class:`~numpy.ndarray`, :class:`.ANM`, :class:`.PCA`, :class:`.GNM`

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg select: a selection string or selection for residues of interest
    :type select: str, :class:`.Selection`

    :arg func_sel: a selection string or selection for functional residues
    :type func_sel: str, :class:`.Selection`

    :arg norm: whether to normalise the covariance to a PRS matrix, default False
        This option is only valid when providing a model.
    :type norm: bool

    .. [AK15] Kumar A, Glembo TJ, Ozkan SB. The Role of Conformational Dynamics and Allostery 
        in the Disease Development of Human Ferritin.
       *Biophys J.* **2015** 109(6):1273-81.

    """
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms should be an Atomic object')

    if not isinstance(select, (str, Selection)):
        raise TypeError('select should be a Selection or selection string')
    
    if not isinstance(func_sel, (str, Selection)):
        raise TypeError('func_sel should be a Selection or selection string')

    if isinstance(matrix, NMA):
        model = matrix

        norm = kwargs.get('norm', False)
        if norm:
            matrix, _, _ = calcPerturbResponse(model, atoms=atoms, **kwargs)
        else:
            matrix = model.getCovariance()

    elif not isinstance(matrix, np.ndarray):
        raise TypeError('matrix must be an array, ANM, GNM or PCA, not {0}'
                        .format(type(model)))

    profiles = sliceAtomicData(matrix, atoms, select, axis=0)
    func_profiles = sliceAtomicData(profiles, atoms, func_sel, axis=1)

    if isinstance(func_sel, str):
        func_sel = atoms.select(func_sel)

    N_functional = func_sel.numAtoms()

    numerator = np.sum(func_profiles, axis=1) / N_functional
    denominator = np.sum(profiles, axis=1) / atoms.numAtoms()
    return numerator/denominator
    
