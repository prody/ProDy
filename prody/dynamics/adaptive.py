# -*- coding: utf-8 -*-
"""This module defines functions for performing adaptive ANM."""

from prody.atomic import Atomic, AtomMap
import time
from numbers import Integral, Number
import numpy as np

from prody import LOGGER
from prody.utilities import getCoords, importLA
from prody.measure import calcRMSD, calcDistance, superpose
from prody.ensemble import Ensemble

from .functions import calcENM
from .modeset import ModeSet

__all__ = ['calcAdaptiveANM', 'AANM_ONEWAY', 'AANM_ALTERNATING', 'AANM_BOTHWAYS', 'AANM_DEFAULT']

AANM_ALTERNATING = 0
AANM_ONEWAY = 1
AANM_BOTHWAYS = 2

AANM_DEFAULT = AANM_ALTERNATING

norm = importLA().norm

def checkInput(a, b, **kwargs):
    coordsA = getCoords(a)
    if isinstance(a, Atomic):
        title = a.getTitle()
        atoms = a
    else:
        title = None
        atoms = None

    coordsB = getCoords(b)
    
    if title is None:
        if isinstance(b, Atomic):
            title = b.getTitle()
            atoms = b
        else:
            title = 'Unknown'
            atoms = None

    maskA = a.getFlags("mapped") if isinstance(a, AtomMap) else 1.
    maskB = b.getFlags("mapped") if isinstance(b, AtomMap) else 1.
    weights = maskA * maskB

    if np.isscalar(weights):
        weights = None
    
    if np.isscalar(maskA):
        maskA = None

    if np.isscalar(maskB):
        maskB = None

    aligned = kwargs.get('aligned', False)
    if not aligned:
        coordsA, _ = superpose(coordsA, coordsB, weights)

    rmsd = calcRMSD(coordsA, coordsB, weights)
    LOGGER.info('Initialized Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd))

    return coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd

def getTitle(structure, def_title='structure'):
    if isinstance(structure, Atomic):
        title = structure.getTitle()
    else:
        title = def_title

    return title

def calcStep(initial, target, n_modes, ensemble, defvecs, rmsds, mask=None, callback_func=None, **kwargs):
    """Runs a single step of adaptive ANM. 
    Modes will be calculated for *initial* with a square cumulative overlap above a threshold defined by 
    *Fmin* and used for transitioning towards *target*.
    """

    Fmin = kwargs.get('Fmin', None)
    f = kwargs.get('f', 0.2)

    Fmin_max = kwargs.get('Fmin_max', 0.6)
    resetFmin = kwargs.get('resetFmin', False)

    weights = ensemble.getWeights()
    if weights is not None:
        weights = weights.flatten()
    #coords_init, _ = superpose(initial, target, weights) # we should keep this off otherwise RMSD calculations are off
    coords_init = initial
    coords_tar = target

    dof = coords_init.shape[0] - 6
    n_max_modes = kwargs.get('n_max_modes', None)
    if n_max_modes is None:
        n_max_modes = dof

    if n_max_modes < 1:
        n_max_modes = int(n_max_modes * dof)
    if n_max_modes > dof:
        n_max_modes = dof

    if n_modes > n_max_modes:
        n_modes = n_max_modes

    model = kwargs.pop('model', 'anm')
    anm, _ = calcENM(coords_init, select=mask, mask=mask, 
                     model=model, trim='trim', n_modes=n_modes, 
                     **kwargs)

    if mask is not None:
        anm.masked = False

    defvec = coords_tar - coords_init
    d = defvec.flatten()
    if weights is not None:
        d *= weights.repeat(3)
    defvecs.append(d)

    if Fmin is None:
        if resetFmin:
            Fmin = 0.  # Select the first mode only
        else:
            Fmin = 1 - np.sqrt(norm(defvecs[-1])/norm(defvecs[0]))

    if Fmin > Fmin_max:
        Fmin = Fmin_max

    overlaps = np.dot(d, anm.getEigvecs())

    normalised_overlaps = overlaps / norm(d)
    c_sq = np.cumsum(np.power(normalised_overlaps, 2), axis=0)

    if Fmin == 0 and resetFmin:
        torf_Fmin = np.zeros(c_sq.shape, dtype=bool)
        argmax_overlap = np.argmax(abs(normalised_overlaps))
        torf_Fmin[argmax_overlap] = True
    else:
        torf_Fmin = c_sq <= Fmin
        if np.any(torf_Fmin) and not np.all(torf_Fmin):
            i = np.where(torf_Fmin)[0].max()
            torf_Fmin[i+1] = True

        if not np.any(torf_Fmin):
            torf_Fmin[0] = True

    selected_mode_indices = np.arange(anm.numModes())[torf_Fmin]

    n_sel_modes = len(selected_mode_indices)

    modes = ModeSet(anm, selected_mode_indices)
    c_sq_crit = c_sq[torf_Fmin].max()

    if n_sel_modes == 1:
        LOGGER.info('Using 1 mode with square overlap {0}'
                    .format('%4.3f'%c_sq_crit))
    else:
        LOGGER.info('Using {0} modes with square cumulative overlap {1}'
                    .format(n_sel_modes, '%4.3f'%c_sq_crit))

    if n_sel_modes > n_modes-5:
        n_modes *= 2

    if n_modes > dof:
        n_modes = dof

    v = modes.getEigvecs().dot(overlaps[torf_Fmin])
    s = f * v.dot(d) / v.dot(v)

    # update coords_init
    coords_init += s * v.reshape(coords_init.shape)
    # initial[:] = coords_init[:] # turn this on in case coords_init is not initial in the future
    rmsd = calcRMSD(coords_init, coords_tar, weights)
    rmsds.append(rmsd)

    if callback_func is not None:
        cbkwargs = {'init': coords_init, 
                    'tar': coords_tar, 
                    'modes': modes, 
                    'defvec': d,
                    'c_sq': c_sq_crit,
                    'rmsd': rmsd}
        callback_func(**cbkwargs)

    # deposit 
    ensemble.addCoordset(coords_init.copy())
    converged = checkConvergence(rmsds, coords_init, **kwargs)

    if converged:
        n_modes = 0

    LOGGER.info('Current RMSD is {:4.3f}\n'.format(rmsd))

    return n_modes

def checkConvergence(rmsds, coords, **kwargs):
    """Check convergence of adaptive ANM. 

    Convergence is reached if one of three conditions is met:
    1. Difference between *rmsds* from previous step to current < *min_rmsd_diff*
    2. Current rmsd < *target_rmsd* for the last five runs
    3. A node in *coords* gets disconnected from another by > *cutoff*
    """
    min_rmsd_diff = kwargs.get('min_rmsd_diff', 0.05)
    target_rmsd = kwargs.get('target_rmsd', 1.0)
    cutoff = kwargs.get('cutoff', 15)

    if len(rmsds) > 4:
        drmsd = np.abs(np.diff(rmsds))

        if np.all(drmsd[-4:] < min_rmsd_diff):
            LOGGER.warn(
                'The RMSD decrease fell below {0}'.format(min_rmsd_diff))
            return True

    if rmsds[-1] < target_rmsd:
        LOGGER.warn('The RMSD fell below target RMSD {0}'.format(target_rmsd))
        return True

    if checkDisconnection(coords, cutoff):
        LOGGER.warn('Disconnections were found in one of the structures {0}')
        return True

    return False


def checkDisconnection(coords, cutoff):
    """Check disconnection of ANM, i.e. a node in *coords* gets 
    disconnected from another by > *cutoff*. This is one of the 
    stopping criteria for adaptive ANM. 
    """

    all_dists = np.array([calcDistance(coords, entry) for entry in coords])
    min_dists = np.array([np.min([np.min(all_dists[i, :i]), np.min(all_dists[i, i+1:])])
                          for i in range(1, coords.shape[0]-1)])
    if max(min_dists) > cutoff:
        LOGGER.warn('A bead has become disconnected. '
                    'Adaptive ANM cannot proceed without unrealistic deformations')
        return True

    return False

def calcAdaptiveANM(a, b, n_steps, mode=AANM_DEFAULT, **kwargs):
    """Runs adaptive ANM analysis of proteins ([ZY09]_) that creates a path that 
    connects two conformations using normal modes.

    This function can be run in three modes:
    
    1. *AANM_ONEWAY*: all steps are run in one direction: from *a* to *b*.

    2. *AANM_ALTERNATING*: steps are run in alternating directions: from *a* to *b*, 
        then *b* to *a*, then back again, and so on.

    3. *AANM_BOTHWAYS*: steps are run in one direction (from *a* to 
        *b*) until convergence is reached and then the other way.

    This also implementation differs from the original one in that it sorts the 
    modes by overlap prior to cumulative overlap calculations for efficiency.

    .. [ZY09] Zheng Yang, Peter Májek, Ivet Bahar. Allosteric Transitions of 
            Supramolecular Systems Explored by Network Models: Application to 
            Chaperonin GroEL. *PLOS Comp Biol* **2009** 40:512-524.

    :arg a: structure A for the transition
    :type a: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg b: structure B for the transition
    :type b: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg n_steps: the maximum number of steps to be calculated. For *AANM_BOTHWAYS*, 
        this means the maximum number of steps from each direction
    :type n_steps: int

    :arg mode: the way of the calculation to be performed, which can be either *AANM_ONEWAY*, 
        *AANM_ALTERNATING*, or *AANM_BOTHWAYS*. Default is *AANM_ALTERNATING*
    :type mode: int

    :kwarg f: step size. Default is 0.2
    :type f: float

    :kwarg Fmin: cutoff for selecting modes based on square cumulative overlaps
        Default is **None**, which automatically determines and adapts *Fmin* on the fly.
    :type Fmin: float

    :kwarg Fmin_max: maximum value for *Fmin* when it is automatically determined
        Default is 0.6
    :type Fmin_max: float

    :arg min_rmsd_diff: cutoff for rmsds converging. Default is 0.05
    :type min_rmsd_diff: float

    :kwarg target_rmsd: target rmsd for stopping. Default is 1.0
    :type target_rmsd: float

    :kwarg n_modes: the number of modes to be calculated for the first run. *n_modes* 
        will be dynamically adjusted later as the calculation progresses. Default is 20
    :type n_modes: int

    :kwarg n_max_modes: the maximum number of modes to be calculated in each run. 
        Default is **None**, which allows as many as degree of freedom
    :type n_max_modes: int

    :kwarg callback_func: a callback function that can be used to collect quantities 
        from each iteration. The function must accept `**kwargs` as its only input. 
        Keywords in `kwargs` are:
        'init': the initial coordinate; 
        'tar': the target coordinate; 
        'modes': a :class:`.ModeSet` of selected modes; 
        'defvec': the deformation vector; 
        'c_sq': the critical square cumulative overlap; 
        'rmsd': the RMSD between the two structures after the deformation.
    :type callback_func: func

    Please see keyword arguments for calculating the modes in :func:`.calcENM`.
    """

    if mode == AANM_ONEWAY:
        return calcOneWayAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == AANM_ALTERNATING:
        return calcAlternatingAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == AANM_BOTHWAYS:
        return calcBothWaysAdaptiveANM(a, b, n_steps, **kwargs)
    else:
        raise ValueError('unknown aANM mode: %d'%mode)

def calcOneWayAdaptiveANM(a, b, n_steps, **kwargs):
    """Runs one-way adaptivate ANM. """

    n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd = checkInput(a, b, **kwargs)
    coordsA = coordsA.copy()

    LOGGER.timeit('_prody_calcAdaptiveANM')
    n = 0
    resetFmin = True
    defvecs = []
    rmsds = [rmsd]
    ensemble = Ensemble(title + '_aANM')
    ensemble.setAtoms(atoms)
    ensemble.setCoords(coordsB)
    ensemble.setWeights(weights)
    ensemble.addCoordset(coordsA.copy())
    while n < n_steps:
        LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(n+1, title))                                                            
        n_modes = calcStep(coordsA, coordsB, n_modes, ensemble, defvecs, rmsds, mask=maskA,
                           resetFmin=resetFmin, **kwargs)
        n += 1
        resetFmin = False
        if n_modes == 0:
            LOGGER.report('One-way Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
            break

    return ensemble


def calcAlternatingAdaptiveANM(a, b, n_steps, **kwargs):
    """Runs alternating adaptivate ANM. """

    n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd = checkInput(a, b, **kwargs)
    coordsA = coordsA.copy()
    coordsB = coordsB.copy()

    LOGGER.timeit('_prody_calcAdaptiveANM')
    n = 0
    resetFmin = True
    defvecs = []
    rmsds = [rmsd]
    ensA = Ensemble('A')
    ensA.setCoords(coordsA)
    ensA.setWeights(weights)
    ensA.addCoordset(coordsA.copy())

    ensB = Ensemble('B')
    ensB.setCoords(coordsB.copy())
    ensB.setWeights(weights)
    ensB.addCoordset(coordsB.copy())

    while n < n_steps:
        LOGGER.info('\nStarting cycle {0} with {1}'.format(n + 1, getTitle(a, 'structure A')))
        n_modes = calcStep(coordsA, coordsB, n_modes, ensA, defvecs, rmsds, mask=maskA,
                           resetFmin=resetFmin, **kwargs)
        resetFmin = False

        if n_modes == 0:
            LOGGER.report('Alternating Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
            break

        LOGGER.info('\nContinuing cycle {0} with structure {1}'.format(n+1, getTitle(b, 'structure B')))
        n_modes = calcStep(coordsB, coordsA, n_modes, ensB, defvecs, rmsds, mask=maskB,
                           resetFmin=resetFmin, **kwargs)
        n += 1

        if n_modes == 0:
            LOGGER.report('Alternating Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
            break
    
    ensemble = ensA + ensB[::-1]
    ensemble.setTitle(title + '_aANM')
    ensemble.setAtoms(atoms)
    ensemble.setCoords(ensB.getCoords())

    return ensemble


def calcBothWaysAdaptiveANM(a, b, n_steps, **kwargs):
    """Runs both-way adaptivate ANM. """

    n_modes0 = n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd = checkInput(a, b, **kwargs)
    coordsA = coordsA.copy()
    coordsB = coordsB.copy()

    LOGGER.timeit('_prody_calcAdaptiveANM')
    n = 0
    resetFmin = True
    defvecs = []
    rmsds = [rmsd]
    ensA = Ensemble('A')
    ensA.setCoords(coordsA)
    ensA.setWeights(weights)
    ensA.addCoordset(coordsA.copy())

    ensB = Ensemble('B')
    ensB.setCoords(coordsB.copy())
    ensB.setWeights(weights)
    ensB.addCoordset(coordsB.copy())
    
    while n < n_steps:
        LOGGER.info('\nStarting cycle {0} with {1}'.format(n + 1, getTitle(a, 'structure A')))
        n_modes = calcStep(coordsA, coordsB, n_modes, ensA, defvecs, rmsds, mask=maskA,
                           resetFmin=resetFmin, **kwargs)
        n += 1
        resetFmin = False

        if n_modes == 0:
            break
    
    n = 0
    n_modes = n_modes0
    resetFmin = True
    while n < n_steps:
        LOGGER.info('\nStarting cycle {0} with structure {1}'.format(n+1, getTitle(b, 'structure B')))
        n_modes = calcStep(coordsB, coordsA, n_modes, ensB, defvecs, rmsds, mask=maskB,
                           resetFmin=resetFmin, **kwargs)
        n += 1
        resetFmin = False

        if n_modes == 0:
            LOGGER.report('Alternating Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
            break

    ensemble = ensA + ensB[::-1]
    ensemble.setTitle(title + '_aANM')
    ensemble.setAtoms(atoms)
    ensemble.setCoords(ensB.getCoords())

    LOGGER.report('Both-way Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')

    return ensemble

