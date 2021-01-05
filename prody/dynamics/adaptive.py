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

AANM_ONEWAY = 0
AANM_ALTERNATING = 1
AANM_BOTHWAYS = 2

AANM_DEFAULT = AANM_ONEWAY

norm = importLA().norm

def checkInput(structureA, structureB, **kwargs):
    weights = kwargs.pop('weights', 1.)

    coordsA = getCoords(structureA)
    if isinstance(structureA, Atomic):
        title = structureA.getTitle()
        atoms = structureA
    else:
        title = None
        atoms = None

    coordsB = getCoords(structureB)
    
    if title is None:
        if isinstance(structureB, Atomic):
            title = structureB.getTitle()
            atoms = structureB
        else:
            title = 'Unknown'
            atoms = None

    weightsA = structureA.getFlags("mapped") if isinstance(structureA, AtomMap) else 1.
    weightsB = structureB.getFlags("mapped") if isinstance(structureB, AtomMap) else 1.
    weights = weights * weightsA * weightsB

    if np.isscalar(weights):
        weights = None

    coordsA, _ = superpose(coordsA, coordsB, weights)
    rmsd = calcRMSD(coordsA, coordsB, weights)
    LOGGER.info('Initialized Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd))

    return coordsA, coordsB, title, atoms, weights, rmsd

def getTitle(structure, def_title='structure'):
    if isinstance(structure, Atomic):
        title = structure.getTitle()
    else:
        title = def_title

    return title

def calcStep(initial, target, n_modes, ensemble, defvecs, rmsds, callback_func=None, **kwargs):
    """Run a single step of adaptive ANM. 
    Modes will be calculated for *initial* and the subset with 
    a cumulative overlap above a threshold defined by *Fmin* 
    is used for transitioning towards *target*.

    By default this function uses values from initialisation but 
    they can be over-ridden if desired. For example, in bi-directional 
    adaptive ANM, we switch *initial* and *target*
    """

    Fmin = kwargs.get('Fmin', None)
    f = kwargs.get('f', 0.2)

    Fmin_max = kwargs.get('Fmin_max', 0.6)
    resetFmin = kwargs.get('resetFmin', False)

    weights = ensemble.getWeights()
    if weights is not None:
        weights = weights.flatten()
    # coords_init, _ = superpose(initial, target, weights) # we should keep this off otherwise RMSD calculations are off
    coords_init = initial
    coords_tar = target

    dof = coords_init.shape[0] - 6
    maxModes = kwargs.get('maxModes', None)
    if maxModes is None:
        maxModes = dof

    if not isinstance(maxModes, (int, float)):
        raise TypeError('maxModes should be an integer or float')
    if maxModes < 1:
        maxModes = int(maxModes * dof)
    if maxModes > dof:
        maxModes = dof

    if n_modes > maxModes:
        n_modes = maxModes

    mask = None if weights is None else weights != 0
    anm, _ = calcENM(coords_init, select=mask, mask=mask, 
                     model='anm', trim='trim', n_modes=n_modes, 
                     **kwargs)

    if mask is not None:
        anm.masked = False

    defvec = coords_tar - coords_init
    d = defvec.flatten()
    defvecs.append(d)

    if Fmin is None:
        if resetFmin:
            Fmin = 0.  # Select the first mode only
        else:
            Fmin = 1 - np.sqrt(norm(defvecs[-1])/norm(defvecs[0]))

    if Fmin > Fmin_max:
        Fmin = Fmin_max

    overlaps = np.abs(np.dot(d, anm.getEigvecs())) 
    sorted_indices = overlaps.argsort()[::-1]
    overlaps = overlaps[sorted_indices]

    sorted_mode_indices = np.arange(anm.numModes())[sorted_indices]

    normalised_overlaps = overlaps / norm(d)
    c_sq = np.cumsum(np.power(normalised_overlaps, 2), axis=0)

    torf_Fmin = c_sq <= Fmin
    if not np.any(torf_Fmin):
        torf_Fmin[0] = True

    selected_mode_indices = sorted_mode_indices[torf_Fmin]

    n_sel_modes = len(selected_mode_indices)

    modes = ModeSet(anm, selected_mode_indices)
    mode_ids = modes.getIndices()

    if n_sel_modes == 1:
        LOGGER.info('Using 1 mode with square overlap {0} (Mode {1})'
                    .format('%4.3f'%c_sq[0], mode_ids[0]+1))
    else:
        LOGGER.info('Using {0} modes with square cumulative overlap {1} (Max mode number {2})'
                    .format(n_sel_modes, '%4.3f'%c_sq[torf_Fmin].max(), np.max(mode_ids)+1))

    if np.max(mode_ids) > n_modes-5:
        n_modes *= 10

    if n_modes > dof:
        n_modes = dof

    v = modes.getEigvecs().dot(overlaps[torf_Fmin])
    s = f * v.dot(d) / v.dot(v)

    # update coords_init
    coords_init += s * v.reshape(coords_init.shape)
    #initial[:] = coords_init[:] # turn this on in case coords_init is not initial in the future
    rmsd = calcRMSD(coords_init, coords_tar, weights)
    rmsds.append(rmsd)

    if callback_func is not None:
        cbkwargs = {'coords_init': coords_init, 
                    'coords_tar': coords_tar, 
                    'modes': modes, 
                    'defvec': d}
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
    1. difference between *rmsds* from previous step to current < *min_rmsd_diff*
    2. Current rmsd < *target_rmsd*
    3. A node in *coords* gets disconnected from another by > *cutoff*

    :arg rmsds: a list of RMSDs from Adaptive ANM
    :type rmsds: list

    :arg coordsA: coordinate set A for checking disconnections
    :type coordsA: :class:`~numpy.ndarray`

    :arg min_rmsd_diff: cutoff for rmsds converging. Default 0.01 A
    :type min_rmsd_diff: float

    :arg target_rmsd: target rmsd for stopping. Default 1.0 A
    :type target_rmsd: float

    :arg cutoff: cutoff for building ANM. Default 15 A
    :type cutoff: float    
    """
    min_rmsd_diff = kwargs.get('min_rmsd_diff', 0.01)
    target_rmsd = kwargs.get('target_rmsd', 1.0)
    cutoff = kwargs.get('cutoff', 15)

    if min_rmsd_diff is not None:
        if rmsds[-2] - rmsds[-1] < min_rmsd_diff:
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

    :arg coords: a coordinate set for checking disconnections
    :type coords: :class:`~numpy.ndarray`

    :arg cutoff: cutoff for building ANM. Default 15 A
    :type cutoff: float    
    """
    all_dists = np.array([calcDistance(coords, entry) for entry in coords])
    min_dists = np.array([np.min([np.min(all_dists[i, :i]), np.min(all_dists[i, i+1:])])
                          for i in range(1, coords.shape[0]-1)])
    if max(min_dists) > cutoff:
        LOGGER.warn('A bead has become disconnected. '
                    'Adaptive ANM cannot proceed without unrealistic deformations')
        return True

    return False

def calcAdaptiveANM(structureA, structureB, n_steps, mode=AANM_DEFAULT, **kwargs):
    if mode == AANM_ONEWAY:
        return calcOneWayAdaptiveANM(structureA, structureB, n_steps, **kwargs)
    elif mode == AANM_ALTERNATING:
        return calcAlternatingAdaptiveANM(structureA, structureB, n_steps, **kwargs)
    elif mode == AANM_BOTHWAYS:
        return calcBothWaysAdaptiveANM(structureA, structureB, n_steps, **kwargs)
    else:
        raise ValueError('unknown aANM mode: %d'%mode)

def calcOneWayAdaptiveANM(structureA, structureB, n_steps, **kwargs):
    """Run a modified version of adaptive ANM analysis of proteins ([ZY09]_) 
    where all steps are run in one direction: from *structureA* to *structureB*.

    This also implementation differs from the original one in that it sorts the 
    modes by overlap prior to cumulative overlap calculations for efficiency.

    .. [ZY09] Zheng Yang, Peter Májek, Ivet Bahar. Allosteric Transitions of 
            Supramolecular Systems Explored by Network Models: Application to 
            Chaperonin GroEL. *PLOS Comp Biol* **2009** 40:512-524.

    :arg structureA: starting structure for the transition
    :type structureA: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg structureB: starting structure for the transition
    :type structureB: :class:`.Atomic`, :class:`~numpy.ndarray`
    """

    n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, rmsd = checkInput(structureA, structureB, **kwargs)
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
        n_modes = calcStep(coordsA, coordsB, n_modes, ensemble, defvecs, rmsds,
                           resetFmin=resetFmin, **kwargs)
        n += 1
        resetFmin = False
        if n_modes == 0:
            LOGGER.report('One-way Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
            break

    return ensemble


def calcAlternatingAdaptiveANM(structureA, structureB, n_steps, **kwargs):
    """Run the traditional version of adaptive ANM analysis of proteins ([ZY09]_) 
    where steps are run in alternating directions: from *structureA* to *structureB*, 
    then *structureB* to *structureA*, then back again, and so on.

    This implementation differs from the original one in that it sorts the 
    modes by overlap prior to cumulative overlap calculations for efficiency.

    .. [ZY09] Zheng Yang, Peter Májek, Ivet Bahar. Allosteric Transitions of 
            Supramolecular Systems Explored by Network Models: Application to 
            Chaperonin GroEL. *PLOS Comp Biol* **2009** 40:512-524.

    :arg structureA: starting structure for the transition
    :type structureA: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg structureB: starting structure for the transition
    :type structureB: :class:`.Atomic`, :class:`~numpy.ndarray`
    """

    n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, rmsd = checkInput(structureA, structureB, **kwargs)
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
        LOGGER.info('\nStarting cycle {0} with {1}'.format(n + 1, getTitle(structureA, 'structure A')))
        n_modes = calcStep(coordsA, coordsB, n_modes, ensA, defvecs, rmsds,
                           resetFmin=resetFmin, **kwargs)
        n += 1
        resetFmin = False

        if n_modes == 0:
            LOGGER.report('Alternating Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
            break

        LOGGER.info('\nStarting cycle {0} with structure {1}'.format(n+1, getTitle(structureB, 'structure B')))
        n_modes = calcStep(coordsB, coordsA, n_modes, ensB, defvecs, rmsds,
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


def calcBothWaysAdaptiveANM(structureA, structureB, n_steps, **kwargs):
    """Run a modified version of adaptive ANM analysis of proteins ([ZY09]_) 
    where all steps are run in one direction (from *structureA* to *structureB*) 
    until convergence is reached and then the other way.

    This also implementation differs from the original one in that it sorts the 
    modes by overlap prior to cumulative overlap calculations for efficiency.

    .. [ZY09] Zheng Yang, Peter Májek, Ivet Bahar. Allosteric Transitions of 
            Supramolecular Systems Explored by Network Models: Application to 
            Chaperonin GroEL. *PLOS Comp Biol* **2009** 40:512-524.

    :arg structureA: starting structure for the transition
    :type structureA: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg structureB: starting structure for the transition
    :type structureB: :class:`.Atomic`, :class:`~numpy.ndarray`
    """

    n_modes0 = n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, rmsd = checkInput(structureA, structureB, **kwargs)
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
        LOGGER.info('\nStarting cycle {0} with {1}'.format(n + 1, getTitle(structureA, 'structure A')))
        n_modes = calcStep(coordsA, coordsB, n_modes, ensA, defvecs, rmsds,
                           resetFmin=resetFmin, **kwargs)
        n += 1
        resetFmin = False

        if n_modes == 0:
            LOGGER.report('Alternating Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
            break
    
    n = 0
    n_modes = n_modes0
    resetFmin = True
    while n < n_steps:
        LOGGER.info('\nStarting cycle {0} with structure {1}'.format(n+1, getTitle(structureB, 'structure B')))
        n_modes = calcStep(coordsB, coordsA, n_modes, ensB, defvecs, rmsds,
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

