# -*- coding: utf-8 -*-
"""This module defines functions for performing adaptive ANM."""

from prody.atomic.atomic import Atomic
import time
from numbers import Integral, Number
import numpy as np

from prody import LOGGER
from prody.utilities import getCoords
from prody.measure import calcRMSD, calcDistance, superpose
from prody.ensemble import Ensemble

from .functions import calcENM
from .modeset import ModeSet


__all__ = ['AdaptiveANM', 'runAdaptiveANM', 'runOneWayAdaptiveANM', 'runBothWaysAdaptiveANM']


class AdaptiveANM(object):
    """Class for initialising adaptive ANM analysis of proteins ([ZY09]_).

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

    def __init__(self, structureA, structureB, **kwargs):

        self.coordsA = getCoords(structureA)
        self.coordsA_0 = self.coordsA
        if isinstance(structureA, Atomic):
            self.titleA = structureA.getTitle()
        else:
            self.titleA = 'structure A'

        self.coordsB = getCoords(structureB)
        self.coordsB, _ = superpose(self.coordsB, self.coordsA)
        self.coordsB_0 = self.coordsB
        if isinstance(structureB, Atomic):
            self.titleB = structureB.getTitle()
        else:
            self.titleB = 'structure B'

        if self.coordsA.shape != self.coordsB.shape:
            raise ValueError('structures should have the same number of atoms')

        self.n_atoms = self.coordsA.shape[0]

        rmsd = calcRMSD(self.coordsA, self.coordsB)
        self.rmsds = [rmsd]
        self.dList = []
        self.numSteps = 0

        self.n_modes = kwargs.get('n_modes', 20)
        if not isinstance(self.n_modes, Integral):
            raise TypeError('n_modes should be an integer')

        self.cutoff = kwargs.get('cutoff', 15)
        if not isinstance(self.cutoff, Number):
            raise TypeError('cutoff should be a number')

        self.anmListA = []
        self.anmListB = []
        self.numModesList = []
        self.whichModesA = []
        self.whichModesB = []

        self.maxModes = kwargs.get('maxModes', None)
        dof = 3*self.n_atoms-6
        if self.maxModes is None:
            self.maxModes = dof

        if not isinstance(self.maxModes, (Integral, float)):
            raise TypeError('maxModes should be an integer or float')
        if self.maxModes < 1:
            self.maxModes = int(self.maxModes * dof)
        if self.maxModes > dof:
            self.maxModes = dof

        self.Fmin = kwargs.get('Fmin', None)
        self.Fmin_max = kwargs.get('Fmin_max', 0.6)

        self.ensembleA = Ensemble(self.titleA)
        self.ensembleA.addCoordset(self.coordsA)

        self.ensembleB = Ensemble(self.titleB)
        self.ensembleB.addCoordset(self.coordsB)

        self.ensemble = Ensemble('combined trajectory')
        self.ensemble.setCoords(self.coordsA)

        self.rmsd_diff_cutoff = kwargs.get('rmsd_diff_cutoff', 0.05)
        self.target_rmsd = kwargs.get('target_rmsd', 1.0)

        LOGGER.info(
            'Initialised Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd))


def _runStep(structureA, structureB, **kwargs):
    """Run a single step of adaptive ANM. 
    Modes will be calculated for *structureA* and the subset with 
    a cumulative overlap above a threshold defined by *Fmin* 
    is used for transitioning towards *structureB*.

    By default this function uses values from initialisation but 
    they can be over-ridden if desired. For example, in bi-directional 
    adaptive ANM, we switch *structureA* and *structureB*
    """
    Fmin = kwargs.get('Fmin', None)
    f = kwargs.get('f', 0.2)

    Fmin_max = kwargs.get('Fmin_max', 0.6)
    resetFmin = kwargs.get('resetFmin', False)

    dList = kwargs.get('dList', [])
    numModesList = kwargs.get('numModesList', [])

    cutoff = kwargs.get('cutoff', 15)
    if not isinstance(cutoff, Number):
        raise TypeError('cutoff should be a number')

    coordsA = getCoords(structureA)
    coordsB = getCoords(structureB)
    coordsB, _ = superpose(coordsB, coordsA)

    dof = coordsA.shape[0] - 6
    maxModes = kwargs.get('maxModes', None)
    if maxModes is None:
        maxModes = dof

    if not isinstance(maxModes, (int, float)):
        raise TypeError('maxModes should be an integer or float')
    if maxModes < 1:
        maxModes = int(maxModes * dof)
    if maxModes > dof:
        maxModes = dof

    n_modes = kwargs.pop('n_modes', 20)
    if not isinstance(n_modes, Integral):
        raise TypeError('n_modes should be an integer')

    if n_modes > maxModes:
        n_modes = maxModes

    anmA, _ = calcENM(coordsA, n_modes=n_modes, **kwargs)

    defvec = coordsB - coordsA
    d = defvec.flatten()
    dList.append(d)

    if Fmin is None:
        if resetFmin:
            Fmin = 0.  # Select the first mode only
        else:
            Fmin = 1 - np.sqrt(np.linalg.norm(dList[-1])
                               / np.linalg.norm(dList[0]))

    if Fmin > Fmin_max:
        Fmin = Fmin_max

    LOGGER.info('Fmin is {:4.3f}, corresponding to a cumulative overlap of {:4.3f}'.format(
        Fmin, np.sqrt(Fmin)))

    overlaps = np.dot(d, anmA.getEigvecs())
    overlap_sorting_indices = list(reversed(list(np.argsort(abs(overlaps)))))
    overlaps = overlaps[overlap_sorting_indices]

    modesetA = ModeSet(anmA, overlap_sorting_indices)

    normalised_overlaps = overlaps / np.linalg.norm(d)
    c_sq = np.cumsum(np.power(normalised_overlaps, 2), axis=0)

    modesCrossingFmin = np.where(c_sq <= Fmin)[0]
    numModes = len(modesCrossingFmin)
    if numModes == 0:
        numModes = 1
        modesCrossingFmin = [0]

    numModesList.append(numModes)

    if numModes == 1:
        LOGGER.info('Using 1 mode with overlap {0} (Mode {1})'
                    .format('{:4.3f}'.format(np.sqrt(c_sq[0])), modesetA.getIndices()[0]+1))
    elif numModes < 11:
        LOGGER.info('Using {0} modes with cumulative overlap {1} (Modes {2} and {3})'
                    .format(numModes, '{:4.3f}'.format(np.sqrt(c_sq[numModes-1])),
                            ', '.join([str(entry)
                                       for entry in modesetA.getIndices()[:numModes-1]+1]),
                            str(modesetA.getIndices()[numModes-1]+1)))
    else:
        LOGGER.info('Using {0} modes with cumulative overlap {1} (Modes {2}, ... and {3}) with max mode number {4} and min mode number {5}'
                    .format(numModes, '{:4.3f}'.format(np.sqrt(c_sq[numModes-1])),
                            ', '.join([str(entry)
                                       for entry in modesetA.getIndices()[:10]+1]),
                            str(modesetA.getIndices()[numModes-1]+1),
                            np.max(modesetA.getIndices()[:numModes]+1),
                            np.min(modesetA.getIndices()[:numModes]+1)))

    if np.max(modesetA.getIndices()[:numModes]) > n_modes-5:
        n_modes *= 10

    if n_modes > dof:
        n_modes = dof

    v = np.sum(np.multiply(overlaps[:numModes], modesetA.getEigvecs()[:, :numModes]),
               axis=1).reshape(coordsA.shape)

    s_min = sum(np.multiply(v.flatten(), d)) / sum(np.power(v.flatten(), 2))

    new_coordsA = coordsA + f * s_min * v

    rmsd = calcRMSD(new_coordsA, coordsB)

    LOGGER.info('Current RMSD is {:4.3f}\n'.format(rmsd))

    return rmsd, new_coordsA, modesetA, modesCrossingFmin, dList, n_modes

def _checkConvergence(rmsds, coordsA, coordsB, **kwargs):
    """Check convergence of adaptive ANM. 

    Convergence is reached if one of three conditions is met:
    1. difference between *rmsds* from previous step to current < *rmsd_diff_cutoff*
    2. Current rmsd < *target_rmsd*
    3. A node in *coordsA* or *coordsB* gets disconnected from another by > *cutoff*

    :arg rmsds: a list of RMSDs from Adaptive ANM
    :type rmsds: list

    :arg coordsA: coordinate set A for checking disconnections
    :type coordsA: :class:`~numpy.ndarray`

    :arg coordsA: coordinate set B for checking disconnections
    :type coordsA: :class:`~numpy.ndarray`

    :arg rmsd_diff_cutoff: cutoff for rmsds converging. Default 0.05 A
    :type rmsd_diff_cutoff: float

    :arg target_rmsd: target rmsd for stopping. Default 1.0 A
    :type target_rmsd: float

    :arg cutoff: cutoff for building ANM. Default 15 A
    :type cutoff: float    
    """
    rmsd_diff_cutoff = kwargs.get('rmsd_diff_cutoff', 0.05)
    target_rmsd = kwargs.get('target_rmsd', 1.0)
    cutoff = kwargs.get('cutoff', 15)

    if rmsds[-2] - rmsds[-1] < rmsd_diff_cutoff:
        LOGGER.warn(
            'The RMSD decrease fell below {0}'.format(rmsd_diff_cutoff))
        return True

    if rmsds[-1] < target_rmsd:
        LOGGER.warn('The RMSD fell below target RMSD {0}'.format(target_rmsd))
        return True

    disconA = _checkDisconnection(coordsA, cutoff)
    disconB = _checkDisconnection(coordsB, cutoff)

    if disconA or disconB:
        return True

    return False


def _checkDisconnection(coords, cutoff):
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


def runOneWayAdaptiveANM(structureA, structureB, n_steps, **kwargs):
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
    self = AdaptiveANM(structureA, structureB, **kwargs)

    LOGGER.timeit('_prody_runManySteps')
    n_start = self.numSteps
    resetFmin = True
    while self.numSteps < n_start + n_steps:
        LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1,
                                                                                self.titleA))                                                                   
        rmsd, new_coordsA, modesetA, modesCrossingFmin, self.dList, self.n_modes = _runStep(self.coordsA,
                                                                                                       self.coordsB,
                                                                                                       n_modes=self.n_modes,
                                                                                                       resetFmin=resetFmin,
                                                                                                       dList=self.dList,
                                                                                                       **kwargs)
        self.numSteps += 1
        resetFmin = False
        self.anmListA.append(modesetA)
        self.coordsA = new_coordsA
        self.ensembleA.addCoordset(new_coordsA)
        self.whichModesA.append(modesetA[modesCrossingFmin])
        self.rmsds.append(rmsd)

        LOGGER.debug('Total time so far is %.2f minutes' % (
            (time.time() - LOGGER._times['_prody_runManySteps'])/60))
        converged = _checkConvergence(
            self.rmsds, self.coordsA, self.coordsB)
        if converged:
            # That way the original object is back to normal
            self.coordsA = self.coordsA_0
            # That way the original object is back to normal
            self.coordsB = self.coordsB_0
            LOGGER.debug('Process completed in %.2f hours' % (
                (time.time() - LOGGER._times['_prody_runManySteps'])/3600))
            break

    return self


def runAlternatingAdaptiveANM(structureA, structureB, n_steps, **kwargs):
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
    self = AdaptiveANM(structureA, structureB, **kwargs)

    LOGGER.timeit('_prody_runManySteps')
    n_start = self.numSteps
    resetFmin = True
    while self.numSteps < n_start + n_steps:
        LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1,
                                                                                self.titleA))
        rmsd, new_coordsA, modesetA, modesCrossingFmin, self.dList, self.n_modes = _runStep(self.coordsA,
                                                                                                       self.coordsB,
                                                                                                       n_modes=self.n_modes,
                                                                                                       resetFmin=resetFmin,
                                                                                                       dList=self.dList,
                                                                                                       **kwargs)
        self.numSteps += 1
        resetFmin = False
        self.anmListA.append(modesetA)
        self.coordsA = new_coordsA
        self.ensembleA.addCoordset(new_coordsA)
        self.whichModesA.append(modesetA[modesCrossingFmin])
        self.rmsds.append(rmsd)

        LOGGER.debug('Total time so far is %.2f minutes' % (
            (time.time() - LOGGER._times['_prody_runManySteps'])/60))

        LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1,
                                                                                self.titleA))
        rmsd, new_coordsA, modesetA, modesCrossingFmin, self.dList, self.n_modes = _runStep(self.coordsB,
                                                                                                       self.coordsA,
                                                                                                       n_modes=self.n_modes,
                                                                                                       resetFmin=resetFmin,
                                                                                                       dList=self.dList,
                                                                                                       **kwargs)
        self.numSteps += 1
        self.anmListB.append(modesetA)
        self.coordsB = new_coordsA
        self.ensembleB.addCoordset(new_coordsA)
        self.whichModesB.append(modesetA[modesCrossingFmin])
        self.rmsds.append(rmsd)

        LOGGER.debug('Total time so far is %.2f minutes' % (
            (time.time() - LOGGER._times['_prody_runManySteps'])/60))

        converged = _checkConvergence(
            self.rmsds, self.coordsA, self.coordsB)
        if converged:
            # That way the original object is back to normal
            self.coordsA = self.coordsA_0
            # That way the original object is back to normal
            self.coordsB = self.coordsB_0
            LOGGER.debug('Process completed in %.2f hours' % (
                (time.time() - LOGGER._times['_prody_runManySteps'])/3600))
            break

    for coordset in self.ensembleA.getCoordsets():
        self.ensemble.addCoordset(coordset)
    for coordset in reversed(self.ensembleB.getCoordsets()):
        self.ensemble.addCoordset(coordset)

    return self

runAdaptiveANM = runAlternatingAdaptiveANM


def runBothWaysAdaptiveANM(structureA, structureB, n_steps, **kwargs):
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
    self = AdaptiveANM(structureA, structureB, **kwargs)

    LOGGER.timeit('_prody_runManySteps')
    n_start = self.numSteps

    LOGGER.info('\n\nStarting from structure A ({0})'.format(self.titleA))
    resetFmin = True
    while self.numSteps < n_start + n_steps:
        LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1,
                                                                                self.titleA))
        rmsd, new_coordsA, modesetA, modesCrossingFmin, self.dList, self.n_modes = _runStep(self.coordsA,
                                                                                                       self.coordsB,
                                                                                                       n_modes=self.n_modes,
                                                                                                       resetFmin=resetFmin,
                                                                                                       dList=self.dList,
                                                                                                       **kwargs)
        self.numSteps += 1
        resetFmin = False
        self.anmListA.append(modesetA)
        self.coordsA = new_coordsA
        self.ensembleA.addCoordset(new_coordsA)
        self.whichModesA.append(modesetA[modesCrossingFmin])
        self.rmsds.append(rmsd)

        LOGGER.debug('Total time so far is %.2f minutes' % (
            (time.time() - LOGGER._times['_prody_runManySteps'])/60))
        converged = _checkConvergence(
            self.rmsds, self.coordsA, self.coordsB)
        if converged:
            # That way the original object is back to normal
            #self.coordsA = self.coordsA_0
            # That way the original object is back to normal
            #self.coordsB = self.coordsB_0
            LOGGER.warn('The part starting from structure A converged.')
            break

    LOGGER.info('\n\nStarting from structure B ({0})'.format(self.titleB))
    resetFmin = True
    while self.numSteps < n_start + n_steps:
        LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1,
                                                                                self.titleB))
        rmsd, new_coordsA, modesetA, modesCrossingFmin, self.dList, self.n_modes = _runStep(self.coordsB,
                                                                                                       self.coordsA,
                                                                                                       n_modes=self.n_modes,
                                                                                                       resetFmin=resetFmin,
                                                                                                       dList=self.dList,
                                                                                                       **kwargs)
        self.numSteps += 1
        resetFmin = False
        self.anmListB.append(modesetA)
        self.coordsB = new_coordsA
        self.ensembleB.addCoordset(new_coordsA)
        self.whichModesB.append(modesetA[modesCrossingFmin])
        self.rmsds.append(rmsd)

        LOGGER.debug('Total time so far is %.2f minutes' % (
            (time.time() - LOGGER._times['_prody_runManySteps'])/60))
        converged = _checkConvergence(
            self.rmsds, self.coordsA, self.coordsB)
        if converged:
            # That way the original object is back to normal
            self.coordsA = self.coordsA_0
            # That way the original object is back to normal
            self.coordsB = self.coordsB_0
            LOGGER.warn('The part starting from structure B converged.')
            break

    for coordset in self.ensembleA.getCoordsets():
        self.ensemble.addCoordset(coordset)
    for coordset in reversed(self.ensembleB.getCoordsets()):
        self.ensemble.addCoordset(coordset)

    LOGGER.debug('Process completed in %.2f hours' %
                    ((time.time() - LOGGER._times['_prody_runManySteps'])/3600))

    return self

