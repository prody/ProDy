# -*- coding: utf-8 -*-
"""This module defines functions for performing adaptive ANM."""

import time
from numbers import Integral, Number
import numpy as np

from prody import LOGGER
from prody.measure import calcRMSD, calcDistance
from prody.ensemble import PDBEnsemble, Ensemble, buildPDBEnsemble

from .functions import calcENM
from .modeset import ModeSet


__all__ = ['AdaptiveANM']


class AdaptiveANM(object):
    """Class for adaptive ANM analysis of proteins ([ZY09]_).

    This implementation differs from the original one in that it sorts the 
    modes by overlap prior to cumulative overlap calculations for efficiency.

    .. [ZY09] Zheng Yang, Peter Májek, Ivet Bahar. Allosteric Transitions of 
              Supramolecular Systems Explored by Network Models: Application to 
              Chaperonin GroEL. *PLOS Comp Biol* **2009** 40:512-524.

    :arg ensemble: ensemble or list containing the two endpoints for the transition
    :type ensemble: :class:`.PDBEnsemble`, list
    """

    def __init__(self, ensemble, **kwargs):
    
        if not isinstance(ensemble, PDBEnsemble):
            try:
                ensemble = buildPDBEnsemble(ensemble)
            except:
                raise TypeError('ensemble should be a PDBEnsemble or list of structures')

        if ensemble.numConfs() != 2:
            raise ValueError('ensemble should have two conformers')

        self.ensemble = ensemble
        
        self.coordsA = ensemble.getCoordsets()[0]
        self.coordsA_0 = self.coordsA
        self.titleA = ensemble.getLabels()[0]

        self.coordsB = ensemble.getCoordsets()[1]
        self.coordsB_0 = self.coordsB
        self.titleB = ensemble.getLabels()[1]

        rmsd = ensemble.getRMSDs(pairwise=True)[0, 1]
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
        dof = 3*self.ensemble.numAtoms()-6
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
        self.resetFmin = kwargs.get('resetFmin', False)
        self.f = kwargs.get('f', 0.2)

        self.ensembleA = Ensemble(self.ensemble.getLabels()[0])
        self.ensembleA.addCoordset(self.coordsA)

        self.ensembleB = Ensemble(self.ensemble.getLabels()[1])
        self.ensembleB.addCoordset(self.coordsB)

        self.rmsd_diff_cutoff = kwargs.get('rmsd_diff_cutoff', 0.05)
        self.target_rmsd = kwargs.get('target_rmsd', 1.0)

        LOGGER.info('Initialised Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd))

    def runStep(self, **kwargs):
        """Run a single step of adaptive ANM. 
        Modes will be calculated for *coordsA* and the subset with 
        a cumulative overlap above a threshold defined by *Fmin* 
        is used for transitioning towards *coordsB*.

        By default this function uses values from initialisation but 
        they can be over-ridden if desired. For example, in bi-directional 
        adaptive ANM, we switch *coordsA* and *coordsA*
        """

        coordsA = kwargs.pop('coordsA', self.coordsA)
        coordsB = kwargs.pop('coordsB', self.coordsB)

        Fmin = kwargs.get('Fmin', self.Fmin)
        f = kwargs.get('f', self.f)

        if np.allclose(coordsA, self.coordsA):
            LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1,
                                                                                 self.titleA))
        else:
            LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1,
                                                                                 self.titleB))            

        dof = 3*self.ensemble.numAtoms()-6
        maxModes = kwargs.get('maxModes', self.maxModes)
        if not isinstance(maxModes, (int, float)):
            raise TypeError('maxModes should be an integer or float')
        if maxModes < 1:
            maxModes = int(maxModes * dof)
        if maxModes > dof:
            maxModes = dof

        if self.n_modes > maxModes:
            self.n_modes = maxModes

        anmA, _ = calcENM(coordsA, n_modes=self.n_modes, cutoff=self.cutoff)

        defvec = coordsB - coordsA
        d = defvec.flatten()
        self.dList.append(d)

        if Fmin is None:
            if self.numSteps == 0 or self.resetFmin:
                Fmin = 0.  # Select the first mode only
            else:
                Fmin = 1 - np.sqrt(np.linalg.norm(self.dList[self.numSteps])
                                   / np.linalg.norm(self.dList[0]))

        if Fmin > self.Fmin_max:
            Fmin = self.Fmin_max

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

        self.numModesList.append(numModes)

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

        if np.max(modesetA.getIndices()[:numModes]) > self.n_modes-5:
            self.n_modes *= 10

        if self.n_modes > dof:
            self.n_modes = dof

        v = np.sum(np.multiply(overlaps[:numModes], modesetA.getEigvecs()[:, :numModes]),
                   axis=1).reshape(coordsA.shape)

        s_min = sum(np.multiply(v.flatten(), d)) / sum(np.power(v.flatten(), 2))

        new_coordsA = coordsA + f * s_min * v

        if np.allclose(coordsA, self.coordsA):
            self.anmA = anmA
            self.anmListA.append(modesetA)
            self.coordsA = new_coordsA
            self.ensembleA.addCoordset(new_coordsA)
            self.whichModesA.append(modesetA[modesCrossingFmin])
        elif np.allclose(coordsA, self.coordsB):
            self.anmB = anmA
            self.anmListB.append(modesetA)
            self.coordsB = new_coordsA
            self.ensembleB.addCoordset(new_coordsA)
            self.whichModesB.append(modesetA[modesCrossingFmin])

        rmsd = calcRMSD(new_coordsA, coordsB)

        LOGGER.info('Current RMSD is {:4.3f}\n'.format(rmsd))

        self.numSteps += 1
        self.rmsds.append(rmsd)

        return

    def checkConvergence(self, **kwargs):
        converged = False
        rmsd_diff_cutoff = kwargs.get('rmsd_diff_cutoff',
                                      self.rmsd_diff_cutoff)
        target_rmsd = kwargs.get('target_rmsd', self.target_rmsd)

        if self.rmsds[-2] - self.rmsds[-1] < rmsd_diff_cutoff:
            LOGGER.warn('The RMSD decrease fell below {0}'.format(
                rmsd_diff_cutoff))
            converged = True

        if self.rmsds[-1] < target_rmsd:
            LOGGER.warn('The RMSD fell below target RMSD {0}'.format(
                target_rmsd))
            converged = True

        all_dists = np.array([calcDistance(self.coordsA, self.coordsA[i])
                              for i in range(self.ensemble.numAtoms())])
        min_dists = np.array([np.min([np.min(all_dists[i, :i]), np.min(all_dists[i, i+1:])])
                              for i in range(1, self.ensemble.numAtoms()-1)])
        if max(min_dists) > self.cutoff:
            LOGGER.warn('A bead has become disconnected. '
                        'Adaptive ANM cannot proceed without unrealistic deformations')
            converged = True

        return converged

    def runManySteps(self, n_steps, **kwargs):
        LOGGER.timeit('_prody_runManySteps')
        n_start = self.numSteps
        while self.numSteps < n_start + n_steps:
            self.runStep(coordsA=self.coordsA, coordsB=self.coordsB, **kwargs)
            LOGGER.debug('Total time so far is %.2f minutes' % (
                (time.time() - LOGGER._times['_prody_runManySteps'])/60))
            converged = self.checkConvergence()
            if converged:
                # That way the original object is back to normal
                self.coordsA = self.coordsA_0
                # That way the original object is back to normal
                self.coordsB = self.coordsB_0
                LOGGER.debug('Process completed in %.2f hours' % (
                    (time.time() - LOGGER._times['_prody_runManySteps'])/3600))
                break

    def runManyStepsAlternating(self, n_steps, **kwargs):
        LOGGER.timeit('_prody_runManySteps')
        n_start = self.numSteps
        while self.numSteps < n_start + n_steps:
            n_modes = self.n_modes

            self.runStep(coordsA=self.coordsA, coordsB=self.coordsB,
                         n_modes=n_modes, **kwargs)
            LOGGER.debug('Total time so far is %.2f minutes' % (
                (time.time() - LOGGER._times['_prody_runManySteps'])/60))

            self.runStep(coordsA=self.coordsB, coordsB=self.coordsA, 
                         n_modes=n_modes, **kwargs)
            LOGGER.debug('Total time so far is %.2f minutes' % (
                (time.time() - LOGGER._times['_prody_runManySteps'])/60))

            converged = self.checkConvergence()
            if converged:
                # That way the original object is back to normal
                self.coordsA = self.coordsA_0
                # That way the original object is back to normal
                self.coordsB = self.coordsB_0
                LOGGER.debug('Process completed in %.2f hours' % (
                    (time.time() - LOGGER._times['_prody_runManySteps'])/3600))
                break

        ensemble = Ensemble('combined trajectory')
        ensemble.setAtoms(self.ensemble.getAtoms())
        ensemble.setCoords(self.coordsA)
        for coordset in self.ensembleA.getCoordsets():
            ensemble.addCoordset(coordset)
        for coordset in reversed(self.ensembleB.getCoordsets()):
            ensemble.addCoordset(coordset)

        return

    def runManyStepsFurthestEachWay(self, n_steps, **kwargs):
        LOGGER.timeit('_prody_runManySteps')
        n_start = self.numSteps
        n_modes = self.n_modes

        LOGGER.info('\n\nStarting from structure A ({0})'.format(self.ensemble.getLabels()[0]))
        while self.numSteps < n_start + n_steps:
            self.runStep(coordsA=self.coordsA, coordsB=self.coordsB, 
                         n_modes=n_modes, **kwargs)
            LOGGER.debug('Total time so far is %.2f minutes' % (
                (time.time() - LOGGER._times['_prody_runManySteps'])/60))
            converged = self.checkConvergence()
            if converged:
                # That way the original object is back to normal
                #self.coordsA = self.coordsA_0
                # That way the original object is back to normal
                #self.coordsB = self.coordsB_0
                LOGGER.warn('The part starting from structure A converged.')
                break

        LOGGER.info('\n\nStarting from structure B ({0})'.format(self.ensemble.getLabels()[1]))
        self.resetFmin = True
        while self.numSteps < n_start + n_steps:
            self.runStep(coordsA=self.coordsB, coordsB=self.coordsA,
                         n_modes=n_modes, **kwargs)
            LOGGER.debug('Total time so far is %.2f minutes' % (
                (time.time() - LOGGER._times['_prody_runManySteps'])/60))
            self.resetFmin = False
            converged = self.checkConvergence()
            if converged:
                # That way the original object is back to normal
                self.coordsA = self.coordsA_0
                # That way the original object is back to normal
                self.coordsB = self.coordsB_0
                LOGGER.warn('The part starting from structure B converged.')
                break

        ensemble = Ensemble('combined trajectory')
        ensemble.setAtoms(self.ensemble.getAtoms())
        ensemble.setCoords(self.coordsA)
        for coordset in self.ensembleA.getCoordsets():
            ensemble.addCoordset(coordset)
        for coordset in reversed(self.ensembleB.getCoordsets()):
            ensemble.addCoordset(coordset)

        LOGGER.debug('Process completed in %.2f hours' %
                     ((time.time() - LOGGER._times['_prody_runManySteps'])/3600))

        return
