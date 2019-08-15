# -*- coding: utf-8 -*-
"""This module defines functions for performing adaptive ANM."""

import time
from numbers import Integral
from numpy import ndarray
import numpy as np
import matplotlib.pyplot as plt

from prody import LOGGER
from prody.measure import calcDeformVector, calcRMSD, superpose, applyTransformation, calcDistance
from prody.atomic import Atomic, AtomGroup, sliceAtomicData
from prody.utilities import getCoords
from .compare import calcOverlap, calcCumulOverlap
from prody.ensemble import PDBEnsemble, Ensemble
from .mode import Vector
from .anm import ANM
from .functions import calcENM, reduceModel, sliceModel
from .compare import matchModes
from .modeset import ModeSet
from prody.trajectory import writeDCD
from prody.proteins import writePDB, mapOntoChain, mapChainByChain

__all__ = ['AdaptiveANM']


class AdaptiveANM(object):
    def __init__(self, structA, structB, **kwargs):
        self.structA = structA
        self.structB = structB
        self.coordsA = structA.getCoords()
        self.coordsB = structB.getCoords()

        self.cutoff = kwargs.pop('cutoff', 15.0) # default cutoff for ANM.buildHessian
        self.gamma = kwargs.pop('gamma', 1.0) # default gamma for ANM.buildHessian

        self.alignSel = kwargs.get('alignSel', None)
        self.alignSelA = kwargs.get('alignSelA', self.alignSel)
        self.alignSelB = kwargs.get('alignSelB', self.alignSel)

        self.reduceSel = kwargs.pop('reduceSel', None)
        self.reduceSelA = kwargs.get('reduceSelA', self.reduceSel)
        self.reduceSelB = kwargs.get('reduceSelB', self.reduceSel)

        if self.alignSelA is None:
            self.alignSelA = self.reduceSelA

        if self.alignSelB is None:
            self.alignSelB = self.reduceSelB

        if self.reduceSelA is None:
            self.reduceSelA = self.alignSelA

        if self.reduceSelB is None:
            self.reduceSelB = self.alignSelB

        if self.alignSelA is None:
            structA_sel = self.structA
        else:
            structA_sel = self.structA.select(self.alignSelA)

        if self.alignSelB is None:
            structB_sel = self.structB
        else:
            structB_sel = self.structB.select(self.alignSelB)

        mapping_func = kwargs.pop('mapping_func', mapChainByChain)
        seqid = kwargs.pop('seqid', 90.) 
        coverage = kwargs.pop('overlap', 70.)
        coverage = kwargs.pop('coverage', coverage)
        pwalign = kwargs.pop('pwalign', None)
        pwalign = kwargs.pop('mapping', pwalign)

        try:
            _, T = superpose(structA_sel, structB_sel)
            structA = applyTransformation(T, structA)
            structB_amap = structB_sel
        except:
            structB_amap = sum(np.array(mapping_func(structB_sel, structA_sel, overlap=coverage, seqid=seqid, pwalign=pwalign))[:,0])
            _, T = superpose(structA_sel, structB_amap)
            structA = applyTransformation(T, structA)

        rmsd = calcRMSD(structA_sel, structB_amap)
        self.rmsds = [rmsd]
        self.dList = []
        self.numSteps = 0

        self.trim = kwargs.pop('trim', 'slice')

        self.anmA = kwargs.pop('anmA', None)
        self.anmB = kwargs.pop('anmB', None)

        if self.anmA is not None:
            self.anmListA = [self.anmA]
        else:
            self.anmListA = []

        if self.anmB is not None:
            self.anmListB = [self.anmB]
        else:
            self.anmListB = []

        self.n_modes = 20
        self.numModesList = []
        self.whichModesA = []
        self.whichModesB = []

        self.plotRMSD = kwargs.get('plotRMSD', False)
        self.plotNumModes = kwargs.get('plotNumModes', False)

        self.maxModes = kwargs.get('maxModes', 0.05)
        if not isinstance(self.maxModes, (int,float)):
            raise TypeError('maxModes should be an integer or float')
        if self.maxModes < 1:
            self.maxModes = int(self.maxModes * 3*self.structA.numAtoms()-6)
        if self.maxModes > 3*self.structA.numAtoms()-6:
            self.maxModes = 3*self.structA.numAtoms()-6

        self.Fmin = kwargs.get('Fmin', None)
        self.Fmin_max = kwargs.get('Fmin_max', 0.6)
        self.resetFmin = kwargs.get('resetFmin', False)
        self.f = kwargs.get('f', 0.2)

        self.ensembleA = Ensemble(structA)
        self.ensembleA.addCoordset(structA)

        self.ensembleB = Ensemble(structB)
        self.ensembleB.addCoordset(structB)

        self.outputDCD = kwargs.get('outputDCD', False)
        self.outputPDB = kwargs.get('outputPDB', False)

        self.filename = kwargs.get('filename', None)
        if self.filename is None:
            self.filename = self.structA.getTitle().replace(' ', '_')

        self.rmsd_diff_cutoff = kwargs.get('rmsd_diff_cutoff', 0.05)
        self.target_rmsd = kwargs.get('target_rmsd', 1.0)

        LOGGER.info('Initialised Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd))

    def runStep(self, structA=None, structB=None, **kwargs):

        if structA is None:
            structA = self.structA

        if structB is None:
            structB = self.structB

        alignSelA = kwargs.pop('alignSel', self.alignSelA)
        alignSelB = kwargs.pop('alignSel', self.alignSelB)

        reduceSel = kwargs.pop('reduceSel', self.reduceSel)
        reduceSelA = kwargs.pop('reduceSelA', self.reduceSelA)
        reduceSelB = kwargs.pop('reduceSelB', self.reduceSelB)

        if reduceSelA is None:
            reduceSelA = reduceSel

        if reduceSelB is None:
            reduceSelB = reduceSel

        if alignSelA is None:
            alignSelA = reduceSelA

        if alignSelB is None:
            alignSelB = reduceSelB

        Fmin = kwargs.get('Fmin', self.Fmin)

        f = kwargs.get('f', self.f)

        outputDCD = kwargs.get('outputDCD', self.outputDCD)
        outputPDB = kwargs.get('outputPDB', self.outputPDB)
        filename = kwargs.get('filename', self.filename)

        plotRMSD = kwargs.get('plotRMSD', self.plotRMSD)
        plotNumModes = kwargs.get('plotNumModes', self.plotNumModes)

        LOGGER.info('\nStarting cycle {0} with initial structure {1}'.format(self.numSteps+1, structA))

        mapping_func = kwargs.get('mapping_func', mapOntoChain)

        if alignSelA is None:
            structA_sel = structA
        else:
            structA_sel = structA.select(alignSelA)

        if alignSelB is None:
            structB_sel = structB
        else:
            structB_sel = structB.select(alignSelB)

        try:
            _, T = superpose(structA_sel, structB_sel)
            structA = applyTransformation(T, structA)
        except:
            structB_amap = sum(np.array(mapping_func(structB_sel, structA_sel))[:,0], **kwargs)
            _, T = superpose(structA_sel, structB_amap)
            structA = applyTransformation(T, structA)

        trim = kwargs.pop('trim', self.trim)
        anmA, _ = calcENM(structA, n_modes=self.n_modes)

        if trim == 'slice':
            trim_anmA, _ = sliceModel(anmA, structA, reduceSelA)
        elif trim == 'reduce':
            trim_anmA, _ = reduceModel(anmA, structA, reduceSelA)
            trim_anmA.calcModes(n_modes=self.n_modes)
        else:
            trim_anmA = anmA

        coordsA = structA.getCoords()
        coordsA_sel = structA_sel.getCoords()
        coordsB_sel = structB_sel.getCoords()

        defvec = coordsB_sel - coordsA_sel
        d = defvec.flatten()
        self.dList.append(d)

        if Fmin is None:
            if self.numSteps == 0 or self.resetFmin:
                Fmin = 0. # Select the first mode only
            else:
                Fmin = 1 - np.sqrt(np.linalg.norm(self.dList[self.numSteps])
                                   / np.linalg.norm(self.dList[0]))

        if Fmin > self.Fmin_max:
            Fmin = self.Fmin_max

        LOGGER.info('Fmin is {:4.3f}, corresponding to a cumulative overlap of {:4.3f}'.format(
            Fmin, np.sqrt(Fmin)))

        trim_d = sliceAtomicData(d, structA_sel, reduceSelA)
        overlaps = np.dot(trim_d, trim_anmA.getEigvecs())
        overlap_sorting_indices = list(reversed(list(np.argsort(abs(overlaps)))))
        overlaps = overlaps[overlap_sorting_indices]

        if trim == 'reduce':
            sliced_anmA, _ = sliceModel(anmA, structA, reduceSelA)
            modesetA = ModeSet(trim_anmA, overlap_sorting_indices)
            _, overlap_sorting_indices = matchModes(modesetA, sliced_anmA, index=True)

        modesetA = ModeSet(anmA, overlap_sorting_indices)

        normalised_overlaps = overlaps / np.linalg.norm(d)
        c_sq = np.cumsum(np.power(normalised_overlaps, 2), axis=0)

        modesCrossingFmin = np.where(c_sq <= Fmin)[0]
        numModes = len(modesCrossingFmin)
        if numModes == 0:
            numModes = 1
            modesCrossingFmin = [0]

        maxModes = kwargs.get('maxModes', self.maxModes)
        if not isinstance(maxModes, (int,float)):
            raise TypeError('maxModes should be an integer or float')
        if maxModes < 1:
            maxModes = int(maxModes * 3*self.structA.numAtoms()-6)
        if maxModes > 3*self.structA.numAtoms()-6:
            maxModes = 3*self.structA.numAtoms()-6

        if numModes > maxModes:
            numModes = maxModes

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
            LOGGER.info('Using {0} modes with cumulative overlap {1} (Modes {2}, ... and {3}) with max mode number {4}'
                        .format(numModes, '{:4.3f}'.format(np.sqrt(c_sq[numModes-1])),
                                ', '.join([str(entry)
                                           for entry in modesetA.getIndices()[:10]+1]),
                                str(modesetA.getIndices()[numModes-1]+1),
                                np.max(modesetA.getIndices()[:numModes]+1)))

        if np.max(modesetA.getIndices()[:numModes]) > self.n_modes-5:
            self.n_modes *= 10
        
        if self.n_modes > 3*self.structA.numAtoms()-6:
            self.n_modes = 3*self.structA.numAtoms()-6

        v = np.sum(np.multiply(overlaps[:numModes], modesetA.getEigvecs()[:, :numModes]),
                   axis=1).reshape(coordsA.shape)

        trim_v = sliceAtomicData(v.reshape(-1), structA, reduceSelA).reshape(-1, 3)
        s_min = sum(np.multiply(trim_v.flatten(), trim_d))/sum(np.power(trim_v.flatten(), 2))

        new_coordsA = coordsA + f * s_min * v

        if structA == self.structA:
            self.anmA = anmA
            self.anmListA.append(modesetA)
            self.structA.setCoords(new_coordsA)
            self.ensembleA.addCoordset(new_coordsA)
            self.whichModesA.append(modesetA[modesCrossingFmin])
        elif structA == self.structB:
            self.anmB = anmA
            self.anmListB.append(modesetA)
            self.structB.setCoords(new_coordsA)
            self.ensembleB.addCoordset(new_coordsA)
            self.whichModesB.append(modesetA[modesCrossingFmin])

        new_coordsA_reduceSel = structA.select(reduceSelA).getCoords()
        coordsB_reduceSel = structB.select(reduceSelB).getCoords()
        rmsd = calcRMSD(new_coordsA_reduceSel, coordsB_reduceSel)

        LOGGER.info('Current RMSD is {:4.3f}\n'.format(rmsd))

        if plotRMSD:
            plt.figure(1); plt.plot(self.numSteps, rmsd)

        if plotNumModes:
            plt.figure(2); plt.bar(self.numSteps, numModes)

        self.numSteps += 1

        self.rmsds.append(rmsd)

        if outputPDB:
            writePDB(filename + '_A', self.ensembleA)
            LOGGER.clear()
            writePDB(filename + '_B', self.ensembleB)
            LOGGER.clear()

        if outputDCD:
            writeDCD(filename + '_A', self.ensembleA)
            LOGGER.clear()
            writeDCD(filename + '_B', self.ensembleB)
            LOGGER.clear()

        return


    def checkConvergence(self, **kwargs):
        converged = False
        rmsd_diff_cutoff = kwargs.get('rmsd_diff_cutoff', self.rmsd_diff_cutoff)
        target_rmsd = kwargs.get('target_rmsd', self.target_rmsd)

        if self.rmsds[-2] - self.rmsds[-1] < rmsd_diff_cutoff:
            LOGGER.warn('The RMSD decrease fell below {0}'.format(rmsd_diff_cutoff))
            converged = True

        if self.rmsds[-1] < target_rmsd:
            LOGGER.warn('The RMSD fell below target RMSD {0}'.format(target_rmsd))
            converged = True

        all_dists = np.array([calcDistance(self.structA, self.structA[i]) 
                              for i in range(self.structA.numAtoms())])
        min_dists = np.array([np.min([np.min(all_dists[i, :i]),np.min(all_dists[i, i+1:])]) 
                              for i in range(1,self.structA.numAtoms()-1)])
        if max(min_dists) > self.cutoff:
            LOGGER.warn('A bead has become disconnected. Adaptive ANM cannot proceed without unrealistic deformations')
            converged = True

        return converged


    def runManySteps(self, n_steps, **kwargs):
        n_start = self.numSteps
        while self.numSteps < n_start + n_steps:
            self.runStep(self.structA, self.structB, **kwargs)
            converged = self.checkConvergence()
            if converged:
                self.structA.setCoords(self.coordsA) # That way the original object is back to normal
                self.structB.setCoords(self.coordsB) # That way the original object is back to normal
                break


    def runManyStepsAlternating(self, n_steps, **kwargs):
        n_start = self.numSteps
        while self.numSteps < n_start + n_steps:
            self.runStep(self.structA, self.structB, **kwargs)
            self.runStep(self.structB, self.structA, **kwargs)

            converged = self.checkConvergence()
            if converged:
                self.structA.setCoords(self.coordsA) # That way the original object is back to normal
                self.structB.setCoords(self.coordsB) # That way the original object is back to normal
                break

        ensemble = Ensemble('combined trajectory')
        ensemble.setAtoms(self.structA)
        for coordset in self.ensembleA.getCoordsets():
            ensemble.addCoordset(coordset)
        for coordset in reversed(self.ensembleB.getCoordsets()):
            ensemble.addCoordset(coordset)

        if self.outputPDB:
            writePDB(self.filename, ensemble)

        if self.outputDCD:
            writeDCD(self.filename, ensemble)

        return


    def runManyStepsFurthestEachWay(self, n_steps, **kwargs):
        n_start = self.numSteps

        LOGGER.info('\n\nStarting from struct A ({0})'.format(self.structA))
        while self.numSteps < n_start + n_steps:
            self.runStep(self.structA, self.structB, **kwargs)
            converged = self.checkConvergence()
            if converged:
                self.structA.setCoords(self.coordsA) # That way the original object is back to normal
                self.structB.setCoords(self.coordsB) # That way the original object is back to normal
                LOGGER.warn('The part starting from structA converged.')
                break

        LOGGER.info('\n\nStarting from structB ({0})'.format(self.structB))
        self.resetFmin = True
        while self.numSteps < n_start + n_steps:
            self.runStep(self.structB, self.structA, **kwargs)
            self.resetFmin = False
            converged = self.checkConvergence()
            if converged:
                self.structA.setCoords(self.coordsA) # That way the original object is back to normal
                self.structB.setCoords(self.coordsB) # That way the original object is back to normal
                LOGGER.warn('The part starting from structB converged.')
                break

        ensemble = Ensemble('combined trajectory')
        ensemble.setAtoms(self.structA)
        for coordset in self.ensembleA.getCoordsets():
            ensemble.addCoordset(coordset)
        for coordset in reversed(self.ensembleB.getCoordsets()):
            ensemble.addCoordset(coordset)

        if self.outputPDB:
            writePDB(self.filename, ensemble)

        if self.outputDCD:
            writeDCD(self.filename, ensemble)

        return
        