# -*- coding: utf-8 -*-
"""This module defines functions for performing adaptive ANM."""

__author__ = 'James Krieger'
__credits__ = ['Hongchun Li', 'She Zhang', 'Burak Kaynak']
__email__ = ['jamesmkrieger@gmail.com', 'hongchun28@gmail.com', 'shz66@pitt.edu', 'burak.kaynak@pitt.edu']

from itertools import product
from multiprocessing import cpu_count, Pool
from collections import OrderedDict
from os import chdir, mkdir
from os.path import isdir
from sys import stdout

import time
from numbers import Integral, Number
from decimal import Decimal, ROUND_HALF_UP
import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomMap
from prody.utilities import getCoords, createStringIO, importLA, checkCoords, copy

from prody.measure.transform import calcTransformation, superpose, applyTransformation, calcRMSD, getRMSD
from prody.measure.measure import calcDeformVector, calcDistance

from prody.ensemble.ensemble import Ensemble
from prody.proteins.pdbfile import writePDBStream, parsePDBStream

from .functions import calcENM
from .modeset import ModeSet
from .nma import NMA
from .sampling import traverseMode
from .hybrid import Hybrid

__all__ = ['calcAdaptiveANM', 'ONEWAY', 'ALTERNATING', 'SERIAL', 'DEFAULT',
           'AdaptiveHybrid']

ALTERNATING = 0
ONEWAY = 1
SERIAL = 2

DEFAULT = ALTERNATING

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
       for the last five runs
    2. Current rmsd < *target_rmsd*
    3. A node in *coords* gets disconnected from another by > *cutoff*
    """
    min_rmsd_diff = kwargs.get('min_rmsd_diff', 0.05)
    target_rmsd = kwargs.get('target_rmsd', 1.0)
    cutoff = kwargs.get('cutoff', 15)

    if len(rmsds) > 4:
        drmsd = np.abs(np.diff(rmsds))

        if np.all(drmsd[-4:] < min_rmsd_diff):
            LOGGER.warn(
                'The RMSD decrease stayed below {0} for five cycles'.format(min_rmsd_diff))
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

def calcAdaptiveANM(a, b, n_steps, mode=DEFAULT, **kwargs):
    """Runs adaptive ANM analysis of proteins ([ZY09]_) that creates a path that 
    connects two conformations using normal modes.

    This function can be run in three modes:
    
    1. *ONEWAY*: all steps are run in one direction: from *a* to *b*.

    2. *ALTERNATING*: steps are run in alternating directions: from *a* to *b*, 
        then *b* to *a*, then back again, and so on.

    3. *SERIAL*: steps are run in one direction (from *a* to 
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

    :arg n_steps: the maximum number of steps to be calculated. For *SERIAL*, 
        this means the maximum number of steps from each direction
    :type n_steps: int

    :arg mode: the way of the calculation to be performed, which can be either *ONEWAY*, 
        *ALTERNATING*, or *SERIAL*. Default is *ALTERNATING*
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

    if mode == ONEWAY:
        return calcOneWayAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == ALTERNATING:
        return calcAlternatingAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == SERIAL:
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

class AdaptiveHybrid(Hybrid):
    '''
    This is a new version of the Adaptive ANM transition sampling algorithm [ZY09]_, 
    that is an ANM-based hybrid algorithm. It requires PDBFixer and OpenMM for performing 
    energy minimization and MD simulations in implicit/explicit solvent.

    Instantiate a ClustENM-like hybrid method object for AdaptiveANM.
    '''
    def __init__(self, title, **kwargs):
        super().__init__(title=title)
        self._atomsB = None
        self._coordsB = None
        self._indicesB = None
        
        self._defvecs = []
        self._rmsds = []
        self._cg_ensA = Ensemble(title=title)
        self._cg_ensB = Ensemble(title=title)
        self._n_modes0 = self._n_modes = kwargs.pop('n_modes', 20)

    def _sample(self, conf, **kwargs):

        tmp = self._atoms.copy()
        tmp.setCoords(conf)
        cg = tmp[self._idx_cg]

        anm_cg = self._buildANM(cg)

        if not self._checkANM(anm_cg):
            return None

        tmpB = self._atomsB.copy()
        cgB = tmpB[self._idx_cg]

        coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd = checkInput(cg, cgB, **kwargs)
        coordsA = coordsA.copy()

        self._direction_mode = kwargs.get('mode', DEFAULT)

        if self._direction_mode == ONEWAY:
            LOGGER.info('\nStarting cycle with structure A')
            self._n_modes = calcStep(coordsA, coordsB, self._n_modes, self._cg_ensA,
                                     self._defvecs, self._rmsds, mask=maskA,
                                     resetFmin=self._resetFmin, **kwargs)
            self._resetFmin = False
            cg_ens = self._cg_ensA

        elif self._direction_mode == ALTERNATING:
            if self._direction == 1:
                LOGGER.info('\nStarting cycle with structure A')
                self._n_modes = calcStep(coordsA, coordsB, self._n_modes, self._cg_ensA,
                                         self._defvecs, self._rmsds, mask=maskA,
                                         resetFmin=self._resetFmin, **kwargs)
                self._resetFmin = False
                cg_ens = self._cg_ensA

            else:
                LOGGER.info('\nStarting cycle with structure B')
                self._n_modes = calcStep(coordsB, coordsA, self._n_modes, self._cg_ensB,
                                         self._defvecs, self._rmsds, mask=maskB,
                                         resetFmin=self._resetFmin, **kwargs)
                cg_ens = self._cg_ensB

        elif self._direction_mode == SERIAL:
            if self._direction == 1:
                LOGGER.info('\nStarting cycle with structure A')
                self._n_modes = calcStep(coordsA, coordsB, self._n_modes, self._cg_ensA,
                                         self._defvecs, self._rmsds, mask=maskA,
                                         resetFmin=self._resetFmin, **kwargs)
                self._resetFmin = False
                cg_ens = self._cg_ensA

            else:
                LOGGER.info('\nStarting cycle with structure B')
                self._n_modes = calcStep(coordsB, coordsA, self._n_modes, self._cg_ensB,
                                         self._defvecs, self._rmsds, mask=maskB,
                                         resetFmin=self._resetFmin, **kwargs)
                self._resetFmin = False
                cg_ens = self._cg_ensB

        else:
            raise ValueError('unknown aANM mode: %d' % self._direction_mode)
        
        if self._direction == 1:
            defvec = calcDeformVector(cg, cg_ens.getCoordsets()[-1])
            model = NMA()
            model.setEigens(defvec.getArray().reshape((defvec.getArray().shape[0], 1)))
            model_ex = self._extendModel(model, cg, tmp)
            def_ens = traverseMode(model_ex[0], tmp, 1, rmsd)
            coordsets = [def_ens.getCoordsets()[-1]]

            if self._direction_mode == ALTERNATING:
                self._direction = 2
        else:
            defvec = calcDeformVector(cgB, cg_ens.getCoordsets()[-1])
            model = NMA()
            model.setEigens(defvec.getArray().reshape((defvec.getArray().shape[0], 1)))
            model_ex = self._extendModel(model, cgB, tmpB)
            def_ens = traverseMode(model_ex[0], tmpB, 1, rmsd)
            coordsets = [def_ens.getCoordsets()[-1]]

            if self._direction_mode == ALTERNATING:
                self._direction = 1

        if self._targeted:
            if self._parallel:
                with Pool(cpu_count()) as p:
                    pot_conf = p.map(self._multi_targeted_sim,
                                     [(conf, coords) for coords in coordsets])
            else:
                pot_conf = [self._multi_targeted_sim((conf, coords)) for coords in coordsets]

            pots, poses = list(zip(*pot_conf))

            idx = np.logical_not(np.isnan(pots))
            coordsets = np.array(poses)[idx]

            LOGGER.debug('%d/%d sets of coordinates were moved to the target' % (len(poses), len(coordsets)))

        return coordsets

    def setAtoms(self, atomsA, atomsB=None, pH=7.0, **kwargs):
        aligned = kwargs.get('aligned', False)
        if not aligned and atomsB is not None:
            T = calcTransformation(atomsA.ca, atomsB.ca, weights=atomsA.ca.getFlags("mapped"))
            _ = applyTransformation(T, atomsA)

        if self._isBuilt():
            super(Hybrid, self).setAtoms(atomsA)
            self._atomsB = atomsB
        else:
            A_B_dict = {0: 'A', 1: 'B'}
            for i, atoms in enumerate([atomsA, atomsB]):
                if i == 1 and atomsB is None:
                    break

                atoms = atoms.select('not hetatm')

                self._nuc = atoms.select('nucleotide')

                if self._nuc is not None:

                    idx_p = []
                    for c in self._nuc.getChids():
                        tmp = self._nuc[c].iterAtoms()
                        for a in tmp:
                            if a.getName() in ['P', 'OP1', 'OP2', 'OP3']:
                                idx_p.append(a.getIndex())

                    if idx_p:
                        nsel = 'not index ' + ' '.join([str(i) for i in idx_p])
                        atoms = atoms.select(nsel)

                LOGGER.info('Fixing structure {0}...'.format(A_B_dict[i]))
                LOGGER.timeit('_clustenm_fix')
                self._ph = pH
                self._fix(atoms, i)
                LOGGER.report('The structure was fixed in %.2fs.',
                            label='_clustenm_fix')

                if self._nuc is None:
                    self._idx_cg = self._atoms.ca.getIndices()
                    self._n_cg = self._atoms.ca.numAtoms()
                else:
                    self._idx_cg = self._atoms.select("name CA C2 C4' P").getIndices()
                    self._n_cg = self._atoms.select("name CA C2 C4' P").numAtoms()

                self._n_atoms = self._atoms.numAtoms()
                self._indices = None

            if i == 1:
                self._cg_ensA.setAtoms(self._atoms[self._idx_cg])
            else:
                self._cg_ensB.setAtoms(self._atomsB[self._idx_cg])

    def _fix(self, atoms, i):
        try:
            from pdbfixer import PDBFixer
            from simtk.openmm.app import PDBFile
        except ImportError:
            raise ImportError('Please install PDBFixer and OpenMM in order to use ClustENM and related hybrid methods.')

        stream = createStringIO()
        title = atoms.getTitle()
        writePDBStream(stream, atoms)
        stream.seek(0)
        fixed = PDBFixer(pdbfile=stream)
        stream.close()

        fixed.missingResidues = {}
        fixed.findNonstandardResidues()
        fixed.replaceNonstandardResidues()
        fixed.removeHeterogens(False)
        fixed.findMissingAtoms()
        fixed.addMissingAtoms()
        fixed.addMissingHydrogens(self._ph)

        stream = createStringIO()
        PDBFile.writeFile(fixed.topology, fixed.positions,
                          stream, keepIds=True)
        stream.seek(0)
        if i == 0:
            self._atoms = parsePDBStream(stream)
            self._atoms.setTitle(title)
        else:
            self._atomsB = parsePDBStream(stream)
            self._atomsB.setTitle(title)            
        stream.close()

        self._topology = fixed.topology
        self._positions = fixed.positions

    def getAtomsA(self, selected=True):
        'Returns atoms for structure A (main atoms).'
        return super(AdaptiveHybrid, self).getAtoms(selected)

    def getAtomsB(self, selected=True):
        'Returns atoms for structure B.'
        if self._atomsB is None:
            return None
        if self._indices is None or not selected:
            return self._atomsB
        return self._atomsB[self._indices]

    def getRMSDsB(self):
        if self._confs is None or self._coords is None:
            return None

        indices = self._indices
        if indices is None:
            indices = np.arange(self._confs.shape[1])
        
        weights = self._weights[indices] if self._weights is not None else None

        return calcRMSD(self._atomsB, self._confs[:, indices], weights)

    def getConvergenceRMSDs(self):
        if self._confs is None or self._coords is None:
            return None

        indices = self._indices
        if indices is None:
            indices = np.arange(self._confs.shape[1])
        
        weights = self._weights[indices] if self._weights is not None else None

        n_confs = self.numConfs()
        n_confsA = int(Decimal(n_confs/2).to_integral(rounding=ROUND_HALF_UP))

        confsA = self._confs[:n_confsA]
        if n_confs % 2:
            confsB = self._confs[n_confsA:]
        else:
            confsB = self._confs[n_confsA:]

        RMSDs = []
        for i in range(n_confsA):
            for j in range(2):
                if i + j > n_confsA - 1:
                    break
                RMSDs.append(getRMSD(confsA[i+j], confsB[n_confsA-(i+1)], weights=weights))

        return np.array(RMSDs)

    def _generate(self, confs, **kwargs):

        LOGGER.info('Sampling conformers in generation %d ...' % self._cycle)
        LOGGER.timeit('_clustenm_gen')

        sample_method = self._sample

        if self._parallel:
            with Pool(cpu_count()) as p:
                tmp = p.map(sample_method, [conf for conf in confs])
        else:
            tmp = [sample_method(conf, **kwargs) for conf in confs]

        tmp = [r for r in tmp if r is not None]

        confs_ex = np.concatenate(tmp)

        return confs_ex, [1]

    def setCoordsB(self, coords):
        """Set *coords* as the ensemble reference coordinate set.  *coords*
        may be an array with suitable data type, shape, and dimensionality, or
        an object with :meth:`getCoords` method."""

        atoms = coords
        try:
            if isinstance(coords, Ensemble):
                coords = copy(coords._coords)
            else:
                coords = coords.getCoords()
        except AttributeError:
            pass
        finally:
            if coords is None:
                raise ValueError('coordinates of {0} are not set'
                                 .format(str(atoms)))

        try:
            checkCoords(coords, natoms=self._n_atoms)
        except TypeError:
            raise TypeError('coords must be a numpy array or an object '
                            'with `getCoords` method')

        if coords.shape == self._coords.shape:
            self._coordsB = coords
            self._n_atomsB = coords.shape[0]

            if isinstance(atoms, Ensemble):
                self._indicesB = atoms._indices
                self._atomsB = atoms._atoms
        else:
            raise ValueError('coordsB must have the same shape as main coords')

    def getCoordsB(self, selected=True):
        """Returns a copy of reference coordinates for selected atoms."""

        if self._coordsB is None:
            return None
        if self._indicesB is None or not selected:
            return self._coordsB.copy()
        return self._coordsB[self._indicesB].copy()
