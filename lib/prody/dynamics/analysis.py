# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines functions for calculating atomic properties from normal
modes.""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import time

import numpy as np

from prody import LOGGER
from prody.atomic import AtomGroup
from prody.ensemble import Ensemble, Conformation
from prody.trajectory import TrajBase

from .nma import NMA
from .modeset import ModeSet
from .mode import VectorBase, Mode, Vector
from .gnm import GNMBase

__all__ = ['calcCollectivity', 'calcCovariance', 'calcCrossCorr',
           'calcFractVariance', 'calcSqFlucts', 'calcTempFactors',
           'calcProjection', 'calcCrossProjection', 'calcPerturbResponse', ]

           
def calcCollectivity(mode, masses=None):
    """Return collectivity of the mode.  This function implements collectivity 
    as defined in equation 5 of [BR95]_.  If *masses* are provided, they will 
    be incorporated in the calculation.  Otherwise, atoms are assumed to have 
    uniform masses.
    
    :arg mode: mode or vector
    :type mode: :class:`.Mode` or :class:`.Vector`
    
    :arg masses: atomic masses
    :type masses: :class:`numpy.ndarray`"""
    
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance')
    
    is3d = mode.is3d()
    if masses is not None:
        if len(masses) != mode.numAtoms(): 
            raise ValueError('length of massesmust be equal to number of atoms')
        if is3d:
            u2in = (mode.getArrayNx3() ** 2).sum(1) / masses
    else:
        if is3d:
            u2in = (mode.getArrayNx3() ** 2 ).sum(1)
        else:
            u2in = (mode.getArrayNx3() ** 2 )
    u2in = u2in * (1 / u2in.sum() ** 0.5)
    coll = np.exp(-(u2in * np.log(u2in)).sum()) / mode.numAtoms()
    return coll
    
    
def calcFractVariance(mode):
    """Return fraction of variance explained by the *mode*.  Fraction of 
    variance is the ratio of the variance along a mode to the trace of the 
    covariance matrix of the model."""
    
    if isinstance(mode, Mode):
        var = mode.getVariance()
        trace = mode.getModel()._trace
    elif isinstance(mode, (ModeSet, NMA)):
        var = mode.getVariances()
        if isinstance(mode, ModeSet):
            trace = mode.getModel()._trace
        else:
            trace = mode._trace
    else:
        raise TypeError('mode must be a Mode instance')
    if trace is None:
        raise ValueError('modes are not calculated')
    
    return var / trace


def calcProjection(ensemble, modes, rmsd=True):
    """Return projection of conformational deviations onto given modes.
    For K conformations and M modes, a (K,M) matrix is returned.
    
    :arg ensemble: an ensemble, trajectory or a conformation for which 
        deviation(s) will be projected, or a deformation vector
    :type ensemble: :class:`.Ensemble`, :class:`.Conformation`, 
        :class:`.Vector`, :class:`.Trajectory`
    :arg modes: up to three normal modes
    :type modes: :class:`.Mode`, :class:`.ModeSet`, :class:`.NMA`
    
    By default root-mean-square deviation (RMSD) along the normal mode is 
    calculated. To calculate the projection pass ``rmsd=True``.
    :class:`.Vector` instances are accepted as *ensemble* argument to allow
    for projecting a deformation vector onto normal modes."""
    
    if not isinstance(ensemble, (Ensemble, Conformation, Vector, TrajBase)):
        raise TypeError('ensemble must be Ensemble, Conformation, Vector, '
                        'or a TrajBase, not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, VectorBase)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    if not modes.is3d(): 
        raise ValueError('modes must be 3-dimensional')
    if isinstance(ensemble, Vector):
        n_atoms = ensemble.numAtoms()
    else:
        n_atoms = ensemble.numSelected()
    if n_atoms != modes.numAtoms():
        raise ValueError('number of atoms are not the same')
    if isinstance(ensemble, Vector):
        if not ensemble.is3d(): 
            raise ValueError('ensemble must be a 3d vector instance')
        deviations = ensemble._getArray()
    elif isinstance(ensemble, (Ensemble, Conformation)):
        deviations = ensemble.getDeviations()
    else:
        nfi = ensemble.nextIndex()
        ensemble.goto(0)
        deviations = np.array([frame.getDeviations() for frame in ensemble])
        ensemble.goto(nfi)
    if deviations.ndim == 3:
        deviations = deviations.reshape((deviations.shape[0], 
                                         deviations.shape[1] * 3))
    elif deviations.ndim == 2:
        deviations = deviations.reshape((1, deviations.shape[0] * 3))
    else:
        deviations = deviations.reshape((1, deviations.shape[0]))
    projection = np.dot(deviations, modes._getArray())
    if rmsd:
        projection =  (1 / (n_atoms ** 0.5)) * projection
    return projection


def calcCrossProjection(ensemble, mode1, mode2, scale=None, **kwargs):
    """Return projection of conformational deviations onto modes from
    different models.
    
    :arg ensemble: ensemble for which deviations will be projected
    :type ensemble: :class:`.Ensemble`
    :arg mode1: normal mode to project conformations onto 
    :type mode1: :class:`.Mode`, :class:`.Vector`
    :arg mode2: normal mode to project conformations onto
    :type mode2: :class:`.Mode`, :class:`.Vector`
    :arg scale: scale width of the projection onto mode ``x`` or ``y``,
        best scaling factor will be calculated and printed on the console,
        absolute value of scalar makes the with of two projection same,
        sign of scalar makes the projections yield a positive correlation"""
    
    if not isinstance(ensemble, (Ensemble, Conformation, Vector, TrajBase)):
        raise TypeError('ensemble must be Ensemble, Conformation, Vector, '
                        'or a Trajectory, not {0:s}'.format(type(ensemble)))
    if not isinstance(mode1, VectorBase):
        raise TypeError('mode1 must be a Mode instance, not {0:s}'
                        .format(type(mode1)))
    if not mode1.is3d():
        raise ValueError('mode1 must be 3-dimensional')
    if not isinstance(mode2, VectorBase):
        raise TypeError('mode2 must be a Mode instance, not {0:s}'
                        .format(type(mode2)))
    if not mode2.is3d():
        raise ValueError('mode2 must be 3-dimensional')
    
    if scale is not None:
        assert isinstance(scale, str), 'scale must be a string'
        scale = scale.lower()
        assert scale in ('x', 'y'), 'scale must be x or y'

    xcoords = calcProjection(ensemble, mode1, kwargs.get('rmsd', True))
    ycoords = calcProjection(ensemble, mode2, kwargs.pop('rmsd', True))
    if scale:
        scalar = kwargs.get('scalar', None)
        if scalar:
            assert isinstance(scalar, (float, int)), 'scalar must be a number'
        else:
            scalar = ((ycoords.max() - ycoords.min()) / 
                      (xcoords.max() - xcoords.min())
                      ) * np.sign(np.dot(xcoords, ycoords))
            if scale == 'x':
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode1, scalar))
            else:
                scalar = 1 / scalar
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode2, scalar))

        if scale == 'x':
            xcoords = xcoords * scalar  
        else:
            ycoords = ycoords * scalar

    return xcoords, ycoords


def calcSqFlucts(modes):
    """Return sum of square-fluctuations for given set of normal *modes*. 
    Square fluctuations for a single mode is obtained by multiplying the 
    square of the mode array with the variance (:meth:`~.Mode.getVariance`) 
    along the mode.  For :class:`.PCA` and :class:`.EDA` models built using
    coordinate data in Å, unit of square-fluctuations is |A2|, for 
    :class:`.ANM` and :class:`.GNM`, on the other hand, it is arbitrary or 
    relative units."""
    
    if not isinstance(modes, (VectorBase, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    is3d = modes.is3d()
    if isinstance(modes, Vector):
        if is3d:
            return (modes._getArrayNx3()**2).sum(axis=1)
        else:
            return (modes._getArray() ** 2)
    else:
        sq_flucts = np.zeros(modes.numAtoms()) 
        if isinstance(modes, VectorBase):
            modes = [modes]
        for mode in modes:
            if is3d:
                sq_flucts += ((mode._getArrayNx3()**2).sum(axis=1) * 
                                mode.getVariance())
            else:
                sq_flucts += (mode._getArray() ** 2)  * mode.getVariance()
        return sq_flucts


def calcCrossCorr(modes, n_cpu=1):
    """Return cross-correlations matrix.  For a 3-d model, cross-correlations 
    matrix is an NxN matrix, where N is the number of atoms.  Each element of 
    this matrix is the trace of the submatrix corresponding to a pair of atoms.
    Covariance matrix may be calculated using all modes or a subset of modes
    of an NMA instance.  For large systems, calculation of cross-correlations 
    matrix may be time consuming.  Optionally, multiple processors may be 
    employed to perform calculations by passing ``n_cpu=2`` or more."""
    
    if not isinstance(n_cpu, int):
        raise TypeError('n_cpu must be an integer')
    elif n_cpu < 1:
        raise ValueError('n_cpu must be equal to or greater than 1')
        
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
        
    if modes.is3d():
        model = modes
        if isinstance(modes, (Mode, ModeSet)):
            model = modes._model
            if isinstance(modes, (Mode)):
                indices = [modes.getIndex()]
                n_modes = 1
            else:
                indices = modes.getIndices()
                n_modes = len(modes)
        else:
            n_modes = len(modes)
            indices = np.arange(n_modes)
        array = model._array
        n_atoms = model._n_atoms
        variances = model._vars
        if n_cpu == 1:
            arvar = (array[:, indices]*variances[indices]).T.reshape(
                                                        (n_modes, n_atoms, 3))
            array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
            covariance = np.tensordot(array.transpose(2, 0, 1),
                                      arvar.transpose(0, 2, 1),
                                      axes=([0, 1], [1, 0]))
        else:
            import multiprocessing
            n_cpu = min(multiprocessing.cpu_count(), n_cpu)
            queue = multiprocessing.Queue()
            size = n_modes / n_cpu
            for i in range(n_cpu):
                if n_cpu - i == 1:
                    indices = modes.indices[i*size:]
                else:
                    indices = modes.indices[i*size:(i+1)*size]
                process = multiprocessing.Process(target=_crossCorrelations, 
                              args=(queue, n_atoms, array, variances, indices))
                process.start()
            while queue.qsize() < n_cpu:
                time.sleep(0.05)
            covariance = queue.get()
            while queue.qsize() > 0:
                covariance += queue.get()
    else:
        covariance = calcCovariance(modes)
    diag = np.power(covariance.diagonal(), 0.5)
    return covariance / np.outer(diag, diag)

def _crossCorrelations(queue, n_atoms, array, variances, indices):
    """Calculate covariance-matrix for a subset of modes."""
    
    n_modes = len(indices)
    arvar = (array[:, indices] * variances[indices]).T.reshape((n_modes,
                                                                n_atoms, 3))
    array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
    covariance = np.tensordot(array.transpose(2, 0, 1),
                              arvar.transpose(0, 2, 1),
                              axes=([0, 1], [1, 0]))
    queue.put(covariance)
    
    
def calcTempFactors(modes, atoms):
    """Return temperature (β) factors calculated using *modes* from a 
    :class:`.ANM` or :class:`.GNM` instance scaled according to the 
    experimental β-factors from *atoms*."""
    
    model = modes.getModel()
    if not isinstance(model, GNMBase):
        raise TypeError('modes must come from GNM or ANM')
    if model.numAtoms() != atoms.numAtoms():
        raise ValueError('modes and atoms must have same number of nodes')
    sqf = calcSqFlucts(modes)
    return sqf / ((sqf**2).sum()**0.5) * (atoms.getBetas()**2).sum()**0.5


def calcCovariance(modes):
    """Return covariance matrix calculated for given *modes*."""
    
    if isinstance(modes, Mode):
        array = modes._getArray()
        return np.outer(array, array) * modes.getVariance()
    elif isinstance(modes, ModeSet):
        array = modes._getArray()
        return np.dot(array, np.dot(np.diag(modes.getVariances()), array.T))
    elif isinstance(modes, NMA):
        return modes.getCovariance()
    else:
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance')


def calcPerturbResponse(model, atoms=None, repeats=100):
    """Return a matrix of profiles from scanning of the response of the 
    structure to random perturbations at specific atom (or node) positions. 
    The function implements the perturbation response scanning (PRS) method 
    described in [CA09]_.  Rows of the matrix are the average magnitude of the 
    responses obtained by perturbing the atom/node position at that row index, 
    i.e. ``prs_profile[i,j]`` will give the response of residue/node *j* to 
    perturbations in residue/node *i*.  PRS is performed using the covariance 
    matrix from *model*, e.t. :class:`.ANM` instance.  Each residue/node is 
    perturbed *repeats* times with a random unit force vector.  When *atoms* 
    instance is given, PRS profile for residues will be added as an attribute 
    which then can be retrieved as ``atoms.getData('prs_profile')``.  *model* 
    and *atoms* must have the same number of atoms. *atoms* must be an :class:`
    ~.AtomGroup` instance.
    
    The RPS matrix can be save as follows::
        
      prs_matrix = calcPerturbationResponse(p38_anm)
      writeArray('prs_matrix.txt', prs_matrix, format='%8.6f', delimiter='\t')
    """
    
    if not isinstance(model, NMA): 
        raise TypeError('model must be an NMA instance')
    elif not model.is3d():
        raise TypeError('model must be a 3-dimensional NMA instance')
    elif len(model) == 0:
        raise ValueError('model must have normal modes calculated')
    if atoms is not None:
        if not isinstance(atoms, AtomGroup):
            raise TypeError('atoms must be an AtomGroup instance')
        elif atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')
            
    assert isinstance(repeats, int), 'repeats must be an integer'
    cov = calcCovariance(model)
    if cov is None:
        raise ValueError('model did not return a covariance matrix')
    
    n_atoms = model.numAtoms()
    response_matrix = np.zeros((n_atoms, n_atoms))
    LOGGER.progress('Calculating perturbation response', n_atoms, '_prody_prs')
    i3 = -3
    i3p3 = 0
    for i in range(n_atoms):
        i3 += 3
        i3p3 += 3
        forces = np.random.rand(repeats * 3).reshape((repeats, 3))
        forces /= ((forces**2).sum(1)**0.5).reshape((repeats, 1))
        for force in forces:
            response_matrix[i] += (np.dot(cov[:, i3:i3p3], force) ** 2
                                            ).reshape((n_atoms, 3)).sum(1)
        LOGGER.update(i, '_prody_prs')

    response_matrix /= repeats
    LOGGER.clear()
    LOGGER.report('Perturbation response scanning completed in %.1fs.', 
                  '_prody_prs')
    if atoms is not None:
        atoms.setData('prs_profile', response_matrix)
    return response_matrix
    
    # save the original PRS matrix
    np.savetxt('orig_PRS_matrix', response_matrix, delimiter='\t', fmt='%8.6f')
    # calculate the normalized PRS matrix
    self_dp = np.diag(response_matrix) # using self displacement (diagonal of
                               # the original matrix) as a
                               # normalization factor     
    self_dp = self_dp.reshape(n_atoms, 1)
    norm_PRS_mat = response_matrix / np.repeat(self_dp, n_atoms, axis=1)
    # suppress the diagonal (self displacement) to facilitate
    # visualizing the response profile
    norm_PRS_mat = norm_PRS_mat - np.diag(np.diag(norm_PRS_mat))
    np.savetxt('norm_PRS_matrix', norm_PRS_mat, delimiter='\t', fmt='%8.6f')
    return response_matrix
