# -*- coding: utf-8 -*-
"""This module defines functions for performing perturbation response scanning
from PCA and normal modes."""

import time

import numpy as np

from prody import LOGGER
from prody.proteins import parsePDB
from prody.atomic import AtomGroup, Selection
from prody.ensemble import Ensemble, Conformation
from prody.trajectory import TrajBase
from prody.utilities import importLA
from numpy import sqrt, arange, log, polyfit, array

from .nma import NMA
from .modeset import ModeSet
from .mode import VectorBase, Mode, Vector
from .gnm import GNMBase
from .analysis import calcCovariance

__all__ = ['calcPerturbResponse']

def calcPerturbResponse(model, atoms=None, **kwargs):

    """Returns a matrix of profiles from scanning the response of the
    structure to random perturbations at specific atom (or node) positions.
    The function implements the perturbation response scanning (PRS) method
    described in [CA09]_.  Rows of the matrix are the average magnitude of the
    responses obtained by perturbing the atom/node position at that row index,
    i.e. ``prs_profile[i,j]`` will give the response of residue/node *j* to
    perturbations in residue/node *i*.  PRS is performed using the covariance
    matrix from *model*, e.g. :class:`.ANM` instance.

    When an *atoms* instance is given, the PRS matrix will be added as data, 
    which can be retrieved with ``atoms.getData('prs_matrix')``.  

    *model* and *atoms* must have the same number of atoms. *atoms* must be an
    :class:`.AtomGroup` instance. 

    .. [CA09] Atilgan C, Atilgan AR, Perturbation-Response Scanning
       Reveals Ligand Entry-Exit Mechanisms of Ferric Binding Protein.
       *PLoS Comput Biol* **2009** 5(10):e1000544.

    """

    if not isinstance(model, (NMA, ModeSet, Mode)):
        raise TypeError('model must be an NMA, ModeSet, or Mode instance')

    if isinstance(model, NMA) and len(model) == 0:
        raise ValueError('model must have normal modes calculated')

    atoms = kwargs.get('atoms',None)
    if atoms is not None:
        if isinstance(atoms, Selection):
            atoms = atoms.copy()
        if not isinstance(atoms, AtomGroup):
            raise TypeError('atoms must be an AtomGroup instance')
        elif atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')

    n_atoms = model.numAtoms()
    LOGGER.timeit('_prody_prs_all')
    LOGGER.info('Calculating covariance matrix')
    LOGGER.timeit('_prody_cov')

    cov = model.getCovariance()

    LOGGER.clear()
    LOGGER.report('Covariance matrix calculated in %.1fs.', '_prody_cov')

    LOGGER.info('Calculating perturbation response')
    LOGGER.timeit('_prody_prs_mat')
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

    LOGGER.clear()
    LOGGER.report('Perturbation response matrix calculated in %.1fs.',
                  '_prody_prs_mat')

    no_diag = kwargs.get('no_diag', False)
    #filename = kwargs.get('filename', None)

    norm_prs_matrix = np.zeros((n_atoms, n_atoms))
    self_dp = np.diag(prs_matrix)  
    self_dp = self_dp.reshape(n_atoms, 1)
    norm_prs_matrix = prs_matrix / np.repeat(self_dp, n_atoms, axis=1)

    effectiveness = np.mean(norm_prs_matrix, axis=1)
    sensitivity = np.mean(norm_prs_matrix, axis=0)

    if no_diag:
       # suppress the diagonal (self displacement) to facilitate
       # visualizing the response profile
       norm_prs_matrix = norm_prs_matrix - np.diag(np.diag(norm_prs_matrix))

    #if filename:
    #    np.savetxt(filename, norm_prs_matrix, delimiter='\t', fmt='%8.6f')

    LOGGER.report('Perturbation response scanning completed in %.1fs.',
                  '_prody_prs_all')

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

    return norm_prs_matrix, effectiveness, sensitivity

