# -*- coding: utf-8 -*-
"""This module defines functions for performing perturbation response scanning
from PCA and normal modes."""

import time

import numpy as np

from prody import LOGGER
from prody.atomic import AtomGroup, Selection, Atomic, sliceAtomicData
from prody.utilities import div0

from .nma import NMA
from .modeset import ModeSet
from .mode import VectorBase, Mode, Vector
from .gnm import GNMBase
from .analysis import calcCovariance

__all__ = ['calcPerturbResponse']

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

    # LOGGER.clear()
    # LOGGER.report('Covariance matrix calculated in %.1fs.', '_prody_cov')

    # LOGGER.info('Calculating perturbation response')
    # LOGGER.timeit('_prody_prs_mat')
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

    # LOGGER.clear()
    # LOGGER.report('Perturbation response matrix calculated in %.1fs.',
    #               '_prody_prs_mat')

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

    return norm_prs_matrix, effectiveness, sensitivity


def calcDynamicFlexibilityIndex(prs_matrix, atoms, select):
    """
    Calculate the dynamic flexibility index for the selected residue(s).
    This function implements the dynamic flexibility index (Dfi) method
    described in [ZNG13]_.

    :arg prs_matrix: a matrix from PRS
    :type prs_matrix: list, tuple, :class:`~numpy.ndarray`

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg select: a selection string or selection for residues of interest
    :type select: str, :class:`.Selection`

    .. [ZNG13] Gerek ZN, Kumar S, Ozkan SB, Structural dynamics flexibility 
       informs function and evolution at a proteome scale.
       *Evol Appl.* **2013** 6(3):423-33.

    """
    profiles = sliceAtomicData(prs_matrix, atoms, select, axis=0)
    return np.sum(profiles, axis=1)/np.sum(prs_matrix)
