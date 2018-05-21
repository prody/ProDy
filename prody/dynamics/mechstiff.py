# -*- coding: utf-8 -*-
"""This module defines functions for Mechanical Stiffness calculations."""

from numbers import Integral

import numpy as np

from prody import LOGGER
from prody.utilities import checkCoords


__all__ = ['calcMechStiff', 'calcStiffnessRange', 'calcMechStiffStatistic', 
           'calcStiffnessRangeSel']

def calcMechStiff(modes, coords, kbt=1.):
    """Calculate stiffness matrix calculated using :class:`.ANM` instance. 
    Method described in [EB08]_. 

    :arg coords: a coordinate set or an object with ``getCoords`` method
    :type coords: :class:`numpy.ndarray`.

    :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
        If **None** is given, all modes will be calculated (3x number of atoms).
    :type n_modes: int or **None**, default is 20.
    
    Author: Mustafa Tekpinar & Karolina Mikulska-Ruminska & Cihan Kaya
    """

    try:
        coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                    coords.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be a Numpy array or an object '
                            'with `getCoords` method')
    try:
        is3d = modes.is3d()
        eigvecs = modes.getArray().T.flatten()
        eigvals = modes.getEigvals()
    except:
        raise TypeError('modes must be either an NMA or ModeSet object')

    if not is3d:
        raise TypeError('modes must be 3-dimensional')

    n_atoms = modes.numAtoms()
    n_modes = modes.numModes()
    
    LOGGER.timeit('_sm')

    sm = np.zeros((n_atoms, n_atoms), np.double)
    from .smtools import calcSM
    LOGGER.info('Calculating stiffness matrix.')

    calcSM(coords, sm, eigvecs, eigvals,
            n_atoms, n_modes, float(kbt))

    LOGGER.report('Stiffness matrix calculated in %.2lfs.', label='_sm')
    
    LOGGER.info('The range of effective force constant is: {0} to {1}.'
                                .format(*calcStiffnessRange(sm)))

    return sm


def calcStiffnessRange(stiffness):
    """ Return the range of effective spring constant."""
    
    return np.min(stiffness[np.nonzero(stiffness)]), np.amax(stiffness)

def calcMechStiffStatistic(stiffness, rangeK, minAA=0, AA='all'):
    """Returns number of effective spring constant with set range of
    amino acids of protein structure.
    ``AA`` can be a list with a range of analysed amino acids as:
    [first_aa, last_aa, first_aa2, last_aa2],
    minAA - eliminate amino acids that are within 20aa and
    ``rangeK`` is a list [minK, maxK]"""
    
    if AA == 'all':
        sm = stiffness
    elif isinstance(AA, Integral): 
        sm = stiffness[0: AA, (-1)*AA-1:-1]
    elif not np.isscalar(AA) and len(AA) == 1:
        sm = stiffness[0: AA, (-1)*AA-1:-1]
    elif not np.isscalar(AA) and len(AA) == 4:
        sm = stiffness[AA[0]:AA[1],AA[2]:AA[3]]
    if minAA > 0:
        sm2 = sm[minAA:-1,0:-1-minAA]  # matrix without close contacts
        sm3 = np.tril(sm2, k=-1)
        #sort_sm2 = np.sort((np.tril(sm2, k=-1)1).flatten())
        a = np.where(np.logical_and(sm3>rangeK[0], sm3<rangeK[1]))
    if minAA == 0:
        sm2 = np.tril(sm, k=-1)
        a = np.where(np.logical_and(sm2>rangeK[0], sm2<rangeK[1]))
    return len(a[0])
    

def calcStiffnessRangeSel(stiffness, value, minAA=20, AA='all'):
    """Returns minimum or maximum value of sping constant from 
    mechanical stiffness calculations for residues that are within 
    more than ``min_aa`` from each other. ``Value`` should be 'minK'
    or 'maxK'. It alow to avoid residues near each other. 
    ``AA`` is a number of residues from both terminus (N and C)
    of protein strcuture, it can be ``all`` or int value (than first 
    and last ``AA`` residues will be analysed. 
    With ``minAA=0`` it can be used to search the highest/lowest 
    values of interactions between N-C terminus if protein structure
    has a shear, zipper or SD1-disconnected mechanical clamp 
    -it is common in FnIII/Ig like domains and determines the maximum 
    unfolding force in AFM or SMD method."""
    
    if AA == 'all':
        sm = stiffness
    elif isinstance(AA, Integral):
        sm = stiffness[0: AA, (-1)*AA-1:-1]
    minK = np.min(sm[np.nonzero(sm)]) 
    maxK = np.amax(sm)
    
    if value == 'minK':
        indices = np.where(sm == minK)
    elif value == 'maxK':
        indices = np.where(sm == maxK)
    try:
        residue_diff = abs(indices[0][0]-indices[0][1])
    except: 
        residue_diff = abs(indices[0]-indices[1])
    sm_mod = sm.copy()
    checking = True
    while checking:
        if residue_diff < minAA:
            sm_mod[indices[0][0],indices[0][1]] = 0
            sm_mod[indices[1][0],indices[1][1]] = 0
            if value == 'minK':
                mK = np.min(sm_mod[np.nonzero(sm_mod)])
                indices = np.where(sm_mod == np.min(sm_mod[np.nonzero(sm_mod)]))
            elif value == 'maxK':
                mK = np.amax(sm_mod)
                indices = np.where(sm_mod == np.amax(sm_mod))
            try:
                residue_diff = abs(indices[0][0]-indices[0][1])
            except: residue_diff = abs(indices[0]-indices[1])
        else: 
            if value == 'minK':
                mK = np.min(sm_mod[np.nonzero(sm_mod)])
                indices = np.where(sm_mod == np.min(sm_mod[np.nonzero(sm_mod)]))
            elif value == 'maxK':
                mK = np.amax(sm_mod)
                indices = np.where(sm_mod == np.amax(sm_mod))
            checking = False

    if len(indices[0]) == 2:
        return mK, list(indices[0])
    return mK, list(indices[0])+list(indices[1])
    