# -*- coding: utf-8 -*-
"""This module defines functions for calculating entropy transfer from normal
modes."""

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


__all__ = ['calcEntropyTransfer', 'calcOverallNetEntropyTransfer']

def calcEntropyTransfer(model, ind1, ind2, tau):
    """This function calculates the entropy transfer from residue indice 
    ind1 to ind2 for a given time constant tau based on GNM.  
    """
    if not isinstance(model, NMA):
        raise TypeError('model must be a NMA instance')
    elif model.is3d():
        raise TypeError('model must be a 1-dimensional NMA instance')

    linalg = importLA()
    n_atoms = model.numAtoms()
    n_modes = model.numModes()

    eigvecs = model.getEigvecs().T
    eigvals = model.getEigvals()

    tau_0 = 1
    T = 0
    dummy1 = 0
    dummy2 = 0
    for k in range(n_modes):
        dummy1 += 1.0 / eigvals[k] * eigvecs[k,ind2] * eigvecs[k,ind2]
        dummy2 += 1.0 / eigvals[k] * eigvecs[k,ind2] * eigvecs[k,ind2] * np.exp(-eigvals[k]*tau/tau_0)
    T += 0.5 * np.log(dummy1**2 - dummy2**2)

    dummy1 = 0
    dummy2 = 0
    dummy3 = 0 
    dummy4 = 0 
    dummy5 = 0 
    dummy6 = 0 
    dummy7 = 0 
    dummy8 = 0 
    dummy9 = 0 
    dummy10 = 0 
    for k in range(n_modes):
        dummy1 += 1.0 / eigvals[k] * eigvecs[k,ind1] * eigvecs[k,ind1]
        dummy2 += 1.0 / eigvals[k] * eigvecs[k,ind2] * eigvecs[k,ind2]
        dummy3 += 1.0 / eigvals[k] * eigvecs[k,ind1] * eigvecs[k,ind2]
        dummy4 += 1.0 / eigvals[k] * eigvecs[k,ind2] * eigvecs[k,ind2] * np.exp(-eigvals[k]*tau/tau_0)
        dummy5 += 1.0 / eigvals[k] * eigvecs[k,ind1] * eigvecs[k,ind2] * np.exp(-eigvals[k]*tau/tau_0)  
        dummy9 += 1.0 / eigvals[k] * eigvecs[k,ind2] * eigvecs[k,ind2] * np.exp(-eigvals[k]*tau/tau_0)       
        
    dummy6 = dummy5
    dummy7 = dummy3 
    dummy8 = dummy2
    dummy10 = dummy1

    T -= 0.5 * np.log(dummy1*dummy2**2+2*dummy3*dummy4*dummy5-(dummy6**2+dummy7**2)*dummy8-dummy9**2*dummy10)
    T -= 0.5 * np.log(dummy2)
    T += 0.5 * np.log(dummy1*dummy2-dummy3**2)
    return T

def calcAllEntropyTransfer(model, tau):
    """This function calculates the net entropy transfer for a whole structure 
    with a given time constant tau based on GNM.  
    """
    if not isinstance(model, NMA):
        raise TypeError('model must be a NMA instance')
    elif model.is3d():
        raise TypeError('model must be a 1-dimensional NMA instance')

    linalg = importLA()
    n_atoms = model.numAtoms()
    n_modes = model.numModes()

    entropyTransfer = np.zeros((n_atoms,n_atoms))
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i != j:
                entropyTransfer[i,j]=calcEntropyTransfer(model,i,j,tau)

    return entropyTransfer

def calcNetEntropyTransfer(entropyTransfer):
    n_atoms = entropyTransfer.shape[0]

    netEntropyTransfer = np.zeros((n_atoms,n_atoms))
    for i in range(n_atoms):
        for j in range(n_atoms):
            netEntropyTransfer[i,j]=entropyTransfer[i,j]-entropyTransfer[j,i]

    return netEntropyTransfer

def calcOverallNetEntropyTransfer(model, turbo=False):
    """This function calculates the net entropy transfer for a whole structure 
    with a given time constant tau based on GNM.  
    """
    if not isinstance(model, NMA):
        raise TypeError('model must be a NMA instance')
    elif model.is3d():
        raise TypeError('model must be a 1-dimensional NMA instance')

    linalg = importLA()
    n_atoms = model.numAtoms()
    n_modes = model.numModes()

    tau_max = 5.0 
    tau_step = 0.1
    taus = np.arange(start=tau_step, stop=tau_max+1e-6, step=tau_step)
    taus = np.insert(taus,0,0.000001)
    numTaus = len(taus)
    netEntropyTransfer = np.zeros((numTaus,n_atoms,n_atoms))
    if turbo:
        try:
            from joblib import Parallel, delayed
            import multiprocessing as mp
        except: 
            LOGGER.report('joblib and multiprocessing is not imported. Running' + 
                'with sequential execution.')
        LOGGER.timeit('_ent_trans')
        n_cpu = mp.cpu_count()
        netEntropyTransfer = Parallel(n_jobs=n_cpu)(delayed(calcAllEntropyTransfer)(model,taus[i]) \
            for i in range(numTaus))
        netEntropyTransfer = np.asarray(netEntropyTransfer)
    
    else:
        LOGGER.timeit('_ent_trans')
        for i in range(len(taus)):
            netEntropyTransfer[i,:,:] = calcAllEntropyTransfer(model,taus[i])

    LOGGER.report('Net Entropy Transfer calculation is completed in %.1fs.',
                  '_ent_trans')

    overallNetEntropyTransfer = np.zeros((n_atoms,n_atoms))

    LOGGER.timeit('_num_int')
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i != j:
                overallNetEntropyTransfer[i,j] = np.trapz(netEntropyTransfer[:,i,j],taus)
    LOGGER.report('Numerical integration is completed in %.1fs.',
                  '_num_int')

    return overallNetEntropyTransfer

def test():
    from prody import parsePDB, GNM
    from prody.dynamics.analysis import calcOverallNetEntropyTransfer
    import matplotlib.pyplot as plt

    pdb = parsePDB('1z83', subset='ca', chain='A')
    gnm = GNM()
    gnm.buildKirchhoff(pdb, cutoff=7.0)
    gnm.calcModes(n_modes=None)
    entTransfer = calcOverallNetEntropyTransfer(gnm,turbo=True)

    # f = open('/data/Manuscript_data/Data/1Z83A/monomer_overallnet_A_cihan2.txt','w')
    # for i in range(gnm.numAtoms()):
    #     for j in range(gnm.numAtoms()):
    #         if i != j:
    #             f.write('%d\t%d\t%f\n' % (i+1,j+1,entTransfer[i,j]))
    # f.close()
    return entTransfer
