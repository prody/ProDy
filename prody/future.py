# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

"""This module contains features that are under development.

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import time

sparse = None
import numpy as np
import prody
from prody import ProDyLogger as LOGGER
from collections import defaultdict

__all__ = ['MarkovModel']

class MarkovModel(object):
    
    """Implementation of the Markov model described in [CC06]_ and [CC07]_.

    .. versionadded:: 0.7.1
    
    """
    
    def __init__(self, name):
        name = str(name)
        if name == '':
            name = 'Unnamed'
        self._name = name 
        
        self._affinity = None
        self._stationary = None
        self._transition = None
        self._n_atoms = 0
        self._n_residues = 0
        self._n_levels = 0
    
    def __repr__(self):
        return '<MarkovModel: {0:s} ({1:d} residues, {2:d} levels)>'.format(
                self._name, self._n_residues, self._n_levels)

    def __str__(self):
        return 'MarkovModel {0:s}'.format(self._name)
    
    def getNumOfHierLevels(self):
        """Return number of levels in the hierarchical coarse graining."""
        
        return self._n_levels
    
    def getNumOfAtoms(self):
        """Return the number of atoms in the system for which model was built.
        
        """
        
        return self._n_atoms
    
    def getNumOfResidues(self):
        """Return the number of residues in the system for which model was 
        built.
        
        """
        
        return self._n_residues
    
    
    def buildAffinityMatrix(self, atoms, cutoff=4):
        """Build the affinity matrix for given *atoms*.

        :arg atoms: atoms for which the affinity matrix will be calculated
        :type atoms: :class:`~prody.atomic.Atomic`
        
        :arg cutoff: pairwise atomic contact cutoff distance, default is 4 Å
        :type cutoff: float

        """
        
        if not isinstance(atoms, prody.Atomic):
            raise TypeError('atoms must be an Atomic instance, '
                            '{0:s} is not valid'.format(type(atoms)))
        cutoff = float(cutoff)
        assert cutoff > 0, 'cutoff distance must be greater than 0'
        
        if prody.dynamics.KDTree is None:
            prody.importBioKDTree()

        start = time.time()

        n_atoms = atoms.getNumOfAtoms()
        hv = prody.HierView(atoms)
        n_res = hv.getNumOfResidues()

        rids = np.zeros(n_atoms, np.int64)
        rlen = np.zeros(n_res)
        for i, res in enumerate(hv.iterResidues()):
            rids[ res.getIndices() ] = i
            rlen[ i ] = len(res)
        
        affinity = defaultdict(int)
        kdtree = prody.dynamics.KDTree(3)
        kdtree.set_coords(atoms.getCoordinates())
        kdtree.all_search(cutoff)
        LOGGER.debug('KDTree was built in {0:.2f}s.'
                     .format(time.time()-start))
        start = time.time()
        for i, j in kdtree.all_get_indices():
            i = rids[i] 
            j = rids[j]
            if i == j:
                affinity[(i,j)] += 0.5
                #affinity[i, j] += 0.5
            else:
                affinity[(i,j)] += 1
                #affinity[i, j] += 1
        
        #affinity += affinity.T  
        rlen = 1 / (rlen ** 0.5)
        length = len(affinity)
        i = np.zeros(length, np.int32) 
        j = np.zeros(length, np.int32)
        v = np.zeros(length, np.float64)
        k = 0
        for key, value in affinity.iteritems():
            i[k] = key[0]
            j[k] = key[1]
            v[k] = value
            k += 1
            
        v = v * rlen[i] * rlen[j]  
        from scipy import sparse
        affinity = sparse.coo_matrix((v, (i,j)), shape=(n_res, n_res))
        self._affinity = affinity + affinity.T 
        self._stationary = None
        self._n_atoms = n_atoms
        self._n_residues = n_res
        
        LOGGER.debug('Affinity matrix was built in {0:.2f}s.'
                     .format(time.time()-start))
    
    def getAffinityMatrix(self):
        """Return a copy of the affinity matrix.
        
        Affinity matrix returned as a :class:`scipy.sparse.csr.csr_matrix`.
        
        """
        
        if self._affinity is not None:
            return self._affinity.copy()

    def getStationaryDistribution(self):
        """Return a copy of the stationary distribution."""
        
        if self._affinity is None:
            return None
        elif self._stationary is None:  
            d = self._affinity.sum(0)
            self._stationary = d * (1 / d.sum())
        return self._stationary.copy()
        


    def buildTransitionMatrix(self):
        """Build Markov transition matrix from the affinity matrix."""
        
        from scipy import sparse
        start = time.time()
        self._transition = np.dot(self._affinity, 
                  sparse.lil_diags([1 / self._affinity.sum(0)], [0], [self._n_residues,self._n_residues]))
        LOGGER.debug('Markov transition matrix was built in {0:.2f}s.'
                     .format(time.time()-start))

    
    def getTransitionMatrix(self):
        """Return a copy of the Markov transition matrix."""
        
        if self._transition is not None:
            return self._transition.copy()
    
    def grainCoarser(self, beta=10):
        
        assert isinstance(beta, int), 'beta must be an integer'
        assert beta > 0, 'beta must be positive'
        
        transition = np.linalg.matrix_power(self._transition, beta) 
        torf = np.ones(self._n_residues, np.bool)
        stationary = self.getStationaryDistribution()
        kernels = []
        whiches = [] 
        #which = (stationary == stationary.max()).nonzeros()[0][0]        
        #kernels = [transition[:, which]]
        #torf[which] = False
        #k = kernels[-1]
        #torf[ k > k.max() / 2 ] = False
        while np.any(torf):
            stationary = stationary * torf
            which = (stationary == stationary.max()).nonzero()[0][0]
            whiches.append(which)
            kernels.append(transition[:, which])
            torf[which] = False
            k = kernels[-1]
            torf[ k > k.max() / 2 ] = False
        return whiches
        temp = zip(self.getStationaryDistribution(), 
                   np.arange(self._n_residues))
        temp.sort(reverse=True)
        kernels = []
        while temp:        
            kernels.append(transition[:,temp[0][0]])
            
        delta = np.ones()
        
if __name__ == '__main__':
    
     #p = parsePDB('1aon')
     m = MarkovModel('')
     m.buildAffinityMatrix(p.select('protein'))
     #m.buildTransitionMatrix()
     #c = m.grainCoarser()
