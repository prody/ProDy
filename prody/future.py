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

"""This module contains features that are under development.

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import time

from scipy import sparse
import numpy as np
import prody
LOGGER = prody.LOGGER
from collections import defaultdict

__all__ = ['MarkovModel']

def parseDCD(filename, indices=None, header=False, first=0, skip=0, last=None):
    """Returns coordinate sets for atoms at *indices* from specified frames."""
    pass

class MarkovModel(object):
    
    """Implementation of the Markov model described in [CC06]_ and [CC07]_.

    .. versionadded:: 0.7.1
    
    """
    
    def __init__(self, title):
        title = str(title)
        if title == '':
            title = 'Unknown'
        self._title = title 
        
        self._affinity = None
        self._stationary = None
        self._transition = None
        self._n_nodes = 0
        self._level = 0
        self._lower = None
        self._higher = None
    
    def __repr__(self):
        return '<MarkovModel: {0:s} ({1:d} nodes at level {2:d})>'.format(
                self._title, self._n_nodes, self._level)

    def __str__(self):
        return 'MarkovModel {0:s}'.format(self._title)
    
    def getTitle(self):
        """Return the title of the model."""
        
        return self._title
    
    def getLevelNumber(self):
        """Return the level number."""
        
        return self._level
    
    def getLowerLevel(self):
        """Return the model at the lower level in the hierarchy."""
        
        return self._lower
    
    def getHigherLevel(self):
        """Return the model at the higher level in the hierarchy."""
        
        return self._higher

    def getNumOfNodes(self):
        """Return the number of nodes in the level."""
        
        return self._n_nodes
    
    
    def buildAffinityMatrix(self, atoms, cutoff=4):
        """Build the affinity matrix for given *atoms*.

        Note that if you do not want to incorporate hydrogen and non-protein
        atoms in calculations, make the selection ``"noh and protein"``.
        
    
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
        if not isinstance(atoms, prody.AtomGroup):
            atoms = atoms.getAtomGroup().copy(atoms)
        n_atoms = atoms.numAtoms()
        hv = prody.HierView(atoms)
        n_res = hv.numResidues()

        rids = np.zeros(n_atoms, int) # residue indices of atoms
        rlen = np.zeros(n_res) # residue lengths
        resmap = {} # used for symmetry purposes
        for i, res in enumerate(hv.iterResidues()):
            rids[ res.getIndices() ] = i
            rlen[ i ] = len(res)
            res = (res.getChid(), res.getNumber(), res.getIcode())
            resmap[i] = res
            resmap[res] = i
        self._resmap = resmap
        LOGGER.debug('Atoms were evalulated in {0:.2f}s.'
                     .format(time.time()-start))

        start = time.time()
        kdtree = prody.dynamics.KDTree(3)
        kdtree.set_coords(atoms.getCoords())
        kdtree.all_search(cutoff)
        LOGGER.debug('KDTree was built in {0:.2f}s.'
                     .format(time.time()-start))

        start = time.time()
        affinity = defaultdict(int)
        for i, j in kdtree.all_get_indices():
            i = rids[i] 
            j = rids[j]
            if i == j:
                affinity[(i,j)] += 0.5
            else:
                affinity[(i,j)] += 1
        
        length = len(affinity)
        i = np.zeros(length, int) 
        j = np.zeros(length, int)
        v = np.zeros(length, float)
        k = 0
        for key, value in affinity.iteritems():
            i[k] = key[0]
            j[k] = key[1]
            v[k] = value
            k += 1
        rlen = rlen ** -0.5
        # = Nij * (1/Ni^0.5) * (1/Nj^0.5)    
        v = v * rlen[i] * rlen[j]  
        affinity = sparse.coo_matrix((v, (i,j)), shape=(n_res, n_res))
        self._affinity = affinity + affinity.T 
        LOGGER.debug('Affinity matrix was built in {0:.2f}s.'
                     .format(time.time()-start))
        self._stationary = None
        self._n_nodes = n_res
    
    def setAffinityMatrix(self, matrix):
        """Set the affinity matrix."""
        
        shape = matrix.shape 
        if shape[0] != shape[1]:# or not np.all(matrix == matrix.T):
            raise ValueError('matrix must be a symmetric matrix')
        if self._n_nodes != 0 and self._n_nodes != shape[0]:
            raise ValueError('matrix shape must match number of nodes')
        self._affinity = matrix
        self._n_nodes = shape[0]


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
            d = np.array(self._affinity.sum(0)).flatten()
            self._stationary = d * (1. / d.sum())
        return self._stationary.copy()
        


    def buildTransitionMatrix(self):
        """Build Markov transition matrix from the affinity matrix."""
        
        start = time.time()
        self._transition = np.dot(self._affinity, 
                  sparse.lil_diags([1 / np.array(self._affinity.sum(0)
                                    ).flatten()], [0], 
                                   [self._n_nodes, self._n_nodes]))
        LOGGER.debug('Markov transition matrix was built in {0:.2f}s.'
                     .format(time.time()-start))

    def setTransitionMatrix(self, matrix):
        """Set the Markov transition matrix."""
        
        shape = matrix.shape 
        if shape[0] != shape[1]:
            raise ValueError('matrix must be a square matrix')
        if not np.all(matrix.sum(0).round(6) == 1):
            raise ValueError('columns of the matrix must add up to 1 '
                             '(sum is rounded to six significant digits)')
        if self._n_nodes != 0 and self._n_nodes != shape[0]:
            raise ValueError('matrix shape must match number of nodes')
        self._transition = matrix
        self._n_nodes = shape[0]

    def getTransitionMatrix(self):
        """Return a copy of the Markov transition matrix."""
        
        if self._transition is not None:
            return self._transition.copy()
    
    def grainCoarser(self, power=4, symmetry=None):
        """
        
        :arg power: the power to which Markov transition matrix is raised
            before initial kernel is estimated 
        :type power: int
        :arg symmetry: symmetrically related chains, e.g. ``"ABC"``, 
            ``["ABC", "DEF"]`` 
        :type symmetry: str, list
        
        """
        
        assert isinstance(power, int), 'power must be an integer'
        assert power > 0, 'power must be positive'
        # Power
        start = time.time()
        transition = np.linalg.matrix_power(self._transition, power) 
        LOGGER.debug('Markov transition matrix was raised to the power {1:d} '
                     'in {0:.2f}s.'.format(time.time()-start, power))
        if symmetry is not None and self._resmap is not None:
            resmap = self._resmap
            if isinstance(symmetry, str):
                symmetry = [symmetry]
            symdict = {}
            for i, chids in enumerate(symmetry):
                if not isinstance(chids, str):
                    raise TypeError('symmetry must be a string or a list of '
                                    'strings')
                chids = set(chids)
                for chid in chids:
                    symdict[chid] = chids 
        n_res = self._n_nodes    
        # Initialize
        start = time.time()
        torf = np.ones(n_res, np.bool)
        stationary = self.getStationaryDistribution()
        mu = stationary.copy()
        columns = []
        
        while np.any(torf):
            mu = mu * torf
            which = (mu == mu.max()).nonzero()[0][0]
            if symmetry:
                chid, rnum, rins = resmap[which]
                chids = symdict.get(chid, None)
                if chids is not None:
                    whiches = [resmap[(chid, rnum, rins)] for chid in chids]
            else:
                whiches = (which, )
            for which in whiches:
                columns.append(which)
                torf[which] = False
                k = transition[:, which].toarray().flatten()
                torf[ k > k.max() / 2 ] = False
        columns.sort()
        length = len(columns)
        kernel = transition[:, columns]
        delta = np.ones(length) * (1. / length)
        LOGGER.debug('Initial kernel matrix was built '
                     'in {0:.2f}s.'.format(time.time()-start))
        start = time.time()
        delta_old = np.ones(length)
        step = 0
        while ((delta - delta_old) ** 2 ).sum() ** 0.5 > 0.000001:
            # E-step
            divisor = sparse.dia_matrix((1/(kernel * delta), 0), 
                                        shape=(n_res, n_res))
            ownership = divisor * kernel * sparse.dia_matrix((delta, 0), 
                                                        shape=(length, length))
            # M-step
            delta_old = delta
            delta = stationary * ownership
            kernel = sparse.dia_matrix((stationary, 0), shape=(n_res, n_res)) * \
                        ownership * sparse.dia_matrix((1 / delta, 0), 
                                                        shape=(length, length))
            step += 1
        LOGGER.debug('EM converged after {1:d} steps '
                     'in {0:.2f}s.'.format(time.time()-start, step))
        self._ownership = ownership 
        self._kernel = kernel
        self._delta = delta
        
        model = MarkovModel(self.getTitle())
        delta_diag = sparse.dia_matrix((delta, 0), shape=(length, length))
        transition = delta_diag * kernel.T * \
            sparse.dia_matrix((1 / (kernel * delta), 0), shape=(n_res, n_res)) * \
            kernel
        
        affinity = transition * delta_diag
        self._lower = model
        model._higher = self
        model.setAffinityMatrix(affinity)
        model.setTransitionMatrix(transition)
        return model
        
if __name__ == '__main__':
    
     #p = parsePDB('1aon')
     m = MarkovModel('')
     m.buildAffinityMatrix(p.select('protein'))
     m.buildTransitionMatrix()
     c = m.grainCoarser(power=8, symmetry=['ABCDEFG', 'HIJKLMN', 'OPQRSTU'])
