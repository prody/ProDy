#!/usr/bin/python
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

"""This module contains unit tests for :mod:`~prody.ensemble` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import os
import os.path
import sys
import unittest
import tempfile
import inspect
import prody
import numpy as np
from prody.ensemble import *
from prody import ensemble

TEMPDIR = tempfile.gettempdir()
TESTS_PATH = os.path.abspath(os.path.split(inspect.getfile(
                                                   inspect.currentframe()))[0])

prody.changeVerbosity('none')

PDB_FILES = {
    'multi_model_truncated': {
        'pdb': '2k39',
        'path': os.path.join(TESTS_PATH, 'data/pdb2k39_truncated.pdb'),
        'atoms': 167,
        'models': 3
    },
}

ATOMS = prody.parsePDB(PDB_FILES['multi_model_truncated']['path'], subset='ca')

ENSEMBLE = Ensemble(ATOMS)
CONF = ENSEMBLE[0]

ENSEMBLEW = Ensemble(ATOMS)
ENSEMBLEW.setWeights(np.ones(len(ATOMS), dtype=float))
CONFW = ENSEMBLEW[0]

PDBENSEMBLE = PDBEnsemble('PDBEnsemble')
PDBENSEMBLE.setCoordinates(ATOMS.getCoordinates())
for i, xyz in enumerate(ATOMS.iterCoordsets()):
    if i > 0:
        weights = np.ones(len(ATOMS), dtype=float)
        weights[i] = 0
        weights[-i] = 0 
        PDBENSEMBLE.addCoordset(ATOMS.getCoordinates(), weights=weights)
    else:
        PDBENSEMBLE.addCoordset(ATOMS.getCoordinates())
        
class TestEnsemble(unittest.TestCase):
    
    def testCoordinates(self):
        
        self.assertTrue(np.all(ATOMS.getCoordinates() == 
                               ENSEMBLE.getCoordinates()),
                        'failed to set reference coordinates for Ensemble')
    
    def testCoordsets(self):

        self.assertTrue(np.all(ATOMS.getCoordsets() == 
                               ENSEMBLE.getCoordsets()),
                        'failed to add coordinate sets for Ensemble')
    
    def testWeights(self):

        self.assertEqual(ENSEMBLEW.getWeights().ndim, 2,
                        'wrong ndim for weights of Ensemble')
        self.assertTupleEqual(ENSEMBLEW.getWeights().shape, 
                              (ENSEMBLEW.getNumOfAtoms(), 1),
                               'wrong shape for weights of Ensemble')
        self.assertTrue(np.all(ENSEMBLEW.getWeights() == 
                               np.ones(ENSEMBLEW.getNumOfAtoms(), float)),
                        'failed to set weights for Ensemble')


    def testSlicingCopy(self):
        
        SLICE = ENSEMBLE[:]
        self.assertTrue(np.all(SLICE.getCoordinates() == 
                               ENSEMBLE.getCoordinates()),
                        'slicing copy failed to set reference coordinates for '
                        'Ensemble')
        self.assertTrue(np.all(SLICE.getCoordsets() == 
                               ENSEMBLE.getCoordsets()),
                        'slicing copy failed to add coordinate sets for '
                        'Ensemble')

    def testSlicing(self):
        
        SLICE = ENSEMBLE[:2]
        self.assertTrue(np.all(SLICE.getCoordinates() == 
                               ENSEMBLE.getCoordinates()),
                        'slicing failed to set reference coordinates for '
                        'Ensemble')
        self.assertTrue(np.all(SLICE.getCoordsets() == 
                               ENSEMBLE.getCoordsets([0,1])),
                        'slicing failed to add coordinate sets for Ensemble')

    def testSlicingWeights(self):
        
        SLICE = ENSEMBLEW[:2]
        self.assertTrue(np.all(SLICE.getWeights() == 
                               ENSEMBLEW.getWeights()),
                        'slicing failed to set weights for ensemble')

    def testIterCoordsets(self):
        
        for i, xyz in enumerate(ENSEMBLE.iterCoordsets()):
            self.assertTrue(np.all(xyz == ATOMS.getCoordsets(i)),
                            'failed iterate coordinate sets for Ensemble')
            
    def testGetNumOfAtoms(self):

        self.assertEqual(ENSEMBLE.getNumOfAtoms(), ATOMS.getNumOfAtoms(),
                         'failed to get correct number of atoms for Ensemble')  
            
    def testGetNumOfCoordsets(self):

        self.assertEqual(ENSEMBLE.getNumOfCoordsets(), 
                         ATOMS.getNumOfCoordsets(),
                         'failed to get correct number of coordinate sets for ' 
                         'Ensemble')  


class TestConformation(unittest.TestCase): 
    
    
    def testCoordinates(self):
        
        self.assertTrue(np.all(ATOMS.getCoordinates() == 
                               CONF.getCoordinates()),
                        'failed to get coordinates for conformation')
                        
    def testWeights(self):
        
        self.assertEqual(CONFW.getWeights().ndim, 2,
                        'wrong ndim for weights of Conformation')
        self.assertTupleEqual(CONFW.getWeights().shape, 
                              (CONFW.getNumOfAtoms(), 1),
                               'wrong shape for weights of Conformation')
        self.assertTrue(np.all(CONFW.getWeights() == 
                               np.ones(ATOMS.getNumOfAtoms(), float)),
                        'failed to set weights for Conformation')
                                                
    def testCoordinatesForAll(self):
        
        for i, conf in enumerate(ENSEMBLE):
            self.assertTrue(np.all(ATOMS.getCoordsets(i) == 
                                   conf.getCoordinates()),
                            'failed to get coordinates for Conformation')
    def testGetIndex(self):

        for i, conf in enumerate(ENSEMBLE):
            self.assertEqual(conf.getIndex(), i,
                             'failed to get correct index for Conformation')        
                        
    def testGetNumOfAtoms(self):

        for i, conf in enumerate(ENSEMBLE):
            self.assertEqual(conf.getNumOfAtoms(), ATOMS.getNumOfAtoms(),
                             'failed to get correct number of atoms for ' 
                             'Conformation')        


if __name__ == '__main__':
    unittest.main()
