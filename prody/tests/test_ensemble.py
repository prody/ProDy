#!/usr/bin/python
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

"""This module contains unit tests for :mod:`~prody.ensemble`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import unittest
import numpy as np
from numpy.testing import *

from prody import *
from prody.tools import *
from test_datafiles import *

prody.setVerbosity('none')

ATOL = 1e-5
RTOL = 0

ATOMS = parseDatafile('multi_model_truncated', subset='ca')
ALLATOMS = parseDatafile('multi_model_truncated')
DCD = parseDatafile('dcd')
COORDS = ATOMS.getCoords()
COORDSETS = ATOMS.getCoordsets()
ENSEMBLE = Ensemble(ATOMS)
CONF = ENSEMBLE[0]
DATA = DATA_FILES['multi_model_truncated']
ENSEMBLE_RMSD = DATA['rmsd_ca']
ENSEMBLE_SUPERPOSE = DATA['rmsd_ca_aligned']

ENSEMBLEW = Ensemble(ATOMS)
ENSEMBLEW.setWeights(np.ones(len(ATOMS), dtype=float))
CONFW = ENSEMBLEW[0]
        
PDBENSEMBLE = PDBEnsemble('PDBEnsemble')
PDBENSEMBLE.setCoords(COORDS)
WEIGHTS = []
for i, xyz in enumerate(ATOMS.iterCoordsets()):
    weights = np.ones((len(ATOMS), 1), dtype=float)
    if i > 0:
        weights[i] = 0
        weights[-i] = 0 
        PDBENSEMBLE.addCoordset(xyz, weights=weights)
    else:
        PDBENSEMBLE.addCoordset(xyz)
    WEIGHTS.append(weights)
PDBCONF = PDBENSEMBLE[0]
WEIGHTS = np.array(WEIGHTS)
WEIGHTS_BOOL = np.tile(WEIGHTS.astype(bool), (1,1,3))  
WEIGHTS_BOOL_INVERSE = np.invert(WEIGHTS_BOOL)


class TestEnsemble(unittest.TestCase):
    
    def testGetCoordinates(self):
        """Test correctness of reference coordinates."""
        
        assert_equal(ENSEMBLE.getCoords(), COORDS,
                     'failed to get correct coordinates')
    
    
    def testGetCoordsets(self):
        """Test correctness of all coordinates."""

        assert_equal(ATOMS.getCoordsets(), ENSEMBLE.getCoordsets(),
                     'failed to add coordinate sets for Ensemble')
    
    
    def testGetWeights(self):
        """Test correctness of all weights."""

        assert_equal(ENSEMBLEW.getWeights().ndim, 2,
                     'failed to get correct weights ndim')
        assert_equal(ENSEMBLEW.getWeights().shape, 
                     (ENSEMBLEW.numAtoms(), 1),
                    'failed to get correct weights shape')
        assert_equal(ENSEMBLEW.getWeights(),
                     np.ones((ENSEMBLEW.numAtoms(), 1), float),
                     'failed to get expected weights')

    def testSlicingCopy(self):
        """Test making a copy by slicing operation."""
        
        SLICE = ENSEMBLE[:]
        assert_equal(SLICE.getCoords(), ENSEMBLE.getCoords(),
                     'slicing copy failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), ENSEMBLE.getCoordsets(),
                     'slicing copy failed to add coordinate sets')

    def testSlicing(self):
        """Test slicing operation."""
        
        SLICE = ENSEMBLE[:2]
        assert_equal(SLICE.getCoords(), ENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), ENSEMBLE.getCoordsets([0,1]),
                     'slicing failed to add coordinate sets')

    def testSlicingList(self):
        """Test slicing with a list."""
        
        SLICE = ENSEMBLE[[0,2]]
        assert_equal(SLICE.getCoords(), ENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), ENSEMBLE.getCoordsets([0,2]),
                     'slicing failed to add coordinate sets for Ensemble')


    def testSlicingWeights(self):
        """Test slicing operation with weights."""
        
        SLICE = ENSEMBLEW[:2]
        assert_equal(SLICE.getWeights(), ENSEMBLEW.getWeights(),
                     'slicing failed to set weights')

    def testIterCoordsets(self):
        """Test coordinate iteration."""
        
        for i, xyz in enumerate(ENSEMBLE.iterCoordsets()):
            assert_equal(xyz, ATOMS.getCoordsets(i),
                         'failed yield correct coordinates')
            
    def testGetNumAtoms(self):

        self.assertEqual(ENSEMBLE.numAtoms(), ATOMS.numAtoms(),
                         'failed to get correct number of atoms')  
            
    def testGetNumCsets(self):

        self.assertEqual(ENSEMBLE.numCoordsets(), ATOMS.numCoordsets(),
                         'failed to get correct number of coordinate sets')  

    def testGetRMSDs(self):
        
        assert_allclose(ENSEMBLE.getRMSDs(), ENSEMBLE_RMSD,
                        rtol=0, atol=1e-3,  
                        err_msg='failed to calculate RMSDs sets')

    def testSuperpose(self):
        
        ensemble = ENSEMBLE[:]
        ensemble.superpose()
        assert_allclose(ensemble.getRMSDs(), ENSEMBLE_SUPERPOSE,
                        rtol=0, atol=1e-3,
                        err_msg='failed to superpose coordinate sets')

    def testGetRMSDsWeights(self):
        
        assert_allclose(ENSEMBLEW.getRMSDs(), ENSEMBLE_RMSD,
                        rtol=0, atol=1e-3,
                        err_msg='failed to calculate RMSDs')

    def testSuperposeWeights(self):
        
        ensemble = ENSEMBLEW[:]
        ensemble.superpose()
        assert_allclose(ensemble.getRMSDs(), ENSEMBLE_SUPERPOSE,
                        rtol=0, atol=1e-3,
                        err_msg='failed to superpose coordinate sets')

    def testDelCoordsetMiddle(self):
        
        ensemble = ENSEMBLE[:]
        ensemble.delCoordset(1)
        assert_equal(ensemble.getCoordsets(), ATOMS.getCoordsets([0,2]),
                     'failed to delete middle coordinate set for Ensemble')
        
    def testDelCoordsetAll(self):
        
        ensemble = ENSEMBLE[:]
        ensemble.delCoordset(range(len(ENSEMBLE)))
        self.assertIsNone(ensemble.getCoordsets(),
                        'failed to delete all coordinate sets for Ensemble')
        assert_equal(ensemble.getCoords(), COORDS,
                     'failed when deleting all coordinate sets')


    def testConcatenation(self):
        """Test concatenation of ensembles without weights."""
        
        ensemble = ENSEMBLE + ENSEMBLE
        assert_equal(ensemble.getCoordsets(range(3)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(range(3,6)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')

    def testConcatenationWeights(self):
        """Test concatenation of ensembles with weights."""
        
        ensemble = ENSEMBLEW + ENSEMBLEW
        assert_equal(ensemble.getCoordsets(range(3)), ATOMS.getCoordsets(), 
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(range(3,6)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')
        assert_equal(ensemble.getWeights(), ENSEMBLEW.getWeights(),
                     'concatenation failed')

    def testConcatenationNoweightsWeights(self):
        """Test concatenation of ensembles without and with weights."""
        
        ensemble = ENSEMBLE + ENSEMBLEW
        assert_equal(ensemble.getCoordsets(range(3)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(range(3,6)), ATOMS.getCoordsets(),
                    'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')
        self.assertIsNone(ensemble.getWeights(), 'concatenation failed')

    def testConcatenationWeightsNoweights(self):
        """Test concatenation of ensembles with and without weights."""
        
        ensemble = ENSEMBLEW + ENSEMBLE 
        assert_equal(ensemble.getCoordsets(range(3)), ATOMS.getCoordsets(),
                     'failed at concatenation for Ensemble')
        assert_equal(ensemble.getCoordsets(range(3,6)), ATOMS.getCoordsets(),
                     'failed at concatenation for Ensemble')
        assert_equal(ensemble.getCoords(), COORDS,
                     'failed at concatenation for Ensemble')
        assert_equal(ensemble.getWeights(), ENSEMBLEW.getWeights(),
                     'failed at concatenation for Ensemble')


class TestConformation(unittest.TestCase): 
    
    
    def testCoordinates(self):
        
        assert_equal(CONF.getCoords(), COORDS,
                     'failed to get coordinates for conformation')
                        
    def testWeights(self):
        
        weights = CONFW.getWeights()
        self.assertEqual(weights.ndim, 2,
                        'wrong ndim for weights of Conformation')
        self.assertTupleEqual(weights.shape, 
                              (CONFW.numAtoms(), 1),
                              'wrong shape for weights of Conformation')
        assert_equal(weights, np.ones((ATOMS.numAtoms(), 1), float),
                     'failed to set weights for Conformation')
                                                
    def testCoordinatesForAll(self):
        
        for i, conf in enumerate(ENSEMBLE):
            assert_equal(conf.getCoords(), ATOMS.getCoordsets(i),
                         'failed to get coordinates for Conformation')
                         
    def testGetIndex(self):
        """Test get index function."""

        for i, conf in enumerate(ENSEMBLE):
            self.assertEqual(conf.getIndex(), i,
                             'failed to get correct index')        
                        
    def testGetNumAtoms(self):
        """Test get index function."""

        for i, conf in enumerate(ENSEMBLE):
            self.assertEqual(conf.numAtoms(), ATOMS.numAtoms(),
                             'failed to get correct number of atoms')        

class TestPDBEnsemble(unittest.TestCase):
    
    def testGetCoordinates(self):
        
        assert_equal(PDBENSEMBLE.getCoords(), COORDS,
                     'failed to set reference coordinates for PDBEnsemble')
    
    def testGetCoordsets(self):

        assert_equal(PDBENSEMBLE.getCoordsets()[WEIGHTS_BOOL], 
                     ATOMS.getCoordsets()[WEIGHTS_BOOL],
                     'failed to add coordinate sets for PDBEnsemble')
    
    def testGetWeights(self):

        self.assertEqual(PDBENSEMBLE.getWeights().ndim, 3,
                        'wrong ndim for weights of PDBEnsemble')
        self.assertTupleEqual(PDBENSEMBLE.getWeights().shape, 
                              (PDBENSEMBLE.numCoordsets(),
                               PDBENSEMBLE.numAtoms(), 1),
                               'wrong shape for weights of PDBEnsemble')
        assert_equal(PDBENSEMBLE.getWeights(), WEIGHTS,
                     'failed to get correct weights')


    def testSlicingCopy(self):
        
        SLICE = PDBENSEMBLE[:]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE.getCoords(),
                     'slicing copy failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE.getCoordsets(),
                     'slicing copy failed to add coordinate sets')

    def testSlicing(self):
        
        SLICE = PDBENSEMBLE[:2]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE.getCoordsets([0,1]),
                     'slicing failed to add coordinate sets')

    def testSlicingList(self):
        
        SLICE = PDBENSEMBLE[[0,2]]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE.getCoordsets([0,2]),
                     'slicing failed to add coordinate sets')

    def testSlicingWeights(self):
        
        SLICE = PDBENSEMBLE[:2]
        assert_equal(SLICE.getWeights(), PDBENSEMBLE.getWeights()[:2],
                     'slicing failed to set weights')

    def testIterCoordsets(self):
        
        for i, xyz in enumerate(ENSEMBLE.iterCoordsets()):
            assert_equal(xyz[WEIGHTS_BOOL[i]], 
                         ATOMS.getCoordsets(i)[WEIGHTS_BOOL[i]],
                         'failed iterate coordinate sets')
            
    def testGetNumAtoms(self):

        self.assertEqual(PDBENSEMBLE.numAtoms(), ATOMS.numAtoms(),
                         'failed to get correct number of atoms')  
            
    def testGetNumCsets(self):

        self.assertEqual(PDBENSEMBLE.numCoordsets(), 
                         ATOMS.numCoordsets(),
                         'failed to get correct number of coordinate sets')  

    def testDelCoordsetMiddle(self):
        
        ensemble = PDBENSEMBLE[:]
        ensemble.delCoordset(1)
        assert_equal(ensemble.getCoordsets()[WEIGHTS_BOOL[[0,2]]],
                     ATOMS.getCoordsets([0,2])[WEIGHTS_BOOL[[0,2]]], 
                    'failed to delete middle coordinate set')
        
    def testDelCoordsetAll(self):
        """Test consequences of deleting all coordinate sets."""
        
        ensemble = PDBENSEMBLE[:]
        ensemble.delCoordset(range(len(PDBENSEMBLE)))
        self.assertIsNone(ensemble.getCoordsets(),
                        'failed to delete all coordinate sets')
        self.assertIsNone(ensemble.getWeights(), 'failed to delete weights '
                          'with all coordinate sets')
        assert_equal(ensemble.getCoords(), COORDS,
                     'failed to delete all coordinate sets')


    def testConcatenation(self):
        """Test concatenation of PDB ensembles."""
        
        ensemble = PDBENSEMBLE + PDBENSEMBLE
        assert_equal(ensemble.getCoordsets(range(3)),
                     PDBENSEMBLE.getCoordsets(), 
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(range(3,6)),
                     PDBENSEMBLE.getCoordsets(), 
                     'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')
        assert_equal(ensemble.getWeights()[range(3)],
                     PDBENSEMBLE.getWeights(),
                     'concatenation failed')
        assert_equal(ensemble.getWeights()[range(3,6)], 
                     PDBENSEMBLE.getWeights(),
                     'concatenation failed')

class TestPDBConformation(unittest.TestCase): 
    
    
    def testCoordinates(self):
        
        assert_equal(PDBCONF.getCoords(), COORDS,
                     'failed to get correct coordinates')
                        
    def testWeights(self):
        
        weights = PDBCONF.getWeights()
        self.assertEqual(weights.ndim, 2,
                        'wrong ndim for weights')
        self.assertTupleEqual(weights.shape, (PDBCONF.numAtoms(), 1),
                              'wrong shape for weights')
        assert_equal(PDBCONF.getWeights(), WEIGHTS[0],
                     'failed to set weights')

    def testWeightsForAll(self):

        for i, conf in enumerate(PDBENSEMBLE):
            assert_equal(conf.getWeights(), WEIGHTS[i],
                         'failed to set weights')

                                                
    def testCoordinatesForAll(self):
        
        for i, conf in enumerate(PDBENSEMBLE):
            assert_equal(conf.getCoords()[WEIGHTS_BOOL[i]], 
                         ATOMS.getCoordsets(i)[WEIGHTS_BOOL[i]],
                         'failed to get coordinates')
                         
    def testGetIndex(self):

        for i, conf in enumerate(PDBENSEMBLE):
            self.assertEqual(conf.getIndex(), i,
                             'failed to get correct index')        
                        
    def testGetNumAtoms(self):

        for i, conf in enumerate(PDBENSEMBLE):
            self.assertEqual(conf.numAtoms(), ATOMS.numAtoms(),
                             'failed to get correct number of atoms')

class TestCalcOccupancies(unittest.TestCase):

    def testResults(self):
        assert_equal(calcOccupancies(PDBENSEMBLE), WEIGHTS.sum(0).flatten(),
                     'calcOccupancies failed')

    def testInvalidType(self):
        
        self.assertRaises(TypeError, calcOccupancies, ENSEMBLE)
        
    def testWeightsNone(self):
        
        self.assertRaises(ValueError, calcOccupancies, PDBEnsemble())


class TestDCDFile(unittest.TestCase):
    
    def setUp(self):
        
        self.dcd = os.path.join(TEMPDIR, 'temp.dcd')
    
    def testWriteDCD(self):
        dcd = writeDCD(self.dcd, ALLATOMS)
        self.assertEqual(dcd, self.dcd, 'failed to write DCD file')
        
    def testParseDCD(self):
        e = parseDCD(writeDCD(self.dcd, ALLATOMS))
        assert_equal(e._getCoordsets(), DCD._getCoordsets(),
                     err_msg='failed to parse DCD file correctly')

    def testWrite(self):
        dcd = DCDFile(self.dcd, 'w')
        dcd.write(ENSEMBLE.getCoordsets())
        dcd.close()
        e = parseDCD(self.dcd)
        assert_allclose(e._getCoordsets(), ENSEMBLE._getCoordsets(),
                        rtol=RTOL, atol=ATOL,
                        err_msg='failed to parse DCD file correctly')

    def testWriteModeAppend(self):
        dcd = DCDFile(writeDCD(self.dcd, ENSEMBLE), 'a')
        dcd.write(ENSEMBLE.getCoordsets())
        dcd.close()
        e = parseDCD(self.dcd)
        n_csets = len(ENSEMBLE)
        coordsets = e._getCoordsets()
        assert_equal(coordsets, coordsets, 
                     'failed to parse DCD file correctly')
        assert_allclose(coordsets[:n_csets], ENSEMBLE._getCoordsets(),
                        rtol=RTOL, atol=ATOL,
                        err_msg='failed to parse DCD file correctly')


    
if __name__ == '__main__':
    prody.changeVerbosity('none')
    unittest.main()
    prody.changeVerbosity('debug')
