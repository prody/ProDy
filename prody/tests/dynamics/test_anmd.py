"""This module contains unit tests for :mod:`~prody.dynamics`."""

import numpy as np
from numpy.testing import *
from prody.utilities import importDec
dec = importDec()

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

DATA = DATA_FILES['anmd']
ENSEMBLE = PDBEnsemble(parseDatafile('anmd'))
ENSEMBLE.setCoords(ENSEMBLE.getCoordsets()[2])

class TestANMD(unittest.TestCase):

    def setUp(self):
        if prody.PY3K:
            from prody import runANMD
            self.runANMD = runANMD
            self.ATOMS = parseDatafile('1ubi')

    def testAnmdAtomsWrongType(self):
        """Test response to wrong type *atoms* argument."""
        if prody.PY3K:
            self.assertRaises(TypeError, self.runANMD, 'nogood')

    def testAnmdNumModesWrongType(self):
        """Test response to wrong type *num_modes* argument."""
        if prody.PY3K:
            self.assertRaises(TypeError, self.runANMD, self.ATOMS, 'nogood', num_steps=2)

    def testAnmdRmsdWrongType(self):
        """Test response to wrong type *max_rmsd* argument."""
        if prody.PY3K:
            self.assertRaises(TypeError, self.runANMD, self.ATOMS, max_rmsd='nogood', num_steps=2)

    def testAnmdStepsWrongType(self):
        """Test response to wrong type *num_steps* argument."""
        if prody.PY3K:
            self.assertRaises(TypeError, self.runANMD, self.ATOMS, num_steps='nogood')

    def testAnmdToleranceWrongType(self):
        """Test response to wrong type *tolerance* argument."""
        if prody.PY3K:
            self.assertRaises(TypeError, self.runANMD, self.ATOMS, tolerance='nogood', num_steps=2)

class TestAnmdResults(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if prody.PY3K:
            from prody import runANMD
            cls.ATOMS = parseDatafile('1ubi')
            # Run calculations once here to speed up tests
            cls.DEFAULT_RESULTS = runANMD(cls.ATOMS, num_steps=2)
            cls.RESULTS_1_MODE = runANMD(cls.ATOMS, num_modes=1, num_steps=2)

    def testResults(self):
        """Test results with default parameters"""
        if prody.PY3K:
            DEFAULT_RESULTS = self.DEFAULT_RESULTS
            ens1 = DEFAULT_RESULTS[0]

            assert_equal(len(DEFAULT_RESULTS), 2,
                        'runANMD with default parameters failed to give 2 ensembles')
            assert_equal(len(ens1), DATA['models'],
                        'runANMD with default parameters failed to give ensembles with 5 conformers')       

            ens1.setCoords(ens1.getCoordsets()[2])
            assert_allclose(ens1.getRMSDs(), ENSEMBLE.getRMSDs(), 
                            rtol=1e-10, atol=0.2, # may not be so close
                            err_msg='runANMD with default parameters failed to give expected RMSDs')

    def testResultsNumModes1(self):
        """Test that num_modes=1 gives 1 ensemble"""
        if prody.PY3K:
            RESULTS = self.RESULTS_1_MODE
            assert_equal(len(RESULTS), 1,
                        'runANMD with num_modes=1 failed to give 1 ensemble')

            ens1 = RESULTS[0]
            assert_equal(len(ens1), DATA['models'],
                        'runANMD with num_modes=1 failed to give ensembles with 5 conformers')       

            ens1.setCoords(ens1.getCoordsets()[2])
            assert_allclose(ens1.getRMSDs(), ENSEMBLE.getRMSDs(), 
                            rtol=1e-10, atol=0.25, # may not be so close
                            err_msg='runANMD with num_modes=1 failed to give expected RMSDs')

if __name__ == '__main__':
    unittest.main()
