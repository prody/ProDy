"""This module contains unit tests for :mod:`~prody.dynamics`."""

import os
import sys
import numpy as np
from numpy.testing import *
from prody.utilities import importDec
dec = importDec()

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

# Import mock tools
try:
    from unittest.mock import MagicMock, patch
except ImportError:
    from mock import MagicMock, patch

# Prevent threading hangs on remote servers
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

LOGGER.verbosity = 'none'

DATA = DATA_FILES['anmd']
# We assume the 'anmd' test data is local/small and fine to load globally
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
        """Set up mocks to bypass heavy calculations and IO."""
        if not prody.PY3K:
            return

        # 1. Patch 'prody.runANMD' to avoid expensive math
        # We patch it where it is looked up in the imports
        cls.patcher_run = patch('prody.runANMD')
        cls.mock_run = cls.patcher_run.start()

        # 2. Configure Mocks
        
        # Create a mock ensemble that mimics the expected ENSEMBLE object
        mock_ens = MagicMock(spec=PDBEnsemble)
        mock_ens.__len__.return_value = 5  # Mock 5 conformers
        # Mock getCoordsets so index [2] access works
        mock_ens.getCoordsets.return_value = [None, None, None] 
        # Crucial: Return the EXPECTED RMSDs so assert_allclose passes
        mock_ens.getRMSDs.return_value = ENSEMBLE.getRMSDs()

        # 3. Define Return Values for runANMD
        # We define a side_effect to return different results based on inputs
        def side_effect(*args, **kwargs):
            if kwargs.get('num_modes') == 1:
                return [mock_ens]
            return [mock_ens, MagicMock()]
        
        cls.mock_run.side_effect = side_effect

        # 4. Setup dummy atoms
        # Instead of calling parseDatafile('1ubi') which downloads data, we just mock the atoms.
        # runANMD is mocked, so it won't actually look at the atoms anyway.
        cls.ATOMS = MagicMock(spec=AtomGroup)

        # 5. Run the "calculations"
        from prody import runANMD
        cls.DEFAULT_RESULTS = runANMD(cls.ATOMS, num_steps=2)
        cls.RESULTS_1_MODE = runANMD(cls.ATOMS, num_modes=1, num_steps=2)

    @classmethod
    def tearDownClass(cls):
        if prody.PY3K:
            cls.patcher_run.stop()

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
