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

class TestLDA(unittest.TestCase):

    def setUp(self):

        from prody import runANMD
        self.runANMD = runANMD
        self.ATOMS = parseDatafile('1ubi')

    def testAnmdAtomsWrongType(self):
        """Test response to wrong type *atoms* argument."""

        self.assertRaises(TypeError, self.runAMND, 'nogood')

    def testAnmdNumModesWrongType(self):
        """Test response to wrong type *num_modes* argument."""

        self.assertRaises(TypeError, self.runAMND, self.ATOMS, 'nogood')

    def testAnmdRmsdWrongType(self):
        """Test response to wrong type *max_rmsd* argument."""

        self.assertRaises(TypeError, self.runAMND, self.ATOMS, max_rmsd='nogood')

    def testAnmdStepsWrongType(self):
        """Test response to wrong type *num_steps* argument."""

        self.assertRaises(TypeError, self.runAMND, self.ATOMS, num_steps='nogood')

    def testAnmdToleranceWrongType(self):
        """Test response to wrong type *tolerance* argument."""

        self.assertRaises(TypeError, self.runAMND, self.ATOMS, tolerance='nogood')

    def testAnmdSolventWrongType(self):
        """Test response to wrong type *solvent* argument."""

        self.assertRaises(TypeError, self.runAMND, self.ATOMS, solvent=1)

    def testAnmdSolventWrongValue(self):
        """Test response to wrong value *solvent* argument."""

        self.assertRaises(ValueError, self.runAMND, self.ATOMS, solvent='nogood')

    def testAnmdFFWrongType(self):
        """Test response to wrong type *force_field* argument."""

        self.assertRaises(TypeError, self.runAMND, self.ATOMS, force_field=1)

    def testAnmdFFWrongValue(self):
        """Test response to wrong value *force_field* argument."""

        self.assertRaises(ValueError, self.runAMND, self.ATOMS, force_field='nogood')

class TestLDAResults(TestLDA):

    def setUp(self):
        self.ATOMS = parseDatafile('1ubi')

    def testResults(self):
        DEFAULT_RESULTS = runANMD(self.ATOMS)
        ens1 = DEFAULT_RESULTS[0]

        assert_equal(len(DEFAULT_RESULTS), 2,
                     'runANMD with default parameters failed to give 2 ensembles')
        assert_equal(len(ens1), DATA['models'],
                     'runANMD with default parameters failed to give ensembles with 5 conformers')       

        ens1.setCoords(ens1.getCoordsets()[2])
        assert_allclose(ens1.getRMSDs(), ENSEMBLE.getRMSDs(), 
                        rtol=1e-10, atol=0.1, # may not be so
                        err_msg='runANMD with default parameters failed to give expected RMSDs')

    def testNumModes3(self):
        RESULTS = runANMD(self.ATOMS, num_modes=1)
        assert_equal(len(RESULTS), 1,
                     'runANMD with num_modes=1 failed to give 1 ensemble')

        ens1 = RESULTS[0]
        assert_equal(len(ens1), DATA['models'],
                     'runANMD with num_modes=1 failed to give ensembles with 5 conformers')       

        ens1.setCoords(ens1.getCoordsets()[2])
        assert_allclose(ens1.getRMSDs(), ENSEMBLE.getRMSDs(), 
                        rtol=1e-10, atol=0.1, # may not be so
                        err_msg='runANMD with num_modes=1 failed to give expected RMSDs')

if __name__ == '__main__':
    unittest.main()
