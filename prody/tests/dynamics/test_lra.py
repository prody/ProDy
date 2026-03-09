"""This module contains unit tests for :mod:`~prody.dynamics`."""

import numpy as np
from numpy.testing import assert_array_equal

from prody.utilities import calcTree, findSubgroups
from prody.dynamics import LRA
from prody.ensemble import PDBEnsemble
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import parseDatafile

LOGGER.verbosity = 'none'

ATOL = 1e-5
RTOL = 0

ENSEMBLE = PDBEnsemble(parseDatafile('2k39_insty_dcd'))

RMSD_MATRIX = ENSEMBLE.getRMSDs(pairwise=True)
RMSD_TREE = calcTree(ENSEMBLE.getLabels(), RMSD_MATRIX)

SUBGROUPS_3P5 = findSubgroups(RMSD_TREE, 3.5)
TWO_CLASSES = [np.nonzero(np.array([i in sg for sg in SUBGROUPS_3P5]))[0][0]
               for i in ENSEMBLE.getLabels()]

SUBGROUPS_3P2 = findSubgroups(RMSD_TREE, 3.2)
THREE_CLASSES = [np.nonzero(np.array([i in sg for sg in SUBGROUPS_3P2]))[0][0]
               for i in ENSEMBLE.getLabels()]

# LRA specific setup
lra_two = None
lra_three = None
try:
    lra_two = LRA()
    lra_two.calcModes(ENSEMBLE, TWO_CLASSES)

    lra_three = LRA()
    lra_three.calcModes(ENSEMBLE, THREE_CLASSES)

    SKLEARN_AVAILABLE = True
except ImportError:    
    SKLEARN_AVAILABLE = False

class TestLRA(unittest.TestCase):

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def setUp(self):

        self.model = LRA()
        self.calcModes = self.model.calcModes

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def testCalcCoordsWrongType(self):
        """Test response to wrong type *coords* argument."""

        self.assertRaises(TypeError, self.calcModes, 'nogood', TWO_CLASSES)

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def testCalcWrongCoords(self):
        """Test response to wrong coords.dtype."""

        array = np.array([['a','a','a'] for _ in range(10)])
        self.assertRaises(ValueError, self.calcModes, array, TWO_CLASSES)

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def testCalcLabelsWrongType(self):
        """Test response to wrong type *coords* argument."""

        self.assertRaises(TypeError, self.calcModes, ENSEMBLE, 'nogood')


class TestLRAResults(TestLRA):

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def testCoefficients(self):
        # Test that coefficients are calculated
        coeffs = lra_two.getEigvecs()
        self.assertIsNotNone(coeffs)
        self.assertEqual(coeffs.shape[1], lra_two._n_modes)

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def testEigenvals(self):
        # Test that eigenvalues are set to ones
        eigvals = lra_two.getEigvals()
        assert_array_equal(eigvals, np.ones(lra_two._n_modes))

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def testNumModesTwoClasses(self):
        # Test number of modes
        expected_modes = len(set(TWO_CLASSES))
        if expected_modes == 2:
            expected_modes = 1
        self.assertEqual(lra_two._n_modes, expected_modes)

    @unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn not available")
    def testNumModesThreeClasses(self):
        # Since n_modes is automatic in LRA from classes, check it sets correctly
        expected_modes = len(set(THREE_CLASSES))
        if expected_modes == 2:
            expected_modes = 1
        self.assertEqual(lra_three._n_modes, expected_modes)

if __name__ == '__main__':
    unittest.main()
