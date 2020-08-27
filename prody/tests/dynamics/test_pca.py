"""This module contains unit tests for :mod:`~prody.dynamics`."""

import numpy as np
from numpy import arange
from numpy.testing import *
try:
    import numpy.testing.decorators as dec
except ImportError:
    from numpy.testing import dec

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

ATOL = 1e-5
RTOL = 0

ATOMS = parseDatafile('2k39_ca')
COORDSETS = ATOMS.getCoordsets()

PCA_COV = parseDatafile('pca2k39_cov', symmetric=True)
PCA_EVALUES = parseDatafile('pca2k39_evalues')
PCA_EVECTORS = parseDatafile('pca2k39_vectors')

pca = PCA()
pca.buildCovariance(ATOMS)
pca.calcModes(n_modes=None)

n_modes = 20
pca_n_modes = PCA()
pca_n_modes.buildCovariance(ATOMS)
pca_n_modes.calcModes(n_modes=n_modes)

class TestPCA(unittest.TestCase):

    def setUp(self):

        self.model = PCA()
        self.buildMatrix = self.model.buildCovariance
        self.setMatrix = self.model.setCovariance
        self.getMatrix = self.model.getCovariance
        self.getExpected = pca.getCovariance

    def testBuildMatrixCoordsWrongType(self):
        """Test response to wrong type *coords* argument."""

        self.assertRaises(TypeError, self.buildMatrix, 'nogood')

    def testBuildMatrixWrongCoords(self):
        """Test response to wrong coords.dtype."""

        array = np.array([['a','a','a'] for i in range(10)])
        self.assertRaises(ValueError, self.buildMatrix, array)

    def testBuildMatrixCoordsArray(self):
        """Test output when *coords* is an array."""

        self.buildMatrix(COORDSETS)
        assert_equal(self.getMatrix(), self.getExpected(),
                     'failed to get correct matrix')

    def testSetMatrixWrongType(self):
        """Test response to wrong matrix type argument."""

        self.assertRaises(TypeError, self.setMatrix, list(np.ones((3,3))))

    def testSetMatrixWrongDim(self):
        """Test response to wrong dim *kirchhoff* argument."""

        self.assertRaises(ValueError, self.setMatrix, np.ones((3,4,3)))

    def testSetMatrixNonSquare(self):
        """Test response to non-square matrix."""

        self.assertRaises(ValueError, self.setMatrix, np.ones((3,4)))

    def testSetMatrixWrongDtype(self):
        """Test response to wrong matrix.dtype."""

        array = np.array([['a','a','a'] for i in range(3)])
        self.assertRaises(ValueError, self.setMatrix, array)

    def testSetMatrixAcceptableDtype(self):
        """Test response to acceptable matrix.dtype."""

        self.assertIsNone(self.setMatrix(np.ones((30, 30), int)),
                          'failed to set an acceptable array')


class TestPCAResults(TestPCA):

    def testEigenvalues(self):
        assert_allclose(pca.getEigvals(), PCA_EVALUES,
                        rtol=RTOL, atol=ATOL*100,
                        err_msg='failed to get correct eigenvalues')

    def testEigenvectors(self):
        _temp = np.abs(np.dot(pca.getEigvecs().T, PCA_EVECTORS))
        assert_allclose(_temp, np.eye(_temp.shape[0]), rtol=RTOL, atol=ATOL*10,
                       err_msg='failed to get correct eigenvectors')

    def testNumModesAll(self):
        assert_equal(pca._n_modes, len(PCA_EVALUES))

    def testNumModesSet(self):
        assert_equal(pca_n_modes._n_modes, n_modes)

    def testCovariance(self):
        assert_allclose(pca.getCovariance(), PCA_COV,
                        rtol=0, atol=ATOL,
                        err_msg='failed to get correct Covariance matrix')

    def testCovarianceSymmetry(self):
        cov = pca.getCovariance()
        assert_equal(cov, cov.T, 'Covariance is not symmetric')

if __name__ == '__main__':
    unittest.main()
