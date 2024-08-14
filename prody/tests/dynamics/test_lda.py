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
ENSEMBLE = refineEnsemble(PDBEnsemble(ATOMS), lower=0., upper=5.)

RMSD_MATRIX = ENSEMBLE.getRMSDs(pairwise=True)
RMSD_TREE = calcTree(ENSEMBLE.getLabels(), RMSD_MATRIX)
SUBGROUPS = findSubgroups(RMSD_TREE, 2.9)
CLASSES = [np.nonzero(np.array([i in sg for sg in SUBGROUPS]))[0][0] 
           for i in ENSEMBLE.getLabels()]

LDA_EVALUES = parseDatafile('lda2k39_evalues')
LDA_EVECTORS = parseDatafile('lda2k39_vectors')

lda = LDA()
lda.calcModes(ENSEMBLE, CLASSES)

n_modes = 1
lda_n_modes = LDA()
lda_n_modes.calcModes(ENSEMBLE, CLASSES, n_modes=n_modes)

class TestLDA(unittest.TestCase):

    def setUp(self):

        self.model = LDA()
        self.calcModes = self.model.calcModes

    def testCalcCoordsWrongType(self):
        """Test response to wrong type *coords* argument."""

        self.assertRaises(TypeError, self.calcModes, 'nogood', CLASSES)

    def testCalcWrongCoords(self):
        """Test response to wrong coords.dtype."""

        array = np.array([['a','a','a'] for _ in range(10)])
        self.assertRaises(ValueError, self.calcModes, array, CLASSES)

    def testCalcLabelsWrongType(self):
        """Test response to wrong type *coords* argument."""

        self.assertRaises(TypeError, self.calcModes, ENSEMBLE, 'nogood')


class TestLDAResults(TestLDA):

    def testEigenvalues(self):
        assert_allclose(lda.getEigvals(), LDA_EVALUES,
                        rtol=RTOL, atol=ATOL*100,
                        err_msg='failed to get correct eigenvalues')

    def testEigenvectors(self):
        _temp = np.abs(np.dot(lda.getEigvecs().T, LDA_EVECTORS))
        assert_allclose(_temp, np.eye(_temp.shape[0]), 
                        rtol=1e-10, atol=0.1, # may not be so close to orthonormal
                       err_msg='failed to get correct eigenvectors')

    def testNumModesAll(self):
        assert_equal(lda._n_modes, len(LDA_EVALUES))

    def testNumModesSet(self):
        assert_equal(lda_n_modes._n_modes, n_modes)

if __name__ == '__main__':
    unittest.main()
