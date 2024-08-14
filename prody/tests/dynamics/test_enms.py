"""This module contains unit tests for :mod:`~prody.dynamics`."""

import numpy as np
from numpy import arange
from numpy.testing import *

from prody.utilities import importDec
dec = importDec()

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

ATOL = 1e-5
RTOL = 0


ATOMS = parseDatafile('1ubi_ca')
COORDS = ATOMS.getCoords()

ANM_HESSIAN = parseDatafile('anm1ubi_hessian', symmetric=True)
ANM_EVALUES = parseDatafile('anm1ubi_evalues')[:,1].flatten()
ANM_EVECTORS = parseDatafile('anm1ubi_vectors')[:,1:]

GNM_KIRCHHOFF = parseDatafile('gnm1ubi_kirchhoff', symmetric=True, skiprows=1)
GNM_EVALUES = parseDatafile('gnm1ubi_evalues')[:,1].flatten()
GNM_EVECTORS = parseDatafile('gnm1ubi_vectors', usecols=arange(3,23))


anm = ANM()
anm.buildHessian(ATOMS)
anm.calcModes(n_modes=None, zeros=True)


gnm = GNM()
gnm.buildKirchhoff(ATOMS)
gnm.calcModes(n_modes=None, zeros=True)

ATOMS2 = parseDatafile('2gb1_truncated')
rtb = RTB()
rtb.buildHessian(ATOMS2, ATOMS2.getBetas().astype(int))
RTB_HESSIAN = parseDatafile('rtb2gb1_hessian', symmetric=True)
RTB_PROJECT = parseDatafile('rtb2gb1_project')

HITTIME = parseDatafile('hit1ubi')
COMMUTETIME = parseDatafile('commute1ubi')


ATOMS3 = parseDatafile('2nwl')
COORDS3 = ATOMS3.ca.getCoords()

EXANM_HESSIAN = parseDatafile('exanm2nwl_hessian.coo', symmetric=True)
EXANM_EVALUES = parseDatafile('exanm2nwl_evalues.dat')[:,1].flatten()
EXANM_EVECTORS = parseDatafile('exanm2nwl_vectors.dat')[:,1:]

exanm = ANM()
exanm.buildHessian(COORDS3)
exanm.calcModes(n_modes=None, zeros=True)

class testGNMBase(unittest.TestCase):

    def setUp(self):

        self.model = gnm


    def testGetCutoff(self):
        """Testing return type of :meth:`~.GNMBase.getCutoff`."""

        self.assertIsInstance(self.model.getCutoff(), float,
                              'getCutoff failed to return a float')

    def testGetGamma(self):
        """Testing return type of :meth:`~.GNMBase.getGamma`."""

        self.assertIsInstance(self.model.getGamma(), float,
                              'getCutoff failed to return a float')


class TestANMResults(testGNMBase):

    def testEigenvalues(self):
        """Test eigenvalues."""

        assert_allclose(anm[:len(ANM_EVALUES)].getEigvals(), ANM_EVALUES,
                        rtol=RTOL, atol=ATOL*10,
                        err_msg='failed to get correct eigenvalues')


    def testEigenvectors(self):
        """Test eigenvectors."""

        _temp = np.abs(np.dot(anm[6:6+ANM_EVECTORS.shape[1]].getEigvecs().T,
                              ANM_EVECTORS))
        assert_allclose(_temp, np.eye(20), rtol=RTOL, atol=ATOL,
                        err_msg='failed to get correct eigenvectors')

    def testHessian(self):
        """Test Hessian matrix."""

        assert_allclose(anm.getHessian(), ANM_HESSIAN, rtol=0, atol=ATOL,
                        err_msg='failed to get correct Hessian matrix')


    def testVariances(self):
        """Test variances."""

        assert_allclose(anm[6:len(ANM_EVALUES)].getVariances(),
                        1/ANM_EVALUES[6:],
                        rtol=0, atol=ATOL*100,
                        err_msg='failed to get correct variances')

    def testGetHessian(self):
        assert_equal(anm.getHessian(), anm._getHessian(),
                     err_msg='failed to _get correct Hessian matrix')

    def testHessianSymmetry(self):
        hessian = anm._getHessian()
        assert_equal(hessian, hessian.T, 'hessian is not symmetric')

    def testHessianSums(self):
        hessian = anm._getHessian()
        zeros = np.zeros(hessian.shape[0])
        assert_allclose(hessian.sum(0), zeros,
                        rtol=0, atol=ATOL,
                        err_msg='hessian columns do not add up to zero')
        assert_allclose(hessian.sum(1), zeros,
                        rtol=0, atol=ATOL,
                        err_msg='hessian rows do not add up to zero')
        
    def testModeSlicing(self):
        indices1 = slice(227, 230)
        assert_equal(anm[indices1].numModes(), 1,
                         'slicing ANM for last mode does not work')
        indices = slice(228, 230)
        self.assertRaises(IndexError, anm.__getitem__, indices)

class TestEXANMResults(testGNMBase):

    def testEigenvalues(self):
        """Test eigenvalues."""

        assert_allclose(exanm[:len(EXANM_EVALUES)].getEigvals(), EXANM_EVALUES,
                        rtol=RTOL, atol=ATOL*10,
                        err_msg='failed to get correct eigenvalues')


    def testEigenvectors(self):
        """Test eigenvectors."""

        _temp = np.abs(np.dot(exanm[6:6+EXANM_EVECTORS.shape[1]].getEigvecs().T,
                              EXANM_EVECTORS))
        assert_allclose(_temp, np.eye(20), rtol=RTOL, atol=ATOL,
                        err_msg='failed to get correct eigenvectors')

    def testHessian(self):
        """Test Hessian matrix."""

        assert_allclose(exanm.getHessian(), EXANM_HESSIAN, rtol=0, atol=ATOL,
                        err_msg='failed to get correct Hessian matrix')


    def testVariances(self):
        """Test variances."""

        assert_allclose(exanm[6:len(EXANM_EVALUES)].getVariances(),
                        1/EXANM_EVALUES[6:],
                        rtol=0, atol=ATOL*100,
                        err_msg='failed to get correct variances')

    def testGetHessian(self):
        assert_equal(exanm.getHessian(), exanm._getHessian(),
                     err_msg='failed to _get correct Hessian matrix')

    def testHessianSymmetry(self):
        hessian = exanm._getHessian()
        assert_equal(hessian, hessian.T, 'hessian is not symmetric')

    def testHessianSums(self):
        hessian = exanm._getHessian()
        zeros = np.zeros(hessian.shape[0])
        assert_allclose(hessian.sum(0), zeros,
                        rtol=0, atol=ATOL,
                        err_msg='hessian columns do not add up to zero')
        assert_allclose(hessian.sum(1), zeros,
                        rtol=0, atol=ATOL,
                        err_msg='hessian rows do not add up to zero')

'''
class TestANMSparse(unittest.TestCase):

    """Test result from using sparse matrices."""

    @dec.slow
    @unittest.skipIf(True, 'not completed')
    def testSparse(self):

        anm = ANM()
        anm.buildHessian(COORDS, sparse=True)
        assert_allclose(anm.getHessian().toarray(), ANM_HESSIAN,
                        rtol=0, atol=ATOL,
                        err_msg='failed to get correct sparse Hessian matrix')
        anm.calcModes(None)
        assert_allclose(anm[:len(ANM_EVALUES)].getEigvals(), ANM_EVALUES,
                        rtol=RTOL, atol=ATOL*10,
                        err_msg='failed to get correct eigenvalues')
        _temp = np.abs(np.dot(anm[6:6+ANM_EVECTORS.shape[1]].getEigvecs().T,
                              ANM_EVECTORS))
        assert_allclose(_temp, np.eye(20), rtol=RTOL, atol=ATOL,
                        err_msg='failed to get correct eigenvectors')
'''

class TestGNMResults(testGNMBase):

    def testEigenvalues(self):
        assert_allclose(gnm[:21].getEigvals(), GNM_EVALUES[:21],
                        rtol=RTOL, atol=ATOL*100,
                        err_msg='failed to get correct slow eigenvalues')

        assert_allclose(gnm[-21:].getEigvals(), GNM_EVALUES[21:],
                        rtol=RTOL, atol=ATOL*100,
                        err_msg='failed to get correct fast eigenvalues')

    def testEigenvectors(self):
        _temp = np.abs(np.dot(gnm[1:21].getEigvecs().T, GNM_EVECTORS))
        assert_allclose(_temp, np.eye(20), rtol=RTOL, atol=ATOL*10,
                       err_msg='failed to get correct eigenvectors')

    def testKirchhoff(self):
        assert_allclose(gnm._getKirchhoff(), GNM_KIRCHHOFF,
                        rtol=0, atol=ATOL,
                        err_msg='failed to get correct Kirchhoff matrix')

    def testGetKirchoff(self):
        assert_equal(gnm.getKirchhoff(), gnm._getKirchhoff(),
                     err_msg='failed to _get correct Kirchhoff matrix')

    def testKirchhoffSymmetry(self):
        kirchhoff = gnm._getKirchhoff()
        assert_equal(kirchhoff, kirchhoff.T, 'kirchhoff is not symmetric')

    def testKirchhoffSums(self):
        kirchhoff = gnm._getKirchhoff()
        zeros = np.zeros(kirchhoff.shape[0])
        assert_equal(kirchhoff.sum(0), zeros,
                     'kirchhoff columns do not add up to zero')
        assert_equal(kirchhoff.sum(1), zeros,
                     'kirchhoff rows do not add up to zero')

    def testBuildKirchoffSlow(self):
        slow = GNM()
        slow.buildKirchhoff(ATOMS)
        assert_equal(slow._getKirchhoff(), gnm._getKirchhoff(),
                     'slow method does not reproduce same Kirchhoff')

    def testCommuteTime(self):
        hitTime, commuteTime = calcHitTime(gnm)

        assert_allclose(hitTime, HITTIME,
                        rtol=0, atol=ATOL,
                        err_msg='failed to get correct hit times')
        assert_allclose(commuteTime, COMMUTETIME,
                        rtol=0, atol=ATOL,
                        err_msg='failed to get correct commute times')             

class TestGNM(unittest.TestCase):

    def setUp(self):

        self.model = GNM()
        self.buildMatrix = self.model.buildKirchhoff
        self.setMatrix = self.model.setKirchhoff
        self.getMatrix = self.model.getKirchhoff
        self.getExpected = gnm.getKirchhoff

    def testBuildMatrixCoordsWrongType(self):
        """Test response to wrong type *coords* argument."""

        self.assertRaises(TypeError, self.buildMatrix, 'nogood')

    def testBuildMatrixWrongCoords(self):
        """Test response to wrong coords.dtype."""

        array = np.array([['a','a','a'] for i in range(10)])
        self.assertRaises(ValueError, self.buildMatrix, array)

    def testBuildMatrixCoordsArray(self):
        """Test output when  *coords* is an array."""

        self.buildMatrix(COORDS)
        assert_equal(self.getMatrix(), self.getExpected(),
                     'failed to get correct matrix')

    def testBuildMatrixCutoffWrongType(self):
        """Test response to cutoff of wrong type."""

        self.assertRaises(TypeError, self.buildMatrix, COORDS, 'none')

    def testBuildMatrixCutoffInvalidValue(self):
        """Test response to wrong *cutoff* argument."""

        self.assertRaises(ValueError, self.buildMatrix, COORDS, -1)

    def testBuildMatrixInvalidGamma(self):
        """Test response to invalid *gamma* argument."""

        self.assertRaises(TypeError, self.buildMatrix, COORDS, gamma='none')

    def testBuildMatrixWrongGamma(self):
        """Test response to wrong *gamma* argument."""

        self.assertRaises(ValueError, self.buildMatrix, COORDS, gamma=0)

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

        self.assertIsNone(self.setMatrix(np.ones((30,30), int)),
                          'failed to set an acceptable array')

class TestANM(TestGNM):

    def setUp(self):

        self.model = ANM()
        self.anm = self.model
        self.buildMatrix = self.model.buildHessian
        self.getMatrix = self.model.getHessian
        self.setMatrix = self.model.setHessian
        self.getExpected = anm.getHessian

    def testSetHessianWrongShape(self):
        """Test response to wrong shape *hessian* argument."""

        self.assertRaises(ValueError, self.model.setHessian, np.ones((5,5)))

    def testBuildHessianSlow(self):
        slow = ANM()
        slow.buildHessian(ATOMS, slow=True)
        assert_allclose(slow._getHessian(), anm._getHessian(),
                        rtol=0, atol=ATOL,
                        err_msg='slow method does not reproduce same Hessian')
        assert_equal(slow._getKirchhoff(), anm._getKirchhoff(),
                     'slow method does not reproduce same Kirchhoff')
    

class TestGNMCalcModes(unittest.TestCase):

    def setUp(self):
        pass

class TestRTB(unittest.TestCase):

    def testHessian(self):

        assert_allclose(RTB_HESSIAN, rtb._getHessian(),
                        rtol=0, atol=ATOL,
                        err_msg='expected Hessian is not produced')

    def testProjection(self):

        assert_allclose(RTB_PROJECT, rtb._getProjection(),
                        rtol=0, atol=ATOL,
                        err_msg='expected projection matrix is not produced')


    def testCalcModes(self):

        rtb.calcModes()

if __name__ == '__main__':
    unittest.main()
