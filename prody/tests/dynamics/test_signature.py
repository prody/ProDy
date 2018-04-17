"""This module contains unit tests for :mod:`~prody.KDTree` module."""

from numpy.testing import assert_array_equal, assert_equal
from numpy.random import rand, randint

from prody.dynamics import sdarray

from prody.tests import unittest
from prody.tests.datafiles import parseDatafile

from prody import _PY3K, LOGGER

LOGGER.verbosity = 'none'


if not _PY3K:
    range = xrange

size = (10, 300, 150)
A = rand(*size)
W = randint(0, 2, size)
labels = [str(i) for i in range(size[0])]

S = sdarray(A, W, labels=labels)

class TestSDArray(unittest.TestCase):

    def testSlicing(self):
        s = S[0]
        assert_equal(s.shape, A[0].shape, 'failed at sdarray slicing')
        assert_array_equal(s.flatten(), A[0].flatten(), 'failed at sdarray slicing')
        assert_array_equal(s.getWeights().flatten(), W[0].flatten(), 'failed at sdarray slicing')

        s = S[:, 0]
        assert_equal(s.shape, A[:, 0].shape, 'failed at sdarray slicing')
        assert_array_equal(s.flatten(), A[:, 0].flatten(), 'failed at sdarray slicing')
        assert_array_equal(s.getWeights().flatten(), W[:, 0].flatten(), 'failed at sdarray slicing')

        s = S[:, :, :]
        assert_equal(s.shape, A.shape, 'failed at sdarray slicing')
        assert_array_equal(s.flatten(), A.flatten(), 'failed at sdarray slicing')
        assert_array_equal(s.getWeights().flatten(), W.flatten(), 'failed at sdarray slicing')

        s = S[0, 0, 0]
        assert_array_equal(s, A[0, 0, 0], 'failed at sdarray slicing')
