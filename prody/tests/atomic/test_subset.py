import numpy
from numpy.testing import *

from prody.tests.datafiles import *

from prody.tests import unittest, TEMPDIR

ATOMS = parseDatafile('1ubi')
ATOMS.setData('ones', numpy.zeros(len(ATOMS)))

class TestUserData(unittest.TestCase):

    def testSetData(self):

        ATOMS[:10].setData('ones', 1)
        ATOMS[10:].setData('ones', 1)
        assert_equal(ATOMS.getData('ones'), numpy.ones(len(ATOMS)))
