from prody.tests import TestCase

from numpy.random import random
from numpy.testing import assert_equal

from prody.utilities import checkCoords, checkTypes

COORDS = random((10, 3))*10
COORDSET = random((2, 10, 3))*10

class TestCheckCoords(TestCase):

    def testInvalidCoords(self):

        self.assertRaises(TypeError, checkCoords, [None])

    def testCoords(self):

        self.assertTrue(checkCoords(COORDS))

    def testCoordset(self):

        self.assertTrue(checkCoords(COORDSET, csets=True))


    def testCoordsetNatoms(self):

        self.assertRaises(ValueError, checkCoords, COORDSET, csets=True,
                          natoms=20)


class testCheckTypes(TestCase):

    def testCorrectMonotypeOneArg(self):

        self.assertTrue(checkTypes({'i': 1}, i=int))

    def testCorrectMonotypeTwoArgs(self):

        self.assertTrue(checkTypes({'i': 1, 'n': 10}, i=int, n=int))

    def testCorrectMultitypeOneArg(self):

        self.assertTrue(checkTypes({'i': 1.}, i=(int, float)))

    def testCorrectMonotypeTwoArgs(self):

        self.assertTrue(checkTypes({'i': 1, 'n': 10.}, i=(int, float),
                                                       n=(int, float)))

    def testWrongMonotypeOneArg(self):

        self.assertRaises(TypeError, checkTypes, {'i': 1.0}, i=int)

    def testWrongMonotypeTwoArgs(self):

        self.assertRaises(TypeError, checkTypes, {'i': 1, 'n': 10},
                                                  i=int, n=int)

    def testWrongMultitypeOneArg(self):

        self.assertRaises(TypeError, checkTypes, {'i': '1.'},
                                      i=(int, float))

    def testWrongMonotypeTwoArgs(self):

        self.assertRaises(TypeError, checkTypes, {'i': 1, 'n': '10.'},
                            i=(int, float), n=(int, float))
