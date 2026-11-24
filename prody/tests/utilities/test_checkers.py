from prody.tests import TestCase

from numpy import arange, array, zeros
from numpy.random import random
from numpy.testing import assert_equal

from prody.utilities import checkBlocks, checkCoords, checkTypes

COORDS = random((10, 3))*10
COORDSET = random((2, 10, 3))*10

NNODES = 4
BLOCKS = [0, 0, 1, 1]

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


class TestCheckBlocks(TestCase):

    def testBlocksList(self):

        self.assertIsNone(checkBlocks(BLOCKS, NNODES))

    def testBlocksTuple(self):

        self.assertIsNone(checkBlocks(tuple(BLOCKS), NNODES))

    def testBlocksArray(self):

        self.assertIsNone(checkBlocks(array(BLOCKS), NNODES))

    def testBlocksNotNumbers(self):

        self.assertIsNone(checkBlocks(['a', 'a', 'b', 'b'], NNODES))

    def testBlocksOnePerNode(self):

        self.assertIsNone(checkBlocks(arange(NNODES), NNODES))

    def testInvalidBlocks(self):

        self.assertRaises(TypeError, checkBlocks, None, NNODES)

    def testScalarBlocks(self):

        self.assertRaises(TypeError, checkBlocks, 0, NNODES)

    def testStringBlocks(self):

        self.assertRaises(TypeError, checkBlocks, 'aabb', NNODES)

    def testNestedBlocks(self):

        self.assertRaises(ValueError, checkBlocks, [[0, 1], [1, 0]], NNODES)

    def testBlocks2D(self):

        self.assertRaises(ValueError, checkBlocks, zeros((NNODES, 2), int),
                          NNODES)

    def testTooFewBlocks(self):

        self.assertRaises(ValueError, checkBlocks, BLOCKS[:-1], NNODES)

    def testTooManyBlocks(self):

        self.assertRaises(ValueError, checkBlocks, BLOCKS + [1], NNODES)

    def testOneBlock(self):

        self.assertRaises(ValueError, checkBlocks, zeros(NNODES, int), NNODES)

    def testTwoBlocksOfOneNode(self):
        # two blocks of one node have 3 degrees of freedom each, which is no
        # more than the 6 of a rigid body

        self.assertRaises(ValueError, checkBlocks, [0, 1], 2)

    def testMixedTypeBlocks(self):

        self.assertRaises(TypeError, checkBlocks, [0, 0, '1', '1'], NNODES)

    def testRaggedBlocks(self):

        self.assertRaises(ValueError, checkBlocks, [[0, 1], [2]], 2)

    def testNodesLabel(self):

        self.assertRaisesRegex(ValueError, 'coarse grained nodes',
                               checkBlocks, BLOCKS[:-1], NNODES,
                               'coarse grained nodes')


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
