"""This module contains unit tests for :mod:`prody.measure.transform` module.
"""

from numpy import zeros, ones, eye, all
from numpy.testing import assert_equal

from prody.tests import unittest
from prody.tests.datafiles import parseDatafile

from prody.measure import moveAtoms, wrapAtoms

UBI = parseDatafile('1ubi')


class TestMoveAtoms(unittest.TestCase):

    def testToArgument(self):

        atoms = UBI.ca

        coords = UBI.getCoords()
        center = atoms._getCoords().mean(0)
        moveAtoms(atoms, to=zeros(3), ag=True)
        assert_equal(UBI._getCoords(), coords - center)

    def testByArgument(self):

        atoms = UBI
        offset = ones(3) * 10.
        coords = atoms.getCoords()
        moveAtoms(atoms, by=offset)
        assert_equal(atoms._getCoords(), coords + offset)

    def testTransformation(self):

        atoms = UBI.ca
        coords = UBI.getCoords()
        matrix = eye(4)
        moveAtoms(atoms, by=matrix, ag=True)
        assert_equal(UBI._getCoords(), coords)


class TestWrapAtoms(unittest.TestCase):

    def testWrap(self):

        xyz = UBI.getCoords()
        xyzmin = xyz.min(0) - 0.5
        unitcell = xyz.max(0) - xyzmin + 0.5
        center = xyzmin - unitcell / 2
        wrapAtoms(UBI, unitcell, center)
        diff = xyz - UBI.getCoords()
        self.assertTrue(all(diff == unitcell))

