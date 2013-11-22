"""This module contains unit tests for :mod:`.frame` module."""

from code import interact
from prody.tests import TestCase

from numpy import array
from numpy.testing import assert_allclose

from prody.trajectory import Trajectory
from prody.tests.datafiles import parseDatafile, pathDatafile

DCD = Trajectory(pathDatafile('dcd'))
PDB = parseDatafile('multi_model_truncated', model=1)
RMSD_ALL = array([0.0, 1.380, 1.745])
RMSD_CARBON = array([0.0, 0.964, 1.148])


class TestSuperpose(TestCase):

    def setUp(self):

        DCD.setCoords(PDB.getCoords())
        DCD.reset()

    def testAll(self):

        rmsd = []
        for frame in DCD:
            frame.superpose()
            rmsd.append(frame.getRMSD())
        assert_allclose(rmsd, RMSD_ALL, atol=0.001)

    def testCarbon(self):

        DCD.setAtoms(PDB.carbon)
        rmsd = []
        for frame in DCD:
            frame.superpose()
            rmsd.append(frame.getRMSD())
        assert_allclose(rmsd, RMSD_CARBON, atol=0.001)
