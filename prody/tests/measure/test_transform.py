"""This module contains unit tests for :mod:`prody.measure.transform` module.
"""

from numpy import array, zeros, ones, eye, all, cos, sin
from numpy.random import RandomState
from numpy.testing import assert_equal, assert_allclose
from scipy.spatial.transform import Rotation

from prody.ensemble import Ensemble
from prody.tests import unittest
from prody.tests.datafiles import parseDatafile, pathDatafile
from prody.trajectory import Trajectory

from prody.measure import (
    applyTransformation,
    calcRMSD,
    calcTransformation,
    moveAtoms,
    wrapAtoms,
)

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


class TestWeightedSuperpose(unittest.TestCase):

    def _make_weighted_problem(self):
        rng = RandomState(42)
        target = rng.normal(size=(10, 3))
        theta = 0.7
        rot_true = array([
            [cos(theta), -sin(theta), 0.0],
            [sin(theta), cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ])
        mobile = target.dot(rot_true.T) + rng.normal(scale=0.5, size=(10, 3))
        weights = array([1., 1., 1., 1., 1., 10., 10., 10., 10., 100.])

        w = weights.reshape(-1, 1)
        mob_com = (mobile * w).sum(axis=0) / w.sum()
        tar_com = (target * w).sum(axis=0) / w.sum()
        rot_opt, _ = Rotation.align_vectors(
            target - tar_com, mobile - mob_com, weights=weights
        )
        aligned_opt = (mobile - mob_com).dot(rot_opt.as_matrix().T) + tar_com
        rmsd_opt = calcRMSD(target, aligned_opt, weights=weights)
        return mobile, target, weights, rot_opt.as_matrix(), rmsd_opt

    def testCalcTransformationNonUniformWeights(self):
        mobile, target, weights, rot_opt, rmsd_opt = self._make_weighted_problem()
        t = calcTransformation(mobile, target, weights=weights)
        aligned = applyTransformation(t, mobile.copy())
        assert_allclose(t.getRotation(), rot_opt, atol=1e-6)
        assert_allclose(calcRMSD(target, aligned, weights=weights), rmsd_opt, atol=1e-6)

    def testEnsembleSuperposeNonUniformWeights(self):
        mobile, target, weights, _, rmsd_opt = self._make_weighted_problem()
        ens = Ensemble('weighted')
        ens.setCoords(target)
        ens.addCoordset(mobile)
        ens.setWeights(weights)
        ens.superpose()
        assert_allclose(ens.getRMSDs()[0], rmsd_opt, atol=1e-6)

    def testFrameSuperposeNonUniformWeights(self):
        traj = Trajectory(pathDatafile('dcd'))
        pdb = parseDatafile('multi_model_truncated', model=1)
        ref = pdb.getCoords()
        weights = ones(len(ref))
        weights[:10] = 50.0
        traj.setCoords(ref)
        traj.setWeights(weights)
        traj.reset()

        w = weights.reshape(-1, 1)
        tar_com = (ref * w).sum(axis=0) / w.sum()
        for frame in traj:
            mob = frame.getCoords().copy()
            mob_com = (mob * w).sum(axis=0) / w.sum()
            rot_opt, _ = Rotation.align_vectors(
                ref - tar_com, mob - mob_com, weights=weights
            )
            aligned_opt = (mob - mob_com).dot(rot_opt.as_matrix().T) + tar_com
            rmsd_opt = calcRMSD(ref, aligned_opt, weights=weights)

            frame.superpose()
            assert_allclose(frame.getRMSD(), rmsd_opt, atol=1e-5)


