"""This module contains unit tests for :mod:`~prody.dynamics.pathtools`."""

import numpy as np
from numpy.testing import assert_allclose

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

# adenylate kinase chain A, open (4AKE, all atoms with hydrogens) and closed (1AKE), matching Calpha atoms
OPEN_ALL = parsePDB(pathDatafile('pdb4ake_fixed'), chain='A')
OPEN = OPEN_ALL.select('calpha')
CLOSED = parsePDB(pathDatafile('pdb1ake'), chain='A').select('calpha')
A, B = OPEN.getCoords(), CLOSED.getCoords()


def rotation(angle, axis):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def line(S, T, n=11):
    """Straight-line path from S to T superposed onto S."""
    T2 = applyTransformation(calcTransformation(T, S), T.copy())
    return np.array([(1 - t) * S + t * T2 for t in np.linspace(0, 1, n)])


class TestMidpoints(unittest.TestCase):

    def testDropsOnlyMidpoints(self):
        L = line(A, B, 5)
        with_mid = np.array([L[0], 0.5 * (L[0] + L[2]), L[2], 0.5 * (L[2] + L[4]), L[4]])
        kept, idx = dropMidpoints(with_mid, return_indices=True)
        self.assertEqual(list(idx), [0, 2, 4])
        assert_allclose(kept, L[[0, 2, 4]])

    def testNeverDropsTwoConsecutiveFrames(self):
        L = line(A, B, 9)                                 # every interior frame is a midpoint of its neighbours
        self.assertEqual(list(dropMidpoints(L, return_indices=True)[1]), [0, 2, 4, 6, 8])

    def testKeepsPathsWithoutMidpoints(self):
        P = np.array([A, A + 1.0, A + np.random.RandomState(0).rand(*A.shape)])
        self.assertEqual(len(dropMidpoints(P)), 3)


class TestPathGeometry(unittest.TestCase):

    def testEndpointsAreAdmissibleAndRMSDIsRight(self):
        g = calcPathGeometry(np.array([A, B]), target=B)
        self.assertTrue(g['admissible'].all())
        self.assertAlmostEqual(g['rmsd_target'][1], 0.0, places=5)
        self.assertAlmostEqual(g['rmsd_start'][0], 0.0, places=5)
        self.assertAlmostEqual(g['rmsd_target'][0], calcRMSD(applyTransformation(calcTransformation(A, B), A.copy()), B),
                               places=5)

    def testDetectsAStretchedBond(self):
        X = A.copy()
        X[100:] += (X[100] - X[99]) / np.linalg.norm(X[100] - X[99]) * 1.0      # bond 99-100 stretched by 1 A
        g = calcPathGeometry(np.array([A, X]), target=B)
        self.assertGreater(g['ceiling_excess'][1], 0.5)
        self.assertFalse(g['admissible'][1])


class TestStepWork(unittest.TestCase):

    def testRigidMotionCostsNothing(self):
        R = rotation(0.4, [1, 2, 3])
        P = np.array([A, A @ R.T + [3.0, -2.0, 1.0]])
        w = calcStepWork(P)
        self.assertLess(w['total_work'], 1e-6)

    def testWorkMatchesHessian(self):
        P = line(A, B, 3)
        w = calcStepWork(P, target=B, midpoints=True)     # a straight line looks like a path with midpoints
        anm = ANM(); anm.buildHessian(P[0], cutoff=15.0, gamma=1.0)
        D = (applyTransformation(calcTransformation(P[1], P[0]), P[1].copy()) - P[0]).ravel()
        self.assertAlmostEqual(w['work'][0], 0.5 * D @ anm.getHessian() @ D, places=6)
        self.assertAlmostEqual(w['stiffness'][0], w['work'][0] / (D @ D), places=8)
        self.assertAlmostEqual(w['cumulative_share'][-1], 1.0, places=6)
        self.assertGreater(w['rmsd_gain'].sum(), 0)

    def testSoftModeStepIsSoft(self):
        """A step along the softest mode has stiffness lambda_1 / 2 and all its motion in the softest mode."""
        anm = ANM(); anm.buildHessian(A, cutoff=15.0, gamma=1.0); anm.calcModes(n_modes=1)
        u = anm.getEigvecs()[:, 0].reshape(-1, 3)
        w = calcStepWork(np.array([A, A + 0.01 * u]))
        self.assertAlmostEqual(w['stiffness'][0], anm.getEigvals()[0] / 2, places=4)
        self.assertGreater(w['share_soft5'][0], 0.999)


class TestDeformationModes(unittest.TestCase):

    def testSharesSumToOne(self):
        for model in ('anm', 'genanm'):
            d = calcDeformationModes(OPEN, CLOSED, model=model)
            self.assertAlmostEqual(d['captured'], 1.0, places=6)
            self.assertTrue(d['n50'] <= d['n80'] <= d['n90'])
            self.assertGreaterEqual(d['n_eff'], 1.0)

    def testSingleModeDeformation(self):
        gen = genANM(); gen.buildHessian(A); gen.calcModes(n_modes=10)
        u = gen.getEigvecs()[:, 6].reshape(-1, 3)
        d = calcDeformationModes(A, A + 0.05 * u, model='genanm')
        self.assertEqual(d['max_mode'], 7)
        self.assertAlmostEqual(d['max_overlap'], 1.0, places=4)
        self.assertAlmostEqual(d['n_eff'], 1.0, places=3)


class TestRoutesAndApproach(unittest.TestCase):

    def testRouteDistance(self):
        L = line(A, B, 11)
        self.assertAlmostEqual(calcRouteDistance(L, L)['distance'], 0.0, places=5)
        self.assertAlmostEqual(calcRouteDistance(L, L[::-1])['distance'], 0.0, places=5)
        self.assertGreater(calcRouteDistance(L, L[:3])['distance'], 0.5)

    def testClosestApproachOfTheLineToItsMidpoint(self):
        L = line(A, B, 21)
        c = calcClosestApproach(L, L[10], target=B)
        self.assertLess(c['d_path'], 1e-3)
        self.assertLess(c['d_line'], 0.05)
        self.assertAlmostEqual(c['position_path'], 0.5, places=2)
        self.assertEqual(c['outcome'], 'tie')

    def testSubsetIndices(self):
        L = line(A, B, 21)
        idx = np.arange(20, 150)
        c = calcClosestApproach(L, L[10][idx], target=B, indices=idx)
        self.assertLess(c['d_path'], 1e-3)

    def testPCAOfAStraightLine(self):
        pca = calcPathPCA(line(A, B, 11), structures=[A, B])
        self.assertGreater(pca['variance'][0], 0.999)


class TestAllAtomPath(unittest.TestCase):

    def testFirstFrameReproducesTheReference(self):
        ref = OPEN_ALL.select('protein')
        for method in ('frame', 'translate'):
            aa = buildAllAtomPath(line(A, B, 3), OPEN_ALL, method=method)
            self.assertEqual(aa.numAtoms(), ref.numAtoms())
            assert_allclose(aa.getCoordsets(0), ref.getCoords(), atol=1e-6)

    def testFrameMethodFollowsRigidRotation(self):
        """Rotating the whole Calpha trace rotates every atom with it under the frame method (but not under
        translation)."""
        R = rotation(0.7, [0.3, -1.0, 0.5])
        c = A.mean(0)
        P = np.array([A, (A - c) @ R.T + c])
        ref = OPEN_ALL.select('protein').getCoords()
        expected = (ref - c) @ R.T + c
        aa = buildAllAtomPath(P, OPEN_ALL, method='frame')
        assert_allclose(aa.getCoordsets(1), expected, atol=1e-6)
        tr = buildAllAtomPath(P, OPEN_ALL, method='translate')
        self.assertGreater(np.abs(tr.getCoordsets(1) - expected).max(), 0.5)

    def testPulchraNeedsTheProgram(self):
        with self.assertRaises(EnvironmentError):
            buildAllAtomPath(np.array([A]), OPEN_ALL, method='pulchra', pulchra='no_such_program_xyz')


class TestPlots(unittest.TestCase):

    def testEveryShowFunctionDraws(self):
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            raise unittest.SkipTest('matplotlib is not installed')
        L = line(A, B, 6)
        showPathGeometry(calcPathGeometry(L, target=B))
        showStepWork(calcStepWork(L, target=B))
        showDeformationModes([calcDeformationModes(OPEN, CLOSED, model='anm'),
                              calcDeformationModes(OPEN, CLOSED, model='genanm')], labels=['ANM', 'genANM'])
        showRouteDistance(calcRouteDistance(L, L[::-1]))
        showClosestApproach(calcClosestApproach(L, L[3], target=B))
        showPathPCA(calcPathPCA([L, L[::-1]], structures=[A, B]))
        plt.close('all')


if __name__ == '__main__':
    unittest.main()
