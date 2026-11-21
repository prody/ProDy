import numpy as np
import pytest

from prody.dynamics.generalized_anm import genANM, calcGenANM

def _sample_coords():
    return np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        [0.0, 0.0, 1.5],
    ], dtype=float)

def test_build_hessian_symmetry_and_shape():
    coords = _sample_coords()
    model = genANM('test')
    H = model.buildHessian(coords, cutoff=5.0, gamma=1.0)
    assert H.shape == (12, 12)
    assert np.allclose(H, H.T, atol=1e-8)

def test_local_assembly_matches_full_outer_products():
    """The angle and torsion terms are assembled on the local blocks of the atoms they involve; this must equal the
    definition H = sum_q k_q J_q J_q^T with the full-length gradient J_q (zero outside those atoms), exactly."""
    from prody.dynamics import generalized_anm as g
    rng = np.random.default_rng(1)
    steps = rng.normal(size=(30, 3)); steps *= 3.8 / np.linalg.norm(steps, axis=1)[:, None]
    X = np.cumsum(steps, axis=0)                                   # 30-residue chain with 3.8 A virtual bonds
    bonded = np.ones(29, dtype=bool); bonded[14] = False           # one chain break
    N = len(X); ref_a = np.zeros((3 * N, 3 * N)); ref_t = np.zeros((3 * N, 3 * N))
    def angle(P):
        b, c = P[0] - P[1], P[2] - P[1]
        return np.arccos(np.dot(b, c) / (np.linalg.norm(b) * np.linalg.norm(c)))
    for j in range(1, N - 1):                                      # angle gradients by central differences, independently
        if not (bonded[j - 1] and bonded[j]): continue
        P = X[j - 1:j + 2].copy(); gr = np.zeros(9); e = 1e-6
        for k in range(9):
            Pp = P.copy(); Pp.flat[k] += e; Pm = P.copy(); Pm.flat[k] -= e
            gr[k] = (angle(Pp) - angle(Pm)) / (2 * e)
        gr = g._project_out_rigid_components(gr, P)
        J = np.zeros(3 * N); J[3 * (j - 1):3 * (j + 2)] = gr; ref_a += np.outer(J, J)
    for i in range(N - 3):
        if not (bonded[i] and bonded[i + 1] and bonded[i + 2]): continue
        g4 = g._torsion_gradient_fd(X[i:i + 4])
        if g4 is None: continue
        J = np.zeros(3 * N); J[3 * i:3 * (i + 4)] = g4.reshape(-1); ref_t += np.outer(J, J)
    np.testing.assert_allclose(g._angle_blocks(X, 1.0, bonded), ref_a, rtol=0, atol=1e-7)
    np.testing.assert_allclose(g._torsion_blocks(X, 1.0, bonded), ref_t, rtol=0, atol=1e-12)

def test_calc_modes_and_calcGenANM_with_array():
    coords = _sample_coords()
    model = genANM('test2')
    H = model.buildHessian(coords, cutoff=5.0)
    model.calcModes(n_modes=6, zeros=False)
    eigvals = model.getEigvals()
    assert eigvals is not None and eigvals.size > 0

    # calcGenANM accepting Hessian array should return an ANM-like model
    m2 = calcGenANM(H, n_modes=6, zeros=False)
    assert hasattr(m2, 'getHessian')
    # calcGenANM sets a symmetric copy via setHessian in our implementation
    assert np.allclose(m2.getHessian(), 0.5*(H + H.T), atol=1e-8)


# ---------------------------------------------------------------------------
# Properties of genANM with the default constants (zero modes, reduction to ANM, calcENM, save/load)
# ---------------------------------------------------------------------------

import os
import tempfile
from numpy.testing import assert_allclose, assert_array_equal
from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

ATOMS = parseDatafile('1ubi_ca')


class TestGenANM(unittest.TestCase):

    def setUp(self):
        self.gen = genANM('1ubi genANM')
        self.gen.buildHessian(ATOMS)
        self.gen.calcModes(n_modes=10, zeros=True)

    def testDefaultConstants(self):
        """The force constants 20, 10, 1 and 0 are the defaults."""
        self.assertEqual(self.gen.getForceConstants(), (20.0, 10.0, 1.0, 0.0))
        self.assertEqual(self.gen.numZeros(), 6)

    def testSixZeroModes(self):
        """Internal-coordinate terms keep the six rigid-body zero modes."""
        vals = self.gen.getEigvals()
        assert_allclose(vals[:6], 0, atol=1e-6)
        self.assertGreater(vals[6], 1e-3)

    def testReducesToANM(self):
        """With zero covalent weights the Hessian is the classical ANM Hessian."""
        anm = ANM()
        anm.buildHessian(ATOMS, cutoff=15.0, gamma=1.0)
        gen0 = genANM()
        gen0.buildHessian(ATOMS, k_bond=0, k_theta=0, k_phi=0)
        assert_allclose(gen0.getHessian(), anm.getHessian(), atol=1e-10)

    def testStifferThanANM(self):
        """Adding springs never softens a mode."""
        anm = ANM()
        anm.buildHessian(ATOMS)
        anm.calcModes(n_modes=10)
        gen = genANM()
        gen.buildHessian(ATOMS)
        gen.calcModes(n_modes=10)
        self.assertTrue(np.all(gen.getEigvals() >= anm.getEigvals() - 1e-8))

    def testCurvatureLeavesThreeZeroModes(self):
        """A nonzero curvature weight keeps only the three translational zero modes."""
        gen = genANM()
        gen.buildHessian(ATOMS, kappa=20.0)
        self.assertEqual(gen.numZeros(), 3)
        gen.calcModes(n_modes=5, zeros=True)
        vals = gen.getEigvals()
        assert_allclose(vals[:3], 0, atol=1e-6)
        self.assertGreater(vals[3], 1e-3)

    def testHessianSymmetricAndMatchesBuilder(self):
        H = buildGeneralizedHessian(ATOMS.getCoords())
        assert_allclose(H, H.T, atol=1e-12)
        assert_allclose(H, self.gen.getHessian(), atol=1e-10)

    def testUnknownKeywordRejected(self):
        """Direct calls are strict, so typos are caught."""
        gen = genANM()
        self.assertRaises(TypeError, gen.buildHessian, ATOMS, k_bnd=20)

    def testCalcGenANM(self):
        model, sel = calcGenANM(ATOMS, selstr='all', n_modes=5)
        self.assertIsInstance(model, genANM)
        # setUp keeps 10 modes including the six zero modes, so four internal modes to compare
        assert_allclose(model.getEigvals()[:4], self.gen.getEigvals()[6:10], rtol=1e-8)


class TestCalcENMGenANM(unittest.TestCase):

    def testModelString(self):
        """calcENM(model='genanm') gives the same modes as genANM directly."""
        enm, _ = calcENM(ATOMS, model='genanm', n_modes=5)
        self.assertIsInstance(enm, genANM)
        gen = genANM()
        gen.buildHessian(ATOMS)
        gen.calcModes(n_modes=5)
        assert_allclose(enm.getEigvals(), gen.getEigvals(), rtol=1e-8)

    def testModelClass(self):
        enm, _ = calcENM(ATOMS, model=genANM, n_modes=5)
        self.assertIsInstance(enm, genANM)

    def testKeywordsForwarded(self):
        """genANM keywords reach the model; unrelated keywords are ignored, as for ANM."""
        enm, _ = calcENM(ATOMS, model='genanm', n_modes=5, k_bond=5.0, unrelated_option=1)
        self.assertEqual(enm.getForceConstants()[0], 5.0)


class TestSaveLoadGenANM(unittest.TestCase):

    def testRoundTrip(self):
        """A saved genANM loads back as a genANM with its constants and modes."""
        gen = genANM('ubi')
        gen.buildHessian(ATOMS, k_bond=15.0, kappa=20.0)
        gen.calcModes(n_modes=5)
        tmp = tempfile.mkdtemp()
        fn = saveModel(gen, os.path.join(tmp, 'ubi'))
        self.assertTrue(fn.endswith('.genanm.npz'))
        loaded = loadModel(fn)
        self.assertIsInstance(loaded, genANM)
        self.assertEqual(loaded.getForceConstants(), gen.getForceConstants())
        self.assertEqual(loaded.numZeros(), 3)
        assert_array_equal(loaded.getEigvals(), gen.getEigvals())


if __name__ == '__main__':
    unittest.main()
