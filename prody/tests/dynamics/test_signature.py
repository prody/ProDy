"""This module contains unit tests for :mod:`~prody.KDTree` module."""

import numpy as np
from numpy.testing import assert_array_equal, assert_equal, assert_allclose
from numpy.random import rand, randint

from prody.dynamics import sdarray
from prody.dynamics.signature import ModeEnsemble

from prody.tests import unittest
from prody.tests.datafiles import parseDatafile

from prody import (_PY3K, LOGGER, ANM, PDBEnsemble, calcEnsembleENMs,
                   calcSignatureModes, calcSignatureOverlaps, calcOverlap)

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
        assert s.shape==A[0].shape, 'failed at sdarray slicing'
        #assert_array_equal(s.flatten(), A[0].flatten(), 'failed at sdarray slicing')
        #assert_array_equal(s.getWeights().flatten(), W[0].flatten(), 'failed at sdarray slicing')

        s = S[:, 0]
        assert s.shape==A[:, 0].shape, 'failed at sdarray slicing'
        #assert_array_equal(s.flatten(), A[:, 0].flatten(), 'failed at sdarray slicing')
        #assert_array_equal(s.getWeights().flatten(), W[:, 0].flatten(), 'failed at sdarray slicing')

        s = S[:, :, :]
        assert s.shape==A.shape, 'failed at sdarray slicing'
        #assert_array_equal(s.flatten(), A.flatten(), 'failed at sdarray slicing')
        #assert_array_equal(s.getWeights().flatten(), W.flatten(), 'failed at sdarray slicing')

        s = S[0, 0, 0]
        #assert_array_equal(s, A[0, 0, 0], 'failed at sdarray slicing')


class TestSignatureOverlapsRefModel(unittest.TestCase):
    """Tests for the *ref_model* option of :func:`.calcSignatureOverlaps`, which
    overlaps each modeset of the ensemble against a single reference model
    rather than within the ensemble."""

    @classmethod
    def setUpClass(cls):
        ca = parseDatafile('pdb1ake').select('calpha')
        cls.ca = ca
        cls.n_atoms = ca.numAtoms()
        cls.n_modes = 5
        cls.n_sets = 3

        # a small matched mode ensemble on a common set of atoms
        me = ModeEnsemble('test')
        for _ in range(cls.n_sets):
            anm = ANM()
            anm.buildHessian(ca)
            anm.calcModes(n_modes=cls.n_modes, zeros=False)
            me.addModeSet(anm[:cls.n_modes])
        cls.me = me

        # a multi-mode reference and a single-mode reference on the same atoms
        ref = ANM()
        ref.buildHessian(ca)
        ref.calcModes(n_modes=cls.n_modes, zeros=False)
        cls.ref_full = ref[:cls.n_modes]
        cls.ref_single = ref[0]

    # -- shapes -----------------------------------------------------------

    def testRefModelFullShape(self):
        """A multi-mode reference gives (n_modes, n_modes_ref, n_sets, 1)."""
        ov = calcSignatureOverlaps(self.me, ref_model=self.ref_full, diag=False)
        self.assertEqual(ov.shape,
                         (self.n_modes, self.n_modes, self.n_sets, 1))

    def testRefModelCollapseShape(self):
        """collapse stacks each set's block: (n_modes*n_sets, n_modes_ref)."""
        ov = calcSignatureOverlaps(self.me, ref_model=self.ref_full,
                                   diag=False, collapse=True)
        self.assertEqual(ov.shape,
                         (self.n_modes * self.n_sets, self.n_modes))
        self.assertTrue(np.all(np.isfinite(ov)))

    def testRefModelSingleModeNonCollapse(self):
        """A single-mode reference gives n_modes_ref == 1 (reshape path)."""
        ov = calcSignatureOverlaps(self.me, ref_model=self.ref_single,
                                   diag=False)
        self.assertEqual(ov.shape, (self.n_modes, 1, self.n_sets, 1))

    def testRefModelSingleModeCollapse(self):
        """Single-mode reference with collapse: (n_modes*n_sets, 1)."""
        ov = calcSignatureOverlaps(self.me, ref_model=self.ref_single,
                                   diag=False, collapse=True)
        self.assertEqual(ov.shape, (self.n_modes * self.n_sets, 1))
        self.assertTrue(np.all(np.isfinite(ov)))

    def testRefModelDiagShape(self):
        """diag against an equal-mode reference: (n_modes, n_sets, 1)."""
        ov = calcSignatureOverlaps(self.me, ref_model=self.ref_full, diag=True)
        self.assertEqual(ov.shape, (self.n_modes, self.n_sets, 1))

    # -- values -----------------------------------------------------------

    def testNonCollapseMatchesCalcOverlap(self):
        """Each block equals abs(calcOverlap(modeset_i, ref))."""
        ov = calcSignatureOverlaps(self.me, ref_model=self.ref_full, diag=False)
        for i in range(self.n_sets):
            expected = np.abs(calcOverlap(self.me[i], self.ref_full))
            assert_allclose(ov[:, :, i, 0], expected, atol=1e-8)

    def testCollapseMatchesNonCollapse(self):
        """The collapsed stack holds the same blocks as the 4-D array."""
        ov4 = calcSignatureOverlaps(self.me, ref_model=self.ref_full, diag=False)
        ov2 = calcSignatureOverlaps(self.me, ref_model=self.ref_full,
                                    diag=False, collapse=True)
        for i in range(self.n_sets):
            block = ov2[i * self.n_modes:(i + 1) * self.n_modes, :]
            assert_allclose(block, ov4[:, :, i, 0], atol=1e-8)

    def testSingleModeCollapseMatchesNonCollapse(self):
        """Single-mode reference: collapse and non-collapse agree."""
        ov4 = calcSignatureOverlaps(self.me, ref_model=self.ref_single,
                                    diag=False)
        ov2 = calcSignatureOverlaps(self.me, ref_model=self.ref_single,
                                    diag=False, collapse=True)
        for i in range(self.n_sets):
            block = ov2[i * self.n_modes:(i + 1) * self.n_modes, :]
            assert_allclose(block, ov4[:, :, i, 0], atol=1e-8)

    def testListInputEqualsEnsembleInput(self):
        """A plain list of modesets is accepted and matches ModeEnsemble input."""
        modesets = [ms for ms in self.me]
        ov_list = calcSignatureOverlaps(modesets, ref_model=self.ref_full,
                                        diag=False)
        ov_ens = calcSignatureOverlaps(self.me, ref_model=self.ref_full,
                                       diag=False)
        assert_allclose(ov_list, ov_ens, atol=1e-8)

    # -- errors -----------------------------------------------------------

    def testRefModelWrongTypeRaises(self):
        self.assertRaises(TypeError, calcSignatureOverlaps, self.me,
                          ref_model='not-a-model')

    def testRefModelWrongNumAtomsRaises(self):
        small = ANM()
        small.buildHessian(self.ca[:self.n_atoms // 2])
        small.calcModes(n_modes=self.n_modes, zeros=False)
        self.assertRaises(ValueError, calcSignatureOverlaps, self.me,
                          ref_model=small[:self.n_modes], diag=False)

    def testRefModelDiagModeMismatchRaises(self):
        """diag needs the reference to have as many modes as the ensemble."""
        self.assertRaises(ValueError, calcSignatureOverlaps, self.me,
                          ref_model=self.ref_single, diag=True)

    # -- regression: within-ensemble path still works ---------------------

    def testWithinEnsembleCollapseUnchanged(self):
        """Without ref_model, collapse gives the full block matrix."""
        ov = calcSignatureOverlaps(self.me, diag=False, collapse=True)
        self.assertEqual(
            ov.shape,
            (self.n_modes * self.n_sets, self.n_modes * self.n_sets))


class TestSignatureOverlapsSignDy(unittest.TestCase):
    """ref_model overlaps on GNM and ANM mode ensembles built the SignDy way,
    i.e. :func:`.calcEnsembleENMs` on a :class:`.PDBEnsemble`, overlapped against
    the mean signature model from :func:`.calcSignatureModes` (as in the SignDy
    tutorial). Exercised on the 2k39 ubiquitin NMR ensemble."""

    @classmethod
    def setUpClass(cls):
        cls.n_modes = 10
        cls.n_sets = 6

        ags = parseDatafile('2k39_ca')
        ens = PDBEnsemble('2k39')
        ens.setAtoms(ags)
        ens.setCoords(ags.getCoords())
        ens.addCoordset(ags.getCoordsets()[:cls.n_sets])
        ens.iterpose()
        cls.n_atoms = ens.numAtoms()

        # one GNM and one ANM mode ensemble, each with its mean-signature model
        cls.mes = {}
        for model in ('gnm', 'anm'):
            me = calcEnsembleENMs(ens, model=model, trim='reduce',
                                  n_modes=cls.n_modes, match=True)
            cls.mes[model] = (me, calcSignatureModes(me))

    def testEnsembleShape(self):
        """Both ensembles are matched with the expected dimensions."""
        for model, (me, ref) in self.mes.items():
            with self.subTest(model=model):
                self.assertEqual(me.numModeSets(), self.n_sets)
                self.assertEqual(me.numModes(), self.n_modes)
                self.assertEqual(me.numAtoms(), self.n_atoms)
                self.assertTrue(me.isMatched())
                self.assertEqual(ref.numModes(), self.n_modes)
                self.assertEqual(ref.numAtoms(), self.n_atoms)

    def testRefModelFullShapeAndValues(self):
        """(n_modes, n_modes_ref, n_sets, 1) and each block == calcOverlap."""
        for model, (me, ref) in self.mes.items():
            with self.subTest(model=model):
                ov = calcSignatureOverlaps(me, ref_model=ref, diag=False)
                self.assertEqual(
                    ov.shape, (self.n_modes, self.n_modes, self.n_sets, 1))
                for i in range(self.n_sets):
                    expected = np.abs(calcOverlap(me[i], ref))
                    assert_allclose(ov[:, :, i, 0], expected, atol=1e-8)

    def testRefModelCollapse(self):
        """Collapsed stack matches the 4-D blocks for both models."""
        for model, (me, ref) in self.mes.items():
            with self.subTest(model=model):
                ov4 = calcSignatureOverlaps(me, ref_model=ref, diag=False)
                ov2 = calcSignatureOverlaps(me, ref_model=ref, diag=False,
                                            collapse=True)
                self.assertEqual(
                    ov2.shape, (self.n_modes * self.n_sets, self.n_modes))
                for i in range(self.n_sets):
                    block = ov2[i * self.n_modes:(i + 1) * self.n_modes, :]
                    assert_allclose(block, ov4[:, :, i, 0], atol=1e-8)

    def testRefModelSingleMode(self):
        """A single signature mode reference exercises the reshape path."""
        for model, (me, ref) in self.mes.items():
            with self.subTest(model=model):
                ov = calcSignatureOverlaps(me, ref_model=ref[0], diag=False)
                self.assertEqual(ov.shape, (self.n_modes, 1, self.n_sets, 1))
                ovc = calcSignatureOverlaps(me, ref_model=ref[0], diag=False,
                                            collapse=True)
                self.assertEqual(ovc.shape, (self.n_modes * self.n_sets, 1))
                for i in range(self.n_sets):
                    block = ovc[i * self.n_modes:(i + 1) * self.n_modes, :]
                    assert_allclose(block, ov[:, :, i, 0], atol=1e-8)

    def testRefModelDiag(self):
        """diag against the equal-mode signature model: (n_modes, n_sets, 1)."""
        for model, (me, ref) in self.mes.items():
            with self.subTest(model=model):
                ov = calcSignatureOverlaps(me, ref_model=ref, diag=True)
                self.assertEqual(ov.shape, (self.n_modes, self.n_sets, 1))

    def testSlicedEnsembleTutorialStyle(self):
        """Mode-sliced ensembles (gnms[:, :5]) overlap against a sliced ref."""
        k = 5
        for model, (me, ref) in self.mes.items():
            with self.subTest(model=model):
                ov = calcSignatureOverlaps(me[:, :k], ref_model=ref[:k],
                                           diag=False)
                self.assertEqual(ov.shape, (k, k, self.n_sets, 1))
