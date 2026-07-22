"""This module contains unit tests for :func:`.calcPerturbResponse` and
:func:`.comparePerturbResponses`, in particular the directional
``return_vectors`` response ensemble and the *force*/*perturb_node* options,
exercised on the adenylate kinase open->closed transition
(4ake chain A -> 1ake chain A)."""

import numpy as np
from numpy.testing import assert_allclose

from prody import (ANM, GNM, matchChains, superpose, calcDeformVector,
                   calcOverlap, calcPerturbResponse, comparePerturbResponses,
                   LOGGER)
from prody.dynamics.signature import ModeEnsemble
from prody.tests import unittest
from prody.tests.datafiles import parseDatafile

LOGGER.verbosity = 'none'


class TestCalcPerturbResponse(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # open (4ake, chain A) and closed (1ake) adenylate kinase; bundled data
        op = parseDatafile('pdb4ake_fixed').select('calpha')
        cl = parseDatafile('pdb1ake').select('calpha')
        matches = matchChains(op, cl, pwalign=True, seqid=90, overlap=90)
        op_ca, cl_ca = matches[0][0], matches[0][1]

        # closure target: open->closed internal deformation (after superposition)
        cl_fit = cl_ca.copy()
        superpose(cl_fit, op_ca)
        cls.target = calcDeformVector(op_ca, cl_fit)

        cls.n_atoms = op_ca.numAtoms()
        cls.anm = ANM('4ake_A')
        cls.anm.buildHessian(op_ca)
        cls.anm.calcModes(n_modes=None, zeros=False)
        cls.gnm = GNM('4ake_A')
        cls.gnm.buildKirchhoff(op_ca)
        cls.gnm.calcModes(n_modes=None, zeros=False)

    def testLegacyReturnsTriple(self):
        """Default call is unchanged: (matrix, effectiveness, sensitivity)."""
        prs, eff, sens = calcPerturbResponse(self.anm)
        self.assertEqual(prs.shape, (self.n_atoms, self.n_atoms))
        self.assertEqual(eff.shape, (self.n_atoms,))
        self.assertEqual(sens.shape, (self.n_atoms,))

    def testReturnVectorsEnsembleShape(self):
        """return_vectors gives a ModeEnsemble: one modeset/node, repeats modes."""
        repeats = 7
        ens = calcPerturbResponse(self.anm, return_vectors=True, repeats=repeats)
        self.assertIsInstance(ens, ModeEnsemble)
        self.assertEqual(ens.numModeSets(), self.n_atoms)
        self.assertEqual(ens.numModes(), repeats)
        self.assertEqual(ens.numAtoms(), self.n_atoms)

    def testResponseVectorEqualsCovDotForce(self):
        """The kept response vector is exactly cov . f (no square/average)."""
        cov = self.anm.getCovariance()
        node = 5
        force = np.array([1.0, -2.0, 0.5])
        ens = calcPerturbResponse(self.anm, return_vectors=True,
                                  perturb_node=node, force=force)
        self.assertEqual(ens.numModeSets(), 1)
        self.assertEqual(ens.numModes(), 1)
        resp = np.asarray(ens[0].getEigvecs()).flatten()
        expected = np.dot(cov[:, node * 3:node * 3 + 3], force)
        assert_allclose(resp, expected, rtol=1e-6, atol=1e-8)

    def testEigvalsAreSquaredResponseMagnitude(self):
        """Eigenvalues store the squared response magnitude of each force."""
        cov = self.anm.getCovariance()
        node = 3
        force = np.array([0.3, 0.4, 0.0])
        ens = calcPerturbResponse(self.anm, return_vectors=True,
                                  perturb_node=node, force=force)
        expected = np.dot(cov[:, node * 3:node * 3 + 3], force)
        eigval = np.atleast_1d(ens[0].getEigvals())[0]
        assert_allclose(eigval, (expected ** 2).sum(), rtol=1e-6, atol=1e-8)

    def testPerturbNodeList(self):
        """perturb_node restricts the scan to the requested nodes."""
        nodes = [0, 4, 9]
        ens = calcPerturbResponse(self.anm, return_vectors=True,
                                  repeats=3, perturb_node=nodes)
        self.assertEqual(ens.numModeSets(), len(nodes))

    def testMultipleForces(self):
        """An (m, 3) force array gives m responses per node, ignoring repeats."""
        forces = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        ens = calcPerturbResponse(self.anm, return_vectors=True, repeats=50,
                                  perturb_node=2, force=forces)
        self.assertEqual(ens.numModeSets(), 1)
        self.assertEqual(ens.numModes(), 2)

    def testBadForceShapeRaises(self):
        """A force that is not 3-component is rejected."""
        self.assertRaises(ValueError, calcPerturbResponse, self.anm,
                          return_vectors=True, force=np.array([1.0, 0.0]))

    def testBadPerturbNodeRaises(self):
        """Out-of-range node indices are rejected."""
        self.assertRaises(ValueError, calcPerturbResponse, self.anm,
                          return_vectors=True, perturb_node=self.n_atoms + 1)

    def testPartialMatrix(self):
        """perturb_node without return_vectors gives a partial PRS matrix: only
        the perturbed rows are filled, the rest stay zero."""
        nodes = [0, 1, 2]
        prs, eff, sens = calcPerturbResponse(self.anm, perturb_node=nodes,
                                             return_vectors=False, repeats=10)
        nonzero = np.where(prs.any(axis=1))[0]
        assert_allclose(nonzero, nodes)

    def testReturnMatrixReturnsBoth(self):
        """return_vectors + return_matrix returns the triple and the ensemble."""
        out = calcPerturbResponse(self.anm, return_vectors=True,
                                  return_matrix=True, repeats=5)
        self.assertEqual(len(out), 4)
        prs, eff, sens, ens = out
        self.assertEqual(prs.shape, (self.n_atoms, self.n_atoms))
        self.assertEqual(eff.shape, (self.n_atoms,))
        self.assertEqual(sens.shape, (self.n_atoms,))
        self.assertIsInstance(ens, ModeEnsemble)

    def testGNMRaisesForVectors(self):
        """Directional responses require a 3d model; GNM has no direction."""
        self.assertRaises(ValueError, calcPerturbResponse, self.gnm,
                          return_vectors=True)

    def testDirectionalPRSRecoversClosure(self):
        """A directed force at some node drives the open->closed closure: the
        best response overlaps the transition better than any single ANM mode."""
        np.random.seed(42)
        ens = calcPerturbResponse(self.anm, return_vectors=True, repeats=100)
        tv = self.target.getArray()
        tv = tv / np.linalg.norm(tv)

        best_per_node = []
        for i in range(ens.numModeSets()):
            R = np.asarray(ens[i].getEigvecs())
            if R.shape[0] != tv.shape[0]:
                R = R.T
            R = R / np.linalg.norm(R, axis=0, keepdims=True)
            best_per_node.append(np.abs(tv @ R).max())
        best_directional = float(np.max(best_per_node))

        best_single_mode = float(np.abs(calcOverlap(self.anm[:20], self.target)).max())

        # closure is genuinely recovered, and beats the best single normal mode
        self.assertGreater(best_directional, 0.6)
        self.assertGreater(best_directional, best_single_mode)


class TestComparePerturbResponses(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        op = parseDatafile('pdb4ake_fixed').select('calpha')
        cl = parseDatafile('pdb1ake').select('calpha')
        matches = matchChains(op, cl, pwalign=True, seqid=90, overlap=90)
        op_ca, cl_ca = matches[0][0], matches[0][1]
        cl_fit = cl_ca.copy()
        superpose(cl_fit, op_ca)
        cls.target = calcDeformVector(op_ca, cl_fit)

        cls.n_atoms = op_ca.numAtoms()
        cls.anm = ANM('4ake_A')
        cls.anm.buildHessian(op_ca)
        cls.anm.calcModes(n_modes=None, zeros=False)

        np.random.seed(7)
        cls.ens = calcPerturbResponse(cls.anm, return_vectors=True, repeats=20)

    def testProfileShape(self):
        """One overlap value per perturbed node."""
        prof = comparePerturbResponses(self.ens, self.target)
        self.assertEqual(prof.shape, (self.n_atoms,))
        self.assertTrue(np.all((prof >= 0) & (prof <= 1 + 1e-9)))

    def testVectorAndModelTargets(self):
        """A Vector target and a multi-mode model target both work; overlapping
        against the whole model is at least as good as against one mode."""
        prof_vec = comparePerturbResponses(self.ens, self.target, stat='max')
        prof_model = comparePerturbResponses(self.ens, self.anm, stat='max')
        prof_mode = comparePerturbResponses(self.ens, self.anm[0], stat='max')
        self.assertEqual(prof_model.shape, (self.n_atoms,))
        # best over all modes >= best over a single mode, per node
        self.assertTrue(np.all(prof_model >= prof_mode - 1e-9))
        self.assertEqual(prof_vec.shape, (self.n_atoms,))

    def testStatOptions(self):
        """max >= mean for the per-node reduction over force directions."""
        pmax = comparePerturbResponses(self.ens, self.target, stat='max')
        pmean = comparePerturbResponses(self.ens, self.target, stat='mean')
        self.assertTrue(np.all(pmax >= pmean - 1e-9))

    def testEnsembleVsSelfIsOne(self):
        """Comparing an ensemble with itself gives 1 (best force pair matches)."""
        prof = comparePerturbResponses(self.ens, self.ens, stat='max')
        assert_allclose(prof, np.ones(self.n_atoms), rtol=0, atol=1e-6)

    def testEnsembleMismatchRaises(self):
        """Two ensembles must have the same number of modesets."""
        small = calcPerturbResponse(self.anm, return_vectors=True, repeats=3,
                                    perturb_node=[0, 1])
        self.assertRaises(ValueError, comparePerturbResponses, self.ens, small)

    def testBadArgumentsRaise(self):
        self.assertRaises(ValueError, comparePerturbResponses, self.ens,
                          self.target, stat='median')
        self.assertRaises(TypeError, comparePerturbResponses, self.ens,
                          'not-a-target')
        self.assertRaises(TypeError, comparePerturbResponses, self.target,
                          self.target)


if __name__ == '__main__':
    unittest.main()
