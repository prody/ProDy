"""This module contains unit tests for :func:`.calcPerturbResponse`, in
particular the directional ``return_vectors`` response ensemble, exercised on
the adenylate kinase open->closed transition (4ake chain A -> 1ake chain A)."""

import numpy as np
from numpy.testing import assert_allclose

from prody import (ANM, GNM, matchChains, superpose,
                   calcDeformVector, calcOverlap, calcPerturbResponse, LOGGER)
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


if __name__ == '__main__':
    unittest.main()
