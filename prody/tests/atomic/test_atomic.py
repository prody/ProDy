"""This module contains unit tests for :mod:`~prody.atomic`."""

import os.path
import pickle

from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.atomic.fields import READONLY
from prody.tests import unittest, TEMPDIR
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

ATOL = 1e-5
RTOL = 0

ATOMS = parseDatafile('multi_model_truncated', subset='ca')


class TestCopying(unittest.TestCase):

    def TestCopyAtomGroup(self):

        atoms = ATOMS.copy()

        assert_equal(atoms.getCoordsets(), ATOMS.getCoordsets())
        for label in ATOMS.getDataLabels():
            if label in READONLY:
                continue
            assert_equal(atoms.getData(label), ATOMS.getData(label))

    def TestCopyChain(self):

        CHAIN = ATOMS['A']
        chain = CHAIN.copy()

        assert_equal(chain.getCoordsets(), CHAIN.getCoordsets())
        for label in ATOMS.getDataLabels():
            if label in READONLY:
                continue
            assert_equal(chain.getData(label), CHAIN.getData(label),
                         'failed to copy ' + label)

    def TestCopyAtom(self):

        ATOM = ATOMS[0]
        atom = ATOM.copy()

        assert_equal(atom[0].getCoordsets(), ATOM.getCoordsets())
        for label in ATOMS.getDataLabels():
            if label in READONLY:
                continue
            assert_equal(atom[0].getData(label), ATOM.getData(label),
                         'failed to copy ' + label)


    def TestCopySelstr(self):

        SELECTION = ATOMS.calpha
        selection = SELECTION.copy()

        assert_equal(selection.getCoordsets(), SELECTION.getCoordsets())
        for label in ATOMS.getDataLabels():
            if label in READONLY:
                continue
            assert_equal(selection.getData(label), SELECTION.getData(label),
                         'failed to copy ' + label)

class TestSaveLoad(unittest.TestCase):

    def testSaveLoad(self):

        atoms = loadAtoms(saveAtoms(ATOMS, os.path.join(TEMPDIR, 'atoms')))
        assert_equal(atoms.getCoordsets(), ATOMS.getCoordsets())
        for label in ATOMS.getDataLabels():
            assert_equal(atoms.getData(label), ATOMS.getData(label),
                         'failed to load ' + label)


class TestPickling(unittest.TestCase):


    def testAtomGroup(self):

        atoms1 = parseDatafile('multi_model_truncated', subset='ca')
        atoms2 = pickle.loads(pickle.dumps(atoms1))
        s1 = atoms1.__getstate__()
        s2 = atoms2.__getstate__()
        for key in s1:
            assert_equal(s1[key], s2[key])
        self.assertEqual(atoms1, atoms2)

    def testSelection(self):

        sel = parseDatafile('multi_model_truncated', subset='ca')[:3]
        self.assertEqual(sel, pickle.loads(pickle.dumps(sel)))

    def testAtom(self):

        atom = parseDatafile('multi_model_truncated', subset='ca')[0]
        self.assertEqual(atom, pickle.loads(pickle.dumps(atom)))
