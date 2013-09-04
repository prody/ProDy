#!/usr/bin/python
# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
#
# Copyright (C) 2010-2012 Ahmet Bakan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module contains unit tests for :mod:`~prody.proteins`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os

import numpy as np
from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR, unittest
from prody.tests.test_datafiles import *

LOGGER.verbosity = 'none'

class TestParsePDB(unittest.TestCase):

    def setUp(self):
        """Set PDB file data and parse the PDB file."""

        self.pdb = DATA_FILES['multi_model_truncated']
        self.one = DATA_FILES['oneatom']
        self.ca = DATA_FILES['1ubi_ca']

    def testUsualCase(self):
        """Test the outcome of a simple parsing scenario."""

        ag = parseDatafile(self.pdb['file'])

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePDB failed to return an AtomGroup instance')

        self.assertEqual(ag.numAtoms(), self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms')

        self.assertEqual(ag.numCoordsets(), self.pdb['models'],
            'parsePDB failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(ag.getTitle(),
             os.path.splitext(self.pdb['file'])[0],
            'failed to set AtomGroup title based on filename')

    def testPDBArgument(self):
        """Test outcome of invalid *pdb* arguments."""

        self.assertRaises(IOError, parsePDB, self.pdb['file'] + '.gz')
        self.assertRaises(TypeError, parsePDB, None)

    def testModelArgument(self):
        """Test outcome of valid and invalid *model* arguments."""

        path = pathDatafile(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, model='0')
        self.assertRaises(ValueError, parsePDB, path, model=-1)
        self.assertRaises(proteins.PDBParseError, parsePDB, path,
                          model=self.pdb['models']+1)
        self.assertIsNone(parsePDB(path, model=0),
            'parsePDB failed to parse no coordinate sets')
        self.assertEqual(parsePDB(path, model=1).numCoordsets(), 1,
            'parsePDB failed to parse the first coordinate set')
        self.assertEqual(parsePDB(path, model=self.pdb['models'])
            .numCoordsets(), 1,
            'parsePDB failed to parse the last coordinate set')

    def testTitleArgument(self):
        """Test outcome of *title* argument."""

        path = pathDatafile(self.pdb['file'])
        title = 'small protein'
        self.assertEqual(parsePDB(path, title=title).getTitle(),
             title, 'parsePDB failed to set user given title')

        name = 1999
        self.assertEqual(parsePDB(path, title=title).getTitle(),
             str(title), 'parsePDB failed to set user given non-string name')

    def testChainArgument(self):
        """Test outcome of valid and invalid *chain* arguments."""

        path = pathDatafile(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, chain=['A'])
        self.assertRaises(ValueError, parsePDB, path, chain='')
        self.assertIsNone(parsePDB(path, chain='$'))
        self.assertEqual(parsePDB(path, chain='A')
            .numAtoms(), self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms when chain is '
            'specified')

    def testSubsetArgument(self):
        """Test outcome of valid and invalid *subset* arguments."""

        path = pathDatafile(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, subset=['A'])
        self.assertEqual(parsePDB(path, subset='ca').numAtoms(), 10,
                        'failed to parse correct number of "ca" atoms')
        self.assertEqual(parsePDB(path, subset='bb').numAtoms(), 40,
                        'failed to parse correct number of "bb" atoms')

    def testAgArgument(self):
        """Test outcome of valid and invalid *ag* arguments."""

        path = pathDatafile(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, ag='AtomGroup')
        ag = prody.AtomGroup('One atom')
        ag.setCoords(np.array([[0., 0., 0.]]))
        self.assertRaises(ValueError, parsePDB, path, ag=ag)
        ag = prody.AtomGroup('Test')
        self.assertEqual(parsePDB(path, ag=ag).numAtoms(),
            self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms')

    def testAgArgMultiModel(self):
        """Test number of coordinate sets when using *ag* arguments."""

        path = pathDatafile(self.pdb['file'])
        ag = parsePDB(path)
        coords = ag.getCoordsets()
        ncsets = ag.numCoordsets()
        ag = parsePDB(path, ag=ag)
        self.assertEqual(ag.numCoordsets(), ncsets*2,
            'parsePDB failed to append coordinate sets to given ag')
        assert_equal(coords, ag.getCoordsets(np.arange(ncsets, ncsets*2)))

    def testAgArgSingleModel(self):
        """Test number of coordinate sets when using *ag* arguments."""

        path = pathDatafile(self.ca['file'])
        ag = parsePDB(path)
        coords = ag.getCoordsets()
        ncsets = ag.numCoordsets()
        ag = parsePDB(path, ag=ag, subset='bb')
        self.assertEqual(ag.numCoordsets(), ncsets*2,
            'parsePDB failed to append coordinate sets to given ag')
        assert_equal(coords, ag.getCoordsets(np.arange(ncsets, ncsets*2)))


    def testBiomolArgument(self):

        self.assertRaises(ValueError, parsePDB, self.one['path'],
                          biomol=True)


    def testSecondaryArgument(self):

        self.assertRaises(ValueError, parsePDB, self.one['path'],
                          secondary=True)

class TestWritePDB(unittest.TestCase):

    @dec.slow
    def setUp(self):
        """Set PDB file data and parse the PDB file."""

        self.pdb = DATA_FILES['multi_model_truncated']
        self.ag = parsePDB(self.pdb['path'])
        self.tmp = os.path.join(TEMPDIR, 'test.pdb')

    msg = 'user does not have write access to temp dir {0:s}'.format(TEMPDIR)

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testParsingOutput(self):
        """Test if parsing output is the same as parsing original file."""

        out = writePDB(self.tmp, self.ag)
        self.assertEqual(self.tmp, out,
            'writePDB failed to return correct output filename')
        self.assertTrue(os.path.isfile(out),
            'writePDB failed to write output')
        out = parsePDB(out)
        self.assertEqual(self.ag.numAtoms(), out.numAtoms(),
            'writePDB failed to write correct number of atoms')
        self.assertEqual(self.ag.numCoordsets(), out.numCoordsets(),
            'writePDB failed to write correct number of atoms')

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testModelArgument(self):
        """Test valid and invalid model arguments."""

        self.assertRaises(IndexError, writePDB, self.tmp, self.ag, csets='s')
        for i in range(self.ag.numCoordsets()):
            out = writePDB(self.tmp, self.ag, csets=i)
            out = parsePDB(out)
            self.assertEqual(out.numCoordsets(), 1,
                'failed to write correct number of models')
            assert_equal(out.getCoords(), self.ag.getCoordsets(i),
                 'failed to write model {0} coordinates correctly'.format(i+1))

    @dec.slow
    def tearDown(self):
        """Remove test file."""

        if os.path.isfile(self.tmp):
            os.remove(self.tmp)


class TestParsePDBHeaderAndAllModels(unittest.TestCase):

    def setUp(self):
        self.atomgroup, self.header = \
            parsePDB(pathDatafile('pdb2k39_truncated.pdb'), header=True)

    def testAtomGroupType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')
        self.assertIsInstance(self.atomgroup, prody.AtomGroup,
            'atom group type is incorrect')

    def testAtomGroupContent(self):

        self.assertEqual(self.atomgroup.numAtoms(), 167,
            'incorrect number of atoms')
        self.assertEqual(self.atomgroup.numCoordsets(), 3,
            'incorrect number of coordinate sets (models)')

    def tearDown(self):

        self.header = None
        self.atomgroup = None


class TestParsePDBAltloc(unittest.TestCase):

    def setUp(self):

        self.pdbfile = pathDatafile('pdb1ejg.pdb')

    def testAltlocNone(self):

        self.assertEqual(len(parsePDB(self.pdbfile)), 637,
            'failed to parse unspecified alternate locations correctly')

    def testAltlocA(self):

        self.assertEqual(len(parsePDB(self.pdbfile, altloc='A')), 637,
            'failed to parse alternate locations A correctly')

    def testAltlocB(self):

        self.assertEqual(len(parsePDB(self.pdbfile, altloc='B')), 634,
            'failed to parse alternate locations B correctly')

    def testAltlocC(self):

        self.assertEqual(len(parsePDB(self.pdbfile, altloc='C')), 496,
            'failed to parse alternate locations C correctly')
