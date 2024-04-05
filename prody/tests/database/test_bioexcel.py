"""This module contains unit tests for :mod:`prody.database.bioexcel` module."""

from prody.tests import unittest
from prody.database.bioexcel import (fetchBioexcelPDB, parseBioexcelPDB, convertXtcToDcd,
                                     fetchBioexcelTrajectory, parseBioexcelTrajectory,
                                     fetchBioexcelTopology, parseBioexcelTopology,
                                     checkSelection, checkQuery, checkConvert, 
                                     checkTimeout, checkFilePath, checkFrames, checkInputs)
import prody

import os
import shutil

from prody import LOGGER
LOGGER.verbosity = 'none'

FULL_N_ATOMS = 52350
SELE_N_ATOMS = 16860

class TestFetchParseBioexcelPDB(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.workdir = 'bioexcel_PDB_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

        cls.query = 'MCV1900370'
        cls.outname = 'outname'

    def testFetchDefault(self):
        """Test the outcome of a simple fetch scenario using 
        default options."""

        a = fetchBioexcelPDB(self.query, folder=self.workdir)

        self.assertIsInstance(a, str,
            'fetchBioexcelPDB failed to return a str instance')
        
        self.assertTrue(os.path.isfile(a),
                        'fetchBioexcelPDB failed to return a file')
        
        self.assertTrue(a.endswith('.pdb'),
                        'fetchBioexcelPDB failed to return a pdb file')
        
        self.assertEqual(a, os.path.join(self.workdir, self.query + '.pdb'),
                         'fetchBioexcelPDB default run did not give the right path')
        
        ag = prody.parsePDB(a)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePDB failed to return an AtomGroup from fetchBioexcelPDB')
        
        self.assertEqual(ag.numAtoms(), FULL_N_ATOMS, 
                         'fetchBioexcelPDB default output does not have correct number of atoms')
        
    def testFetchSelection(self):
        """Test the outcome of a simple fetch scenario 
        using selection='_C'."""

        a = fetchBioexcelPDB(self.query, folder=self.workdir,
                             selection='_C')
        
        ag = prody.parsePDB(a)
        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePDB failed to return an AtomGroup from fetchBioexcelPDB')
        self.assertEqual(ag.numAtoms(), SELE_N_ATOMS, 
                         'fetchBioexcelPDB selection _C output does not have correct number of atoms')
        
    def testFetchOutname(self):
        """Test the outcome of a simple fetch scenario using 
        outname='outname'."""

        a = fetchBioexcelPDB(self.query, folder=self.workdir,
                             outname=self.outname)

        self.assertEqual(a, os.path.join(self.workdir, self.outname + '.pdb'),
                         'fetchBioexcelPDB default run did not give the right path')

    def testParseDefault(self):
        """Test the outcome of a simple fetch and parse scenario 
        with default parameters."""

        ag = parseBioexcelPDB(self.query, folder=self.workdir)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parseBioexcelPDB failed to return an AtomGroup instance')
        
        self.assertEqual(ag.numAtoms(), FULL_N_ATOMS, 
                         'parseBioexcelPDB default output does not have correct number of atoms')
        
    def testParseSelection(self):
        """Test the outcome of a simple fetch and parse scenario 
        using selection='_C'."""

        ag = parseBioexcelPDB(self.query, folder=self.workdir,
                              selection='_C')
        
        self.assertIsInstance(ag, prody.AtomGroup,
            'parseBioexcelPDB with selection failed to return an AtomGroup')
        
        self.assertEqual(ag.numAtoms(), SELE_N_ATOMS, 
                         'parseBioexcelPDB selection _C output does not have correct number of atoms')
    
    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)


class TestFetchConvertParseBioexcelTop(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.workdir = 'bioexcel_top_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

        cls.query = 'MCV1900370'
        cls.outname = 'outname'

    def testFetchDefault(self):
        """Test the outcome of a simple fetch scenario using 
        default options."""

        a = fetchBioexcelTopology(self.query, folder=self.workdir)

        self.assertIsInstance(a, str,
            'fetchBioexcelTopology failed to return a str instance')
        
        self.assertTrue(os.path.isfile(a),
                        'fetchBioexcelTopology failed to return a file')
        
        self.assertTrue(a.endswith('.psf'),
                        'fetchBioexcelTopology default failed to return a psf file')
        
        self.assertEqual(a, os.path.join(self.workdir, self.query + '.psf'),
                         'fetchBioexcelTopology default run did not give the right path')
        
        ag = prody.parsePSF(a)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePSF failed to return an AtomGroup from fetchBioexcelTopology default')
        
        self.assertEqual(ag.numAtoms(), FULL_N_ATOMS, 
                         'fetchBioexcelTopology default output does not have correct number of atoms')
        
    def testFetchSelection(self):
        """Test the outcome of a simple fetch scenario 
        using selection='_C'."""

        a = fetchBioexcelTopology(self.query, folder=self.workdir,
                                  selection='_C')
        
        ag = prody.parsePSF(a)
        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePSF failed to return an AtomGroup from fetchBioexcelTopology')
        self.assertEqual(ag.numAtoms(), SELE_N_ATOMS, 
                         'fetchBioexcelTopology selection _C output does not have correct number of atoms')
        
    def testFetchOutname(self):
        """Test the outcome of a simple fetch scenario using 
        outname='outname'."""

        a = fetchBioexcelTopology(self.query, folder=self.workdir,
                             outname=self.outname)

        self.assertEqual(a, os.path.join(self.workdir, self.outname + '.psf'),
                         'fetchBioexcelPDB default run did not give the right path')

    def testFetchConvertFalse(self):
        """Test the outcome of a simple fetch scenario using 
        convert=False."""

        a = fetchBioexcelTopology(self.query, folder=self.workdir, convert=False)

        self.assertIsInstance(a, str,
            'fetchBioexcelTopology failed to return a str instance')
        
        self.assertTrue(os.path.isfile(a),
                        'fetchBioexcelTopology failed to return a file')
        
        self.assertTrue(a.endswith('.json'),
                        'fetchBioexcelTopology default failed to return a json file')
        
        self.assertEqual(a, os.path.join(self.workdir, self.query + '.json'),
                         'fetchBioexcelTopology default run did not give the right path')

    def testParseDefault(self):
        """Test the outcome of a simple fetch and parse scenario 
        with default parameters."""

        ag = parseBioexcelTopology(self.query, folder=self.workdir)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parseBioexcelTopology failed to return an AtomGroup instance')
        
        self.assertEqual(ag.numAtoms(), FULL_N_ATOMS, 
                         'parseBioexcelTopology default output does not have correct number of atoms')
        
    def testParseSelection(self):
        """Test the outcome of a simple fetch and parse scenario 
        using selection='_C'."""

        ag = parseBioexcelTopology(self.query, folder=self.workdir,
                                   selection='_C')
        
        self.assertIsInstance(ag, prody.AtomGroup,
            'parseBioexcelTopology with selection failed to return an AtomGroup')
        
        self.assertEqual(ag.numAtoms(), SELE_N_ATOMS, 
                         'parseBioexcelTopology selection _C output does not have correct number of atoms')

    def testFetchAndParse(self):
        """Test the outcome of a simple fetch and parse scenario"""

        a = fetchBioexcelTopology(self.query, folder=self.workdir)
        
        ag = parseBioexcelTopology(a, folder=self.workdir)
        
        self.assertIsInstance(ag, prody.AtomGroup,
            'fetch then parseBioexcelTopology failed to return an AtomGroup')
        
        self.assertEqual(ag.numAtoms(), FULL_N_ATOMS, 
                         'fetch then parseBioexcelTopology output does not have correct number of atoms')

    def testFetchConvParse(self):
        """Test the outcome of a simple fetch, convert and parse scenario."""

        a = fetchBioexcelTopology(self.query, folder=self.workdir, convert=False)
        
        ag = parseBioexcelTopology(a, folder=self.workdir)
        
        self.assertIsInstance(ag, prody.AtomGroup,
            'fetch, then convert & parseBioexcelTopology failed to return an AtomGroup')
        
        self.assertEqual(ag.numAtoms(), FULL_N_ATOMS, 
                         'fetch, then convert & parseBioexcelTopology output does not have correct number of atoms')

    def testConvertWrongType(self):
        with self.assertRaises(TypeError):
            fetchBioexcelTopology(self.query, folder=self.workdir, convert='False')

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)


class TestCheckSelection(unittest.TestCase):
    """Test that checkSelection gives the right errors and outputs."""
    
    def testWrongType(self):
        with self.assertRaises(TypeError):
            checkSelection(**{"selection": 1})

    def testWrongValue(self):
        with self.assertRaises(ValueError):
            checkSelection(**{"selection": '1'})

    def testNothing(self):
        self.assertIsNone(checkSelection(**{}))

    def testCarbon(self):
        self.assertEqual(checkSelection(**{'selection': '_C'}), '_C')

    def testBB(self):
        self.assertEqual(checkSelection(**{'selection': 'backbone'}), 'backbone')
        
    def testBoth(self):
        self.assertEqual(checkSelection(**{'selection': 'backbone and _C'}), 'backbone and _C')


class TestCheckQuery(unittest.TestCase):
    """Test that checkQuery gives the right errors and outputs."""

    @classmethod
    def setUpClass(cls):
        cls.query = 'MCV1900370'

    def testWrongType(self):
        with self.assertRaises(TypeError):
            checkQuery(1)

    def testCorrect(self):
        self.assertEqual(checkQuery(self.query), self.query)


class TestCheckConvert(unittest.TestCase):
    """Test that checkConvert gives the right errors and outputs."""

    def testWrongType(self):
        with self.assertRaises(TypeError):
            checkConvert(**{'convert': '1'})

    def testTrue(self):
        self.assertTrue(checkConvert(**{'convert': True}))

    def testFalse(self):
        self.assertFalse(checkConvert(**{'convert': False}))

    def testNothing(self):
        self.assertTrue(checkConvert(**{}))


class TestCheckTimeout(unittest.TestCase):
    """Test that checkTimeout gives the right errors and outputs."""

    def testWrongType(self):
        with self.assertRaises(TypeError):
            checkTimeout(**{'timeout': '1'})

    def testReplace(self):
        self.assertEqual(checkTimeout(**{'timeout': 50}), 50)

    def testDefault(self):
        self.assertEqual(checkTimeout(**{}), 60)

class TestCheckFilePath(unittest.TestCase):
    """Test that checkFilePath gives the right errors and outputs."""

    @classmethod
    def setUpClass(cls):
        cls.query = 'MCV1900370'
        cls.workdir = 'bioexcel_PDB_tests'
        cls.outname = 'outname'

        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)

    def testDefault(self):
        self.assertEqual(checkFilePath(self.query, **{}), 
                         os.path.join('.', self.query))

    def testWrongTypeFolder(self):
        self.assertEqual(checkFilePath(self.query, **{'folder': 1}), 
                         os.path.join('1', self.query))

    def testReplaceFolder(self):
        self.assertEqual(checkFilePath(self.query, **{'folder': self.workdir}), 
                         os.path.join(self.workdir, self.query))

    def testWrongTypeOutname(self):
        with self.assertRaises(TypeError):
            checkFilePath(self.query, **{'outname': 1})

    def testReplaceOutname(self):
        self.assertEqual(checkFilePath(self.query, **{'outname': self.outname}), 
                         os.path.join('.', self.outname))

    def testReplaceBoth(self):
        kwargs = {'outname': self.outname, 'folder': self.workdir}
        self.assertEqual(checkFilePath(self.query, **kwargs), 
                         os.path.join(self.workdir, self.outname))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.workdir)
