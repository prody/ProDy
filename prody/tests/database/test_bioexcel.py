"""This module contains unit tests for :mod:`prody.database.bioexcel` module."""

from prody.tests import unittest
from prody.database.bioexcel import (fetchBioexcelPDB, parseBioexcelPDB,
                                     fetchBioexcelTrajectory, parseBioexcelTrajectory,
                                     convertXtcToDcd,
                                     fetchBioexcelTopology, parseBioexcelTopology)
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
    

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)

