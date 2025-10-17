# """This module contains unit tests for :mod:`prody.database.pfam` module."""

from prody.tests import unittest
from prody.database.pfam import searchPfam
from prody.database.pfam import fetchPfamMSA
from prody.database.pfam import parsePfamPDBs

from prody.atomic.selection import Selection

import os
import shutil

from prody import LOGGER
LOGGER.verbosity = 'none'

class TestSearchPfam(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.workdir = 'pfam_search_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

        cls.queries = ['P19491', '6qkcB', '6qkcI', 'PF00047',
                       'hellow', 'hello']

    def testUniprotAccMulti(self):
        """Test the outcome of a simple search scenario using a Uniprot Accession
        for a multi-domain protein, AMPAR GluA2."""

        a = searchPfam(self.queries[0])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())),
                         ['PF00060', 'PF01094', 'PF10613'],
                         'searchPfam failed to return the right domain family IDs')
        
    def testPdbIdChMulti(self):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID for the same multi-domain protein from specifying chain B."""

        a = searchPfam(self.queries[1])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), ['PF00060', 'PF01094', 'PF10613'],
                         'searchPfam failed to return the right domain family IDs for AMPAR')

    def testPdbIdChSingle(self):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID to get the single domain protein TARP g8 from chain I."""

        a = searchPfam(self.queries[2])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')

        self.assertEqual(sorted(list(a.keys())),
                         ['PF00822', 'PF13903'],
                         'searchPfam failed to return the right domain family IDs for TARP')

    def testPfamInput(self):
        """Test the outcome of a search scenario where a Pfam ID is
        provided as input."""

        a = searchPfam(self.queries[3])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return None for Pfam ID input {0}'.format(self.queries[3]))

    def testWrongInput1(self):
        """Test the outcome of a search scenario where a 6-char text is
        provided as input."""

        with self.assertRaises(OSError):
            searchPfam(self.queries[4])

    def testWrongInput2(self):
        """Test the outcome of a search scenario where a 5-char text is
        provided as input."""

        with self.assertRaises(ValueError):
            searchPfam(self.queries[5])

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)


class TestFetchPfamMSA(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.query = 'PF00822'

        cls.workdir = 'pfam_msa_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

    def testDefault(self):
        """Test the outcome of fetching the domain MSA for claudins
        with default parameters."""

        b = fetchPfamMSA(self.query)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))


    def testSeed(self):
        """Test the outcome of fetching the domain MSA for claudins
        with the alignment type argument set to seed"""

        b = fetchPfamMSA(self.query, "seed")

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))

    def testFolder(self):
        """Test the outcome of fetching the domain MSA for claudins
        with keyword folder set to a folder that is made especially."""

        folder = "new_folder"
        os.mkdir(folder)
        b = fetchPfamMSA(self.query, folder=folder)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'new_folder/PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))
    
    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)
            

class TestParsePfamPDBs(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.queries = ['PF20446', 'Q57ZF2', 'P40682']

        cls.workdir = 'pfam_pdb_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

    def testPfamIdDefault(self):
        """Test the outcome of parsing PDBs for a tiny family
        of ABC class ATPase N-terminal domains (5 members)
        with the Pfam ID and default parameters."""

        b = parsePfamPDBs(self.queries[0])

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 5,
            'parsePfamPDBs failed to return a list of length 5')


    def testUniprotDefault(self):
        """Test the outcome of parsing PDBs for a tiny family
        of ABC class ATPase N-terminal domains (5 members)
        with the Uniprot long ID and default parameters."""

        b = parsePfamPDBs(self.queries[1])

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 5,
            'parsePfamPDBs failed to return a list of length 5')

        
    def testMultiDomainDefault(self):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1,
        which has two domains but few relatives. Default parameters should
        return Selection objects containing the first domain."""

        b = parsePfamPDBs(self.queries[2])

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(b[0].getResnums()[0], 264,
            'parsePfamPDBs failed to return a first Selection with first resnum 264')

    def testMultiDomainStart1(self):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1,
        which has two domains but few relatives. Using start=1 should be like default and
        return Selection objects containing the first domain."""

        b = parsePfamPDBs(self.queries[2], start=1)

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(b[0].getResnums()[0], 264,
            'parsePfamPDBs failed to return a first Selection with first resnum 264')
        
    def testMultiDomainStart2(self):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1,
        which has two domains but few relatives. Setting start to 418 should
        return Selection objects containing the second domain."""

        b = parsePfamPDBs(self.queries[2], start=418)

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(b[0].getResnums()[0], 217,
            'parsePfamPDBs failed to return a first Selection with first resnum 217')

    def testPfamIdNumPdbs(self):
        """Test the outcome of parsing PDBs for a tiny family
        of ABC class ATPase N-terminal domains (5 members)
        with the Pfam ID and default parameters."""

        b = parsePfamPDBs(self.queries[0], num_pdbs=2)

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 2,
            'parsePfamPDBs failed to return a list of length 2 with num_pdbs=2')

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)

