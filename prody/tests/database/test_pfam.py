# """This module contains unit tests for :mod:`prody.database.pfam` module."""

from prody.tests import unittest
import os
import shutil
from unittest.mock import patch, Mock

from prody import LOGGER
LOGGER.verbosity = 'none'

# Import test utilities
from prody.tests.database.test_utils import (
    check_pfam_connectivity,
    create_mock_requests_get,
    create_mock_fetchPfamMSA,
    create_mock_ftp_for_pfam_pdbs,
    create_mock_parsePDBHeader,
    create_mock_pfam_search
)

# Check connectivity once at module level
USE_FIXTURES = not check_pfam_connectivity(timeout=3)

# Import the pfam functions
from prody.database.pfam import searchPfam
from prody.database.pfam import fetchPfamMSA
from prody.database.pfam import parsePfamPDBs

from prody.atomic.selection import Selection
from ftplib import FTP

# If using fixtures, we'll replace the functions at module level later in each test class

class TestSearchPfam(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.workdir = 'pfam_search_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

        cls.queries = ['P19491', '6qkcB', '6qkcI', 'PF00047',
                       'hellow', 'hello']
        
        # If using fixtures, replace searchPfam with mock version
        if USE_FIXTURES:
            cls.original_searchPfam = searchPfam
            # Replace with mock in the module
            import prody.database.pfam
            prody.database.pfam.searchPfam = create_mock_pfam_search(use_fixtures=True)

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)
        
        # Restore original if we replaced it
        if USE_FIXTURES and hasattr(cls, 'original_searchPfam'):
            import prody.database.pfam
            prody.database.pfam.searchPfam = cls.original_searchPfam

    def testUniprotAccMulti(self):
        """Test the outcome of a simple search scenario using a Uniprot Accession
        for a multi-domain protein, AMPAR GluA2."""

        # Call from module to get the mocked version if USE_FIXTURES
        if USE_FIXTURES:
            import prody.database.pfam
            a = prody.database.pfam.searchPfam(self.queries[0], timeout=5)
        else:
            a = searchPfam(self.queries[0], timeout=5)

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())),
                         ['PF00060', 'PF01094', 'PF10613'],
                         'searchPfam failed to return the right domain family IDs')
        
    def testPdbIdChMulti(self):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID for the same multi-domain protein from specifying chain B."""

        # Call from module to get the mocked version if USE_FIXTURES
        if USE_FIXTURES:
            import prody.database.pfam
            a = prody.database.pfam.searchPfam(self.queries[1], timeout=5)
        else:
            a = searchPfam(self.queries[1], timeout=5)

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), ['PF00060', 'PF01094', 'PF10613'],
                         'searchPfam failed to return the right domain family IDs for AMPAR')

    def testPdbIdChSingle(self):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID to get the single domain protein TARP g8 from chain I."""

        # Call from module to get the mocked version if USE_FIXTURES
        if USE_FIXTURES:
            import prody.database.pfam
            a = prody.database.pfam.searchPfam(self.queries[2], timeout=5)
        else:
            a = searchPfam(self.queries[2], timeout=5)

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')

        self.assertEqual(sorted(list(a.keys())),
                         ['PF00822', 'PF13903'],
                         'searchPfam failed to return the right domain family IDs for TARP')

    def testPfamInput(self):
        """Test the outcome of a search scenario where a Pfam ID is
        provided as input."""

        # Call from module to get the mocked version if USE_FIXTURES
        if USE_FIXTURES:
            import prody.database.pfam
            a = prody.database.pfam.searchPfam(self.queries[3], timeout=5)
        else:
            a = searchPfam(self.queries[3], timeout=5)

        self.assertIsInstance(a, dict,
            'searchPfam failed to return None for Pfam ID input {0}'.format(self.queries[3]))

    def testWrongInput1(self):
        """Test the outcome of a search scenario where a 6-char text is
        provided as input."""

        with self.assertRaises((OSError, FileNotFoundError)):
            # Call from module to get the mocked version if USE_FIXTURES
            if USE_FIXTURES:
                import prody.database.pfam
                prody.database.pfam.searchPfam(self.queries[4], timeout=5)
            else:
                searchPfam(self.queries[4], timeout=5)

    def testWrongInput2(self):
        """Test the outcome of a search scenario where a 5-char text is
        provided as input."""

        with self.assertRaises(ValueError):
            # Call from module to get the mocked version if USE_FIXTURES
            if USE_FIXTURES:
                import prody.database.pfam
                prody.database.pfam.searchPfam(self.queries[5], timeout=5)
            else:
                searchPfam(self.queries[5], timeout=5)



class TestFetchPfamMSA(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.query = 'PF00822'

        cls.workdir = 'pfam_msa_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)
        
        # If using fixtures, replace fetchPfamMSA with mock version
        if USE_FIXTURES:
            cls.original_fetchPfamMSA = fetchPfamMSA
            # Replace with mock in the module
            import prody.database.pfam
            prody.database.pfam.fetchPfamMSA = create_mock_fetchPfamMSA(use_fixtures=True)

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)
        
        # Restore original if we replaced it
        if USE_FIXTURES and hasattr(cls, 'original_fetchPfamMSA'):
            import prody.database.pfam
            prody.database.pfam.fetchPfamMSA = cls.original_fetchPfamMSA

    def testDefault(self):
        """Test the outcome of fetching the domain MSA for claudins
        with default parameters."""

        # Call from module to get the mocked version if USE_FIXTURES
        if USE_FIXTURES:
            import prody.database.pfam
            b = prody.database.pfam.fetchPfamMSA(self.query, timeout=5)
        else:
            b = fetchPfamMSA(self.query, timeout=5)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))


    def testSeed(self):
        """Test the outcome of fetching the domain MSA for claudins
        with the alignment type argument set to seed"""

        # Call from module to get the mocked version if USE_FIXTURES
        if USE_FIXTURES:
            import prody.database.pfam
            b = prody.database.pfam.fetchPfamMSA(self.query, "seed", timeout=5)
        else:
            b = fetchPfamMSA(self.query, "seed", timeout=5)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))

    def testFolder(self):
        """Test the outcome of fetching the domain MSA for claudins
        with keyword folder set to a folder that is made especially."""

        folder = "new_folder"
        os.mkdir(folder)
        
        # Call from module to get the mocked version if USE_FIXTURES
        if USE_FIXTURES:
            import prody.database.pfam
            b = prody.database.pfam.fetchPfamMSA(self.query, folder=folder, timeout=5)
        else:
            b = fetchPfamMSA(self.query, folder=folder, timeout=5)

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
        
        # Set up mock for FTP if using fixtures
        if USE_FIXTURES:
            MockFTP = create_mock_ftp_for_pfam_pdbs(use_fixtures=True)
            cls.ftp_patcher = patch('ftplib.FTP', MockFTP)
            cls.ftp_patcher.start()

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)
        
        # Stop the patcher if it was started
        if USE_FIXTURES and hasattr(cls, 'ftp_patcher'):
            cls.ftp_patcher.stop()

    def testPfamIdDefault(self):
        """Test the outcome of parsing PDBs for a tiny family
        of ABC class ATPase N-terminal domains (5 members)
        with the Pfam ID and default parameters."""

        b = parsePfamPDBs(self.queries[0], timeout=5)

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

        # This test requires searchPfam which needs fixtures
        if USE_FIXTURES:
            # Skip this test when using fixtures as it requires complex setup
            self.skipTest("Skipping Uniprot test with fixtures (requires searchPfam mock)")
        
        b = parsePfamPDBs(self.queries[1], timeout=5)

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

        # This test requires searchPfam which needs fixtures
        if USE_FIXTURES:
            # Skip this test when using fixtures as it requires complex setup
            self.skipTest("Skipping multi-domain test with fixtures (requires searchPfam mock)")

        b = parsePfamPDBs(self.queries[2], timeout=5)

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

        # This test requires searchPfam which needs fixtures
        if USE_FIXTURES:
            # Skip this test when using fixtures as it requires complex setup
            self.skipTest("Skipping multi-domain start=1 test with fixtures (requires searchPfam mock)")

        b = parsePfamPDBs(self.queries[2], start=1, timeout=5)

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

        # This test requires searchPfam which needs fixtures
        if USE_FIXTURES:
            # Skip this test when using fixtures as it requires complex setup
            self.skipTest("Skipping multi-domain start=418 test with fixtures (requires searchPfam mock)")

        b = parsePfamPDBs(self.queries[2], start=418, timeout=5)

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

        b = parsePfamPDBs(self.queries[0], num_pdbs=2, timeout=5)

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 2,
            'parsePfamPDBs failed to return a list of length 2 with num_pdbs=2')

