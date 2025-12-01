# """This module contains unit tests for :mod:`prody.database.pfam` module."""

from prody.tests import unittest
from unittest.mock import patch, MagicMock
from io import BytesIO

from prody.database.pfam import searchPfam
from prody.database.pfam import fetchPfamMSA
from prody.database.pfam import parsePfamPDBs

from prody.atomic.selection import Selection

import os
import shutil
import numpy as np

from prody import LOGGER
LOGGER.verbosity = 'none'

# --- MOCK DATA ---
# Updated to support multiple domains for P40682
# PF_DOM1 -> Corresponds to the first domain (starts ~1 in sequence, PDB res 264)
# PF_DOM2 -> Corresponds to the second domain (starts ~418 in sequence, PDB res 217)
MOCK_PFAM_MAPPING = b"""Pfam Mapping File Generated
PDB\tCHAIN\tPDB_START\tPDB_END\tPFAM_ACC\tPFAM_NAME\tPFAM_DESC\teValue
1MK1\tA\t10\t50\tPF20446\tTinyFam\tTinyDesc\t1.0
1MK2\tB\t10\t50\tPF20446\tTinyFam\tTinyDesc\t1.0
1MK3\tC\t10\t50\tPF20446\tTinyFam\tTinyDesc\t1.0
1MK4\tD\t10\t50\tPF20446\tTinyFam\tTinyDesc\t1.0
1MK5\tE\t10\t50\tPF20446\tTinyFam\tTinyDesc\t1.0
2XYZ\tA\t264\t300\tPF_DOM1\tFirstDom\tDesc\t1.0
3XYZ\tA\t217\t300\tPF_DOM2\tSecondDom\tDesc\t1.0
"""

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
        a = searchPfam(self.queries[0])
        self.assertIsInstance(a, dict, 'searchPfam failed to return a dict instance')
        self.assertEqual(sorted(list(a.keys())),
                         ['PF00060', 'PF01094', 'PF10613'])
        
    def testPdbIdChMulti(self):
        a = searchPfam(self.queries[1])
        self.assertIsInstance(a, dict, 'searchPfam failed to return a dict instance')
        self.assertEqual(sorted(list(a.keys())), ['PF00060', 'PF01094', 'PF10613'])

    def testPdbIdChSingle(self):
        a = searchPfam(self.queries[2])
        self.assertIsInstance(a, dict, 'searchPfam failed to return a dict instance')
        self.assertEqual(sorted(list(a.keys())), ['PF00822', 'PF13903'])

    def testPfamInput(self):
        a = searchPfam(self.queries[3])
        self.assertIsInstance(a, dict)

    def testWrongInput1(self):
        with self.assertRaises(OSError):
            searchPfam(self.queries[4])

    def testWrongInput2(self):
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
        b = fetchPfamMSA(self.query)
        self.assertIsInstance(b, str)
        self.assertEqual(b, 'PF00822_seed.sth')
        self.assertTrue(os.path.exists(b))

    def testSeed(self):
        b = fetchPfamMSA(self.query, "seed")
        self.assertIsInstance(b, str)
        self.assertEqual(b, 'PF00822_seed.sth')
        self.assertTrue(os.path.exists(b))

    def testFolder(self):
        folder = "new_folder"
        os.mkdir(folder)
        b = fetchPfamMSA(self.query, folder=folder)
        self.assertIsInstance(b, str)
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

    def setUp(self):
        # 1. Patch FTP (The mapping file download)
        self.patcher_ftp = patch('ftplib.FTP')
        self.mock_ftp = self.patcher_ftp.start()
        
        self.mock_ftp_instance = self.mock_ftp.return_value
        def write_fake_data(command, callback):
            callback(MOCK_PFAM_MAPPING)
        self.mock_ftp_instance.retrbinary.side_effect = write_fake_data

        # 2. Patch parsePDB (The structure download)
        self.patcher_pdb = patch('prody.database.pfam.parsePDB')
        self.mock_parsePDB = self.patcher_pdb.start()

        self.mock_ag = MagicMock(spec=Selection)
        self.mock_selection = MagicMock(spec=Selection)
        self.mock_ag.select.return_value = self.mock_selection

        # --- KEY FIX: Dynamic return values ---
        # Instead of returning a hardcoded list of 5, we look at how many IDs were requested
        # This prevents the IndexError when the test only asks for 1 or 2 items
        def fake_parse_pdb(pdb_ids, **kwargs):
            if isinstance(pdb_ids, list):
                count = len(pdb_ids)
            else:
                count = 1
            return ([self.mock_ag] * count, ['H'] * count)

        self.mock_parsePDB.side_effect = fake_parse_pdb

        # 3. Patch searchPfam (The ID lookup)
        self.patcher_search = patch('prody.database.pfam.searchPfam')
        self.mock_search = self.patcher_search.start()

        def fake_search(query, **kwargs):
            if query == 'P40682': 
                # Simulate a protein with two domains:
                # PF_DOM1: starts at residue 1
                # PF_DOM2: starts at residue 418
                return {
                    'PF_DOM1': {'locations': [{'start': 1, 'end': 100}]},
                    'PF_DOM2': {'locations': [{'start': 418, 'end': 500}]}
                }
            else:
                return {'PF20446': {'locations': [{'start': 10, 'end': 50}]}}
        
        self.mock_search.side_effect = fake_search

    def tearDown(self):
        self.patcher_ftp.stop()
        self.patcher_pdb.stop()
        self.patcher_search.stop()

    def testPfamIdDefault(self):
        b = parsePfamPDBs(self.queries[0])
        self.assertIsInstance(b, list)
        self.assertEqual(len(b), 5) 

    def testUniprotDefault(self):
        b = parsePfamPDBs(self.queries[1])
        self.assertIsInstance(b, list)
        self.assertEqual(len(b), 5)

    def testMultiDomainDefault(self):
        """Expects PF_DOM1 (starts at 264 in PDB)"""
        self.mock_selection.getResnums.return_value = [264]
        
        b = parsePfamPDBs(self.queries[2]) # start defaults to 1 -> matches PF_DOM1
        
        self.assertIsInstance(b, list)
        self.assertEqual(b[0].getResnums()[0], 264)

    def testMultiDomainStart1(self):
        """Expects PF_DOM1 (starts at 264 in PDB)"""
        self.mock_selection.getResnums.return_value = [264]
        
        b = parsePfamPDBs(self.queries[2], start=1) # Matches PF_DOM1
        
        self.assertIsInstance(b, list)
        self.assertEqual(b[0].getResnums()[0], 264)
        
    def testMultiDomainStart2(self):
        """Expects PF_DOM2 (starts at 217 in PDB)"""
        self.mock_selection.getResnums.return_value = [217]

        b = parsePfamPDBs(self.queries[2], start=418) # Matches PF_DOM2
        
        self.assertIsInstance(b, list)
        self.assertEqual(b[0].getResnums()[0], 217)

    def testPfamIdNumPdbs(self):
        # We don't need to override the mock here anymore, 
        # because the smart 'side_effect' in setUp will see num_pdbs=2 
        # (which limits the pdb_ids list) and return 2 items automatically.
        b = parsePfamPDBs(self.queries[0], num_pdbs=2)

        self.assertIsInstance(b, list)
        self.assertEqual(len(b), 2)

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)
