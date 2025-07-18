"""This module contains unit tests for :mod:`prody.database.bioexcel` module."""
import prody
if prody.PY3K:
    from prody.tests import unittest
    from prody.tests.datafiles import pathDatafile
    from prody.database.bioexcel import (fetchBioexcelPDB, parseBioexcelPDB, convertXtcToDcd,
                                         fetchBioexcelTrajectory, parseBioexcelTrajectory,
                                         fetchBioexcelTopology, parseBioexcelTopology,
                                         checkSelection, checkQuery, checkConvert,
                                         checkTimeout, checkFilePath, checkFrames)

    import os
    import shutil

    from prody import LOGGER
    LOGGER.verbosity = 'none'

    FULL_N_ATOMS = 12152
    SELE_N_ATOMS = 3908
    FULL_N_ATOMS_CV = 52350
    N_FRAMES_1 = 10
    N_FRAMES_2 = 6

    class TestFetchParseBioexcelPDB(unittest.TestCase):
        
        @classmethod
        def setUpClass(cls):
            cls.workdir = 'bioexcel_PDB_tests'
            if not os.path.exists(cls.workdir):
                os.mkdir(cls.workdir)
            os.chdir(cls.workdir)

            cls.query = 'A01Z9'
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
            """Test the outcome of a simple fetch scenario
            using outname='outname'."""

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

            cls.query = 'A01Z9'
            cls.outname = 'outname'

        def testFetchDefault(self):
            """Test the outcome of a simple fetch scenario
            using default options."""

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
            """Test the outcome of a simple fetch scenario
            using outname='outname'."""

            a = fetchBioexcelTopology(self.query, folder=self.workdir,
                                outname=self.outname)

            self.assertEqual(a, os.path.join(self.workdir, self.outname + '.psf'),
                            'fetchBioexcelPDB default run did not give the right path')

        def testFetchConvertFalse(self):
            """Test the outcome of a simple fetch scenario
            using convert=False."""

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
            """Test the outcome of a simple parse from file scenario
            with default parameters."""

            ag = parseBioexcelTopology(self.query, folder=self.workdir)

            self.assertIsInstance(ag, prody.AtomGroup,
                'parseBioexcelTopology failed to return an AtomGroup instance')
            
            self.assertEqual(ag.numAtoms(), FULL_N_ATOMS,
                            'parseBioexcelTopology default output does not have correct number of atoms')
            
        def testParseSelection(self):
            """Test the outcome of a simple parse from file scenario
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
            cls.query = 'A01Z9'

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
            self.assertEqual(checkTimeout(**{}), 200)

    class TestCheckFrames(unittest.TestCase):
        """Test that checkFrames gives the right errors and outputs."""

        def testWrongType(self):
            with self.assertRaises(TypeError):
                checkFrames(**{'frames': 1})

        def testDefault(self):
            self.assertIsNone(checkFrames(**{}))

        def testReplace1(self):
            self.assertEqual(checkFrames(**{'frames': '1-5,11-15'}), 
                            '1-5,11-15')

        def testReplace2(self):
            self.assertEqual(checkFrames(**{'frames': '10:20:2'}), 
                            '10:20:2')
            
        def testBadValue1(self):
            with self.assertRaises(ValueError):
                checkFrames(**{'frames': '1 '})

        def testBadValue2(self):
            with self.assertRaises(ValueError):
                checkFrames(**{'frames': '1, 1:2'})

        def testBadValue3(self):
            with self.assertRaises(ValueError):
                checkFrames(**{'frames': '1,1:2:3:4'})

        def testBadValue4(self):
            with self.assertRaises(ValueError):
                checkFrames(**{'frames': '1,1-2-3'})


    class TestCheckFilePath(unittest.TestCase):
        """Test that checkFilePath gives the right errors and outputs."""

        @classmethod
        def setUpClass(cls):
            cls.query = 'A01Z9'
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

    class TestFetchConvertParseBioexcelTraj(unittest.TestCase):
        
        @classmethod
        def setUpClass(cls):
            cls.workdir = 'bioexcel_traj_tests'
            if not os.path.exists(cls.workdir):
                os.mkdir(cls.workdir)
            os.chdir(cls.workdir)

            cls.query = 'A01Z9'
            cls.outname = 'outname'
            cls.frames1 = '1-5,11-15'
            cls.frames2 = '10:20:2'

        def testFetchFrames1(self):
            """Test the outcome of a simple fetch scenario
            using default options."""

            try:
                a = fetchBioexcelTrajectory(self.query, folder=self.workdir,
                                            frames=self.frames1)
            except OSError:
                pass
            else:
                self.assertIsInstance(a, str,
                    'fetchBioexcelTrajectory failed to return a str instance')
                
                self.assertTrue(os.path.isfile(a),
                                'fetchBioexcelTrajectory failed to return a file')
                
                self.assertTrue(a.endswith('.dcd'),
                                'fetchBioexcelTrajectory default failed to return a dcd file')
                
                self.assertEqual(a, os.path.join(self.workdir, self.query + '.dcd'),
                                'fetchBioexcelTrajectory default run did not give the right path')

            ens = prody.parseDCD(a)

            self.assertIsInstance(ens, prody.Ensemble,
                'parseDCD failed to return an Ensemble from fetchBioexcelTrajectory default')
            self.assertEqual(ens.numAtoms(), FULL_N_ATOMS,
                            'fetchBioexcelTrajectory default output does not have correct number of atoms')
            self.assertEqual(ens.numCoordsets(), N_FRAMES_1,
                            'fetchBioexcelTrajectory output with example frames 1 does not have correct number of frames')

        def testFetchSelectionFrames2(self):
            """Test the outcome of a simple fetch scenario
            using selection='_C'."""

            try:
                a = fetchBioexcelTrajectory(self.query, folder=self.workdir,
                                            selection='_C', frames=self.frames2)
            except OSError:
                pass
            else:
                ens = prody.parseDCD(a)
                self.assertIsInstance(ens, prody.Ensemble,
                    'parseDCD failed to return an Ensemble from fetchBioexcelTrajectory')
                self.assertEqual(ens.numAtoms(), SELE_N_ATOMS,
                                'fetchBioexcelTrajectory selection _C output does not have correct number of atoms')
                self.assertEqual(ens.numCoordsets(), N_FRAMES_2,
                                'fetchBioexcelTrajectory output with example frames 2 does not have correct number of frames')

        def testFetchConvertFalse(self):
            """Test the outcome of a simple fetch scenario
            using convert=False."""

            try:
                a = fetchBioexcelTrajectory(self.query, folder=self.workdir,
                                            convert=False, frames=self.frames1)
            except OSError:
                pass
            else:
                self.assertIsInstance(a, str,
                    'fetchBioexcelTrajectory failed to return a str instance')
                
                self.assertTrue(os.path.isfile(a),
                                'fetchBioexcelTrajectory failed to return a file')
                
                self.assertTrue(a.endswith('.xtc'),
                                'fetchBioexcelTrajectory default failed to return a xtc file')
                
                self.assertEqual(a, os.path.join(self.workdir, self.query + '.xtc'),
                                'fetchBioexcelTrajectory default run did not give the right path')

        def testParseFrames1(self):
            """Test the outcome of a simple parse from file scenario
            with default parameters."""

            try:
                ens = parseBioexcelTrajectory(self.query, folder=self.workdir,
                                              frames=self.frames1)
            except OSError:
                pass
            else:
                self.assertIsInstance(ens, prody.Ensemble,
                    'parseBioexcelTrajectory failed to return an Ensemble instance')
                self.assertEqual(ens.numAtoms(), FULL_N_ATOMS,
                                'parseBioexcelTrajectory default output does not have correct number of atoms')
                self.assertEqual(ens.numCoordsets(), N_FRAMES_1,
                                'parseBioexcelTrajectory output with example frames 1 does not have correct number of frames')

        def testParseSelectionFrames2(self):
            """Test the outcome of a simple parse from file scenario
            using selection='_C'."""
            try:
                ens = parseBioexcelTrajectory(self.query, folder=self.workdir,
                                              selection='_C', frames=self.frames2)
            except OSError:
                pass
            else:
                self.assertIsInstance(ens, prody.Ensemble,
                    'parseBioexcelTrajectory with selection failed to return an Ensemble')
                self.assertEqual(ens.numAtoms(), SELE_N_ATOMS,
                                'parseBioexcelTrajectory selection _C output does not have correct number of atoms')
                self.assertEqual(ens.numCoordsets(), N_FRAMES_2,
                                'parseBioexcelTrajectory output with example frames 2 does not have correct number of frames')

        def testFetchAndParse(self):
            """Test the outcome of a simple fetch and parse scenario"""
            try:
                a = fetchBioexcelTrajectory(self.query, folder=self.workdir,
                                            frames=self.frames1)
            except OSError:
                pass
            else:
                ens = parseBioexcelTrajectory(a, folder=self.workdir)
                
                self.assertIsInstance(ens, prody.Ensemble,
                    'parseBioexcelTrajectory failed to return an Ensemble instance')
                self.assertEqual(ens.numAtoms(), FULL_N_ATOMS,
                                'parseBioexcelTrajectory default output does not have correct number of atoms')
                self.assertEqual(ens.numCoordsets(), N_FRAMES_1,
                                'parseBioexcelTrajectory output with example frames 1 does not have correct number of frames')

        def testFetchNoConvParse(self):
            """Test the outcome of a simple fetch, then internally convert and parse scenario."""
            try:
                a = fetchBioexcelTrajectory(self.query, folder=self.workdir,
                                            convert=False, frames=self.frames1)
            except OSError:
                pass
            else:
                ens = parseBioexcelTrajectory(a)
                
                self.assertIsInstance(ens, prody.Ensemble,
                    'parseBioexcelTrajectory failed to return an Ensemble instance')
                self.assertEqual(ens.numAtoms(), FULL_N_ATOMS,
                                'parseBioexcelTrajectory default output does not have correct number of atoms')
                self.assertEqual(ens.numCoordsets(), N_FRAMES_1,
                                'parseBioexcelTrajectory output with example frames 1 does not have correct number of frames')

        def testFetchConvParse(self):
            """Test the outcome of a simple fetch, externally convert and then parse scenario."""
            try:
                a = fetchBioexcelTrajectory(self.query, folder=self.workdir,
                                            convert=False, frames=self.frames1)
            except OSError:
                pass
            else:
                b = convertXtcToDcd(a)
                ens = parseBioexcelTrajectory(b)
                
                self.assertIsInstance(ens, prody.Ensemble,
                    'parseBioexcelTrajectory failed to return an Ensemble instance')
                self.assertEqual(ens.numAtoms(), FULL_N_ATOMS,
                                'parseBioexcelTrajectory default output does not have correct number of atoms')
                self.assertEqual(ens.numCoordsets(), N_FRAMES_1,
                                'parseBioexcelTrajectory output with example frames 1 does not have correct number of frames')

        def testConvertWrongType(self):
            with self.assertRaises(TypeError):
                fetchBioexcelTrajectory(self.query, folder=self.workdir,
                                        convert='False')

        @classmethod
        def tearDownClass(cls):
            os.chdir('..')
            shutil.rmtree(cls.workdir)

    class TestOnlyConvertParseBioexcel(unittest.TestCase):
        
        @classmethod
        def setUpClass(cls):
            cls.query = 'MCV1900370'
            cls.psfPath = pathDatafile(cls.query + '.psf')
            cls.xtcPath = pathDatafile(cls.query + '.xtc')
            cls.dcdPath = pathDatafile(cls.query + '.dcd')

            cls.jsonPath = pathDatafile('MCV1900193.json')
            cls.PROTEIN_GLYCAN_N_ATOMS = 72759
            cls.CA_N_ATOMS = 3768

        def testParseBioexcelTop(self):
            ag = parseBioexcelTopology(self.psfPath)
            self.assertIsInstance(ag, prody.AtomGroup,
                'parseBioexcelTopology failed to return an AtomGroup from data files')
            self.assertEqual(ag.numAtoms(), FULL_N_ATOMS_CV,
                            'parseBioexcelTopology data files output does not have correct number of atoms')

        def testParseBioexcelTopJsonGlycan(self):
            ag = parseBioexcelTopology(self.jsonPath)
            self.assertIsInstance(ag, prody.AtomGroup,
                'parseBioexcelTopology failed to return an AtomGroup from data files')
            self.assertEqual(ag.numAtoms(), self.PROTEIN_GLYCAN_N_ATOMS, 
                            'parseBioexcelTopology data files output using MCV1900193 with glycans does not have correct number of atoms')
            self.assertEqual(ag.ca.numAtoms(), self.CA_N_ATOMS, 
                            'parseBioexcelTopology data files output using MCV1900193 with glycans does not have correct number of CA atoms')
            
        def testConvertToDCD(self):
            a = convertXtcToDcd(self.xtcPath, top=self.psfPath)
            self.assertTrue(os.path.isfile(a),
                            'convertXtcToDcd failed to return a file')
            self.assertTrue(a.endswith('.dcd'),
                            'convertXtcToDcd output file does not end with .dcd')

        def testParseConvertBioexcelTraj(self):
            ens = parseBioexcelTrajectory(self.xtcPath, top=self.psfPath)
            self.assertIsInstance(ens, prody.Ensemble,
                'parseBioexcelTrajectory failed to return an Ensemble from xtc and psf data files')
            self.assertEqual(ens.numAtoms(), FULL_N_ATOMS_CV,
                            'parseBioexcelTrajectory output from xtc and psf data files does not have correct number of atoms')
            self.assertEqual(ens.numCoordsets(), N_FRAMES_2,
                            'parseBioexcelTrajectory output from xtc and psf data files does not have correct number of frames')

        def testOnlyParseBioexcelTraj(self):
            ens = parseBioexcelTrajectory(self.dcdPath, top=self.psfPath)
            self.assertIsInstance(ens, prody.Ensemble,
                'parseBioexcelTrajectory failed to return an Ensemble from xtc and psf data files')
            self.assertEqual(ens.numAtoms(), FULL_N_ATOMS_CV,
                            'parseBioexcelTrajectory output from xtc and psf data files does not have correct number of atoms')
            self.assertEqual(ens.numCoordsets(), N_FRAMES_2,
                            'parseBioexcelTrajectory output from xtc and psf data files does not have correct number of frames')
