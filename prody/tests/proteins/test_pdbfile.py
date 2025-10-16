"""This module contains unit tests for :mod:`~prody.proteins`."""

import os

import numpy as np
from numpy.testing import assert_equal, assert_allclose

from prody.utilities import importDec
dec = importDec()

import prody
from prody import LOGGER
from prody.proteins import parsePDB, writePDB, PDBParseError
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import parseDatafile, pathDatafile, DATA_FILES

LOGGER.verbosity = 'none'

class TestParsePDB(unittest.TestCase):

    def setUp(self):
        """Set PDB file data."""

        self.pdb = DATA_FILES['multi_model_truncated']
        self.one = DATA_FILES['oneatom']
        self.ca = DATA_FILES['1ubi_ca']

        self.ficlve_dig = DATA_FILES['five_digits']
        self.hex = DATA_FILES['hex']
        self.h36 = DATA_FILES['h36']

        self.probes = DATA_FILES['probes']

        self.altlocs = DATA_FILES['6flr']
        self.his_selstr = 'resname HIS and chain B and resnum 234 and name CA'

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
        self.assertRaises(PDBParseError, parsePDB, path,
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
        self.assertEqual(parsePDB(path, title=name).getTitle(),
             str(name), 'parsePDB failed to set user given non-string name')

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
        """Test outcome of 2 invalid and 2 valid *ag* arguments."""

        path = pathDatafile(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, ag='AtomGroup')

        ag = prody.AtomGroup('One atom')
        ag.setCoords(np.array([[0., 0., 0.]]))
        self.assertRaises(ValueError, parsePDB, path, ag=ag)

        ag = prody.AtomGroup('Test')
        self.assertEqual(parsePDB(path, ag=ag).numAtoms(),
            self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms')
        
        ag = prody.AtomGroup('Test')
        ag.setCoords(np.array([[0., 0., 0.]]*self.pdb['atoms']))
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

    def testFiveDigitResnum(self):
        """Test parsing PDB files with five digit resnums
        (using integer insertion code)."""

        path = pathDatafile(self.five_dig['file'])
        resnum = '10000'
        self.assertEqual(str(parsePDB(path).getResnums()[33107]),
             resnum, 'parsePDB failed to parse five digit resnum')

    def testHexResnum(self):
        """Test parsing PDB files with hexadecimal resnums
        ('2710' is 10000)."""

        path = pathDatafile(self.hex['file'])
        resnum = '10000'
        self.assertEqual(str(parsePDB(path).getResnums()[33107]),
             resnum, 'parsePDB failed to parse hexadecimal resnum')

    def testHybrid36Resnum(self):
        """Test parsing PDB files with Hybrid36 resnums
        ('A000' is 10000)."""

        path = pathDatafile(self.h36['file'])
        resnum = '10000'
        self.assertEqual(str(parsePDB(path).getResnums()[33107]),
             resnum, 'parsePDB failed to parse Hybrid36 resnum')

    def testHexSerial(self):
        """Test parsing PDB files with hexadecimal serial numbers
        ('2710' is 10000)."""

        path = pathDatafile(self.hex['file'])
        serial = '100000'
        self.assertEqual(str(parsePDB(path).getSerials()[100000-1]),
             serial, 'parsePDB failed to hexadecimal serial number')

    def testHybrid36Serial(self):
        """Test parsing PDB files with Hybrid36 serial numbers
        ('A0000' is 100000)."""

        path = pathDatafile(self.h36['file'])
        serial = '100000'
        self.assertEqual(str(parsePDB(path).getSerials()[100000-1]),
             serial, 'parsePDB failed to parse Hybrid36 serial number')
        
    def testAltlocAllToMultiAtoms(self):
        """Test number of coordinate sets and atoms with altloc='all'."""

        path = pathDatafile(self.altlocs['file'])

        ag = parsePDB(path, altloc="all")
        self.assertEqual(ag.numAtoms(), self.altlocs['atoms_altloc'],
            'parsePDB failed to parse correct number of atoms with altloc "all"')
        self.assertEqual(ag.numCoordsets(), 1,
            'parsePDB failed to parse correct number of coordsets (1) with altloc "all"')

        hisB234 = ag.select(self.his_selstr)
        self.assertEqual(hisB234.numAtoms(), self.altlocs['num_altlocs'],
            'parsePDB failed to parse correct number of His B234 CA atoms (2) with altloc "all"')

        self.assertEqual(hisB234.getAnisous().shape, (self.altlocs['num_altlocs'], 6),
            'parsePDB failed to have right shape for His B234 CA atoms getAnisous (2, 6) with altloc "all"')

        assert_allclose(hisB234.getAnisous()[0], self.altlocs['anisousA'][0],
            err_msg='parsePDB failed to have right His B234 CA atoms getAnisous A with altloc "all"')

        assert_allclose(hisB234.getAnisous()[1], self.altlocs['anisousB'][0],
            err_msg='parsePDB failed to have right His B234 CA atoms getAnisous B with altloc "all"')
        
    def testAltlocNoneToMultiCoordsets(self):
        """Test number of coordinate sets and atoms with altloc=None."""

        path = pathDatafile(self.altlocs['file'])

        ag = parsePDB(path, altloc=None)
        self.assertEqual(ag.numAtoms(), self.altlocs['atoms_single'],
            'parsePDB failed to parse correct number of atoms with altloc None')
        self.assertEqual(ag.numCoordsets(), self.altlocs['num_altlocs'],
            'parsePDB failed to parse correct number of coordsets (2) with altloc None')

        hisB234 = ag.select(self.his_selstr)
        self.assertEqual(hisB234.numAtoms(), 1,
            'parsePDB failed to parse correct number of His B234 CA atoms (1) with altloc None')

        self.assertEqual(hisB234.getAnisous().shape, (1, 6),
            'parsePDB failed to have right shape for His B234 CA atoms getAnisous (1, 6) with altloc None')

        assert_allclose(hisB234.getAnisous(), self.altlocs['anisousA'],
            err_msg='parsePDB failed to have right His B234 CA atoms getAnisous A with altloc None')

        hisB234.setACSIndex(1)

        assert_allclose(hisB234.getAnisous(), self.altlocs['anisousB'],
            err_msg='parsePDB failed to have right His B234 CA atoms getAnisous B with altloc None')

    def testParseLongResnameFalse(self):
        """Test the outcome of parsing with clipping long resnames."""

        ag = parseDatafile(self.probes['file'])

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePDB failed to return an AtomGroup instance')

        self.assertEqual(ag.getResnames()[0],
                         self.probes['short_resname'])

    def testParseLongResnameTrue(self):
        """Test the outcome of parsing with clipping long resnames."""

        ag = parseDatafile(self.probes['file'], long_resname=True)

        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePDB failed to return an AtomGroup instance')

        self.assertEqual(ag.getResnames()[0],
                         self.probes['long_resname'])


class TestWritePDB(unittest.TestCase):

    @dec.slow
    def setUp(self):
        """Set PDB file data and parse the PDB file."""

        self.pdb = DATA_FILES['multi_model_truncated']
        self.ag = parsePDB(self.pdb['path'])
        self.tmp = os.path.join(TEMPDIR, 'test.pdb')

        self.ens = prody.PDBEnsemble()
        self.ens.setAtoms(self.ag)
        self.ens.setCoords(self.ag.getCoords())
        self.ens.addCoordset(self.ag.getCoordsets())

        self.ubi = parsePDB(DATA_FILES['1ubi']['path'], secondary=True)
        self.probes = parsePDB(DATA_FILES['probes']['path'], long_resname=True)

        self.hex = parsePDB(DATA_FILES['hex']['path'])
        self.h36 = parsePDB(DATA_FILES['h36']['path'])

        self.hex_ter = parsePDB(DATA_FILES['hex_ter']['path'])
        self.h36_ter = parsePDB(DATA_FILES['h36_ter']['path'])

        self.altlocs = DATA_FILES['6flr']
        self.altloc_full = parsePDB(self.altlocs['path'], altloc=None,
                                    secondary=False)
        self.altloc_sel = DATA_FILES['6flr_sel']['path']

        self.sort_sel = DATA_FILES['6zu5_sel']
        self.sort_sel_ag = parsePDB(self.sort_sel['path'])

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
            'writePDB failed to write correct number of coordsets')

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingHex(self):
        """Test if output from writing hexadecimal is as expected."""

        out = writePDB(self.tmp, self.hex)
        fi = open(out, 'r')
        lines = fi.readlines()
        fi.close()

        resnum_9999_line = lines[33107]
        self.assertEqual(resnum_9999_line[22:26], '9999',
            'writePDB failed to write correct pre-hex resnum')

        resnum_2710_line = lines[33108]
        self.assertEqual(resnum_2710_line[22:26], '2710',
            'writePDB failed to write correct hex resnum')

        serial_99999_line = lines[99999]
        self.assertEqual(serial_99999_line[6:11], '99999',
            'writePDB failed to write correct pre-hex serial')        

        serial_186a0_line = lines[100000]
        self.assertEqual(serial_186a0_line[6:11], '186a0',
            'writePDB failed to write correct hex serial')  

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingHybrid36(self):
        """Test if output from writing Hybrid36 is as expected."""

        out = writePDB(self.tmp, self.h36, hybrid36=True)
        fi = open(out, 'r')
        lines = fi.readlines()
        fi.close()

        resnum_9999_line = lines[33107]
        self.assertEqual(resnum_9999_line[22:26], '9999',
            'writePDB failed to write correct pre-h36 resnum')

        resnum_A000_line = lines[33108]
        self.assertEqual(resnum_A000_line[22:26], 'A000',
            'writePDB failed to write correct h36 resnum')

        serial_99999_line = lines[99999]
        self.assertEqual(serial_99999_line[6:11], '99999',
            'writePDB failed to write correct pre-h36 serial')

        serial_A0000_line = lines[100000]
        self.assertEqual(serial_A0000_line[6:11], 'A0000',
            'writePDB failed to write correct h36 serial')
        
    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingSecstrs(self):
        """Test if output from writing secstrs is as expected."""

        out = writePDB(self.tmp, self.ubi)
        ubi_new = parsePDB(out, secondary=True)
        self.assertListEqual(list(self.ubi.getSecstrs()), list(ubi_new.getSecstrs()))

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingLongResnums(self):
        """Test if output from writing secstrs is as expected."""

        out = writePDB(self.tmp, self.ubi)
        ubi_new = parsePDB(out, secondary=True)
        self.assertListEqual(list(self.ubi.getSecstrs()), list(ubi_new.getSecstrs()))

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
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingHexTer(self):
        """Test if output from writing hexadecimal with TER lines is as expected."""

        out = writePDB(self.tmp, self.hex_ter)
        fi = open(out, 'r')
        lines = fi.readlines()
        fi.close()

        pre_ter_line = lines[4060]
        self.assertEqual(pre_ter_line[22:26], ' 550',
            'writePDB failed to write correct pre-ter resnum')
        self.assertEqual(pre_ter_line[6:11], ' 4060',
            'writePDB failed to write correct pre-ter serial')

        post_ter_line = lines[4062]
        self.assertEqual(post_ter_line[22:26], '   4',
            'writePDB failed to write correct post-ter resnum')
        self.assertEqual(post_ter_line[6:11], ' 4062',
            'writePDB failed to write correct post-ter serial')
        
        serial_99999_line = lines[99999]
        self.assertEqual(serial_99999_line[6:11], '99999',
            'writePDB failed to write correct pre-hex serial')

        serial_186a0_line = lines[100000]
        self.assertEqual(serial_186a0_line[6:11], '186a0',
            'writePDB failed to write correct hex serial')

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingHybrid36Ter(self):
        """Test if output from writing Hybrid36 with TER lines is as expected."""

        out = writePDB(self.tmp, self.h36_ter, hybrid36=True)
        fi = open(out, 'r')
        lines = fi.readlines()
        fi.close()

        pre_ter_line = lines[4060]
        self.assertEqual(pre_ter_line[22:26], ' 550',
            'writePDB failed to write correct pre-ter resnum')
        self.assertEqual(pre_ter_line[6:11], ' 4060',
            'writePDB failed to write correct pre-ter serial')

        post_ter_line = lines[4062]
        self.assertEqual(post_ter_line[22:26], '   4',
            'writePDB failed to write correct post-ter resnum')
        self.assertEqual(post_ter_line[6:11], ' 4062',
            'writePDB failed to write correct post-ter serial')

        serial_99999_line = lines[99999]
        self.assertEqual(serial_99999_line[6:11], '99999',
            'writePDB failed to write correct pre-h36 serial')        

        serial_A0000_line = lines[100000]
        self.assertEqual(serial_A0000_line[6:11], 'A0000',
            'writePDB failed to write correct h36 serial')
        
    def testWritingAltlocModels(self):
        """Test if output from writing hexadecimal with TER lines is as expected."""

        hisB234 = self.altloc_full.select('resname HIS and chain B and resnum 234 and name CA')
        out = writePDB(self.tmp, hisB234)

        fi = open(out, 'r')
        lines1 = fi.readlines()
        fi.close()

        fi = open(self.altloc_sel, 'r')
        lines2 = fi.readlines()
        fi.close()
        
        self.assertEqual(lines1[3], lines2[3],
            'writePDB failed to write correct ANISOU line 3 for 6flr selection with altloc None')
        
        self.assertEqual(lines1[7], lines2[7],
            'writePDB failed to write correct ANISOU line 7 for 6flr selection with altloc None')
        
    def testWriteEnsembleToPDB(self):
        """Test that writePDB can handle ensembles."""

        out = writePDB(self.tmp, self.ens)
        out = parsePDB(out)
        self.assertEqual(out.numCoordsets(), self.ens.numCoordsets(),
            'failed to write correct number of models from ensemble')
        assert_equal(out.getCoords(), self.ag.getCoordsets(0),
                'failed to write ensemble model 1 coordinates correctly')

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingAtomMap(self):
        """Test if output from writing a sorted AtomMap works and is as expected."""

        sorted_sel = prody.sortAtoms(self.sort_sel_ag, 'chain')
        out = writePDB(self.tmp, sorted_sel)
        new = parsePDB(out)
        self.assertListEqual(list(new.getChids()), self.sort_sel['sorted_order'])

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingParsingLongResnames(self):
        """Test if parsing output is the same as parsing original file."""

        out = writePDB(self.tmp, self.probes, long_resname=True)
        self.assertEqual(self.tmp, out,
            'writePDB failed to return correct output filename')
        self.assertTrue(os.path.isfile(out),
            'writePDB failed to write output')

        probes = parsePDB(out, long_resname=True)
        self.assertEqual(self.probes.numAtoms(), probes.numAtoms(),
            'writePDB failed to write correct number of atoms')
        self.assertEqual(str(self.probes.getResnames()[0]), str(probes.getResnames()[0]),
            'writePDB failed to write long resnames')

    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testWritingParsingLongResnamesFalse(self):
        """Test if parsing output is the same as parsing original file."""

        out = writePDB(self.tmp, self.probes)
        self.assertEqual(self.tmp, out,
            'writePDB failed to return correct output filename')
        self.assertTrue(os.path.isfile(out),
            'writePDB failed to write output')

        probes = parsePDB(out, long_resname=True)
        self.assertEqual(self.probes.numAtoms(), probes.numAtoms(),
            'writePDB failed to write correct number of atoms')
        self.assertEqual(probes.getResnames()[0], self.probes.getResnames()[0][:3],
            'writePDB failed to write short resnames')

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
