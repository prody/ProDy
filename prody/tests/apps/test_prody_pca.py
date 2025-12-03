from os import remove
import shlex
from os.path import isfile, join, split, splitext
from prody.tests import TestCase, skipIf, skipUnless

from numpy.testing import *

from prody.utilities import importDec
dec = importDec()

from prody.tests.datafiles import TEMPDIR, pathDatafile

from prody.apps import prody_parser

from prody.tests import MATPLOTLIB, NOPRODYCMD, WINDOWS

# NEW: Import ProDy functions to generate valid test data
from prody import parsePDB, writePDB, writeDCD

class TestPCACommand(TestCase):

    def setUp(self):
        # We define the base command. Note that {pdb} is the REFERENCE structure.
        self.command = ('pca --pdb {pdb} '
                        '-e -r -o {outdir} -v -z -t all -j '
                        '-f %8g -d , -x .dat '
                        '-R -Q -J 1,2 '
                        '-F png -D 120 -W 5 -H 4 ').format(outdir=TEMPDIR,
                        pdb=pathDatafile('multi_model_truncated'))
        self.suffixes = [
            '_pca_cc.png',
            '_pca.pca.npz',
            '_pca_covariance.dat',
            '_pca_cross-correlations.dat',
            '_pca_proj_1_2.png',
            '_pca_evalues.dat',
            '_pca_proj.dat',
            '_pca_evectors.dat',
            '_pca_sf.png',
            '_pca_extended_all.nmd',
            '_pca.nmd',
        ]

        self.tearDown()

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipUnless(MATPLOTLIB, 'matplotlib not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testPCACommandDCD(self):

        # 1. Load the original data
        pdb_path = pathDatafile('multi_model_truncated')
        atoms = parsePDB(pdb_path)
        
        # 2. Duplicate coordsets to ensure we have > 3 frames (ProDy requirement)
        while atoms.numCoordsets() <= 3:
            atoms.addCoordset(atoms.getCoordsets())

        # 3. Write a single temporary DCD file
        temp_dcd = join(TEMPDIR, 'sufficient_frames.dcd')
        writeDCD(temp_dcd, atoms)

        # 4. Pass the SINGLE new file to the command
        command = self.command + temp_dcd
        prefix = splitext(split(temp_dcd)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipUnless(MATPLOTLIB, 'matplotlib not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testPCACommandPDB(self):

        # 1. Load original data
        pdb_path = pathDatafile('multi_model_truncated')
        atoms = parsePDB(pdb_path)

        # 2. Duplicate coordsets to ensure > 3 frames
        while atoms.numCoordsets() <= 3:
            atoms.addCoordset(atoms.getCoordsets())

        # 3. Write a single temporary PDB file
        temp_pdb = join(TEMPDIR, 'sufficient_frames.pdb')
        writePDB(temp_pdb, atoms)

        # 4. Pass the SINGLE new file to the command
        command = self.command + temp_pdb
        prefix = splitext(split(temp_pdb)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')


    def tearDown(self):
        # Clean up the generated temp files as well as outputs
        files_to_clean = [
            'multi_model_truncated',
            'dcd',
            'sufficient_frames.pdb', # Clean up our temp input
            'sufficient_frames.dcd'  # Clean up our temp input
        ]

        for fname in files_to_clean:
            # Handle full paths vs datafiles
            if 'sufficient' in fname:
                prefix = splitext(fname)[0]
                base_path = join(TEMPDIR, fname) # The input file itself
                if isfile(base_path): remove(base_path)
            else:
                base_path = pathDatafile(fname)
                prefix = splitext(split(base_path)[1])[0]

            for suffix in self.suffixes:
                fn = join(TEMPDIR, prefix + suffix)
                if isfile(fn): remove(fn)
