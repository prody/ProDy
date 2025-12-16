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
from prody import parsePDB, Trajectory


class TestPCACommand(TestCase):

    def setUp(self):

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

        dcd = pathDatafile('dcd')

        # Skip if there are not enough coordinate sets for PCA
        try:
            traj = Trajectory(dcd)
            if traj.numFrames() <= 3:
                self.skipTest('coordsets must have more than 3 coordinate sets')
        except Exception:
            # If trajectory cannot be read for any reason, skip to avoid false failures
            self.skipTest('unable to read trajectory for PCA test')

        command = self.command + dcd
        prefix = splitext(split(dcd)[1])[0]

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

        pdb = pathDatafile('multi_model_truncated')

        # Skip if there are not enough coordinate sets for PCA
        try:
            atoms = parsePDB(pdb, model='all')
            if atoms is None or atoms.numCoordsets() <= 3:
                self.skipTest('coordsets must have more than 3 coordinate sets')
        except Exception:
            # If PDB cannot be parsed for any reason, skip to avoid false failures
            self.skipTest('unable to parse PDB for PCA test')

        command = self.command + pdb
        prefix = splitext(split(pdb)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')


    def tearDown(self):

        for dcd in [pathDatafile('multi_model_truncated'),
                    pathDatafile('dcd')]:
            prefix = splitext(split(dcd)[1])[0]
            for suffix in self.suffixes:
                fn = join(TEMPDIR, prefix + suffix)
                if isfile(fn): remove(fn)
