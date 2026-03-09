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

class TestLRACommand(TestCase):

    def setUp(self):
        # We define the base command. Note that {pdb} is the REFERENCE structure.
        self.command = ('lra --pdb {pdb} '
                        '-e -r -o {outdir} -v -z -t all -j '
                        '-f %8g -d , -x .dat '
                        '-R -Q -J 1 -S '
                        '-F png -D 120 -W 5 -H 4 ').format(outdir=TEMPDIR,
                        pdb=pathDatafile('2k39_pdb_model0'))
        self.suffixes = [
            '_lra_cc.png',
            '_lra.lra.npz',
            '_lra_covariance.dat',
            '_lra_cross-correlations.dat',
            '_lra_proj_1.png',
            '_lra_evalues.dat',
            '_lra_proj.dat',
            '_lra_evectors.dat',
            '_lra_sf.png',
            '_lra_extended_all.nmd',
            '_lra.nmd',
        ]

        self.tearDown()

    # NOTE: I have kept the skip decorators commented out so you can verify the fix works.
    # In a real PR, you might want to uncomment them.
    # @dec.slow
    # @skipIf(NOPRODYCMD, 'prody command not found')
    # @skipUnless(MATPLOTLIB, 'matplotlib not found')
    # @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testLRACommandDCD(self):
        dcd = pathDatafile('2k39_insty_dcd')
        command = self.command + dcd
        prefix = splitext(split(dcd)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')

    # @dec.slow
    # @skipIf(NOPRODYCMD, 'prody command not found')
    # @skipUnless(MATPLOTLIB, 'matplotlib not found')
    # @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testLRACommandLabels(self):
        dcd = pathDatafile('2k39_insty_dcd')
        command = self.command + dcd + ' -l "[1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0]"'
        prefix = splitext(split(dcd)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')


    def tearDown(self):
        # Clean up the generated temp files as well as outputs
        files_to_clean = ['2k39_insty_dcd']

        for fname in files_to_clean:
            # Handle full paths vs datafiles
            base_path = pathDatafile(fname)
            prefix = splitext(split(base_path)[1])[0]

            for suffix in self.suffixes:
                fn = join(TEMPDIR, prefix + suffix)
                if isfile(fn):
                    remove(fn)
