from os import remove
import shlex
from os.path import isfile, join, split, splitext
from prody.tests import TestCase, skipIf, skipUnless

from numpy.testing import *

from prody.tests.datafiles import TEMPDIR, pathDatafile

from prody.apps import prody_parser

from prody.tests import MATPLOTLIB, NOPRODYCMD, WINDOWS


class TestANMCommand(TestCase):

    def setUp(self):

        self.command = ('anm -e -r -o {outdir} -v -z -t all '
                        '-f %8g -d , -x .dat '
                        '-R -Q '
                        '-F png -D 120 -W 5 -H 4 ').format(outdir=TEMPDIR)
        self.suffixes = [
            '_anm_cc.png',
            '_anm.anm.npz',
            '_anm_covariance.dat',
            '_anm_cross-correlations.dat',
            '_anm_evalues.dat',
            '_anm_evectors.dat',
            '_anm_sf.png',
            '_anm_extended_all.nmd',
            '_anm.nmd',
        ]

        self.tearDown()

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipUnless(MATPLOTLIB, 'matplotlib not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testANMCommand(self):

        pdb = pathDatafile('multi_model_truncated')
        command = self.command + pdb
        prefix = splitext(split(pdb)[1])[0]

        namespace = prody_parser.parse_args(shlex.split(command))
        namespace.func(namespace)

        for suffix in self.suffixes:
            fn = join(TEMPDIR, prefix + suffix)
            self.assertTrue(isfile(fn), msg=fn+' not found')

    def tearDown(self):

        for pdb in [pathDatafile('multi_model_truncated')]:
            prefix = splitext(split(pdb)[1])[0]
            for suffix in self.suffixes:
                fn = join(TEMPDIR, prefix + suffix)
                if isfile(fn):
                    remove(fn)
