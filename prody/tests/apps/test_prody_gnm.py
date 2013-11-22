from os import remove
import shlex
from os.path import isfile, join, split, splitext
from prody.tests import TestCase, skipIf, skipUnless

from numpy.testing import *

from prody.tests.datafiles import TEMPDIR, pathDatafile

from prody.apps import prody_parser

from prody.tests import MATPLOTLIB, NOPRODYCMD, WINDOWS

class TestGNMCommand(TestCase):

    def setUp(self):

        self.command = ('gnm -e -r -o {outdir} -v -z -t all '
                        '-f %8g -d , -x .dat '
                        '-R -Q '
                        '-F png -D 120 -W 5 -H 4 ').format(outdir=TEMPDIR)
        self.suffixes = [
            '_gnm_cc.png',
            '_gnm.gnm.npz',
            '_gnm_covariance.dat',
            '_gnm_cross-correlations.dat',
            '_gnm_evalues.dat',
            '_gnm_evectors.dat',
            '_gnm_sf.png',
            '_gnm_extended_all.nmd',
            '_gnm.nmd',
        ]

        self.tearDown()

    @dec.slow
    @skipIf(NOPRODYCMD, 'prody command not found')
    @skipUnless(MATPLOTLIB, 'matplotlib not found')
    @skipIf(WINDOWS, 'command tests are not run on Windows')
    def testGNMCommand(self):

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
                if isfile(fn): remove(fn)
