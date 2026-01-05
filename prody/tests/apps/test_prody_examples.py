import os
import glob
import shlex
import os.path
from prody.tests import TestCase, skipIf, skipUnless
from subprocess import Popen, PIPE

from numpy.testing import *

from prody.utilities import importDec
dec = importDec()

from prody import LOGGER
LOGGER.verbosity = None

from prody.apps.prody_apps import prody_commands
from prody.tests.datafiles import TEMPDIR

from prody.tests import MATPLOTLIB, NOPRODYCMD, WINDOWS

TESTDIR = os.path.join(TEMPDIR, 'prody_tests')


class TestCommandExamples(TestCase):

    def setUp(self):
        self.cwd = os.getcwd()
        if not os.path.isdir(TESTDIR):
            os.mkdir(TESTDIR)
        os.chdir(TESTDIR)

    def tearDown(self):
        if os.path.isdir(TESTDIR):
            for fn in glob.glob(os.path.join(TESTDIR, '*')):
                try:
                    os.remove(fn)
                except OSError:
                    pass
        os.chdir(self.cwd)


count = 0
for cmd in prody_commands.choices:

    parser = prody_commands.choices[cmd]

    test_examples = parser._defaults.get('test_examples', [])

    usage = parser._defaults['usage_example']

    examples = []
    for line in usage.split('\n'):
        line = line.strip()
        if line.startswith('$'):
            examples.append(line[1:].strip())

    for i in test_examples:
        try:
            egs = [examples[i]]
        except TypeError:
            egs = [examples[x] for x in i]
        count += 1

        @dec.slow
        @skipIf(NOPRODYCMD, 'prody command not found')
        @skipIf(WINDOWS, 'examples are not tested in Windows')
        def func(self, examples=egs):

            for eg in examples:
                # --- FIX APPLIED ---
                # 1. stdin=PIPE allows us to send input to the subprocess
                # 2. stdout/stderr=PIPE captures output
                pipe = Popen(shlex.split(eg + ' --quiet'),
                             stdout=PIPE, stderr=PIPE, stdin=PIPE)
                
                # communicate() reads both stdout and stderr buffers simultaneously
                # preventing the buffer fill deadlock.
                # input=b'n\n' sends "n" + Enter to the process. 
                # If the app asks "Overwrite? [y/n]", this answers "no" and unblocks it.
                stdout, stderr = pipe.communicate(input=b'n\n')

        func.__name__ = 'testCommandExample{0:d}'.format(count)
        func.__doc__ = 'Test example: $ {0:s}'.format(' $ '.join(egs))
        setattr(TestCommandExamples, func.__name__, func)

del func
