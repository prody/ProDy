# ProDy: A Python Package for Protein Dynamics Analysis
#
# Copyright (C) 2010-2012 Ahmet Bakan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'


import os
import glob
import shlex
import os.path
from prody.tests import TestCase, skipIf, skipUnless
from subprocess import Popen, PIPE

from numpy.testing import *

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
                os.remove(fn)
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
                pipe = Popen(shlex.split(eg + ' --quiet'),
                             stdout=PIPE, stderr=PIPE)
                stdout = pipe.stdout.read()
                pipe.stdout.close()
                stderr = pipe.stderr.read()
                pipe.stderr.close()
                self.assertTrue(not stderr, msg=stderr)

        func.__name__ = 'testCommandExample{0:d}'.format(count)
        func.__doc__ = 'Test example: $ {0:s}'.format(' $ '.join(egs))
        setattr(TestCommandExamples, func.__name__, func)

del func
