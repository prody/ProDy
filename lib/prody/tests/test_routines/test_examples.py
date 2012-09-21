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
import shlex
from unittest import TestCase
from subprocess import Popen, PIPE, call 

from numpy.testing import *

from prody.routines import prody_commands
from prody.tests.test_datafiles import TEMPDIR



class TestCommandsMeta(type):
    
    def __init__(cls, name, bases, dict): 
        
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
                def testFunction(self, examples=egs):
                    
                    for eg in examples:
                        pipe = Popen(shlex.split(eg + ' --quiet'),
                                     stdout=PIPE, stderr=PIPE)
                        stdout = pipe.stdout.read()
                        pipe.stdout.close() 
                        stderr = pipe.stderr.read()
                        pipe.stderr.close()
                        self.assertTrue(stderr == '', msg=stderr)
                    
                testFunction.__name__ = 'testCommandExample{0:d}'.format(
                                                                        count)
                testFunction.__doc__ = 'Test example: $ {0:s}'.format(
                                            ' $ '.join(egs))
                setattr(cls, testFunction.__name__, testFunction)
            

class TestCommandExamples(TestCase):    
    
    __metaclass__ = TestCommandsMeta


    def setUp(self):
        
        self.cwd = os.getcwd()
        os.chdir(TEMPDIR)
        
    def tearDown(self):
        
        os.chdir(self.cwd)


