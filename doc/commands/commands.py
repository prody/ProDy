#!/usr/bin/python
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

from os import path
from glob import glob
from subprocess import Popen, PIPE

import prody
from prody import routines

pipe = Popen(['prody', '-h'], stdout=PIPE, stderr=PIPE)
with open('prody.txt', 'w') as rst:
    rst.write(pipe.stdout.read())

for cmd in routines.PRODY_COMMANDS: 
    
    pkg = getattr(getattr(prody, 'routines'), 'prody_' + cmd)
    
    with open('prody_' + cmd + '.rst', 'w') as rst:
        rst.write(""".. _prody-{cmd:s}:

*******************************************************************************
prody {cmd:s}
*******************************************************************************

{doc:s}

Usage
===============================================================================

Running :command:`prody {cmd} -h` displays::

""".format(cmd=cmd, doc=pkg.__doc__))
    
        pipe = Popen(['prody', cmd, '-h'], stdout=PIPE, stderr=PIPE)
        for line in pipe.stdout.readlines():    
            rst.write('  ' + line)

        rst.write("""
Examples
===============================================================================

Running :command:`prody {cmd:s} --examples` displays::

""".format(cmd=cmd, doc='deneme'))
    

        pipe = Popen(['prody', cmd, '--examples'], stdout=PIPE, stderr=PIPE)
        for line in pipe.stdout.readlines():    
            rst.write('  ' + line)
