#!/usr/bin/python
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import imp
from subprocess import Popen, PIPE
from os.path import join

path = [imp.find_module('prody')[1]]
routines = imp.load_module('prody.routines', 
                           *imp.find_module('routines', path))


pipe = Popen(['prody', '-h'], stdout=PIPE, stderr=PIPE)
with open('prody.txt', 'w') as rst:
    rst.write(pipe.stdout.read())

for cmd in routines.PRODY_COMMANDS: 
    
    #pkg = getattr('routines', 'prody_' + cmd)
    
    with open(join('prody_commands', cmd + '.rst'), 'w') as rst:
        rst.write(""".. _prody-{cmd:s}:

*******************************************************************************
prody {cmd:s}
*******************************************************************************

Usage
===============================================================================

Running :command:`prody {cmd} -h` displays::

""".format(cmd=cmd))
    
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
