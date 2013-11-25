#!/usr/bin/python

import os
import imp
from subprocess import Popen, PIPE

path = [imp.find_module('prody')[1]]
apps = imp.load_module('prody.apps',
                           *imp.find_module('apps', path))

for cmd, subcmds in [('prody', apps.PRODY_APPS), ('evol', apps.EVOL_APPS)]:

    pipe = Popen([cmd, '-h'], stdout=PIPE, stderr=PIPE)
    with open(os.path.join(cmd, cmd + '.txt'), 'w') as rst:
        rst.write(pipe.stdout.read())

    for sub in subcmds:

        with open(os.path.join(cmd, sub + '.rst'), 'w') as rst:
            rst.write(""".. _{cmd:s}-{sub:s}:

{cmd:s} {sub:s}
====================

Usage
--------------------

Running :command:`{cmd:s} {sub:s} -h` displays::

""".format(cmd=cmd, sub=sub))

            pipe = Popen([cmd , sub, '-h'], stdout=PIPE, stderr=PIPE)
            for line in pipe.stdout.readlines():
                rst.write('  ' + line)

            rst.write("""
Examples
--------------------

Running :command:`{cmd:s} {sub:s} --examples` displays::

""".format(cmd=cmd, sub=sub))


            pipe = Popen([cmd, sub, '--examples'], stdout=PIPE, stderr=PIPE)
            for line in pipe.stdout.readlines():
                rst.write('  ' + line)
