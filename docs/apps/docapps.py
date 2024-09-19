#!/usr/bin/python

import os
import sys
import importlib
from subprocess import Popen, PIPE

from prody.utilities.misctools import impLoadModule

path_apps = importlib.util.find_spec("prody.apps").submodule_search_locations[0]

prody_apps = impLoadModule('prody.apps.prody_apps', path_apps + '/prody_apps/', '__init__')
evol_apps = impLoadModule('prody.apps.evol_apps', path_apps + '/evol_apps/', '__init__')

for cmd, subcmds in [('prody', prody_apps.PRODY_APPS), ('evol', evol_apps.EVOL_APPS)]:

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
