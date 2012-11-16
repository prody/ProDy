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

"""This module defines some routines used as command line programs."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import imp
import sys

try:
    import argparse
except ImportError:
    from .. import argparse

from ..actions import *

path_prody = imp.find_module('prody')[1]
path_routines = imp.find_module('routines', [path_prody])[1]
path_commands = imp.find_module('evol_routines', [path_routines])[1]

EVOL_COMMANDS = ['search', 'fetch'] 

#if sys.version_info[:2] > (2,6):
#    EVOL_COMMANDS.append('test')

__all__ = ['evol_main']

evol_parser = argparse.ArgumentParser(
    description="Evol: A Python Package for Sequence Evolution",
    epilog="See 'evol <command> -h' for more information on a specific "
           "command."
    )

evol_parser.add_argument('-c', '--cite', 
    help="print citation info and exit",
    action=ProDyCitation, nargs=0)

evol_parser.add_argument('-v', '--version', 
    help="print ProDy version and exit",
    action=ProDyVersion, nargs=0)

evol_commands = evol_parser.add_subparsers(
    title='subcommands')
    

for cmd in EVOL_COMMANDS:
    cmd = 'evol_' + cmd
    mod = imp.load_module('prody.routines.evol_routines.' + cmd, 
                          *imp.find_module(cmd, [path_commands]))
    mod.APP.addApplication(evol_commands)


def evol_main():
    
    if len(sys.argv) == 1:    
        evol_parser.print_help()
    else:
        namespace = evol_parser.parse_args()
        namespace.func(namespace)

    
if __name__ == '__main__':
    evol_main()
