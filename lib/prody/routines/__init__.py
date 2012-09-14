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

"""This module defines ProDy routines that are used as command line programs:

  * :func:`.prody_align`
  * :func:`.prody_contacts`
  * :func:`.prody_fetch`
  * :func:`.prody_select`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import imp
import sys
import argparse

from .actions import *

path_prody = imp.find_module('prody')[1]
path_routines = imp.find_module('routines', [path_prody])[1]

PRODY_COMMANDS = ['anm', 'gnm', 'pca', 'eda', 'align', 'blast', 'biomol', 
                  'catdcd', 'contacts', 'fetch', 'select', ] 

__all__ = ['prody_main']

prody_parser = argparse.ArgumentParser(
    description="ProDy: A Python Package for Protein Dynamics Analysis",
    epilog="See 'prody <command> -h' for more information on a specific "
           "command."
    )

prody_parser.add_argument('-c', '--cite', 
    help="print citation info and exit",
    action=ProDyCitation, nargs=0)

prody_parser.add_argument('-v', '--version', 
    help="print ProDy version and exit",
    action=ProDyVersion, nargs=0)

prody_commands = prody_parser.add_subparsers(
    title='subcommands')
    

for cmd in PRODY_COMMANDS:    
    cmd = 'prody_' + cmd
    mod = imp.load_module('prody.routines.' + cmd, 
                          *imp.find_module(cmd, [path_routines]))
    mod.addCommand(prody_commands)  


def prody_main():
    
    if len(sys.argv) == 1:    
        prody_parser.print_help()
    else:
        args = prody_parser.parse_args()
        args.func(args)

    
if __name__ == '__main__':
    prody_main()
