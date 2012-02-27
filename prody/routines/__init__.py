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

"""This module contains functions which are used as command line programs."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import sys
import os.path
import textwrap

import argparse

from actions import *

__all__ = ['main']

parser = argparse.ArgumentParser(
    description="ProDy: A Python Package for Protein Dynamics Analysis",
    epilog="See 'prody <command> -h' for more information on a specific "
           "command."
    )

parser.add_argument('-c', '--cite', help="print citation info and exit",
    action=ProDyCitation, nargs=0)

parser.add_argument('-v', '--version', help="print ProDy version and exit",
    action=ProDyVersion, nargs=0)

commands = parser.add_subparsers(
    title='subcommands')
    
    
for cmd in ['prody_anm', 'prody_gnm', 'prody_pca', 'prody_align',
           'prody_blast', 'prody_biomol', 'prody_catdcd', 'prody_fetch',  
           'prody_select', ]:    
    pkg = __import__(cmd, globals(), locals(), [], -1)
    pkg.addCommand(commands)

def main():
    
    if len(sys.argv) == 1:    
        parser.print_help()
    else:
        args = parser.parse_args()
        args.func(args)

    
if __name__ == '__main__':
    main()
