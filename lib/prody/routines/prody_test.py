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

"""Extract a selection of atoms from a PDB file."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from .actions import *

__all__ = ['prody_test']

def prody_test(label='fast'):
    """Write selected atoms from a PDB file in PDB format."""

    import prody
    prody.test(label=label)


def addCommand(commands):

    subparser = commands.add_parser('test', 
    help='run ProDy unittests')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.set_defaults(usage_example=
    """This command selects specified atoms and writes them in a PDB file.

Fetch PDB files 1p38 and 1r39 and write backbone atoms in a file:
        
  $ prody select backbone 1p38 1r39""",
    test_examples=[0])


    group = subparser.add_argument_group('output options')
    
    group.add_argument('-l', '--label', dest='label', metavar='STR', 
        type=str, help='output PDB filename (default: pdb_selected.pdb)',
        default='fast')
        
    subparser.set_defaults(func=lambda ns: prody_test(label=ns.label))
    subparser.set_defaults(subparser=subparser)
