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

from actions import *

def prody_select(opt):
    """Write selected atoms from a PDB file in PDB format."""

    import prody
    LOGGER = prody.LOGGER
    pdb = prody.parsePDB(opt.pdb)
    prefix = opt.output
    if not prefix:
        prefix = pdb.getTitle() + '_selected'
    pdbselect = pdb.select(opt.selstr)
    if pdbselect is None:
        opt.subparser('Selection "{0:s}" do not match any atoms.'
                      .format(opt.selstr))
    LOGGER.info('Writing ' + prefix + '.pdb')
    prody.writePDB(prefix + '.pdb', pdbselect)


def addCommand(commands):

    subparser = commands.add_parser('select', 
    help='select atoms and write a PDB file')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """This command selects specified atoms and writes them in a PDB file.

Fetch PDB 2bfu and write backbone atoms in a file:
        
  $ prody select 2bfu "backbone" """
    )


    subparser.add_argument('-o', '--output', dest='output', type=str, 
        metavar='STR', help="output filanem (default: 'pdb_selected.pdb')")
        
    subparser.add_argument('pdb', help='PDB identifier or filename')
    subparser.add_argument('selstr', help='atom selection string')

    subparser.set_defaults(func=prody_select)
    subparser.set_defaults(subparser=subparser)
