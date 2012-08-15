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

def prody_select(**kwargs):
    """Write selected atoms from a PDB file in PDB format."""

    from prody import LOGGER, parsePDB, writePDB
    
    pdb = kwargs.get('pdb', None)
    if not pdb:
        raise ValueError('pdb argument must be provided')
    pdb = parsePDB(pdb)
        
    selstr = kwargs.get('selstr')
    if not selstr:
        raise ValueError('selstr argument must be provided')
        
    pdbselect = pdb.select(selstr)
    if pdbselect is None:
        LOGGER.warn('Selection {0:s} did not match any atoms.'
                    .format(repr(selstr)))
        return
    LOGGER.info('Selection {0:s} matched {1:d} atoms.'
                .format(repr(selstr), len(pdbselect)))

    prefix = kwargs.get('output', '')
    if not prefix:
        prefix = pdb.getTitle() + '_selected'

    ext = '.pdb'
    outname = prefix + (ext if not ext in prefix else '')
    LOGGER.info('Writing selection into ' + outname)
    writePDB(outname, pdbselect)


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
        
  $ prody select 2bfu "backbone" """)


    subparser.add_argument('-o', '--output', dest='output', metavar='STR', 
        type=str, help="output PDB filename (default: 'pdb_selected.pdb')")
        
    subparser.add_argument('pdb', help='PDB identifier or filename')
    subparser.add_argument('selstr', help='atom selection string')

    subparser.set_defaults(func=lambda ns: prody_select(**ns.__dict__))
    subparser.set_defaults(subparser=subparser)
