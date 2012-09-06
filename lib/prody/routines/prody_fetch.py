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

"""Download PDB files for given identifiers."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

from actions import *

def prody_fetch(opt):
    """Fetch PDB files from PDB FTP server."""
    
    import prody
    pdblist = opt.pdb
    if len(pdblist) == 1 and os.path.isfile(pdblist[0]):
        from prody.utilities import openFile
        with openFile(pdblist[0]) as inp:
            for item in inp.read().strip().split():
                for pdb in item.split(','):
                    if len(pdb) == 4 and pdb.isalnum(): 
                        pdblist.append(pdb)

    prody.fetchPDB(pdblist, opt.folder, compressed=opt.gzip, copy=True)

def addCommand(commands):
    
    subparser = commands.add_parser('fetch', 
        help='fetch a PDB file')
        
    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Download PDB file(s) by specifying identifiers:
        
$ prody fetch 1mkp 1p38"""
    )

    subparser.add_argument('-d', '--dir', dest='folder', type=str,
        default='.', metavar='PATH', 
        help=('target directory for saving PDB file(s)'))
    subparser.add_argument('-z', '--gzip', dest='gzip', action='store_true', 
        default=False, help='write compressed PDB file(S)')

    subparser.add_argument('pdb', nargs='+', 
        help='PDB identifier(s) or a file that contains them')

    subparser.set_defaults(func=prody_fetch)
    subparser.set_defaults(subparser=subparser)
