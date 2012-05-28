# -*- coding: utf-8 -*-
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

"""This module defined a routine for contact identification."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from actions import *

__all__ = ['prody_contacts']

def prody_contacts(**kwargs):
    """Identify contacts of a target structure with one or more other 
    structures, such as ligands.
    
    :arg target: target PDB identifier or filename
    :type target: str
    
    :arg target: ligand PDB identifier(s) or filename(s)
    :type target: list

    :arg select: atom selection string for target structure (first PDB file)
    :type select: str
    
    :arg radius: contact radius, default is 4.0 Å
    :type radius: float
    
    :arg extend: output same residue, chain, or segment as contacting atoms
    :type extend: str
    
    :arg prefix: prefix for output file, default is PDB filename itself
    :type prefix: str
    
    :arg suffix: output filename suffix, default is ``'_suffix'``
    :type suffix: str"""
            
    import prody
    LOGGER = prody.LOGGER

    target = prody.parsePDB(kwargs.get('target'))
    title = kwargs.get('prefix') or target.getTitle()
    selstr = kwargs.get('select')
    if selstr:
        target = target.select(selstr)
    contacts = prody.Contacts(target)
    suffix = kwargs.get('suffix', '_contacts')
    extend = kwargs.get('extend')
    radius = float(kwargs.get('radius', 4.0))
    ligands = kwargs.get('ligand')
    if len(ligands) > 1:
        outfn = lambda fn: title + suffix + '_' + fn + '.pdb'
    else:
        outfn = lambda fn: title + suffix + '.pdb'
    for pdb in ligands:
        ligand = prody.parsePDB(pdb)
        sel = contacts(radius, ligand)
        if sel:
            if extend:
                sel = target.select('same ' + extend + ' as sel', sel=sel)
            pdbfn = outfn(ligand.getTitle())
            LOGGER.info('Writing file: ' + pdbfn)
            prody.writePDB(pdbfn, sel)
   
                
def addCommand(commands):

    subparser = commands.add_parser('contacts', 
        help='identify contacts between a target and ligand(s)')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Identify contacts of a target structure with one or more ligands.

Fetch PDB structure 1zz2, save PDB files for individual ligands, and identify 
contacting residues of the target protein:
    
    $ prody select -o B11 1zz2 "resname B11"
    $ prody select -o BOG 1zz2 "resname BOG"
    $ prody contacts -r 4.0 -t residue -s protein 1zz2 B11.pdb BOG.pdb
    """
    )
        
    subparser.add_argument('-s', '--select', dest='select', type=str, 
        metavar='SELSTR', help='selection string for target')
    
    subparser.add_argument('-r', '--radius', dest='radius', type=float, 
        default=4.0, metavar='FLOAT',
        help='contact radius (default: %(default)s)')
    
    subparser.add_argument('-t', '--extend', dest='extend', type=str, 
        default='_aligned', metavar='STR', 
        choices=set(['residue', 'chain', 'segment']),
        help=('output same residue residue, chain, or segment as the '
              'contacts'))

    subparser.add_argument('-p', '--prefix', dest='prefix', type=str, 
        default='', metavar='STR', 
        help=('output filename prefix (default: target filename)'))
        
    subparser.add_argument('-x', '--suffix', dest='suffix', type=str, 
        default='_contacts', metavar='STR', 
        help=('output filename suffix (default: %(default)s)'))


    subparser.add_argument('target',
        help='target PDB identifier or filename')

    subparser.add_argument('ligand', nargs='+',
        help='ligand PDB identifier(s) or filename(s)')
            

    subparser.set_defaults(func=lambda opt: prody_contacts(**opt.__dict__))
    subparser.set_defaults(subparser=subparser)
