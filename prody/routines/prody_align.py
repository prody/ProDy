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

"""Align models in a PDB file or multiple structures in separate PDB files."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from actions import *

def prody_align(opt):
    """Align models in a PDB file or a PDB file onto others."""
            
    import prody
    LOGGER = prody.LOGGER

    args = opt.pdb
    if len(args) == 1:
        pdb = args[0]
        LOGGER.info('Aligning multiple models in: ' + pdb)
        selstr, prefix, model = opt.select, opt.prefix, opt.model
        pdb = prody.parsePDB(pdb)
        pdbselect = pdb.select(selstr)
        if pdbselect is None:
            opt.subparser.error('Selection {0:s} do not match any atoms.'
                               .format(repr(selstr)))
        LOGGER.info('{0:d} atoms will be used for alignment.'
                    .format(len(pdbselect)))
        pdbselect.setACSIndex(model-1)
        prody.printRMSD(pdbselect, msg='Before alignment ')
        prody.alignCoordsets(pdbselect)
        prody.printRMSD(pdbselect, msg='After alignment  ')
        if prefix == '':
            prefix = pdb.getTitle() + '_aligned'
        outfn = prefix + '.pdb'
        LOGGER.info('Writing file: ' + outfn)
        prody.writePDB(outfn, pdb)
    else:
        reffn = args.pop(0)
        seqid=opt.seqid
        overlap=opt.overlap
        LOGGER.info('Aligning structures onto: ' + reffn)
        ref = prody.parsePDB(reffn)
        for arg in args:
            if arg == reffn:
                continue
            if '_aligned.pdb' in arg:
                continue
            pdb = prody.parsePDB(arg)
            result = prody.matchAlign(pdb, ref, seqid=seqid, overlap=overlap, 
                                      tarsel=opt.select, allcsets=True,
                                      cslabel='Model', csincr=1) 
            if result:
                outfn = pdb.getTitle() + '_aligned.pdb'
                LOGGER.info('Writing file: ' + outfn)
                prody.writePDB(outfn, pdb)
            else:
                LOGGER.warning('Failed to align ' + arg)
                
def addCommand(commands):

    subparser = commands.add_parser('align', 
        help='align models or structures')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Align models in PDB structure or multiple PDB structures and save \
aligned coordinate sets.  When multiple structures are aligned, ProDy will \
match chains based on sequence alignment and use best match for aligning the \
structures.

Fetch PDB structure 2k39 and align models (reference model is the first model):
    
    $ prody align 2k39
    
Fetch PDB structure 2k39 and align models using backbone of residues with \
number less than 71:

    $ prody align 2k39 --select "backbone and resnum < 71" 
    
Align 1r39 and 1zz2 onto 1p38 using residues with number less than 300:

    $ prody align --select "resnum < 300" 1p38 1r39 1zz2
    
Align all models of 2k39 onto 1aar using residues 1 to 70 (inclusive):

    $ prody align --select "resnum 1 to 70" 1aar 2k39 
    """
    )
        
    subparser.add_argument('-p', '--prefix', dest='prefix', type=str, 
        default='', metavar='STR', 
        help=('prefix for output files, default is PDB_aligned'))
    subparser.add_argument('-s', '--select', dest='select', type=str, 
        default='calpha', metavar='SELSTR',
        help='selection string (default: "%(default)s")')
    subparser.add_argument('-m', '--model', dest='model', type=int, 
        default=1, metavar='INT',
        help=('for NMR files, reference model index (default: %(default)s)'))
    subparser.add_argument('-i', '--seqid', dest='seqid', type=int, 
        default=90, metavar='INT',
        help=('percent sequence identity (default: %(default)s)'))
    subparser.add_argument('-o', '--overlap', dest='overlap', type=int, 
        default=90, metavar='INT',
        help=('percent sequence overlap (default: %(default)s)'))

    subparser.add_argument('pdb', nargs='+', 
        help='PDB identifier(s) or filename(s)')
            
    subparser.set_defaults(func=prody_align)
    subparser.set_defaults(subparser=subparser)
