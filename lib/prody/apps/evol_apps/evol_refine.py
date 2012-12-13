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

"""Refine MSA application."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..apptools import DevelApp

__all__ = ['evol_refine']

APP = DevelApp('refine', 'refine an MSA by removing gapped rows/colums')

APP.setExample(
"""Refines MSA by removing gapped columns (residue positions)
and rows (sequences).  

In the following, columns that are gaps in sequence with label 
GTHB2_ONCKE will be removed from MSA:

  $ evol refine piwi.slx -l GTHB2_ONCKE""", [])


APP.addArgument('msa', 
    help='MSA filename to be refined')

APP.addGroup('refine', 'refinement options')
APP.addArgument('-l', '--label',
    dest='label',
    help='sequence label, UniProt ID code, or PDB and chain identifier',
    type=str,
    metavar='STR',
    group='refine')

APP.addArgument('-s', '--seqid',
    dest='seqid',
    help='identity threshold for selecting unique sequences',
    default=None,
    type=float,
    metavar='FLOAT',
    group='refine')

APP.addArgument('-c', '--colocc',
    dest='colocc',
    help='column (residue position) occupancy',
    default=None,
    type=float,
    metavar='FLOAT',
    group='refine')

APP.addArgument('-r', '--rowocc',
    dest='rowocc',
    help='row (sequence) occupancy',
    default=None,
    type=float,
    metavar='FLOAT',
    group='refine')

APP.addGroup('output', 'output options')
APP.addArgument('-o', '--outname',
    dest='outname',
    help='output filename, default is msa filename with _refined suffix',
    type=str,
    metavar='STR',
    group='output')
    
APP.addArgument('-f', '--format', 
    dest='format', 
    type=str,
    metavar='STR', 
    help='output MSA file format, default is same as input', 
    group='output')
    
APP.addArgument('-z', '--compressed', 
    dest='compressed', 
    action='store_true', 
    help='gzip refined MSA output',
    group='output')

def evol_refine(msa, **kwargs):
    
    import prody
    from prody import parseMSA, refineMSA, writeMSA, LOGGER
    from os.path import splitext

    outname = kwargs.get('outname')
    if outname is None: 
        outname, ext = splitext(msa)
        if ext.lower() == '.gz': 
            outname, _ = splitext(msa)
        outname += '_refined' + ext 
    
    writeMSA(outname, refineMSA(parseMSA(msa), **kwargs), **kwargs)
    LOGGER.info('Refined MSA is written in file: ' + outname)     
    

APP.setFunction(evol_refine)
