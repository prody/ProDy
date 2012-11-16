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

"""Pfam MSA download application."""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..actions import *
from ..apptools import *

__all__ = ['evol_fetch']

APP = DevelApp('fetch', help='fetch MSA files from Pfam')


APP.setExample(
"""Given a Pfam ID or accession code, this program fetches the MSA of \
that family. Supported alignment options are full, seed, ncbi or metagenomics \
and alignment formats are selex, stockholm or fasta. The output MSA is \
downloaded and saved in the specified or default '.' directory.

Fetch PFAM ID Cys_knot:
    
    $ evol fetch Cys_knot
    
Fetch PFAM accession with specific parameters:

    $ evol fetch PF00007 --compressed --format fasta --outname mymsa""",
    [0, 1])


APP.addArgument('acc', help='Pfam ID or accession code')

APP.addGroup('download', 'download options')
APP.addArgument('-a', '--alignment', 
    dest='alignment', 
    type=str, 
    default='full', 
    metavar='STR', 
    help='alignment type, one of full, seed, ncbi or metagenomics',
    group='download')
APP.addArgument('-f', '--format', 
    dest='format', 
    type=str, 
    default='selex', 
    metavar='STR', 
    help='Pfam supported MSA format, one of selex, fasta or stockholm',
    group='download')
APP.addArgument('-o', '--order', 
    dest='order', 
    type=str, 
    default='tree', 
    metavar='STR', 
    help='ordering of sequences, tree or alphabetical',
    group='download')
APP.addArgument('-i', '--inserts', 
    dest='inserts', 
    type=str, 
    default='upper', 
    metavar='STR', 
    help='letter case for inserts, upper or lower',
    group='download')
APP.addArgument('-g', '--gaps', 
    dest='gaps', 
    type=str, 
    default='dashes', 
    metavar='STR', 
    help='gap character, one of dashes, dots or mixed',
    group='download')

APP.addGroup('output', 'output options')
APP.addArgument('-d', '--outdir', 
    dest='folder', 
    type=str, 
    default='.', 
    metavar='PATH', 
    help='output directory',
    group='output')
APP.addArgument('-p', '--outname', 
    dest='outname', 
    type=str, 
    default=None, 
    metavar='STR', 
    help='out filename',
    group='output')
APP.addArgument('-z', '--compressed', 
    dest='compressed', 
    action='store_true', 
    help='gzip downloaded MSA file',
    group='output')
    
def evol_fetch(acc, **kwargs):
   
    import prody
    
    alignment = kwargs.pop('alignment', DEFAULTS['alignment'])
    compressed = kwargs.pop('compressed', DEFAULTS['compressed'])
    
    prody.fetchPfamMSA(acc, alignment=alignment,
                        compressed=compressed, **kwargs)

        
APP.setFunction(evol_fetch)
