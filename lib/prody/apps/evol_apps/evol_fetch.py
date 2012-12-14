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

from ..apptools import DevelApp

__all__ = ['evol_fetch']

APP = DevelApp('fetch', help='fetch MSA files from Pfam')


APP.setExample(
"""Fetch MSA files from Pfam database for given Pfam ID code or 
accession.

Fetch PFAM ID Cys_knot:
    
  $ evol fetch Cys_knot
    
Fetch PFAM accession with specific parameters:

  $ evol fetch PF00007 --compressed --format fasta --outname mymsa""", [])


APP.addArgument('acc', 
    help='Pfam accession or ID',
    type=str)

APP.addGroup('download', 'download options')
APP.addArgument('-a', '--alignment', 
    dest='alignment', 
    type=str, 
    default='full', 
    metavar='STR', 
    help='alignment type',
    choices='full seed ncbi metagenomics'.split(),
    group='download')
APP.addArgument('-f', '--format', 
    dest='format', 
    type=str, 
    default='selex', 
    metavar='STR', 
    help='Pfam supported MSA format',
    choices='selex fasta stockholm'.split(),
    group='download')
APP.addArgument('-o', '--order', 
    dest='order', 
    type=str, 
    default='tree', 
    metavar='STR', 
    help='ordering of sequences',
    choices='tree alphabetical'.split(),
    group='download')
APP.addArgument('-i', '--inserts', 
    dest='inserts', 
    type=str, 
    default='upper', 
    metavar='STR', 
    help='letter case for inserts',
    choices='upper lower'.split(),
    group='download')
APP.addArgument('-g', '--gaps', 
    dest='gaps', 
    type=str, 
    default='dashes', 
    metavar='STR', 
    help='gap character',
    choices='dashes dots mixed'.split(),
    group='download')
APP.addArgument('-t', '--timeout', 
    dest='timeout', 
    type=int, 
    default=5, 
    metavar='INT', 
    help='timeout for blocking connection attempts',
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
    help='output filename, default is accession and alignment type',
    group='output')
APP.addArgument('-z', '--compressed', 
    dest='compressed', 
    action='store_true', 
    help='gzip downloaded MSA file',
    group='output')
    
def evol_fetch(acc, **kwargs):
   
    import prody
    prody.fetchPfamMSA(acc, **kwargs)

        
APP.setFunction(evol_fetch)
