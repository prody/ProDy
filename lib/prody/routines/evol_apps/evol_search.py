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

"""Pfam search application."""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from collections import defaultdict

from ..actions import *
from ..apptools import *

__all__ = ['evol_search']

APP = DevelApp('search', 'search Pfam with given query')

APP.addGroup('sequence', 'sequence search options')
APP.addGroup('output', 'output options')

APP.addArgument('query', 
    help='protein UniProt ID or sequence, or a sequence file, where '
         'sequence have no gaps and more than 11 characters')
APP.addArgument('-b', '--searchBs',
    dest='search_b',
    help='search Pfam-B families',
    default=False,
    action='store_true',
    group='sequence')

APP.addArgument('-s', '--skipAs',
    dest='skip_a',
    help='do not search Pfam-A families',
    default=False,
    action='store_true',
    group='sequence')
APP.addArgument('-g', '--ga',
    dest='ga',
    help='use gathering threshold',
    default=False,
    action='store_true',
    group='sequence')
APP.addArgument('-e', '--evalue',
    dest='evalue',
    help='e-value cutoff, must be less than 10.0',
    default=None,
    type=float,
    metavar='FLOAT',
    group='sequence')
APP.addArgument('-t', '--timeout',
    dest='timeout',
    help='timeout in seconds for blocking connection attempt',
    default=30,
    type=int,
    metavar='INT',
    group='sequence')
APP.addArgument('-o', '--outname',
    dest='outname',
    help='name for output file, default is standard output',
    type=str,
    metavar='STR',
    group='output')
APP.addArgument('-d', '--delimiter', 
    dest='delimiter', 
    type=str, 
    default='\t', 
    metavar='STR', 
    help='delimiter for output data columns',
    group='output')

APP.setExample(
"""This application searches Pfam database with a UniProt ID or protein \
sequence or filename containing the sequence.  Some specific search options \
are included for sequence search.  Minimum length of query sequence should \
be 12 and should not contain gaps.  If outname is specified it will output \
the results obtained in a file or the output will be directed to standard \
output.

Search Pfam with UniProt ID and write output into a file:
    
    $ evol search P08581 --outname families.txt
    
Search Pfam with a sequence with search options:

    $ evol search PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGSSEPCALCS\
LHSIGKIGGAQNRSYSKLLCGLLAERLRISPDRVYINYYDMNAANVGWNNSTFA --evalue 2 \
--searchBs
    """, [0, 1])

def evol_search(query, **kwargs):
    
    import prody
    from os.path import join, split
    
    pfam_results =  prody.searchPfam(query, **kwargs)
    
    outname = kwargs.get('outname', None)
    delimiter = kwargs.get('delimiter', DEFAULTS['delimiter'])
    if outname:
        folder, outname = split(outname)
        filepath = join(prody.utilities.makePath(folder), outname)
        out = open(filepath, 'wb')
    else:
        from sys import stdout as out
    title = delimiter.join(['acc', 'id', 'type', 'e-value']) + '\n'
    out.write(title)
    for key in pfam_results:
        val = pfam_results[key]
        evalue = ''
        for i, location in enumerate(val.get('locations', [])):
            temp = location.get('evalue', None)
            if temp:
                if i==0:
                    evalue = float(temp)
                else:
                    if float(temp) < evalue:
                        evalue = float(temp)        
        output = delimiter.join([val.get('accession', '    '),
                                 val.get('id', '    '),
                                 val.get('type', '    '),
                                 str(evalue)]) + '\n'               
        out.write(output)
    if outname:
        prody.LOGGER.info('Search results written in {0:s}.'.format(filepath))    
        out.close()
        
APP.setFunction(evol_search)
