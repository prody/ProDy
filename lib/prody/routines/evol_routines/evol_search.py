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

""" Search Pfam with UniProt ID or protein sequence as input for corresponding
    accession code"""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..actions import *

__all__ = ['evol_search']

DEFAULTS = {}
HELPTEXT = {} 
for key, txt, val in [
    ('search_b', 'search Pfam-B families', False),
    ('skip_a', 'do not search Pfam-A families', False),
    ('ga', 'use gathering threshold', True),   
    ('evalue', 'E-value cutoff: must be < 10.0', None),
    ('timeout', 'timeout for blocking connection attempt in seconds', 30),
    ('outname', 'name for output file', None),
    ('folder', 'path where the output should be saved', '.'),
    ('delimiter', 'delimiter for output file', '\t')]:
    
    DEFAULTS[key] = val
    HELPTEXT[key] = txt

def evol_search(pfam_query, **kwargs):
    """Search Pfam with **query**.
    
    """
    import prody
    import sys
    from os.path import join
    
    search_b = kwargs.pop('search_b', DEFAULTS['search_b'])
    skip_a = kwargs.pop('skip_a', DEFAULTS['skip_a'])
    
    pfam_results =  prody.searchPfam(pfam_query, search_b=search_b,
                                     skip_a=skip_a, **kwargs)
    
    outname = kwargs.get('outname', None)
    delimiter = kwargs.get('delimiter', DEFAULTS['delimiter'])
    folder = kwargs.get('folder', DEFAULTS['folder'])
    if outname:
        filepath = join(prody.utilities.makePath(folder), outname)
        out = open(filepath, 'wb')
    else:
        out = sys.stdout
    title = ('accession' + delimiter + 'identifier' + delimiter + 'type' +
             delimiter + 'smallest-evalue' + '\n')
    out.write(title)
    for key in pfam_results:
        val = pfam_results[key]
        output = (val.get('accession', None) + delimiter + val.get('id', None) + delimiter +
                  val.get('type', None) + delimiter)
        locations = val.get('locations', None)
        evalue = ''
        if locations:
            for i, location in enumerate(locations):
                temp = location.get('evalue', None)
                if temp:
                    if i==0:
                        evalue = float(temp)
                    else:
                        if float(temp) < evalue:
                            evalue = float(temp)        
        output =  output + str(evalue) + '\n'               
        out.write(output)
        
    out.close()
        
        
def addCommand(commands):

    subparser = commands.add_parser('search', 
        help='search Pfam with given query')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """This program searches Pfam database with a UniProt ID or protein \
sequence or filename containing the sequence. Some specific search options \
are included for sequence search. Minimum length of query sequence should \
be 12 and should not contain gaps. If outname is specified it will output \
the results obtained in a file or the output will be directed to standard output.

Search Pfam with UniProt ID:
    
    $ prody search P08581
    
Search Pfam with a sequence with search options:

    $ prody search PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGSSEPCALCSLHSIGKIGGAQNRSYSKLLC\
GLLAERLRISPDRVYINYYDMNAANVGWNNSTFA' --evalue=2 --searchBs""")

    group = subparser.add_argument_group('search options for sequence input')
    
    group.add_argument('-b', '--searchBs', dest='search_b', 
        action='store_true', 
        default=DEFAULTS['search_b'], help=HELPTEXT['search_b'])
    group.add_argument('-s', '--skipAs', dest='skip_a', 
        action='store_true', 
        default=DEFAULTS['skip_a'], help=HELPTEXT['skip_a'])
    group.add_argument('-g', '--ga', dest='ga', action='store_true', 
        default=DEFAULTS['ga'], help=HELPTEXT['ga'])
    group.add_argument('-e', '--evalue', dest='evalue', type=float, 
        default=DEFAULTS['evalue'], metavar='FLOAT', 
        help=HELPTEXT['evalue'] + ' (default: %(default)s)')
    group.add_argument('-t', '--timeout', dest='timeout', type=int, 
        default=DEFAULTS['timeout'], metavar='INT', 
        help=HELPTEXT['timeout'] + ' (default: %(default)s)')
    
    group = subparser.add_argument_group('output options')
    
    group.add_argument('-o', '--outdir', dest='folder', type=str, 
        default=DEFAULTS['folder'], metavar='PATH', 
        help=HELPTEXT['folder'] + ' (default: %(default)s)')
    group.add_argument('-p', '--outname', dest='outname', type=str, 
        default=DEFAULTS['outname'], metavar='STR', 
        help=HELPTEXT['outname'] + ' (default: %(default)s)')
    group.add_argument('-d', '--delimiter', dest='delimiter', type=str, 
        default=DEFAULTS['delimiter'], metavar='STR', 
        help=HELPTEXT['delimiter'] + ' (default: %(default)s)')
        
    subparser.add_argument('query', help=('UniProt ID or protein sequence or'
                                          ' filename containing sequence '
                                          'with no gaps and >11 characters'))
            
    subparser.set_defaults(func=lambda ns: evol_search(ns.__dict__.pop('query'), 
                                                        **ns.__dict__))
    subparser.set_defaults(subparser=subparser)

