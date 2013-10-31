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

from ..apptools import DevelApp

__all__ = ['evol_search']

APP = DevelApp('search', 'search Pfam with given query')

APP.addGroup('sequence', 'sequence search options')
APP.addGroup('output', 'output options')

APP.addArgument('query',
    help='protein UniProt ID or sequence, a PDB identifier, or a sequence '
         'file, where sequence have no gaps and 12 or more characters')
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
    default=60,
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
"""Search Pfam database with a UniProt ID or protein sequence.
Sequence input can be parsed from a FASTA or Selex.  Minimum length
of query sequence should be 16 characters and should not contain gaps.
If output name is specified, results will be written in a file.
Otherwise, or the output will be directed to standard output.

Search Pfam with PDB and chain identifier and output results to screen:

  $ evol search 1mkpA

Search Pfam with UniProt ID and write output into a file:

  $ evol search P08581 --outname families.txt

Search Pfam with a sequence and with some search options:

  $ evol search PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGSSEPCALCS\
LHSIGKIGGAQNRSYSKLLCGLLAERLRISPDRVYINYYDMNAANVGWNNSTFA --evalue 2 \
--searchBs""", [0, 1])

def evol_search(query, **kwargs):

    import prody
    from os.path import join, split

    pfam_results =  prody.searchPfam(query, **kwargs)
    if pfam_results is None:
        return
    outname = kwargs.get('outname', None)
    delimiter = kwargs.get('delimiter', '\t')
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
        prody.LOGGER.info('Search results written in {0}.'.format(filepath))
        out.close()

APP.setFunction(evol_search)
