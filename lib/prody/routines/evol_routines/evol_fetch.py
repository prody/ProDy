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

""" Download MSA for given Pfam ID or accession code and print
    the path to the downloaded file"""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..actions import *

__all__ = ['evol_fetch']

DEFAULTS = {}
HELPTEXT = {} 
for key, txt, val in [
    ('alignment', 'alignment type: one of full, seed, ncbi or metagenomics',
     'full'),
    ('format', 'Pfam supported MSA file format: selex, fasta or stockholm',
     'selex'),
    ('order', 'ordering of sequences: tree or alphabetical', 'tree'),   
    ('inserts', 'letter case for inserts: upper or lower', 'upper'),
    ('gaps', 'gap character: dashed, dots or mixed', 'dashes'),
    ('folder', 'output directory', '.'),
    ('compressed', 'gzip downloaded MSA file', False),
    ('outname', 'out filename', None)]:
    
    DEFAULTS[key] = val
    HELPTEXT[key] = txt

def evol_fetch(acc, **kwargs):
    """Download Pfam MSA for **acc**.
    
    """
    import prody
    
    alignment = kwargs.pop('alignment', DEFAULTS['alignment'])
    compressed = kwargs.pop('compressed', DEFAULTS['compressed'])
    
    prody.fetchPfamMSA(acc, alignment=alignment,
                        compressed=compressed, **kwargs)

                
def addCommand(commands):

    subparser = commands.add_parser('fetch', 
        help='fetch MSA from Pfam')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Given a Pfam ID or accession code, this program fetches the MSA of \
that family. Supported alignment options are full, seed, ncbi or metagenomics \
and alignment formats are selex, stockholm or fasta. The output MSA is \
downloaded and saved in the specified or default '.' directory.

Fetch PFAM ID Cys_knot:
    
    $ evol fetch Cys_knot
    
Fetch PFAM accession with specific parameters:

    $ evol fetch PF00007 --compressed --format fasta --outname mymsa""",
    )

    group = subparser.add_argument_group('download options')
    
    group.add_argument('-a', '--alignment', dest='alignment', type=str, 
        default=DEFAULTS['alignment'], metavar='STR', 
        help=HELPTEXT['alignment'] + ' (default: %(default)s)')
    group.add_argument('-f', '--format', dest='format', type=str, 
        default=DEFAULTS['format'], metavar='STR', 
        help=HELPTEXT['format'] + ' (default: %(default)s)')
    group.add_argument('-o', '--order', dest='order', type=str, 
        default=DEFAULTS['order'], metavar='STR', 
        help=HELPTEXT['order'] + ' (default: %(default)s)')
    group.add_argument('-i', '--inserts', dest='inserts', type=str, 
        default=DEFAULTS['inserts'], metavar='STR', 
        help=HELPTEXT['inserts'] + ' (default: %(default)s)')
    group.add_argument('-g', '--gaps', dest='gaps', type=str, 
        default=DEFAULTS['gaps'], metavar='STR', 
        help=HELPTEXT['gaps'] + ' (default: %(default)s)')
    
    group = subparser.add_argument_group('output options')
    
    group.add_argument('-d', '--outdir', dest='folder', type=str, 
        default=DEFAULTS['folder'], metavar='PATH', 
        help=HELPTEXT['folder'] + ' (default: %(default)s)')
    group.add_argument('-p', '--outname', dest='outname', type=str, 
        default=DEFAULTS['outname'], metavar='STR', 
        help=HELPTEXT['outname'] + ' (default: %(default)s)')
    group.add_argument('-z', '--compressed', dest='compressed', 
        action='store_true', 
        default=DEFAULTS['compressed'], help=HELPTEXT['compressed'])
    
    subparser.add_argument('acc', help='Pfam ID or accession code')
            
    subparser.set_defaults(func=lambda ns: evol_fetch(ns.__dict__.pop('acc'), 
                                                        **ns.__dict__))
    subparser.set_defaults(subparser=subparser)

