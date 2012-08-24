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

"""Blast Protein Data Bank for structures matching a user given sequence."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

from actions import *

def readFirstSequenceFasta(filename):
    """Return first sequence from a file."""
    
    fasta = open(filename)
    seq = []
    title = '' 
    first = True
    for line in fasta:
        if line[0] == '>': 
            if first:
                title = line[1:].strip()
                first = False
            else:    
                break
        else:
            seq.append( line.strip() )
    fasta.close()
    return title, ''.join(seq)

def prody_blast(opt):
    """Blast search PDB based on command line arguments."""
    
    import prody
    LOGGER = prody.LOGGER
    seq = opt.seq
    title = None
    if os.path.isfile(seq):
        title, seq = readFirstSequenceFasta(seq)
        LOGGER.info("First sequence ({0:s}) is parsed from {1:s}."
                    .format(title, repr(seq)))
    if not seq.isalpha() or not seq.isupper():
        opt.subparser.error("{0:s} is not a valid sequence or a file"
                            .format(repr(seq)))
        
    folder, identity, coverage = opt.folder, opt.identity, opt.coverage
    if not 0 < identity < 100: 
        opt.subparser.error('identity must be between 0 and 100')
    if not 0 < coverage < 100:
        opt.subparser.error('overlap must be between 0 and 100')
    
    blast_results = prody.blastPDB(seq)
    hits = blast_results.getHits(percent_identity=identity, 
                                 percent_coverage=coverage)
    
    #sort hits by decreasing percent identity
    hits2 = []
    for pdb in hits:
        hits2.append( (-hits[pdb]['percent_identity'], pdb) )
    hits2.sort()
    
    for identity, pdb in hits2:
        chain = hits[pdb]['chain_id']
        percent_identity = hits[pdb]['percent_identity']
        title = hits[pdb]['title']
        print(pdb + ' ' + chain + ' ' + ('%5.1f%%' % (percent_identity)) + 
              ' ' + title)
    
    # download hits if --folder is given
    if opt.folder:
        LOGGER.info('Downloading hits to ' + opt.folder)
        pdblist = [ pdb for identity, pdb in hits2 ]
        pdblist2 = prody.fetchPDB(pdblist, opt.folder, 
                                  compressed=opt.gzip, copy=True)
def addCommand(commands):
    
    subparser = commands.add_parser('blast', 
        help='blast search Protein Data Bank')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Blast search PDB for the first sequence in a fasta file:
    
  $ prody blast seq.fasta -i 70

Blast search PDB for the sequence argument:

  $ prody blast MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ\
KESTLHLVLRLRGG

Blast search PDB for avidin structures, download files, and align all files \
onto the 2avi structure:
    
  $ prody blast -d . ARKCSLTGKWTNDLGSNMTIGAVNSRGEFTGTYITAVTATSNEIKESPLHGTQNTIN\
KRTQPTFGFTVNWKFSESTTVFT
  $ prody align 2avi.pdb *pdb """)

    subparser.add_argument('-i', '--identity', dest='identity', type=float, 
        default=90.0, metavar='FLOAT', 
        help='percent sequence identity (default: %(default)s)')
    subparser.add_argument('-o', '--overlap', dest='coverage', type=float, 
        default=90.0, metavar='FLOAT', 
        help='percent sequence overlap (default: %(default)s)')
    subparser.add_argument('-d', '--dir', dest='folder', type=str,
        default=None, metavar='PATH', 
        help=('download uncompressed PDB files to given path'))

    subparser.add_argument('-z', '--gzip', dest='gzip', action='store_true', 
                     default=False, help='write compressed PDB file')

    subparser.add_argument('seq', type=str,  
        help=('sequence or file in fasta format'))

    subparser.set_defaults(func=prody_blast)
    subparser.set_defaults(subparser=subparser)
