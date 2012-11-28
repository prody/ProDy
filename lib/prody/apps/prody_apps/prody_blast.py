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

from ..apptools import *

__all__ = ['prody_blast']

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

def prody_blast(sequence, **kwargs):
    """Blast search PDB and download hits.
    
    :arg sequence: sequence or file in fasta format
    
    :arg identity: percent sequence identity for blast search, default is 90.0
    :type identity: float
    
    :arg overlap: percent sequence overlap between sequences, default is 90.0
    :type overlap: float
    
    :arg outdir: download uncompressed PDB files to given directory
    :type outdir: str
    
    :arg gzip: write compressed PDB file
    
    *Blast Parameters*
    
    :arg filename: a *filename* to save the results in XML format 
    :type filename: str
    
    :arg hitlist_size: search parameters, default is 250
    :type hitlist_size: int
    
    :arg expect: search parameters, default is 1e-10 
    :type expect: float
    
    :arg sleep: how long to wait to reconnect for results, default is 2
                sleep time is doubled when results are not ready.
    :type sleep: int        
    
    :arg timeout: when to give up waiting for results. default is 30
    :type timeout: int"""  
    
    import prody
    LOGGER = prody.LOGGER
    title = None
    if os.path.isfile(sequence):
        title, sequence = readFirstSequenceFasta(sequence)
        LOGGER.info("First sequence ({0:s}) is parsed from {1:s}."
                    .format(title, repr(sequence)))
    if not sequence.isalpha() or not sequence.isupper():
        raise ValueError("{0:s} is not a valid sequence or a file"
                        .format(repr(sequence)))
        
    outdir = kwargs.get('outdir')
    identity, overlap = kwargs.get('identity', 90), kwargs.get('overlap', 90)
    if not 0 < identity < 100: 
        raise ValueError('identity must be between 0 and 100')
    if not 0 < overlap < 100:
        raise ValueError('overlap must be between 0 and 100')
    


    filename = kwargs.get('filename', None)
    hitlist_size = kwargs.get('hitlist_size', 250)
    expect = kwargs.get('expect', 1e-10)
    sleep, timeout = kwargs.get('sleep', 2), kwargs.get('timeout', 30)
    
    blast_results = prody.blastPDB(sequence,filename=filename,
                                   hitlist_size=hitlist_size, expect=expect,
                                   sleep=sleep, timeout=timeout)

    if blast_results is None:
        raise IOError('blast search timed out, please try again')
        
    hits = blast_results.getHits(percent_identity=identity, 
                                 percent_overlap=overlap)
       
    #sort hits by decreasing percent identity
    hits2 = []
    for pdb in hits:
        hits2.append( (-hits[pdb]['percent_identity'], pdb) )
    hits2.sort()
    
    stdout = kwargs.get('stdout', False)
    
    if not stdout:
        finalHits = []
    else:
        from sys import stdout
        
    for identity, pdb in hits2:
        chain = hits[pdb]['chain_id']
        percent_identity = hits[pdb]['percent_identity']
        title = hits[pdb]['title']
        if stdout:
            stdout.write(pdb + ' ' + chain + ' ' + 
                         ('%5.1f%%' % (percent_identity)) + ' ' + title)
        else:
            finalHits.append((pdb, chain, ('%5.1f%%' % (percent_identity)),
                               title))
            
    
    # download hits if --output-dir is given
    if outdir:
        LOGGER.info('Downloading hits to ' + outdir)
        pdblist = [ pdb for identity, pdb in hits2 ]
        pdblist2 = prody.fetchPDB(pdblist, outdir, 
                                  compressed=kwargs.get('gzip'), copy=True)
        
    if not stdout:
        return finalHits
    
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
    subparser.add_argument('-o', '--overlap', dest='overlap', type=float, 
        default=90.0, metavar='FLOAT', 
        help='percent sequence overlap (default: %(default)s)')
    subparser.add_argument('-d', '--output-dir', dest='outdir', type=str,
        default=None, metavar='PATH',  
        help='download uncompressed PDB files to given directory')

    subparser.add_argument('-z', '--gzip', dest='gzip', action='store_true', 
                     default=False, help='write compressed PDB file')

    subparser.add_argument('sequence', type=str,  
        help='sequence or file in fasta format')
    
    group = subparser.add_argument_group('Blast Parameters')
    
    group.add_argument('-f', '--filename', dest='filename', type=str,
        default=None, metavar='STR',
        help='a filename to save the results in XML format')
        
    group.add_argument('-e', '--expect', dest='expect', type=float,
        default=1e-10, metavar='FLOAT', help='blast search parameter')
        
    group.add_argument('-l', '--hit-list-size', dest='hitlist_size', type=int,
        default=250, metavar='INT', help='blast search parameter')
    
    
    group.add_argument('-s', '--sleep-time', dest='sleep', type=int,
        default=2, metavar='INT',
        help='how long to wait to reconnect for results '
              '(sleep time is doubled when results are not ready)')
    
    group.add_argument('-t', '--timeout', dest='timeout', type=int,
        default=30, metavar='INT',
        help='when to give up waiting for results')
    
    subparser.set_defaults(func=lambda ns: prody_blast(
                        ns.__dict__.pop('sequence'), **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
    subparser.set_defaults(stdout=True)
    
    
