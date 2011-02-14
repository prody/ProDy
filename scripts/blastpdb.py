#!/usr/bin/python
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan
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

__author__ = 'Lidio Meireles'
__copyright__ = 'Copyright (C) 2010 Lidio Meireles, Ahmet Bakan'

import sys
from optparse import OptionParser
from prody import *

def read_fasta_sequence (seqfn):
    # returns first sequence on file
    f = open(seqfn)
    seqlines = []
    n = 0
    for line in f.xreadlines():
        n += 1
        if n == 1 and line.startswith('>'): continue
        if n  > 1 and line.startswith('>'): break
        seqlines.append( line.strip() )
    f.close()
    
    return ''.join(seqlines)       
        
def main():
    usage = '%prog <sequence.fasta> <identity>\n\nBlast sequence and return pdb hits with given identity threshold [0-100].\n'
    parser = OptionParser(usage=usage)
    parser.add_option('','--folder',dest="folder",type="string",default="",metavar="PATH",help="if given, download pdb files to the folder")
    parser.add_option('','--hitlist_size',dest='hitlist_size',type='int',default=250,metavar='INT',help='max hit list size (%default)')
    parser.add_option('-e','--expect',dest='expect',type='float',default=1E-10,metavar='FLOAT',help='Blast e-score threshold (1e-10)')

    opt, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit(-1)
        
    seqfn, identity = args[0:2]
    try:
        identity = float(identity)
        if identity < 0 or identity > 100: raise Exception()
    except Exception:
        parser.print_help()
        print '\nError: <identity> must be float in range [0-100]'
        sys.exit(-1)

    try:
        seq = read_fasta_sequence(seqfn)
    except Exception:
        print '\nError: parsing fasta file',seqfn
        sys.exit(-1)
    
    # blast
    blast_results = blastPDB(seq,hitlist_size=opt.hitlist_size,expect=opt.expect)
    
    hits = blast_results.getHits(percent_identity=identity)
    
    #sort hits by decreasing percent identity
    hits2 = []
    for pdb in hits:
        hits2.append( (-hits[pdb]['percent_identity'],pdb) )
    hits2.sort()
    
    # print hits
    for identity,pdb in hits2:
        chain = hits[pdb]['chain_id']
        percent_identity = hits[pdb]['percent_identity']
        title = hits[pdb]['pdb_title']
        print pdb,chain,'%5.1f%%'%(percent_identity),title
    
    # download hits if --folder is given
    if opt.folder != '':
        print 'Downloading hits to',opt.folder
        pdblist = [ pdb for identity,pdb in hits2 ]
        pdblist2 = fetchPDB(pdblist,opt.folder)

if __name__ == '__main__':
    main()
