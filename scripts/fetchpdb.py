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

def main():
    usage = '%prog 1p38 1r39...\n\nDownload pdb files 1p38 and 1r39 from RCSB PDB.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('','--folder',dest="folder",type="string",default=".",metavar="PATH",help="destination folder (%default)")
    parser.add_option('','--list',dest="listfn",type="string",default="",metavar="FILE",help="file with pdb identifiers")
    
    opt, args = parser.parse_args()
    if len(args) == 0 and opt.listfn == '':
        parser.print_help()
        sys.exit(-1)

    pdblist = []
    pdblist += args
    if opt.listfn != '':
        f = open(opt.listfn)
        for line in f.xreadlines():
            line = line.strip()       # remove trailing spaces
            for s in line.split(','): # split by comma
                for pdb in s.split(): # split by space
                    if len(pdb) == 4: pdblist.append(pdb)
        f.close()
    
    pdblist2 = fetchPDB(pdblist, opt.folder)
    

if __name__ == '__main__':
    main()
