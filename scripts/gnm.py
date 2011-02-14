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
import os.path

# TO DO LIST
# improve output formatting
# give more usage examples
    
def main():
    usage = '%prog [options] <pdb>\n\nProDy v{0:s} - Gaussian Network Model'.format(prody.__version__)
    parser = OptionParser(usage=usage)
    parser.add_option('-n','--nmodes',dest="nmodes",type="int",default=20,metavar="INT",help="number of non-zero eigenvalues/vectors to calculate (%default)")    
    parser.add_option('-c','--cutoff',dest="cutoff",type="float",default=10.0,metavar="FLOAT",help="cutoff distance (%default)")
    parser.add_option('-g','--gamma',dest="gamma",type="float",default=1.0,metavar="FLOAT",help="spring constant (%default)")
    parser.add_option('-p','--prefix',dest="prefix",type="string",default="gnm",metavar="STRING",help="prefix for output files (%default)")
    parser.add_option('-s','--select',dest="select",type="string",default="protein and name CA",metavar="STRING",help="selection string (%default)")
    parser.add_option('','--silent',dest='silent',action='store_true',default=False,help='omit verbose information (%default)')
    parser.add_option('-e','--examples',dest='examples',action='store_true',default=False,help='show usage examples')
    
    usage_examples="""
$> gnm 1p38
Fetch pdb 1p38 and run GNM using default parameters.
    
$> gnm 1p38 --select "protein and name CA and chain A'
Fetch pdb 1p38 and run GNM using default parameters on all CA atoms of protein chain A.
"""
    
    opt, args = parser.parse_args()
    
    if opt.examples:
        print usage_examples
        sys.exit(-1)
    
    if len(args) != 1:
        parser.print_help()
        print "\nError: <pdb> missing\n"
        sys.exit(-1)

    # parameters
    pdbfn = args[0]
    nmodes, select, cutoff, gamma, prefix, silent = opt.nmodes, opt.select, opt.cutoff, opt.gamma, opt.prefix, opt.silent

    if silent:
        ProDySetVerbosity('warning')
    
    # parsing pdb
    pdb = parsePDB(pdbfn)
    pdbselect = pdb.select(select)
        
    # gnm calculation
    gnm = GNM(os.path.splitext(os.path.split(pdbfn)[1])[0])
    gnm.buildKirchhoff(pdbselect,cutoff,gamma)
    gnm.calcModes(nmodes)
    
    # output
    eigenvaluesfn = "{0:s}_evalues.csv".format(prefix)
    eigenvectorsfn = "{0:s}_evectors.csv".format(prefix)
    
    writeArray(eigenvectorsfn, gnm.getArray(), delimiter=',')
    writeArray(eigenvaluesfn, gnm.getEigenvalues(), delimiter=',')
    
if __name__ == '__main__':
    main()
