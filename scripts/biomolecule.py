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
from prody import *
import os.path

def main():
    if len(sys.argv) != 2:
        print 'Usage: %s <pdb>\n\nGenerate biomolecular transformation recorded on pdb header.' % (sys.argv[0])
        sys.exit(-1)
    
    pdbfn = sys.argv[1]
    
    pdb, header = parsePDB(pdbfn, header=True)
    biomols = applyBiomolecularTransformations(header, pdb)
    if not isinstance(biomols, list):
        biomols = [biomols]
    fn = os.path.splitext(os.path.split(pdbfn)[1])[0]
    for i, biomol in enumerate(biomols):
        if isinstance(biomol, Atomic):
            outfn = '{0:s}_biomol_{1:d}.pdb'.format(fn, i+1)
            print('Writing {0:s}'.format(outfn))
            writePDB(outfn, biomol)
        elif isinstance(biomol, tuple):
            for j, part in enumerate(biomol):
                outfn = '{0:s}_biomol_{1:d}_part_{2:d}.pdb'.format(fn, i+1, j+1)
                print('Writing {0:s}'.format(outfn))
                writePDB(outfn, part)
                
if __name__ == '__main__':
    main()
