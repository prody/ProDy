# ProDy: A Python Package for Protein Structural Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
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

""":mod:`hierview` module defines a class for a hieararchical view of an atom 
group and atom selections

Classes:
    
  * :class:`HierView`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

from collections import defaultdict

import numpy as np

import prody
from prody import ProDyLogger as LOGGER
from . import Chain, Residue, Selection, DTYPES


__all__ = ['HierView']


class HierView(object):
    
    
    __slots__ = ['_atoms', '_chains']
    
    def __init__(self, atoms):
        """Instantiate a hierarchical view for atoms in an :class:`AtomGroup` 
        or :class:`Selection`."""
        self._atoms = atoms
        self._chains = dict()
        self.build()
        
    def build(self):
        """Build hierarchical view of the atom group.
        
        This method is called at instantiation, but can be used to rebuild
        the hierarchical view when attributes of atoms change.
        
        """
        
        acsi = self._atoms.getActiveCoordsetIndex()
        if isinstance(self._atoms, prody.proteins.AtomGroup):
            atomgroup = self._atoms
            _indices = np.arange(atomgroup._n_atoms)
            if atomgroup._chainids is None:
                #LOGGER.warning('chain identifiers of {0:s} were not assigned'.format(str()))
                atomgroup._chainids = np.zeros(atomgroup._n_atoms, dtype=DTYPES['chainids'])
            chainids = atomgroup._chainids
        else:
            atomgroup = self._atoms._ag
            _indices = self._atoms._indices
            if atomgroup._chainids is None:
                atomgroup._chainids = np.zeros(atomgroup._n_atoms, dtype=DTYPES['chainids'])
            chainids = atomgroup._chainids[_indices]


        for chid in np.unique(chainids):
            ch = Chain(atomgroup, _indices[chainids == chid], chid, acsi)
            self._chains[chid] = ch
        
        if atomgroup._resnums is None:
            atomgroup._resnums = np.zeros(atomgroup._n_atoms, dtype=DTYPES['resnums'])
        if atomgroup._resnames is None:
            atomgroup._resnames = np.zeros(atomgroup._n_atoms, dtype=DTYPES['resnames'])
        if atomgroup._icodes is None:
            atomgroup._icodes = np.zeros(atomgroup._n_atoms, dtype=DTYPES['icodes'])

        icodes = atomgroup._icodes
        

        for chain in self.iterChains():
            chid = chain._chid
            rd = defaultdict(list)
            indices = chain._indices
            resnums = chain.getResidueNumbers()
            for i in xrange(len(resnums)):
                rd[resnums[i]].append(indices[i])
            
            resnums = rd.keys()
            resnums.sort()
            for resnum in resnums:
                resindices = np.array(rd[resnum])
                res_icodes = icodes[resindices]
                
                for ic in np.unique(res_icodes): 
                    subindices = resindices[res_icodes == ic]
                    temp = subindices[0]
                    res = Residue(atomgroup, subindices, resnum, 
                                 atomgroup._resnames[temp], 
                                 atomgroup._icodes[temp], chain, acsi)   
                    chain._dict[(resnum, ic)] = res
        
    def __repr__(self):
        return '<HierView: {0:s}>'.format(str(self._atoms))
    
    def __str__(self):
        return 'HierView of {0:s}'.format(str(self._atoms))
    
    def __iter__(self):
        """Iterate over chains."""
        return self.iterChains()
    
    def iterResidues(self):
        """Iterate over residues."""
        chids = self._chains.keys()
        chids.sort()
        for chid in chids:
            chain = self._chains[chid]
            for res in chain.iterResidues():
                yield res

                
    def getResidue(self, chainid, resnum, icode=''):
        """Return residue with number *resnum* and insertion code *icode* from 
        the chain with identifier *chainid*, if it exists."""
        ch = self._chains.get(chainid, None)
        if ch is not None:
            return ch.getResidue(resnum, icode)
        return None

    def getNumOfResidues(self):
        """Returns number of residues."""
        return sum([ch.getNumOfResidues() for ch in self._chains.itervalues()])    

    def iterChains(self):
        """Iterate over chains."""
        chids = self._chains.keys()
        chids.sort()
        for chid in chids:
            yield self._chains[chid]
    
    def getChain(self, chainid):
        """Return chain with identifier *chainid*, if it exists."""
        return self._chains.get(chainid, None)

    def getNumOfChains(self):
        """Return number of chains."""
        return len(self._chains)    
