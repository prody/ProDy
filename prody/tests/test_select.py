# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

"""This module contains unit tests for :mod:`~prody.proteins` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import os
import os.path
import unittest
import prody

prody.changeVerbosity('none')

SELECTION_TESTS = {'data/pdb3mht.pdb':
    {'n_atoms': 3211,
     'keyword':     [('none', 0),
                     ('all', 3211),
                     ('acidic', 334),
                     ('acyclic', 2040),
                     ('aliphatic', 730),
                     ('aromatic', 475),
                     ('at', 0),
                     ('basic', 450),
                     ('buried', 944),
                     ('cg', 0),
                     ('charged', 784),
                     ('cyclic', 566),
                     ('heme', 0),
                     ('hydrophobic', 999),
                     ('ion', 0),
                     ('large', 1629),
                     ('lipid', 0),
                     ('medium', 689),
                     ('neutral', 1822),
                     ('nucleic', 509),
                     ('polar', 1607),
                     ('protein', 2606),
                     ('purine', 0),
                     ('pyrimidine', 0),
                     ('small', 288),
                     ('sugar', 0),
                     ('surface', 1662),
                     ('water', 70),],
     'string':      [('name P', 24),
                     ('name P CA', 352),
                     ('chain C', 248),
                     ('chain C D', 521),
                     ('chain CD', 0),
                     ('chain CD', 0),
                     ('resname DG', 132),
                     ('resname DG ALA', 212),
                     ('altloc A', 0),
                     ('altloc _', 3211),
                     #('altloc ``', 3211),
                     ('secondary H', 763),
                     ('secondary H E', 1266),
                     ('secondary _', 605),
                     ('segment _', 3211),],
     'integer':     [('index 0', 1),
                     ('index 10 20 30', 3),
                     ('serial 0', 0),
                     ('serial 1 2', 2),
                     ('resnum 0', 0),
                     ('resnum 100 105', 13),
                     ('resid 0', 0),
                     ('resid 100 105', 13),],
     'range':       [('index 0:10', 10),
                     ('index 0to10', 11),
                     ('index 0 to 10', 11),
                     ('serial 0:10:2', 4),
                     ('serial 0:10:10', 0),
                     ('resnum 10to15', 49),
                     ('resnum 10:16:1', 49),
                     ('resid 10to15', 49),
                     ('resid 10:16:1', 49),],
     'float':       [('none', 0),
                     ('none', 0),
                     ('charge 0', None),
                     ('mass 1', None),],
     'comparison':  [('none', 0),
                     ('none', 0),
                     ('none', 0),
                     ('none', 0),],
     'operation':   [('none', 0),
                     ('none', 0),
                     ('none', 0),
                     ('none', 0),],
     'function':    [('none', 0),
                     ('none', 0),
                     ('none', 0),
                     ('none', 0),],
     'composite':   [('same residue as within 4 of resname SAH', 177),
                     ('name CA and same residue as within 4 of resname SAH', 20),
                     ('none', 0),
                     ('none', 0),
                     ('none', 0),],
     'within':      [('within 10 of index 0', 72),
                     ('exwithin 100 of index 0', 3210),
                     ('exwithin 4 of resname SAH', 61),],
     'sameas':      [('same residue as index 0', 22),
                     ('same chain as index 0', 248),   
                     ('none', 0),
                     ('none', 0),],
    }

}

class TestSelectMeta(type):
    
    def __init__(cls, name, bases, dict):
        
        test_types = set()
        for case in SELECTION_TESTS.itervalues():
            for key, item in case.iteritems():
                if isinstance(item, list):
                    test_types.add(key.lower())
        
        for test in test_types:
            
            def testFunction(self, test=test):
                
                for key, testsets in SELECTION_TESTS.iteritems():
                    atoms = self.atomgroups[key]
                    
                    tests = testsets.get(test, [])
                    for selstr, natoms in tests:
                        if natoms is None:
                            self.assertRaises(prody.select.SelectionError,
                                self.select.getIndices, atoms, selstr)
                        else:
                            sel = self.select.getIndices(atoms, selstr)
                            self.assertEqual(len(sel), natoms,
                                'selection "{0:s}" for {1:s} failed, expected '
                                '{2:d}, selected {3:d}'
                                .format(selstr, str(atoms), natoms, len(sel)))
            testFunction.__name__ = 'test' + test.title() + 'Selections'
            testFunction.__doc__ = 'Test {0:s} selections.'.format(test)
            setattr(cls, testFunction.__name__, testFunction)

class TestSelect(unittest.TestCase):
    
    """Test :func:`~prody.proteins.fetchPDB` function."""
    __metaclass__ = TestSelectMeta
    
    def setUp(self):
        """Instantiate a list for storing downloaded file names."""
        
        self.select = prody.Select()
        self.atomgroups = {}
        for pdb in SELECTION_TESTS.iterkeys(): 
            self.atomgroups[pdb] = prody.parsePDB(pdb, secondary=True)
            
    def testAtomGroups(self):    
        
        for key, atoms in self.atomgroups.iteritems():
            self.assertEqual(
                atoms.getNumOfAtoms(), SELECTION_TESTS[key]['n_atoms'],
                'parsePDB failed to parse correct number of atoms from {0:s}'
                .format(key))

if __name__ == '__main__':
    unittest.main()
