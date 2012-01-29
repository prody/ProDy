#!/usr/bin/python
# -*- coding: utf-8 -*-
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

"""This module contains unit tests for :mod:`~prody.select`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os
import os.path
import unittest
import inspect
import sys
import numpy as np

from prody import *
from test_datafiles import *


prody.select.DEBUG = False
prody.setVerbosity('none')

TESTS_PATH = os.path.abspath(os.path.split(inspect.getfile(
                                                   inspect.currentframe()))[0])

# If a selection string is paired with None, SelectionError is expected
# If two selection strings are paired, they must select exactly same of atoms
# Else, number must be the number atoms that the string is expected to select 

SELECTION_TESTS = {'pdb3mht':
    {'n_atoms': 3211,
     'path': getDatafilePath('pdb3mht.pdb'),
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
                     ('water', 70),
                     ('hetero', 96),
                     ('calpha', 327),
                     ('backbone', 1308),
                     ('bb', 1308),
                     ('backbonefull', 1309),
                     ('bbfull', 1309),
                     ('sidechain', 1297),
                     ('sc', 1297),
                     ('carbon', 1920),
                     ('hydrogen', 0),
                     ('noh', 3211),
                     ('nitrogen', 542),
                     ('oxygen', 711),
                     ('sulfur', 14),
                     ('extended', 503),
                     ('helix', 763),
                     ('helix_3_10', 0),
                     ('turn', 0),
                     ('bridge', 0),
                     ('bend', 0),
                     ('coil', 1222),],
     'string':      [('name P', 24),
                     ('name P CA', 352),
                     ('name `A 1`', 0), 
                     ('name `A *`', 0),
                     ('chain C', 248),
                     ('chain C D', 521),
                     ('chain CD', 0),
                     ('chain CD', 0),
                     ('resname DG', 132),
                     ('resname DG ALA', 212),
                     ('altloc A', 0),
                     ('altloc _', 3211),
                     ('secondary H', 763),
                     ('secondary H', 763, 'helix'),
                     ('secondary H E', 1266),
                     ('secondary _', 605),
                     ('segment _', 3211),],
     'integer':     [('index 0', 1),
                     ('index 10 20 30', 3),
                     ('index 10 20 10000', 2),
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
                     ('resnum `-3:16:1`', 125),
                     ('resid 10to15', 49),
                     ('resid 10:16:1', 49),
                     ('x `-10:20`', 673),
                     ('x `-10 to 20`', 673),
                     ('x 0:20:1', None),],
     'float':       [('beta 5.0 41.15 11.85', 2),
                     ('occupancy 1.0', 3211),
                     ('x 6.665', 1),
                     ('y 69.99 13.314', 2),
                     ('z 115.246 45.784', 2),
                     ('charge 0', 0),
                     ('mass 1', 0),
                     ('radius 0', None),],
     'comparison':  [('x = -51.659', 1),
                     ('x != -51.659', 3210),
                     ('z >= 82.813', 1670),
                     ('z < 82.813', 1541),
                     ('beta > 10', 2874),
                     ('beta < 10', 336),
                     ('occupancy > 0.999999', 3211),
                     ('radius > 10', None),
                     ('chain = A', None),],
     'operation':   [('x ** 2 < 10', 238),
                     ('x ** 2 ** 2 ** 2 < 10', 99),
                     ('x ** (+2 ** (+2 ** +2)) < 10', 87),
                     ('occupancy % 2 == 1', 3211),
                     ('x**2 + y**2 + z**2 < 10000', 1975),],
     'function':    [('sqrt(x**2 + y**2 + z**2) < 100', 1975,
                      'x**2 + y**2 + z**2 < 10000'),
                     ('sqrt(x**2 + y**2 + z**2) == '
                      '(x**2 + y**2 + z**2) ** 0.5', 3211),
                     ('beta % 3 < 1', 1070),
                     ('beta % 4 % 3 < 1', 1530),
                     ('ceil(beta) == 10', 60),
                     ('floor(beta) == 10', 58),
                     ('abs(x) == sqrt(sq(x))', 3211), 
                     ('sq(x-5)+sq(y+4)+sq(z) > sq(100)', 1444),
                     ],
     'composite':   [('same residue as within 4 of resname SAH', 177),
                     ('name CA and same residue as within 4 of resname SAH', 
                      20),
                     ('water and within 5 of not protein', 70),
                     ('backbone and sqrt((x - 25)**2 + (y - 74)**2 + '
                      '(z - 13)**2) <= 500', 1308),
                     ('not resname SAH and (protein and name CA) or '
                      '(nucleic and name P)', 351,
                      '(protein and name CA) or (nucleic and name P)'), 
                      ('protein and (backbone or name H)', 1308),
                      ('same residue as within 4 of and resname SAH', None),],
     'within':      [('within 10 of index 0', 72),
                     ('exwithin 100 of index 0', 3210),
                     ('exwithin 4 of resname SAH', 61),
                     ('within 4 of water and not water', 534, 
                      'exwithin 4 of water'),],
     'sameas':      [('same residue as index 0', 22),
                     ('same chain as index 0', 248),   
                     ('same segment as index 0', 3211),
                     ('same residue as resname DG ALA', 212),
                     ('same chain as chain C', 248),],
     'regexp':      [('resname "S.."', 122),
                     ('name "C.*"', 1920),
                     ('name ".*\'"', 208),
                     ('name "C(A|B)"', 628),],
     'specialchar': [('altloc ``', 3211),
                     ('altloc ` `', 3211),
                     ('z `+100.291`', 1),],
     'logical':     [('name or name', None),
                     ('name and name', None),
                     ('name CA and name CA', 328),
                     ('name CA or name CA', 328),],
      'kwargs':     [('within 100 of origin', 1975, None, 
                      {'origin': np.zeros(3)}),
                     ('within 100 of origin', 1975, None, 
                      {'origin': np.zeros((1, 3))}),                      
                     ('within 100 of origin', 1975, None, 
                      {'origin': np.zeros((10, 3))}),],
      'equivalent': [('chain C', 248, 'not not chain C'),
                     ('chain C', 248, 'not not not not chain C'),],
      'invalid':    [('chain C and and chain C', None),
                     ('chain C or or chain D', None),
                     ('chain C or not or chain D', None),
                     ('chain C + 3', None),
                     ('sqr(x-5)+sqr(y+4)+sqr(z) > sqr(100)', None),
                     ('x > sq(calpha)', None),
                     ('x > sq(name CA and resname ALA)', None),
                     ('resname ALA and +1', None)],
      'attribute':  [('temp < 10', 336, 'beta < 10'),
                     ('temp < 10 and chain D', 37, 'temp < 10 and chain D'),
                     ('oc10 - 9 == 1', 3211, 'occupancy 1'),
                     ('temp < 10', 336, 'temp + oc10 < 20'),
                     ('occ', 3211, 'occupancy != 0'),
                     ('occ and occupancy == 1', 3211, 'occupancy != 0'),
                     ('occ and occupancy == 1 and oc10 - 9 == 1', 3211),
                     ('occ and occupancy == 1 and temp < 10', 336),],
      'synonyms':   [('chain C', 248, 'chid C'),
                     ('chain C D', 521, 'chid C D'),],
      'docexamples':[('serial 1 2 3', 3),
                     ('serial 1 to 10', 10),
                     ('serial 1:10:2', 5),
                     ('serial < 10', 9),
                     ('beta 555.55', 0),
                     ('beta 1 to 500', 3211),
                     ('beta 1:500', 3211),
                     ('beta < 500', 3211),
                     ('resnum 120A 120B', 0),
                     ('icode A', 0),
                     ('icode _', 3211),
                     ('charge 1', 3211),
                     ('abs(charge) == 1', 3211),
                     ('charge < 0', 0),
                     ('0 < mass < 500', 3211),
                     ('abs(mass) <= mass <= 10', 337),],
    }

}
for key, item in SELECTION_TESTS.iteritems():
    ag = prody.parsePDB(item['path'], secondary=True)
    ag.setCharges(ag.getOccupancies())
    ag.setMasses(ag.getBetas())
    ag.setData('temp', ag.getBetas())
    ag.setData('oc10', ag.getOccupancies() * 10)
    ag.setData('occ', ag.getOccupancies().astype(bool))
    SELECTION_TESTS[key]['ag'] = ag
    
SELECT = prody.Select()

EMPTYDICT = {}

class TestSelectMeta(type):
    
    def __init__(cls, name, bases, dict):
        
        test_types = set()
        for case in SELECTION_TESTS.itervalues():
            for key, item in case.iteritems():
                if isinstance(item, list):
                    test_types.add(key.lower())

        for type_ in test_types:
            count = 0        
            for key, testsets in SELECTION_TESTS.iteritems():
                tests = testsets.get(type_, [])
                for test in tests:
                    def testFunction(self, pdb=key, test=test, type_=type_, 
                                     **kwargs):
                
                        atoms = SELECTION_TESTS[key]['ag']
                    
                        selstr = test[0]
                        natoms = test[1]
                        selstr2 = None
                        kwargs = EMPTYDICT
                        if len(test) == 3:
                            selstr2 = test[2]
                        if len(test) == 4:
                            kwargs = test[3]

                            
                        if natoms is None:
                            self.assertRaises(prody.select.SelectionError,
                                SELECT.getIndices, atoms, selstr, **kwargs)
                        elif selstr2 is None:
                            sel = SELECT.getIndices(atoms, selstr, **kwargs)
                            self.assertEqual(len(sel), natoms,
                                'selection "{0:s}" for {1:s} failed, expected '
                                '{2:d}, selected {3:d}'
                                .format(selstr, str(atoms), natoms, len(sel)))
                        else:
                            sel = SELECT.getIndices(atoms, selstr, **kwargs)
                            sel2 = SELECT.getIndices(atoms, selstr2, **kwargs)
                            self.assertTrue(len(sel) == len(sel2) == natoms and
                                    np.all(sel == sel2),
                                'selection strings "{0:s}" and "{1:s}" for '
                                '{2:s} failed to select same number of atoms, '
                                'expected ({3:d})'
                                .format(selstr, selstr2, str(atoms), natoms))
                                
                    #testFunction.__name__ = 'test' + type_.title() + 'Selections'
                    count += 1
                    testFunction.__name__ = 'test{0:s}Selection{1:d}'.format(
                                                        type_.title(), count)
                    testFunction.__doc__ = 'Test {0:s} selections "{1:s}"'\
                                            .format(type_, test[0])
                    setattr(cls, testFunction.__name__, testFunction)

class TestSelect(unittest.TestCase):
    
    """Test :class:`~prody.select.select`."""
    
    __metaclass__ = TestSelectMeta


class TestGetSetFunctions(unittest.TestCase):
    
    def testGetBackboneAtomNames(self):
        bban = list(prody.select.BACKBONE_ATOM_NAMES)
        bban.sort()
        self.assertListEqual(bban, prody.getBackboneAtomNames())
        bban = list(prody.select.BACKBONE_FULL_ATOM_NAMES)
        bban.sort()
        self.assertListEqual(bban, prody.getBackboneAtomNames(True))
    
        
    def testSetBackboneAtomNames(self):

        for full, torf, defn in [
                  ('', False, list(prody.select.BACKBONE_ATOM_NAMES)), 
                  ('full', True, list(prody.select.BACKBONE_FULL_ATOM_NAMES))]:
            for key, case in SELECTION_TESTS.iteritems():
                atoms = case['ag']
                sel1 = SELECT.getIndices(atoms, 'backbone' + full)
                sel2 = SELECT.getIndices(atoms, 'protein and name CB or '
                                                'backbone' + full)        
                prody.setBackboneAtomNames(defn + ['CB'], torf)
                sel3 = SELECT.getIndices(atoms, 'backbone' + full)
                self.assertListEqual(list(sel2), list(sel3),
                                     'failed to change "backbone' + full + '" '
                                     'atom names definition')
                prody.setBackboneAtomNames(defn, torf)
                sel4 = SELECT.getIndices(atoms, 'backbone' + full)
                self.assertListEqual(list(sel1), list(sel4),
                                     'failed to reset "backbone' + full + '" '
                                     'atom names definition')
    
MACROS = [('cacb', 'name CA CB')]

class TestSelectionMacros(unittest.TestCase):

    def testMacroFunctions(self):

        for name, macro in MACROS:
            prody.defSelectionMacro(name, macro)            
            self.assertEqual(prody.getSelectionMacro(name), macro,
                             'failed to get correct macro definition')        
            prody.delSelectionMacro(name)            

    def testSelections(self):
        
        for name, macro in MACROS:
            prody.defSelectionMacro(name, macro)
            for key, case in SELECTION_TESTS.iteritems():
                atoms = case['ag']
                sel1 = SELECT.getIndices(atoms, macro)
                sel2 = SELECT.getIndices(atoms, name)
                self.assertListEqual(list(sel1), list(sel2),
                                     'failed to select correct selection '
                                     'using macro')        
            prody.delSelectionMacro(name)


