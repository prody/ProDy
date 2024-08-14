"""This module contains unit tests for :mod:`~prody.select`."""

import os
import os.path
import inspect
import numpy as np
from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *
from prody.atomic.atommap import DUMMY

from prody.utilities import importDec
dec = importDec()

prody.atomic.select.DEBUG = False
LOGGER.verbosity = 'none'

TESTS_PATH = os.path.abspath(os.path.split(inspect.getfile(
                                                   inspect.currentframe()))[0])

# If a selection string is paired with None, SelectionError is expected
# If two selection strings are paired, they must select exactly same of atoms
# Else, number must be the number atoms that the string is expected to select

pdb3mht = prody.parsePDB(pathDatafile('pdb3mht.pdb'), secondary=True)

SELECTION_TESTS = {'pdb3mht':
    {'n_atoms': len(pdb3mht),
     'ag': pdb3mht,
     'all': pdb3mht.all,
     'atommap': AtomMap(pdb3mht, [DUMMY] + list(range(1500)) + [DUMMY] +
                        list(range(1500, len(pdb3mht))) + [DUMMY]),

     'test_flags':  [('none', 0),
        ('all', 3211),
        ('acidic', 334),
        ('acyclic', 2040),
        ('aliphatic', 821),
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
        ('nucleotide', 509),
        ('nucleobase', 0),
        ('nucleoside', 0),
        ('polar', 1607),
        ('protein', 2606, 'aminoacid'),
        ('stdaa', 2606),
        ('nonstdaa', 0),
        ('purine', 0),
        ('pyrimidine', 0),
        ('small', 288),
        ('sugar', 0),
        ('surface', 1662),
        ('water', 70),
        ('hetero', 96),
        ('hetatm', 96),
        ('calpha', 327, 'ca'),
        ('backbone', 1308, 'bb'),
        ('backbonefull', 1309, 'bbfull'),
        ('sidechain', 1297, 'sc'),
        ('carbon', 1920),
        ('hydrogen', 0),
        ('noh', 3211),
        ('nitrogen', 542),
        ('oxygen', 711),
        ('sulfur', 14),
        ('extended', 503),
        ('helix', 763),
        ('helix310', 118),
        ('turn', 0),
        ('bridge', 0),
        ('bend', 0),
        ('coil', 1222),],

     'test_without_and': [
        ('coil protein', 1222),
        ('sidechain sc protein', 1297),
        ('bbfull bb', 1308),
        ('(charged basic)', 450),
        ('(protein nucleic)', 0),
        ('noh hetero water', 70, 'water hetero noh'),
        ('ca occupancy > 0', 327, 'occupancy > 0 ca'),
        ('ca occupancy - 0 > 0', 327, 'occupancy - 0 > 0 ca'),
        ('ca occupancy - 0 > 0 + 0', 327,
        'occupancy - 0 > 0 + 0 ca'),
        ('occupancy > ca 0', None),
        ('noh hetero (water)', 70),
        ('noh hetero not (water)', 26),
        ('(water) hetero', 70),
        ('ca abs(beta) = beta + abs(0)', 327,
        'abs(beta) = beta + abs(0) ca'),],

     'test_unintended': [
        ('abs beta = beta', 3211, 'abs (beta) = beta')],

     'test_string': [
        ('name P', 24),
        ('name P CA', 352),
        ('name `A 1`', 0),
        ('chain C', 248),
        ('chain x', 0),
        ('chain x y', 0),
        ('chain x y z', 0),
        ('chain x y z C', 248),
        ('chain C D', 521),
        ('chain CD', 0),
        ('resname DG', 132),
        ('resname DG ALA', 212),
        ('altloc A', 0),
        ('altloc _', 3211),
        ('secondary H', 763, 'helix'),
        ('secondary H E', 1266),
        ('secondary _', 605),
        ('segment _', 3211),],

     'test_integer': [
        ('index 10 20 10000', 2),
        ('serial 0', 0),
        ('serial 1 2', 2),
        ('resnum 0', 0),
        ('resnum 100 105', 13),
        ('resid 0', 0),
        ('resid 100 105', 13),
        ('resid 100 A 105', 13),
        ('fragindex 0', None),
        ('fragment 0', None),],

     'test_range':  [
        ('index 0:10', 10),
        ('index 0to10', 11, 'index 0 to 10'),
        ('serial 0:10:2', 4),
        ('serial 0:10:10', 0),
        ('resnum 10to15', 49),
        ('resnum 10:16:1', 49),
        ('resnum `-3:16:1`', 125),
        ('resid 10to15', 49),
        ('resid 10:16:1', 49),
        ('x `-10:20`', 673, 'x `-10 to 20`'),
        ('x 0:20:1', 0),
        ('beta 13.02:13.01', None)],

     'test_float':  [
        ('beta 5.0 41.15 11.85', 2),
        ('occupancy 1.0', 3211),
        ('x 6.665', 1),
        ('y 69.99 13.314', 2),
        ('z 115.246 45.784', 2),
        ('charge 0', 0),
        ('mass 1', 0),
        ('radius 0', None),
        ('beta "1."', 0),
        ('beta = "1."', None),],

     'test_comparisons': [
        ('x = -51.659', 1),
        ('x != -51.659', 3210),
        ('z >= 82.813', 1670),
        ('z < 82.813', 1541),
        ('beta > 10', 2874),
        ('beta < 10', 336),
        ('occupancy > 0.999999', 3211),
        ('-10 <= x < 0', 557, '-10 <= x and x < 0'),
        ('11 > 10', None),
        ('radius > 10', None),
        ('chain = A', None),
        ('x x < 1', None),
        ('name < 1', None),],

     'test_operation': [
        ('x ** 2 < 10', 238),
        ('x ** 2 ** 2 ** 2 < 10', 87),
        ('x ** (+2 ** (+2 ** +2)) < 10', 87),
        ('occupancy % 2 == 1', 3211),
        ('x**2 + y**2 + z**2 < 10000', 1975),],

     'test_function': [
        ('sqrt(x**2 + y**2 + z**2) < 100', 1975,
        'x**2 + y**2 + z**2 < 10000'),
        ('sqrt(x**2 + y**2 + z**2) == '
        '(x**2 + y**2 + z**2) ** 0.5', 3211),
        ('beta % 3 < 1', 1070),
        ('beta % 4 % 3 < 1', 1530),
        ('ceil(beta) == 10', 60),
        ('floor(beta) == 10', 58),
        ('abs(x) == sqrt(sq(x))', 3211),
        ('sq(x-5)+sq(y+4)+sq(z) > sq(100)', 1444),
        ('1 > sq(occ)', None),
        ('sq(x x) > 1', None),],

     'test_composite': [
        ('chain x y z C and x 10', 0),
        ('resnum `1` `2`', 16, 'resnum 1 2'),
        ('same residue as within 4 of resname SAH', 177),
        ('name CA and same residue as within 4 of resname SAH', 20),
        ('water and within 5 of not protein', 70),
        ('backbone and sqrt((x - 25)**2 + (y - 74)**2 + '
         '(z - 13)**2) <= 500', 1308),
        ('(not resname SAH) and (protein and name CA) or '
         '(nucleic and name P)', 351,
         '(protein and name CA) or (nucleic and name P)'),
        ('protein and (backbone or name H)', 1308),
        ('same residue as within 4 of and resname SAH', None),
        ('protein and name CA CB and same residue as '
         '((x+21.2)**2 + (y-35.9)**2 + (z-80.0)**2)**0.5 < 10',
         78, 'protein and name CA CB and same residue as '
         'within 10 of center', {'center': np.array([21.2, 35.9, 80.0])})],

     'test_within': [
        ('within 10 of index 0', 72),
        ('exwithin 100 of index 0', 3210),
        ('exwithin 4 of resname SAH', 61),
        ('(within 4 of water) and not water', 534, 'exwithin 4 of water'),
        ('within 5 of within 5 of within 5 of index 0', 135),
        ('exwithin 5 of exwithin 5 of exwithin 5 of index 0', 99),
        ('within 1 of pdb', 3211, None, {'pdb': pdb3mht}),
        ('exwithin 1 of pdb', 0, None, {'pdb': pdb3mht}),
        ('exwithin 1 of ag', None, None, {'ag': AtomGroup()}),
        ('within 100 of index 10000', 0),],

     'test_sameas': [
         ('same residue as index 0', 22),
         ('same chain as index 0', 248),
         ('same segment as index 0', 3211),
         ('same residue as resname DG ALA', 212),
         ('same chain as chain C', 248),
         ('same residue as chain X', 0),
         ('same none as chain C', None),
         ('same residue as same residue as same residue as index 0', 22,
          'resindex 0')],

     'test_regexp': [
        ('resname "S.."', 122),
        ('name "C.*"', 1920),
        ('name ".*\'"', 208),
        ('name "C(A|B)"', 628),
        ('name "C((A|B)"', None),],

     'test_specialchar': [
        ('altloc ` `', 3211),
        ('name A` `CA`', 328),
        ('name `A``', 0),
        ('z `+100.291`', 1),],

     'test_logical': [
        ('name or name', None),
        ('name and name', None),
        ('name CA and name CA', 328),
        ('name CA or name CA', 328),
        ('index 0 or index 1 ', 2),
        ('not not not not index 1', 1),
        ('index 0 or index 1 or index 2', 3, 'index 0 1 2'),
        ('index 0 or index 1 or index 2 or index 4', 4, 'index 0 1 2 4'),
        ('index 0 and index 1 ', 0),
        ('index < 50 and index < 5', 5, 'index < 5'),
        ('index < 50 and index < 25 and index < 5', 5),
        ('index < 5 and index < 25 and index < 50', 5),
        ('index 0 to 5 and index 0 to 25 and index 0 to 50', 6),
        ('index < 5 and index < 25 and index < 50 or index < 50 or index < 5',
         50),],

     'test_kwargs': [
        ('within 100 of origin', 1975, None, {'origin': np.zeros(3)}),
        ('within 100 of origin', 1975, None, {'origin': np.zeros((1, 3))}),
        ('within 100 of origin', 1975, None, {'origin': np.zeros((10, 3))}),
        ('within 100 of origin', 1975, None, {'origin': np.zeros((50, 3))}),
        ('within 100 of none', None, None, {'none': np.zeros((50, 3))}),],

     'test_equivalent':[
        ('chain C', 248, 'not not chain C'),
        ('chain C', 248, 'not not not not chain C'),
        ('nucleic', 509, 'nucleoside or nucleotide or nucleobase'),],

     'test_invalid':[
        ('chain C and and chain C', None),
        ('chain C or or chain D', None),
        ('chain C or not or chain D', None),
        ('chain C + 3', None),
        ('sqr(x-5)+sqr(y+4)+sqr(z) > sqr(100)', None),
        ('x > sq(calpha)', None),
        ('x > sq(name CA and resname ALA)', None),
        ('resname ALA and +1', None)],

     'test_userdata':[
        ('temp < 10', 336, 'beta < 10'),
        ('temp < 10 and chain D', 37, 'temp < 10 and chain D'),
        ('oc10 - 9 == 1', 3211, 'occupancy 1'),
        ('temp < 10', 336, 'temp + oc10 < 20'),
        ('occ', 3211, 'occupancy != 0'),
        ('occ and occupancy == 1', 3211, 'occupancy != 0'),
        ('occ and occupancy == 1 and oc10 - 9 == 1', 3211),
        ('occ and occupancy == 1 and temp < 10', 336),
        ('occ > 0', None),],

     'test_synonyms':    [
        ('chain C', 248, 'chid C'),
        ('chain C D', 521, 'chid C D'),],

     'test_sequence': [
        ('sequence al', 0),
        ('sequence A', 80, 'resname ALA'),
        ('sequence MIEIK', 42, 'resindex 25 to 29'),
        ('sequence VLNAL', 36, 'resindex 175 to 179'),
        ('sequence FKPY', 40, 'resindex 348 to 351'),
        ('sequence "SS."', 20, 'resindex 344 to 346'),
        ('sequence "S[A-Z]{2}G"', 79,
         'resindex 109 to 112 267 to 270 276 to 279'),
        ('sequence "S.S.S"', 0),
        ('sequence "."', 2606),],

     'test_docexamples': [
        ('serial 1 2 3', 3),
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

subsets = []
for ch in pdb3mht.iterChains():
    subsets.append((ch.getSelstr(), ch.numAtoms()))

for i, res in enumerate(pdb3mht.iterResidues()):
    if i % 80 == 0:
        subsets.append((res.getSelstr(), res.numAtoms()))

SELECTION_TESTS['pdb3mht']['subsets'] = subsets

ligand = fetchPDBLigand(pathDatafile('sti'))['ideal']

SELECTION_TESTS['imatinib'] = {
    'n_atoms': len(ligand),
    'ag': ligand,
    'all': ligand.all,
    'atommap': AtomMap(ligand, [DUMMY] + list(range(10)) + [DUMMY] +
                       list(range(10, len(ligand))) + [DUMMY]),

    'test_bondedto': [
        ('bonded to index 0', ligand[0].numBonds() + 1),
        ('exbonded to index 0', ligand[0].numBonds()),
        ('bonded to index 67', ligand[67].numBonds() + 1),
        ('exbonded to index 67', ligand[67].numBonds()),
        ('bonded 2 to index 0', 8, 'bonded to bonded to index 0'),
        ('bonded 0 to index 0', 0),
        ('bonded 3 to index 0', 10,
         'bonded to bonded to bonded to index 0'),
        ('bonded 4 to index 0', 13,
         'bonded to bonded to bonded to bonded to index 0'),
        ('bonded 4 to index 0', 13,
         'bonded 2 to bonded 2 to index 0'),
        ('bonded 4 to index 0', 13,
         'bonded to bonded 3 to index 0'),
        ('exbonded 1 to index 0', 3, 'exbonded to index 0'),
        ('exbonded 2 to index 0', 5,
         'exbonded to exbonded to index 0'),
        ('exbonded 3 to index 0', 5,
         'exbonded to exbonded to exbonded to index 0'),
        ('bonded 20 to index 0', 64),
        ('bonded 20 to index 10', 68),
        ('bonded to index 1000', 0),
        ('fragment 0', len(ligand)),
        ('fragment 1', 0),
        ('fragment 0 1 2', len(ligand)),
        ('fragindex 0 1 2', len(ligand)),
        ('fragindex 0:2', len(ligand)),],
}

pdb_gromacs = prody.parsePDB(pathDatafile('pdb6fpj_Bb_fixed_solv_ions.pdb'), secondary=True)
SELECTION_TESTS['gromacs'] = {
    'n_atoms': len(pdb_gromacs),
    'ag': pdb_gromacs,
    'all': pdb_gromacs.all,
    'test_flags':  [('protein', 12194, 'aminoacid')],
}

data_nucleoside = DATA_FILES['pdb7pbl']
pdb_nucleoside = prody.parsePDB(pathDatafile(data_nucleoside['file']), secondary=True)
SELECTION_TESTS['nucleoside'] = {
    'n_atoms': data_nucleoside['atoms'],
    'ag': pdb_nucleoside,
    'all': pdb_nucleoside.all,
    'test_flags':  [('nucleoside', data_nucleoside['nucleoside'], 'nucleoside')],
}

pdb3mht = SELECTION_TESTS['pdb3mht']['ag']
pdb3mht.setCharges(pdb3mht.getOccupancies())
pdb3mht.setMasses(pdb3mht.getBetas())
pdb3mht.setData('temp', pdb3mht.getBetas())
pdb3mht.setData('oc10', pdb3mht.getOccupancies() * 10)
pdb3mht.setFlags('occ', pdb3mht.getOccupancies().astype(bool))

SELECT = prody.Select()

EMPTYDICT = {}


class TestSelect(unittest.TestCase):

    """Test :class:`.Select`."""

    pass

count = 0
for case, items in SELECTION_TESTS.items():

    for key, tests in items.items():

        if not key.startswith('test_'):
            continue
        type_ = key[5:]

        for test in tests:
            def func(self, pdb=case, test=test, type_=type_, **kwargs):

                atoms = SELECTION_TESTS[pdb]['ag']

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
                        'selection {0} for {1} failed, expected '
                        '{2}, selected {3}'.format(repr(selstr),
                        str(atoms), natoms, len(sel)))
                else:
                    sel = SELECT.getIndices(atoms, selstr, **kwargs)
                    sel2 = SELECT.getIndices(atoms, selstr2, **kwargs)
                    self.assertTrue(len(sel) == len(sel2) == natoms and
                            np.all(sel == sel2),
                        'selection strings {0} and {1} for '
                        '{2} failed to select same number of atoms, '
                        'expected ({3})'.format(repr(selstr),
                        repr(selstr2), str(atoms), natoms))

            count += 1
            func.__name__ = 'test{0}Selection{1}'.format(
                                                type_.title(), count)
            func.__doc__ = ('Test {0} selections {1} for '
                                    '{2}').format(type_,
                                    repr(test[0]), case)
            setattr(TestSelect, func.__name__, func)

            @dec.slow
            def func(self, pdb=case, test=test, type_=type_, **kwargs):

                atoms = SELECTION_TESTS[pdb]['all']

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
                        'selection {0} for {1} failed, expected '
                        '{2}, selected {3}'.format(repr(selstr),
                        str(atoms), natoms, len(sel)))
                else:
                    sel = SELECT.getIndices(atoms, selstr, **kwargs)
                    sel2 = SELECT.getIndices(atoms, selstr2, **kwargs)
                    self.assertTrue(len(sel) == len(sel2) == natoms and
                            np.all(sel == sel2),
                        'selection strings {0} and {1} for '
                        '{2} failed to select same number of atoms, '
                        'expected ({3})'.format(repr(selstr),
                        repr(selstr2), str(atoms), natoms))

            count += 1
            func.__name__ = 'test{0}Selection{1}'.format(type_.title(),
                                                             count)
            func.__doc__ = 'Test {0} selections "{1}"'.format(type_,
                                                                  test[0])
            setattr(TestSelect, func.__name__, func)
del func

MACROS = [('cacb', 'name CA CB'),
          ('donors', '(protein) and (name N NE NH2 ND2 NE2 ND1 OG OH NH1 '
                                         'SG OG1 NE2 NZ NE1 ND1 NE2)')]


class TestMacros(unittest.TestCase):

    """Test selection macros."""


    def testMacroFunctions(self):

        for name, macro in MACROS:
            prody.defSelectionMacro(name, macro)
            self.assertEqual(prody.getSelectionMacro(name), macro,
                             'failed to get correct macro definition')
            prody.delSelectionMacro(name)

count = 0
for name, macro in MACROS:

    def func(self, name=name, macro=macro):

        prody.defSelectionMacro(name, macro)
        for key, case in SELECTION_TESTS.items():
            atoms = case['ag']
            assert_equal(
                SELECT.getBoolArray(atoms, macro),
                SELECT.getBoolArray(atoms, name),
                'failed to select correct selection using macro')
        prody.delSelectionMacro(name)
    count += 1

    func.__name__ = 'testMacro{0}'.format(count)
    func.__doc__ = 'Test macro *{0}*: {1}'.format(name, repr(macro))
    setattr(TestMacros, func.__name__, func)
del func


def testGetBoolArray():

    ca = pdb3mht.ca
    assert_equal(len(ca), len(SELECT.getBoolArray(ca, 'index 510')))

