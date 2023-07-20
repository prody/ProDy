# -*- coding: utf-8 -*-

"""This module detects, predicts and analyzes water bridges.
"""

__author__ = 'Karolina Mikulska-Ruminska'
__credits__ = ['Frane Doljanin', 'Karolina Mikulska-Ruminska']
__email__ = ['karolamik@fizyka.umk.pl', 'fdoljanin@pmfst.hr']


from collections import deque
from enum import Enum, auto
from numpy import *

from prody import LOGGER
from prody.measure import calcAngle
from prody.measure.contacts import findNeighbors

from timeit import default_timer as timer



__all__ = ['calcWaterBridges']


class ResType(Enum):
    WATER = auto()
    PROTEIN = auto()


class HydrophilicGroup:
    def __init__(self, hydrophilicAtom, hydrogens, type):
        self.hydrophilic = hydrophilicAtom
        self.hydrogens = hydrogens if hydrogens else []
        self.type = type


def getResidueAtoms(atom, atomGroup):
    return atomGroup.select(f'resnum {atom.getResnum()}')


def checkHBondAngle(donorGroup, acceptorGroup, angleRange):
    for hydrogen in donorGroup.hydrogens:
        bondAngle = calcAngle(donorGroup.hydrophilic,
                              hydrogen, acceptorGroup.hydrophilic)
        if angleRange[0] < bondAngle < angleRange[1]:
            return True

    return False


def getNeighborsForAtom(atom, selection, dist):
    return list(map(lambda p: p[1], findNeighbors(atom, dist, selection)))


def calcBridges_cluster(relations, hydrophilicList):
    waterBridges = []

    for atom in hydrophilicList:
      observedGroup = relations[atom.getIndex()][0]
      if not observedGroup:
        continue

      waterBridges += getBridgeChains_BFS(observedGroup, relations)

    return waterBridges

def getRelationsList(hydrogenBonds, numAtoms):
    relations = [[None, [], [], False] for _ in range(numAtoms)]
    for bond in hydrogenBonds:
        index_0, index_1 = bond[0].hydrophilic.getIndex(), bond[1].hydrophilic.getIndex()
        relations[index_0][0] = bond[0]
        relations[index_1][0] = bond[1]

        relations[index_0][1].append(bond)
        relations[index_1][1].append(bond)
        
        relations[index_0][2].append(index_1)
        relations[index_1][2].append(index_0)

    return relations

def resetRelationsList(relationList):
    for el in relationList:
      el[3] = False

def getBridges_BFS(watchedAtom, relationsList, maxDepth):
    queue = deque([(watchedAtom.hydrophilic.getIndex(), [])])
    waterBridges = []

    while queue:
        currentAtomIndex, currentPath = queue.popleft()
        currentAtom = relationsList[currentAtomIndex][0]

        if relationsList[currentAtomIndex][3]:
            continue
        
        relationsList[currentAtomIndex][3] = True
        
        if currentAtom.type == ResType.PROTEIN and len(currentPath):
            waterBridges.append(currentPath + [currentAtom])
            continue
        
        if len(currentPath) == maxDepth:
            continue
        
        for index in relationsList[currentAtomIndex][2]:
            queue.append((index, currentPath + [currentAtom]))
    
    return waterBridges

def getBridgeChains_BFS(watchedAtom, relationsList):
    waterBridges = []
    residueWaterIndices = relationsList[watchedAtom.hydrophilic.getIndex()][2]

    for waterIndex in residueWaterIndices:
        bridgeAtoms = []
        queue = deque([waterIndex])
        while queue:
            currentAtomIndex = queue.popleft()
            currentAtom = relationsList[currentAtomIndex][0]

            if relationsList[currentAtomIndex][3]:
                continue
            
            bridgeAtoms.append(currentAtomIndex)

            if currentAtom.type == ResType.WATER:
              relationsList[currentAtomIndex][3] = True

            if currentAtom.type == ResType.PROTEIN:
                continue
            
            for index in relationsList[currentAtomIndex][2]:
                queue.append(index)
        
        bridgeAtoms = list(set(bridgeAtoms))
        if sum(relationsList[atomIndex][0].type == ResType.PROTEIN for atomIndex in bridgeAtoms) >= 2:
          waterBridges.append(bridgeAtoms)
    
    return waterBridges

def calcWaterBridges(atoms, **kwargs):
    """Compute water bridges for a protein that has water molecules.

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`

    :arg distDA: maximal distance between water/protein donor and acceptor
        default is 3.5
    :type distDA: int, float

    :arg distWR: maximal distance between water and residue
        default is 4
    :type distWR: int, float

    :arg anglePDWA: angle range where protein is donor and water is acceptor
        default is (100, 200)
    :type anglePDWA: (int, int)

    :arg anglePAWD: angle range where protein is acceptor and water is donor
        default is (100, 140)
    :type anglePDWA: (int, int)

    :arg angleWW: angle between water donor/acceptor
        default is (140, 180)
    :type angleWW: (int, int)

    :arg waterDepth: maximum number of waters in chain
      default is 2
    :type waterDepth: int

    :arg donors: which atoms to count as donors 
        default is ['N', 'O', 'S', 'F']
    :type donors: list

    :arg acceptors: which atoms to count as acceptors 
        default is ['N', 'O', 'S', 'F']
    :type acceptors: list
    """

    distDA = kwargs.pop('distDA', 3.5)
    distWR = kwargs.pop('distWR', 4)
    anglePDWA = kwargs.pop('anglePDWA', (100, 200))
    anglePAWD = kwargs.pop('anglePAWD', (100, 140))
    angleWW = kwargs.pop('angleWW', (140, 180))
    waterDepth = kwargs.pop('waterDepth', 3)
    donors = kwargs.pop('donors', ['nitrogen', 'oxygen', 'sulfur'])
    acceptors = kwargs.pop('acceptors', ['nitrogen', 'oxygen', 'sulfur'])

    consideredAtoms = ~atoms.select(f'water and not within {distWR} of protein')
    waterHydrogens = consideredAtoms.select('water and hydrogen')
    waterOxygens = consideredAtoms.select('water and oxygen')
    neighborWaterOxygens = findNeighbors(waterOxygens, distDA)
    contactingWaters = []

    for nbPair in neighborWaterOxygens:
        hydrogens_0, hydrogens_1 = getResidueAtoms(
            nbPair[0], waterHydrogens), getResidueAtoms(nbPair[1], waterHydrogens)
        water_0, water_1 = HydrophilicGroup(
            nbPair[0], hydrogens_0, ResType.WATER), HydrophilicGroup(nbPair[1], hydrogens_1, ResType.WATER)
        contactingWaters.append((water_0, water_1))
    LOGGER.info(f'Water contacts ({len(contactingWaters)}) saved.')

    proteinHydrogens = atoms.select('protein and hydrogen')
    proteinHydrophilic = atoms.select(
        f'protein and ({" or ".join(donors + acceptors)})')
    neighborWaterProtein = findNeighbors(
        waterOxygens, distDA, proteinHydrophilic)
    contactingWaterProtein = []

    for nbPair in neighborWaterProtein:
        hydrogens_w, hydrogens_p = getResidueAtoms(
            nbPair[0], waterHydrogens), getNeighborsForAtom(nbPair[1], proteinHydrogens, 1.4)
        water, protein = HydrophilicGroup(
            nbPair[0], hydrogens_w, ResType.WATER), HydrophilicGroup(nbPair[1], hydrogens_p, ResType.PROTEIN)
        contactingWaterProtein.append((water, protein))
    LOGGER.info(
        f'Water-protein contacts ({len(contactingWaterProtein)}) saved.')

    hydrogenBonds = []

    for nbPair in contactingWaters:
        if checkHBondAngle(nbPair[0], nbPair[1], angleWW):
            hydrogenBonds.append((nbPair[0], nbPair[1]))

        if checkHBondAngle(nbPair[1], nbPair[0], angleWW):
            hydrogenBonds.append((nbPair[1], nbPair[0]))
    LOGGER.info(f'Hydrogen bonds ({len(hydrogenBonds)}) calculated.')

    for nbPair in contactingWaterProtein:
        if nbPair[1].hydrophilic.getName()[0] in ['N', 'S', 'F', 'O'] and checkHBondAngle(nbPair[0], nbPair[1], anglePDWA):
            hydrogenBonds.append((nbPair[0], nbPair[1]))

        if nbPair[1].hydrophilic.getName()[0] in ['N', 'S', 'F', 'O'] and checkHBondAngle(nbPair[1], nbPair[0], anglePAWD):
            hydrogenBonds.append((nbPair[1], nbPair[0]))
    LOGGER.info(f'Hydrogen bonds ({len(hydrogenBonds)}) calculated.')

    start = timer()
    relations = getRelationsList(hydrogenBonds, len(atoms))
    waterBridges_2 = calcBridges_cluster(relations, proteinHydrophilic)
    end = timer()
    LOGGER.info(end-start)
    #waterBridges = []
    #for atom in proteinHydrophilic:
    #  observedGroup = relations[atom.getIndex()][0]
    #  if not observedGroup:
    #    continue
#
    #  waterBridges += getBridges_BFS(relations[atom.getIndex()][0], relations, waterDepth+1)
    #  resetRelationsList(relations)
    LOGGER.info(f'Water bridges (2) calculated.')
    return waterBridges_2
