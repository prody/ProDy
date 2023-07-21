# -*- coding: utf-8 -*-

"""This module detects, predicts and analyzes water bridges.
"""

__author__ = 'Karolina Mikulska-Ruminska'
__credits__ = ['Frane Doljanin', 'Karolina Mikulska-Ruminska']
__email__ = ['karolamik@fizyka.umk.pl', 'fdoljanin@pmfst.hr']


from collections import deque
from enum import Enum, auto

from prody import LOGGER, Atom
from prody.measure import calcAngle
from prody.measure.contacts import findNeighbors

from timeit import default_timer as timer


__all__ = ['calcWaterBridges']


class ResType(Enum):
    WATER = auto()
    PROTEIN = auto()


class AtomNode:
    def __init__(self, hydrophilicAtom, hydrogens, type):
        self.atom = hydrophilicAtom
        self.hydrogens = hydrogens if hydrogens else []
        self.type = type


class AtomNodeWithRelations(AtomNode):
    def __init__(self, atomNode):
        self.atom = atomNode.atom
        self.hydrogens = atomNode.hydrogens
        self.type = atomNode.type
        self.bonds = []
        self.bondedAtoms = []
        self.isVisited = False


class RelationList:
    def __init__(self, numNodes):
        self.nodes = [None for _ in range(numNodes)]

    def resetVisits(self):
        for node in self.nodes:
            if node:
                node.isVisited = False

    def createIfNotExists(self, node):
        if not self[node.atom]:
            self[node.atom.getIndex()] = AtomNodeWithRelations(node)

        return self[node.atom]

    def __getitem__(self, key):
        if isinstance(key, Atom):
            key = key.getIndex()

        return self.nodes[key]

    def __setitem__(self, key, value):
        self.nodes[key] = value


class HydrogenBond:
    def __init__(self, donor, acceptor):
        self.donor = donor
        self.acceptor = acceptor


def getSameResidueAtoms(observedAtom, atoms):
    return atoms.select(f'resnum {observedAtom.getResnum()}')


def checkHBondAngle(donor, acceptor, angleRange):
    for hydrogen in donor.hydrogens:
        bondAngle = calcAngle(donor.atom, hydrogen, acceptor.atom)
        if angleRange[0] < bondAngle < angleRange[1]:
            return True

    return False


def getNeighborsForAtom(observedAtom, atoms, dist):
    return list(map(lambda p: p[1], findNeighbors(observedAtom, dist, atoms)))


def getRelationsList(hydrogenBonds, numAtoms):
    relations = RelationList(numAtoms)
    for bond in hydrogenBonds:
        atomNode_0 = relations.createIfNotExists(bond.donor)
        atomNode_1 = relations.createIfNotExists(bond.acceptor)

        atomNode_0.bonds.append(bond)
        atomNode_1.bonds.append(bond)

        atomNode_0.bondedAtoms.append(atomNode_1.atom)
        atomNode_1.bondedAtoms.append(atomNode_0.atom)

    return relations


def calcBridges(relations, hydrophilicList, waterDepth=None):
    waterBridges = []

    for atom in hydrophilicList:
        observedNode = relations[atom]
        if not observedNode:
            continue

        newBridges = []
        if waterDepth:
            newBridges = getBridgeChain_BFS(
                observedNode, relations, waterDepth+1)
            relations.resetVisits()
        else:
            newBridges = getBridgeForResidue_BFS(observedNode, relations)

        waterBridges += newBridges

    return waterBridges


def getChainBridgeTuple(bridge):
    i_1, i_2 = bridge[0][0].getIndex(), bridge[0][-1].getIndex()
    return (min(i_1, i_2), max(i_1, i_2))


def getUniqueElements(list, key):
    unique, uniqueKeys = [], []
    for element in list:
        elementKey = key(element)
        if elementKey in uniqueKeys:
            continue

        unique.append(element)
        uniqueKeys.append(elementKey)

    return unique


def getBridgeChain_BFS(observedNode, relationsList, maxDepth):
    queue = deque([(observedNode.atom, [], [])])
    waterBridges = []

    while queue:
        currentAtom, currentPath, currentBonds = queue.popleft()
        currentNode = relationsList[currentAtom]

        if currentNode.isVisited:
            continue

        currentNode.isVisited = True

        if currentNode.type == ResType.PROTEIN and len(currentPath):
            waterBridges.append((currentPath + [currentAtom], currentBonds))
            continue

        if len(currentPath) == maxDepth:
            continue

        for atom, bond in zip(currentNode.bondedAtoms, currentNode.bonds):
            queue.append(
                (atom, currentPath + [currentAtom], currentBonds + [bond]))

    return waterBridges


def getBridgeForResidue_BFS(observedNode, relationsList):
    waterBridges = []
    observedRelNode = relationsList[observedNode.atom]

    for waterAtom in observedRelNode.bondedAtoms:
        bridgeAtomIndices = []
        queue = deque([waterAtom])
        while queue:
            currentAtom = queue.popleft()
            currentNode = relationsList[currentAtom]

            if currentNode.isVisited:
                continue

            bridgeAtomIndices.append(currentAtom)

            if currentNode.type == ResType.PROTEIN:
                continue
            if currentNode.type == ResType.WATER:
                currentNode.isVisited = True

            for atom in currentNode.bondedAtoms:
                queue.append(atom)

        bridgeAtomIndices = list(
            set(map(lambda a: a.getIndex(), bridgeAtomIndices)))
        if sum(relationsList[atomIndex].type == ResType.PROTEIN for atomIndex in bridgeAtomIndices) >= 2:
            waterBridges.append(bridgeAtomIndices)

    return waterBridges


def calcWaterBridges(atoms, **kwargs):
    """Compute water bridges for a protein that has water molecules.

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`

    :arg method: 'cluster' | 'chain'
      default is 'chain'
    :type method: string

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

    method = kwargs.pop('method', 'chain')
    distDA = kwargs.pop('distDA', 3.5)
    distWR = kwargs.pop('distWR', 4)
    anglePDWA = kwargs.pop('anglePDWA', (100, 200))
    anglePAWD = kwargs.pop('anglePAWD', (100, 140))
    angleWW = kwargs.pop('angleWW', (140, 180))
    waterDepth = kwargs.pop('waterDepth', 3)
    donors = kwargs.pop('donors', ['nitrogen', 'oxygen', 'sulfur'])
    acceptors = kwargs.pop('acceptors', ['nitrogen', 'oxygen', 'sulfur'])

    if method not in ['chain', 'cluster']:
        raise TypeError('Method should be chain or cluster.')

    consideredAtoms = ~atoms.select(
        f'water and not within {distWR} of protein')
    waterHydrogens = consideredAtoms.select('water and hydrogen')
    waterOxygens = consideredAtoms.select('water and oxygen')
    neighborWaterOxygens = findNeighbors(waterOxygens, distDA)
    contactingWaters = []

    for nbPair in neighborWaterOxygens:
        hydrogens_0, hydrogens_1 = getSameResidueAtoms(
            nbPair[0], waterHydrogens), getSameResidueAtoms(nbPair[1], waterHydrogens)
        water_0, water_1 = AtomNode(
            nbPair[0], hydrogens_0, ResType.WATER), AtomNode(nbPair[1], hydrogens_1, ResType.WATER)
        contactingWaters.append((water_0, water_1))
    LOGGER.info(f'Water contacts ({len(contactingWaters)}) saved.')

    proteinHydrogens = atoms.select('protein and hydrogen')
    proteinHydrophilic = atoms.select(
        f'protein and ({" or ".join(donors + acceptors)})')
    neighborWaterProtein = findNeighbors(
        waterOxygens, distDA, proteinHydrophilic)
    contactingWaterProtein = []

    for nbPair in neighborWaterProtein:
        hydrogens_w, hydrogens_p = getSameResidueAtoms(
            nbPair[0], waterHydrogens), getNeighborsForAtom(nbPair[1], proteinHydrogens, 1.4)
        water, protein = AtomNode(
            nbPair[0], hydrogens_w, ResType.WATER), AtomNode(nbPair[1], hydrogens_p, ResType.PROTEIN)
        contactingWaterProtein.append((water, protein))
    LOGGER.info(
        f'Water-protein contacts ({len(contactingWaterProtein)}) saved.')

    hydrogenBonds = []

    for nbPair in contactingWaters:
        if checkHBondAngle(nbPair[0], nbPair[1], angleWW):
            newHBond = HydrogenBond(nbPair[0], nbPair[1])
            hydrogenBonds.append(newHBond)

        if checkHBondAngle(nbPair[1], nbPair[0], angleWW):
            newHBond = HydrogenBond(nbPair[1], nbPair[0])
            hydrogenBonds.append(newHBond)
    LOGGER.info(f'Hydrogen bonds ({len(hydrogenBonds)}) calculated.')

    for nbPair in contactingWaterProtein:
        if nbPair[1].atom.getName()[0] in ['N', 'S', 'F', 'O'] and checkHBondAngle(nbPair[0], nbPair[1], anglePDWA):
            newHBond = HydrogenBond(nbPair[0], nbPair[1])
            hydrogenBonds.append(newHBond)

        if nbPair[1].atom.getName()[0] in ['N', 'S', 'F', 'O'] and checkHBondAngle(nbPair[1], nbPair[0], anglePAWD):
            newHBond = HydrogenBond(nbPair[1], nbPair[0])
            hydrogenBonds.append(newHBond)
    LOGGER.info(f'Hydrogen bonds ({len(hydrogenBonds)}) calculated.')

    relations = getRelationsList(hydrogenBonds, len(atoms))
    waterBridges = []

    if method == 'chain':
        waterBridges = calcBridges(relations, proteinHydrophilic, waterDepth)
        waterBridges = getUniqueElements(waterBridges, getChainBridgeTuple)
    elif method == 'cluster':
        waterBridges = calcBridges(relations, proteinHydrophilic)

    LOGGER.info(f'Water bridges calculated.')
    return waterBridges
