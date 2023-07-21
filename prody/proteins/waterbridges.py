# -*- coding: utf-8 -*-

"""This module detects, predicts and analyzes water bridges.
"""

__author__ = 'Karolina Mikulska-Ruminska'
__credits__ = ['Frane Doljanin', 'Karolina Mikulska-Ruminska']
__email__ = ['karolamik@fizyka.umk.pl', 'fdoljanin@pmfst.hr']


from collections import deque
from enum import Enum, auto

from prody import LOGGER, Atom, Atomic
from prody.measure import calcAngle
from prody.measure.contacts import findNeighbors
from prody.ensemble import Ensemble

from timeit import default_timer as timer


__all__ = ['calcWaterBridges', 'calcWaterBridgesTrajectory']


class ResType(Enum):
    WATER = auto()
    PROTEIN = auto()


class AtomNode:
    def __init__(self, atom, type):
        self.atom = atom
        self.hydrogens = []
        self.type = type
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

    def addNode(self, atom, type):
        self[atom.getIndex()] = AtomNode(atom, type)

        return self[atom]

    def addHBond(self, bond):
        donor, acceptor = bond.donor, bond.acceptor
        self[donor].bonds.append(bond)
        self[acceptor].bonds.append(bond)

        if acceptor.getIndex() not in map(lambda a: a.getIndex(), self[donor].bondedAtoms):
            self[donor].bondedAtoms.append(acceptor)
        if donor.getIndex() not in map(lambda a: a.getIndex(), self[acceptor].bondedAtoms):
            self[acceptor].bondedAtoms.append(donor)

    def removeUnnecessary(self):
        for i in range(len(self.nodes)):
            if self.nodes[i] and not self.nodes[i].bonds:
                self[i] = None

    def __getitem__(self, key):
        if isinstance(key, Atom):
            key = key.getIndex()

        return self.nodes[key]

    def __setitem__(self, key, value):
        self.nodes[key] = value


class HBondConstraints:
    def __init__(self, donors, acceptors, angleWW, anglePAWD, anglePDWA):
        self.donors = donors
        self.acceptors = acceptors
        self.angleWW = angleWW
        self.anglePAWD = anglePAWD
        self.anglePDWA = anglePDWA


class HydrogenBond:
    def __init__(self, donor, acceptor):
        self.donor = donor
        self.acceptor = acceptor

    @staticmethod
    def checkIsHBond(donor, acceptor, constraints):
        if donor.type != ResType.WATER and donor.atom.getName()[0] not in constraints.donors:
            return False
        if acceptor.type != ResType.WATER and acceptor.atom.getName()[0] not in constraints.acceptors:
            return False

        angleRange = ()
        if donor.type == ResType.WATER and acceptor.type == ResType.WATER:
            angleRange = constraints.angleWW
        elif donor.type == ResType.WATER and acceptor.type == ResType.PROTEIN:
            angleRange = constraints.anglePAWD
        else:
            angleRange = constraints.anglePDWA

        return HydrogenBond.checkIsHBondAngle(donor, acceptor, angleRange)

    @staticmethod
    def checkIsHBondAngle(donor, acceptor, angleRange):
        for hydrogen in donor.hydrogens:
            bondAngle = calcAngle(donor.atom, hydrogen, acceptor.atom)
            if angleRange[0] < bondAngle < angleRange[1]:
                return True

        return False


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
        bridgeAtoms = []
        queue = deque([waterAtom])
        while queue:
            currentAtom = queue.popleft()
            currentNode = relationsList[currentAtom]

            if currentNode.isVisited:
                continue

            bridgeAtoms.append(currentAtom)

            if currentNode.type == ResType.PROTEIN:
                continue
            if currentNode.type == ResType.WATER:
                currentNode.isVisited = True

            for atom in currentNode.bondedAtoms:
                queue.append(atom)

        bridgeAtomIndices = list(
            set(map(lambda a: a.getIndex(), bridgeAtoms)))
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
    waterDepth = kwargs.pop('waterDepth', 2)
    donors = kwargs.pop('donors', ['nitrogen', 'oxygen', 'sulfur'])
    acceptors = kwargs.pop('acceptors', ['nitrogen', 'oxygen', 'sulfur'])

    constraints = HBondConstraints(
        ['N', 'O', 'S', 'F'], ['N', 'O', 'S', 'F'], angleWW, anglePAWD, anglePDWA)

    if method not in ['chain', 'cluster']:
        raise TypeError('Method should be chain or cluster.')

    relations = RelationList(len(atoms))
    consideredAtoms = ~atoms.select(
        f'water and not within {distWR} of protein')

    waterHydrogens = consideredAtoms.select('water and hydrogen')
    waterOxygens = consideredAtoms.select('water and oxygen')
    waterHydroOxyPairs = findNeighbors(waterOxygens, 1.4, waterHydrogens)
    for oxygen in waterOxygens:
        relations.addNode(oxygen, ResType.WATER)
    for pair in waterHydroOxyPairs:
        oxygen, hydrogen = pair[0], pair[1]
        relations[oxygen].hydrogens.append(hydrogen)

    proteinHydrophilic = consideredAtoms.select(
        f'protein and ({" or ".join(donors + acceptors)}) and within {distWR} of water')
    proteinHydrogens = consideredAtoms.select(f'protein and hydrogen')
    proteinHydroPairs = findNeighbors(
        proteinHydrophilic, 1.4, proteinHydrogens)
    for hydrophilic in proteinHydrophilic:
        relations.addNode(hydrophilic, ResType.PROTEIN)
    for pair in proteinHydroPairs:
        hydrophilic, hydrogen = pair[0], pair[1]
        relations[hydrophilic].hydrogens.append(hydrogen)

    contactingWaters = findNeighbors(waterOxygens, distDA)
    contactingWaterProtein = findNeighbors(
        waterOxygens, distDA, proteinHydrophilic)

    contactingWaters = list(
        map(lambda ww: (relations[ww[0]], relations[ww[1]]), contactingWaters))
    contactingWaterProtein = list(
        map(lambda wp: (relations[wp[0]], relations[wp[1]]), contactingWaterProtein))

    for pair in contactingWaters + contactingWaterProtein:
        for a, b in [(0, 1), (1, 0)]:
            if HydrogenBond.checkIsHBond(pair[a], pair[b], constraints):
                newHBond = HydrogenBond(pair[a].atom, pair[b].atom)
                relations.addHBond(newHBond)
                x_cnt += 1

    relations.removeUnnecessary()
    waterBridges = []

    if method == 'chain':
        waterBridges = calcBridges(relations, proteinHydrophilic, waterDepth)
        waterBridges = getUniqueElements(waterBridges, getChainBridgeTuple)
    elif method == 'cluster':
        waterBridges = calcBridges(relations, proteinHydrophilic)

    LOGGER.info(f'{len(waterBridges)} water bridges calculated.')
    return waterBridges


# took from interactions.py
def calcWaterBridgesTrajectory(atoms, trajectory, **kwargs):
    """Compute selected type interactions for DCD trajectory or multi-model PDB 
    using default parameters."""

    interactions_all = []
    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)

    if trajectory is not None:
        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)

        # nfi = trajectory._nfi
        # trajectory.reset()
        numFrames = trajectory._n_csets

        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame+1]

        atoms_copy = atoms.copy()
        for j0, frame0 in enumerate(traj, start=start_frame):
            LOGGER.info('Frame: {0}'.format(j0))
            atoms_copy.setCoords(frame0.getCoords())
            interactions = calcWaterBridges(atoms_copy, **kwargs)
            interactions_all.append(interactions)
        # trajectory._nfi = nfi

    else:
        if atoms.numCoordsets() > 1:
            for i in range(len(atoms.getCoordsets()[start_frame:stop_frame])):
                LOGGER.info('Model: {0}'.format(i+start_frame))
                atoms.setACSIndex(i+start_frame)
                interactions = calcWaterBridges(atoms, **kwargs)
                interactions_all.append(interactions)
        else:
            LOGGER.info('Include trajectory or use multi-model PDB file.')

    return interactions_all
