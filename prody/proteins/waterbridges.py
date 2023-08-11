# -*- coding: utf-8 -*-

"""This module detects, predicts and analyzes water bridges.
"""

__author__ = 'Karolina Mikulska-Ruminska'
__credits__ = ['Frane Doljanin', 'Karolina Mikulska-Ruminska']
__email__ = ['karolamik@fizyka.umk.pl', 'fdoljanin@pmfst.hr']

import numpy as np

from itertools import combinations
from collections import deque
from enum import Enum, auto
from copy import copy

from prody import LOGGER
from prody.atomic import Atom, Atomic, AtomGroup
from prody.dynamics.plotting import showAtomicMatrix
from prody.ensemble import Ensemble
from prody.measure import calcAngle, calcDistance
from prody.measure.contacts import findNeighbors
from prody.proteins import writePDB


__all__ = ['calcWaterBridges', 'calcWaterBridgesTrajectory', 'getWaterBridgesInfoOutput',
           'calcWaterBridgesStatistics', 'getWaterBridgeStatInfo', 'calcWaterBridgeMatrix', 'showWaterBridgeMatrix',
           'calcBridgingResiduesHistogram', 'calcWaterBridgesDistribution',
           'savePDBWaterBridges', 'savePDBWaterBridgesTrajectory',
           'saveWaterBridges', 'parseWaterBridges']


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

    def resetVisits(self, atomsToReset=[]):
        if atomsToReset:
            for atom in atomsToReset:
                self[atom].isVisited = False

            return

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
            self[acceptor].bondedAtoms.append(donor)

    def removeUnbonded(self):
        for i, node in enumerate(self.nodes):
            if node and not node.bonds:
                self[i] = None

    def __getitem__(self, key):
        if isinstance(key, Atom):
            key = key.getIndex()

        return self.nodes[key]

    def __setitem__(self, key, value):
        self.nodes[key] = value


class HBondConstraints:
    def __init__(self, donors, acceptors, angleWW=None, anglePAWD=None, anglePDWA=None):
        self.donors = donors
        self.acceptors = acceptors
        self.angleWW = angleWW
        self.anglePAWD = anglePAWD
        self.anglePDWA = anglePDWA

    def addAngles(self, angleWW, anglePAWD, anglePDWA):
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

        anglesForComb = {
            (ResType.WATER, ResType.WATER): constraints.angleWW,
            (ResType.WATER, ResType.PROTEIN): constraints.anglePAWD,
            (ResType.PROTEIN, ResType.WATER): constraints.anglePDWA,
        }

        angleRange = anglesForComb[(donor.type, acceptor.type)]
        if not angleRange:
            return True

        return HydrogenBond.checkIsHBondAngle(donor, acceptor, angleRange)

    @staticmethod
    def checkIsHBondAngle(donor, acceptor, angleRange):
        for hydrogen in donor.hydrogens:
            bondAngle = calcAngle(donor.atom, hydrogen, acceptor.atom)
            if angleRange[0] < bondAngle < angleRange[1]:
                return True

        return False


def calcBridges(relations, hydrophilicList, method, maxDepth, maxNumResidues):
    waterBridges = []

    for atom in hydrophilicList:
        observedNode = relations[atom]
        if not observedNode:
            continue

        newBridges = []
        if method == 'chain':
            newBridges = getBridgeChain_BFS(
                observedNode, relations, maxDepth+1)
            relations.resetVisits()
        else:
            newBridges = getBridgeForResidue_BFS(
                observedNode, relations, maxDepth, maxNumResidues)

        waterBridges += newBridges

    return waterBridges


def getChainBridgeTuple(bridge):
    i_1, i_2 = bridge[0].getIndex(), bridge[-1].getIndex()
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
    queue = deque([(observedNode.atom, [])])
    waterBridges = []

    while queue:
        currentAtom, currentPath = queue.popleft()
        currentNode = relationsList[currentAtom]

        if currentNode.isVisited:
            continue

        currentNode.isVisited = True

        if currentNode.type == ResType.PROTEIN and len(currentPath):
            waterBridges.append(currentPath + [currentAtom])
            continue

        if len(currentPath) == maxDepth:
            continue

        for atom in currentNode.bondedAtoms:
            queue.append(
                (atom, currentPath + [currentAtom]))

    return waterBridges


def getBridgeForResidue_BFS(observedNode, relationsList, maxDepth, maxNumResidues):
    waterBridges = []

    for waterAtom in observedNode.bondedAtoms:
        bridgeAtoms = []
        newWaters = []
        queue = deque([(waterAtom, 0)])
        while queue:
            currentAtom, depth = queue.popleft()
            currentNode = relationsList[currentAtom]

            if currentNode.isVisited or (maxDepth and depth > maxDepth):
                continue

            bridgeAtoms.append(currentAtom)

            if currentNode.type == ResType.PROTEIN:
                continue
            if currentNode.type == ResType.WATER:
                currentNode.isVisited = True
                newWaters.append(currentNode.atom)

            for atom in currentNode.bondedAtoms:
                queue.append((atom, depth+1))

        bridgeAtomIndices = list(
            set(map(lambda a: a.getIndex(), bridgeAtoms)))

        numProteinResidues = sum(
            relationsList[atomIndex].type == ResType.PROTEIN for atomIndex in bridgeAtomIndices)
        numResidues = len(bridgeAtomIndices)

        if (not maxNumResidues or numResidues <= maxNumResidues) and numProteinResidues >= 2:
            waterBridges.append(bridgeAtomIndices)
        else:
            relationsList.resetVisits(newWaters)

    return waterBridges


class AtomicOutput:
    def __init__(self, proteins, waters):
        self.proteins = proteins
        self.waters = waters


def getInfoOutput(waterBridgesAtomic):
    output = []
    for bridge in waterBridgesAtomic:
        bridgeOutput = []

        for atom in bridge.proteins:
            residueInfo = f"{atom.getResname()}{atom.getResnum()}"
            atomInfo = f"{atom.getName()}_{atom.getIndex()}"
            chainInfo = atom.getChid()
            bridgeOutput += [residueInfo, atomInfo, chainInfo]

        for atom_1, atom_2 in combinations(bridge.proteins, r=2):
            bridgeOutput += [calcDistance(atom_1, atom_2)]

        bridgeOutput += [len(bridge.waters)]
        bridgeOutput += [
            list(map(lambda w: f"{w.getChid()}_{w.getIndex()}", bridge.waters))]

        output.append(bridgeOutput)

    return output


def getWaterBridgesInfoOutput(waterBridgesAtomic):
    """Converts single frame/trajectory atomic output from calcWaterBridges/Trajectory to info output.

    :arg waterBridgesAtomic: water bridges from calcWaterBridges/Trajectory
    :type waterBridgesAtomic: list
    """
    isSingleFrame = isinstance(waterBridgesAtomic[0], AtomicOutput)
    if isSingleFrame:
        return getInfoOutput(waterBridgesAtomic)

    output = []
    for frame in waterBridgesAtomic:
        currentFrameInfo = getInfoOutput(frame)
        output.append(currentFrameInfo)

    return output


def getAtomicOutput(waterBridges, relations):
    output = []
    for bridge in waterBridges:
        proteinAtoms, waterAtoms = [], []
        for atomIndex in bridge:
            if relations[atomIndex].type == ResType.PROTEIN:
                proteinAtoms.append(relations[atomIndex].atom)
            else:
                waterAtoms.append(relations[atomIndex].atom)

        output.append(AtomicOutput(proteinAtoms, waterAtoms))

    return output


def getElementsRegex(elements):
    return f'[{"|".join(elements)}].*'


def calcWaterBridges(atoms, **kwargs):
    """Compute water bridges for a protein that has water molecules.

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`

    :arg method: cluster or chain, where chain find shortest water bridging path between two protein atoms
        default is 'chain'
    :type method: string 'cluster' | 'chain'

    :arg distDA: maximal distance between water/protein donor and acceptor
        default is 3.5
    :type distDA: int, float

    :arg distWR: maximal distance between considered water and any residue
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

    :arg maxDepth: maximum number of waters in chain/depth of residues in cluster
        default is 2
    :type maxDepth: int, None

    :arg maxNumRes: maximum number of water+protein residues in cluster
        default is None
    :type maxNumRes: int, None

    :arg donors: which atoms to count as donors 
        default is ['N', 'O', 'S', 'F']
    :type donors: list

    :arg acceptors: which atoms to count as acceptors 
        default is ['N', 'O', 'S', 'F']
    :type acceptors: list

    :arg output: return information arrays, (protein atoms, water atoms), or just atom indices per bridge
        default is 'atomic'
    :type output: 'info' | 'atomic' | 'indices'

    :arg isInfoLog: should log information
        default is True
    :type output: bool
    """

    method = kwargs.pop('method', 'chain')
    distDA = kwargs.pop('distDA', 3.5)
    distWR = kwargs.pop('distWR', 4)
    anglePDWA = kwargs.pop('anglePDWA', (100, 200))
    anglePAWD = kwargs.pop('anglePAWD', (100, 140))
    angleWW = kwargs.pop('angleWW', (140, 180))
    maxDepth = kwargs.pop('maxDepth', 2)
    maxNumResidues = kwargs.pop('maxNumRes', None)
    donors = kwargs.pop('donors', ['N', 'O', 'S', 'F'])
    acceptors = kwargs.pop('acceptors', ['N', 'O', 'S', 'F'])
    outputType = kwargs.pop('output', 'atomic')
    isInfoLog = kwargs.pop('isInfoLog', True)
    DIST_COVALENT_H = 1.4

    if method not in ['chain', 'cluster']:
        raise TypeError('Method should be chain or cluster.')
    if outputType not in ['info', 'atomic', 'indices']:
        raise TypeError('Output can be info, atomic or indices.')

    relations = RelationList(len(atoms))
    consideredAtoms = ~atoms.select(
        f'water and not within {distWR} of protein')

    waterHydrogens = consideredAtoms.select('water and hydrogen') or []
    waterOxygens = consideredAtoms.select('water and oxygen')
    waterHydroOxyPairs = findNeighbors(
        waterOxygens, DIST_COVALENT_H, waterHydrogens)
    for oxygen in waterOxygens:
        relations.addNode(oxygen, ResType.WATER)
    for pair in waterHydroOxyPairs:
        oxygen, hydrogen, _ = pair
        relations[oxygen].hydrogens.append(hydrogen)

    proteinHydrophilic = consideredAtoms.select(
        f'protein and name "{getElementsRegex(set(donors+acceptors))}" and within {distWR} of water')

    proteinHydrogens = consideredAtoms.select(f'protein and hydrogen') or []
    proteinHydroPairs = findNeighbors(
        proteinHydrophilic, DIST_COVALENT_H, proteinHydrogens)
    for hydrophilic in proteinHydrophilic:
        relations.addNode(hydrophilic, ResType.PROTEIN)
    for pair in proteinHydroPairs:
        hydrophilic, hydrogen, _ = pair
        relations[hydrophilic].hydrogens.append(hydrogen)

    contactingWaters = findNeighbors(waterOxygens, distDA)
    contactingWaterProtein = findNeighbors(
        waterOxygens, distDA, proteinHydrophilic)

    contactingWaterNodes = list(
        map(lambda ww: (relations[ww[0]], relations[ww[1]]), contactingWaters))
    contactingWaterProteinNodes = list(
        map(lambda wp: (relations[wp[0]], relations[wp[1]]), contactingWaterProtein))

    constraints = HBondConstraints(acceptors, donors)
    if len(waterHydrogens) + len(proteinHydrogens):
        constraints.addAngles(angleWW, anglePAWD, anglePDWA)
    else:
        LOGGER.info('No hydrogens detected, angle criteria will not be used.')

    for pair in contactingWaterNodes + contactingWaterProteinNodes:
        for a, b in [(0, 1), (1, 0)]:
            if HydrogenBond.checkIsHBond(pair[a], pair[b], constraints):
                newHBond = HydrogenBond(pair[a].atom, pair[b].atom)
                relations.addHBond(newHBond)

    relations.removeUnbonded()

    waterBridgesWithIndices = calcBridges(
        relations, proteinHydrophilic, method, maxDepth, maxNumResidues)
    if method == 'chain':
        waterBridgesWithIndices = getUniqueElements(
            waterBridgesWithIndices, getChainBridgeTuple)

    LOGGER.info(
        f'{len(waterBridgesWithIndices)} water bridges detected.')
    if method == 'atomic':
        LOGGER.info('Call getInfoOutput to convert atomic to info output.')

    atomicOutput = getAtomicOutput(waterBridgesWithIndices, relations)
    infoOutput = getInfoOutput(atomicOutput)

    if isInfoLog:
        for bridge in infoOutput:
            LOGGER.info(' '.join(map(str, bridge)))

    if outputType == 'info':
        return infoOutput
    if outputType == 'atomic':
        return atomicOutput

    return waterBridgesWithIndices


# took from interactions.py
def calcWaterBridgesTrajectory(atoms, trajectory, **kwargs):
    """Computes water bridges for a given trajectory. Kwargs for options are the same as in calcWaterBridges.

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`

    :arg trajectory: trajectory object, DCD or multimodal PDB

    :arg start_frame: frame to start from
    :type start_frame: int

    :arg stop_frame: frame to stop
    :type stop_frame: int
    """

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
            interactions = calcWaterBridges(
                atoms_copy, isInfoLog=False, **kwargs)
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


def getResidueName(atom):
    return f'{atom.getResname()}{atom.getResnum()}{atom.getChid()}'


class DictionaryList:
    def __init__(self, initialValue):
        self.values = {}
        self.initialValue = initialValue

    def __getitem__(self, key):
        if key in self.values:
            return self.values[key]

        return copy(self.initialValue)

    def __setitem__(self, key, value):
        self.values[key] = value

    def keys(self):
        return self.values.keys()


def getResInfo(atoms):
    dict = {}
    nums = atoms.select('protein').getResnums()
    names = atoms.select('protein').getResnames()
    chids = atoms.select('protein').getChids()

    for i, num in enumerate(nums):
        dict[num] = f"{names[i]}{num}{chids[i]}"

    return dict


def getWaterBridgeStatInfo(stats, atoms):
    """Converts calcWaterBridgesStatistic indices output to info output from stat.

    :arg stats: statistics returned by calcWaterBridgesStatistics, output='indices'
    :type stats: dictionary

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`
    """
    residueInfo = getResInfo(atoms)
    infoOutput = {}
    for key, value in stats.items():
        x_id, y_id = key
        x_info, y_info = residueInfo[x_id], residueInfo[y_id]
        newKey = (x_info, y_info)

        infoOutput[newKey] = value

    return infoOutput


def calcWaterBridgesStatistics(frames, trajectory, **kwargs):
    """Returns statistics. Value is percentage of bridge appearance of frames for each residue.

    :arg frames: list of water bridges from calcWaterBridgesTrajectory(), output='atomic'
    :type frames: list

    :arg output: return dictorinary whose keys are tuples of resnames or resids
        default is 'resid'
    :type output: 'info' | 'indices'

    :arg filename: name of file to save statistic information if wanted
        default is None
    :type filename: string
    """
    output = kwargs.pop('output', 'indices')
    filename = kwargs.pop('filename', None)
    if output not in ['info', 'indices']:
        raise TypeError('Output should be info or indices!')

    allCoordinates = trajectory.getCoordsets()
    interactionCount = DictionaryList(0)
    distances = DictionaryList([])
    resNames = {}

    for frameIndex, frame in enumerate(frames):
        frameCombinations = []

        for bridge in frame:
            proteinAtoms = bridge.proteins
            for atom_1, atom_2 in combinations(proteinAtoms, r=2):
                ind_1, ind_2 = atom_1.getIndex(), atom_2.getIndex()
                res_1, res_2 = atom_1.getResnum(), atom_2.getResnum()

                if res_1 == res_2:
                    continue

                if (res_1, res_2) not in frameCombinations:
                    interactionCount[(res_1, res_2)] += 1
                    interactionCount[(res_2, res_1)] += 1
                    frameCombinations += [(res_1, res_2), (res_2, res_1)]

                coords = allCoordinates[frameIndex]
                atom_1_coords, atom_2_coords = coords[ind_1], coords[ind_2]
                distance = calcDistance(atom_1_coords, atom_2_coords)
                distances[(res_1, res_2)] += [distance]
                distances[(res_2, res_1)] += [distance]

                res_1_name, res_2_name = getResidueName(
                    atom_1),  getResidueName(atom_2)
                resNames[res_1] = res_1_name
                resNames[res_2] = res_2_name

    tableHeader = f'{"RES1":<15}{"RES2":<15}{"PERC":<10}{"DIST_AVG":<10}{"DIST_STD":<10}'
    LOGGER.info(tableHeader)
    info = {}
    file = open(filename, 'w') if filename else None
    if file:
        file.write(tableHeader + '\n')

    for key in interactionCount.keys():
        percentage = 100 * interactionCount[key]/len(frames)
        distAvg = np.average(distances[key])
        distStd = np.std(distances[key])
        pairInfo = {"percentage": percentage,
                    "distAvg": distAvg, "distStd": distStd}

        outputKey = key
        x, y = key
        if output == 'resname':
            outputKey = (resNames[x], resNames[y])

        info[outputKey] = pairInfo

        tableRow = f'{resNames[x]:<15}{resNames[y]:<15}{percentage:<10.3f}{distAvg:<10.3f}{distStd:<10.3f}'
        LOGGER.info(tableRow)
        if file:
            file.write(tableRow + '\n')

    if file:
        file.close()

    return info


def calcWaterBridgeMatrix(data, metric):
    """Returns matrix which has metric as value and residue ids as ax indices.

    :arg data: dictionary returned by calcWaterBridgesStatistics, output='indices' 
    :type data: dict

    :arg metric: dict key from data
    :type metric: 'percentage' | 'distAvg' | 'distStd'
    """
    maxResnum = max(max(key) for key in data.keys()) + 1
    resMatrix = np.zeros((maxResnum, maxResnum), dtype=float)

    for key, value in data.items():
        resMatrix[key] = value[metric]

    return resMatrix


def showWaterBridgeMatrix(data, metric):
    """Shows matrix which has percentage/avg distance as value and residue ids as ax indices.

    :arg data: dictionary returned by calcWaterBridgesStatistics, output='indices' 
    :type data: dict

    :arg metric: dict key from data
    :type metric: 'percentage' | 'distAvg' | 'distStd'
    """
    import matplotlib.pyplot as plt
    matrix = calcWaterBridgeMatrix(data, metric)
    titles = {
        'percentage': 'Interaction percentage',
        'distAvg': 'Average distance',
        'distStd': 'Distance standard deviation'
    }

    showAtomicMatrix(matrix)
    plt.title(titles[metric])


def reduceTo1D(list, elementSel=lambda x: x, sublistSel=lambda x: x):
    return [elementSel(element) for sublist in list for element in sublistSel(sublist)]


def mofifyBeta(bridgeFrames, atoms):
    residueOccurances = {}

    for frame in bridgeFrames:
        frameResnums = set(reduceTo1D(
            frame, lambda a: a.getResnum(), lambda b: b.proteins))
        for res in frameResnums:
            residueOccurances[res] = residueOccurances.get(res, 0) + 1

    atoms.setBetas(0)
    for resnum, value in residueOccurances.items():
        residueAtoms = atoms.select(
            f'resnum {resnum}')
        beta = value/len(bridgeFrames)

        residueAtoms.setBetas(beta)


def calcBridgingResiduesHistogram(frames, **kwargs):
    """Calculates, plots and returns number of frames that each residue is involved in making water bridges, sorted by value. 

    :arg frames: list of water bridges from calcWaterBridgesTrajectory(), output='atomic'
    :type frames: list

    :arg clip: maximal number of residues on graph; to represent all set None
        default is 20
    :type clip: int
    """
    import matplotlib.pyplot as plt

    clip = kwargs.pop('clip', 20)
    if clip == None:
        clip = 0

    residuesWithCount = {}
    for frame in frames:
        frameResidues = set(reduceTo1D(
            frame, getResidueName, lambda wb: wb.proteins))

        for res in frameResidues:
            residuesWithCount[res] = residuesWithCount.get(res, 0) + 1

    sortedResidues = sorted(residuesWithCount.items(),
                            key=lambda r: r[1])

    labels, values = zip(*sortedResidues[-clip:])

    plt.figure(figsize=(5, 3 + 0.11 * len(labels)))
    plt.barh(labels, values)
    plt.xlabel('Number of frame appearances')
    plt.ylabel('Residue')
    plt.title('Water bridging residues')
    plt.tight_layout()
    plt.margins(y=0.01)
    plt.gca().xaxis.set_label_position('top')
    plt.gca().xaxis.tick_top()
    plt.show()

    return sortedResidues


def getBridgingResidues(frames, residue):
    residuesWithCount = {}

    for frame in frames:
        frameResidues = []
        for bridge in frame:
            resNames = list(map(getResidueName, bridge.proteins))
            if residue in resNames:
                frameResidues += resNames

        if not frameResidues:
            continue

        frameResidues = set(frameResidues)
        frameResidues.remove(residue)
        for res in frameResidues:
            residuesWithCount[res] = residuesWithCount.get(res, 0) + 1

    sortedResidues = sorted(residuesWithCount.items(),
                            key=lambda r: r[1], reverse=True)
    return sortedResidues


def getWaterCountDistribution(frames, res_a, res_b):
    waters = []

    for frame in frames:
        for bridge in frame:
            resNames = list(map(getResidueName, bridge.proteins))
            if not (res_a in resNames and res_b in resNames):
                continue

            waterInvolvedCount = len(bridge.waters)
            waters.append(waterInvolvedCount)

    return waters


def getDistanceDistribution(frames, res_a, res_b, trajectory):
    distances = []
    allCoordinates = trajectory.getCoordsets()

    for frameIndex, frame in enumerate(frames):
        for bridge in frame:
            resNames = list(map(getResidueName, bridge.proteins))
            if not (res_a in resNames and res_b in resNames):
                continue

            for atom_1, atom_2 in combinations(bridge.proteins, r=2):
                coords = allCoordinates[frameIndex]
                atom_1_coords, atom_2_coords = coords[atom_1.getIndex(
                )], coords[atom_2.getIndex()]

                dist = calcDistance(atom_1_coords, atom_2_coords)
                distances.append(dist)

    return distances


def getResidueLocationDistrubtion(frames, res_a, res_b):
    locationInfo = {"backbone": 0, "side": 0}
    result = {res_a: locationInfo.copy(), res_b: locationInfo.copy()}

    for frame in frames:
        for bridge in frame:
            atoms_a = filter(lambda a: getResidueName(a)
                             == res_a, bridge.proteins)
            atoms_b = filter(lambda a: getResidueName(a)
                             == res_b, bridge.proteins)

            if not atoms_a or not atoms_b:
                continue

            def atomType(a): return 'backbone' if a.getName() in [
                'CA', 'C', 'N', 'O'] else 'side'
            for atom in atoms_a:
                result[res_a][atomType(atom)] += 1
            for atom in atoms_b:
                result[res_b][atomType(atom)] += 1

    return result


def calcWaterBridgesDistribution(frames, res_a, res_b=None, **kwargs):
    """Returns distribution for certain metric and plots if possible.

    :arg res_a: name of first residue
    :type frames: str

    :arg res_b: name of second residue
        default is None
    :type frames: str

    :arg metric: 'residues' returns names and frame count of residues interacting with res_a,
                'waters' returns water count for each bridge between res_a and res_b
                'distance' returns distance between each pair of protein atoms involved in bridge between res_a and res_b
                'location' returns dictionary with backbone/sidechain count information
    :type metric: 'residues' | 'waters' | 'distance' | 'location'

    :trajectory: DCD file - necessary for distance distribution

    :arg output: return 2D matrices or dictionary where key is residue info
        default is 'dict'
    :type output: 'dict' | 'indices'
    """
    import matplotlib.pyplot as plt
    
    metric = kwargs.pop('metric', 'residues')
    trajectory = kwargs.pop('trajectory', None)

    if metric == 'distance' and not trajectory:
        raise TypeError(
            'Distance distribution measurement needs trajectory argument!')

    methods = {
        'residues': lambda: getBridgingResidues(frames, res_a),
        'waters': lambda: getWaterCountDistribution(frames, res_a, res_b),
        'distance': lambda: getDistanceDistribution(frames, res_a, res_b, trajectory),
        'location': lambda: getResidueLocationDistrubtion(frames, res_a, res_b)
    }

    result = methods[metric]()

    if metric in ['waters', 'distance']:
        plt.hist(result, rwidth=0.95, density=True)
        plt.xlabel('Value')
        plt.ylabel('Probability')
        plt.title(f'Distribution: {metric}')
        plt.show()

    return methods[metric]()


def savePDBWaterBridges(bridges, atoms, filename):
    """Saves single PDB with occupancy on protein atoms and waters involved bridges.

    :arg bridges: atomic output from calcWaterBridges
    :type bridges: list

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`

    :arg filename: name of file to be saved
    :type filename: string
    """
    atoms = atoms.copy()

    mofifyBeta([bridges], atoms)
    atoms.setOccupancies(0)
    atoms.select('beta = 1').setOccupancies(1)

    proteinAtoms = atoms.select('protein')
    waterOxygens = reduceTo1D(
        bridges, lambda w: w.getIndex(), lambda b: b.waters)
    waterAtoms = atoms.select(
        f'same residue as water within 1.6 of index {" ".join(map(str, waterOxygens))}')

    atomsToSave = proteinAtoms.toAtomGroup() + waterAtoms.toAtomGroup()
    return writePDB(filename, atomsToSave)


def savePDBWaterBridgesTrajectory(bridgeFrames, atoms, filename, trajectory=None):
    """Saves one PDB per frame with occupancy and beta on protein atoms and waters forming bridges in frame.

    :arg bridgeFrames: atomic output from calcWaterBridgesTrajectory
    :type bridgeFrames: list

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`

    :arg filename: name of file to be saved; must end in .pdb
    :type filename: string

    :arg trajectory: DCD trajectory (not needed for multimodal PDB)
    """
    if not trajectory and atoms.numCoordsets() < len(bridgeFrames):
        raise TypeError('Provide parsed trajectory!')

    filename = filename[:filename.rfind('.')]

    atoms = atoms.copy()
    mofifyBeta(bridgeFrames, atoms)

    for frameIndex, frame in enumerate(bridgeFrames):
        if trajectory:
            coords = trajectory[frameIndex].getCoords()
            atoms.setCoords(coords)
        else:
            atoms.setACSIndex(frameIndex)

        waterAtoms = reduceTo1D(frame, sublistSel=lambda b: b.waters)
        waterResidues = atoms.select(
            f'same residue as water within 1.6 of index {" ".join(map(lambda a: str(a.getIndex()), waterAtoms))}')

        bridgeProteinAtoms = reduceTo1D(
            frame, lambda p: p.getResnum(), lambda b: b.proteins)
        atoms.setOccupancies(0)
        atoms.select(
            f'resid {" ".join(map(str, bridgeProteinAtoms))}').setOccupancies(1)

        atomsToSave = atoms.select(
            'protein').toAtomGroup() + waterResidues.toAtomGroup()

        if trajectory:
            writePDB(f'{filename}_{frameIndex}.pdb', atomsToSave)
        else:
            writePDB(f'{filename}_{frameIndex}.pdb',
                     atomsToSave, csets=frameIndex)


def getBridgeIndicesString(bridge):
    return ' '.join(map(lambda a: str(a.getIndex()), bridge.proteins))\
        + '|'\
        + ' '.join(map(lambda a: str(a.getIndex()), bridge.waters))


def saveWaterBridges(atomicBridges, filename):
    """Save water bridges as information (.txt) or WaterBridges (.wb) parsable file.

    :arg atomicBridges: atomic output from calcWaterBridges/Trajectory
    :type atomicBridges: list

    :arg filename: path where file should be saved
    :type filename: string
    """
    isSingleBridge = isinstance(atomicBridges[0], AtomicOutput)
    isInfoOutput = filename.split('.')[-1] == 'txt'
    file = open(filename, 'w')

    if isSingleBridge:
        atomicBridges = [atomicBridges]

    if isInfoOutput:
        info = getWaterBridgesInfoOutput(atomicBridges)
        for frameIndex, frame in enumerate(info):
            file.write(f'FRAME {frameIndex}\n')
            for bridge in frame:
                file.write(' '.join(map(str, bridge)) + '\n')

        return file.close()

    for frame in atomicBridges:
        for bridge in frame:
            line = getBridgeIndicesString(bridge)
            file.write(line + '\n')
        file.write('ENDFRAME\n')

    return file.close()


def getBridgeFromIndices(proteinIndices, waterIndices, atoms):
    proteins = list(map(lambda i: atoms[i], proteinIndices))
    waters = list(map(lambda i: atoms[i], waterIndices))
    atomicBridge = AtomicOutput(proteins, waters)

    return atomicBridge


def parseWaterBridges(filename, atoms):
    """Parse water bridges from .wb file saved by saveWaterBridges, returns atomic type.

    :arg filename: path of file where bridges are stored
    :type filename: string

    :arg atoms: Atomic object on which calcWaterBridges was performed
    :type atoms: :class:`.Atomic`
    """
    bridgesFrames = []
    currentBridges = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line == 'ENDFRAME':
                bridgesFrames.append(currentBridges)
                currentBridges = []
                continue

            proteinString, waterString = line.split('|')
            proteinIndices = map(int, proteinString.split())
            waterIndices = map(int, waterString.split())

            bridge = getBridgeFromIndices(proteinIndices, waterIndices, atoms)
            currentBridges.append(bridge)

    if len(bridgesFrames) == 1:
        bridgesFrames = bridgesFrames[0]

    return bridgesFrames
