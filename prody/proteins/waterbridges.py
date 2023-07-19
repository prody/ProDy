# -*- coding: utf-8 -*-

"""This module detect, predicts and analyzes water bridges.
"""

__author__ = 'Karolina Mikulska-Ruminska'
__credits__ = ['Frane Doljanin', 'Karolina Mikulska-Ruminska']
__email__ = ['karolamik@fizyka.umk.pl', 'fdoljanin@pmfst.hr']


import numpy as np

from numpy import *
from prody import LOGGER, SETTINGS, PY3K
from prody.atomic import AtomGroup, Atom, Atomic, Selection, Select
from prody.atomic import flags
from prody.utilities import importLA, checkCoords, showFigure, getCoords
from prody.measure import calcDistance, calcAngle, calcCenter
from prody.measure.contacts import findNeighbors
from prody.proteins import writePDB, parsePDB
from collections import Counter

from prody.trajectory import TrajBase, Trajectory
from prody.ensemble import Ensemble

from enum import Enum, auto

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


def checkDoesGroupFormBridge(waterGroup):
    resAtomIndices = []
    for pair in waterGroup:
        if pair[0].type == ResType.PROTEIN:
            resAtomIndices.append(pair[0].hydrophilic.getIndex())
        if pair[1].type == ResType.PROTEIN:
            resAtomIndices.append(pair[1].hydrophilic.getIndex())

    return len(set(resAtomIndices)) >= 2


def getWaterGroupAtoms(group):
    currentGroupAtoms = []
    for bond in group:
        currentGroupAtoms += [bond[0].hydrophilic,
                              bond[1].hydrophilic]

    uniqueAtomIndices = set([cga.getIndex() for cga in currentGroupAtoms])

    return uniqueAtomIndices


def calcWaterBridges(atoms, **kwargs):
    """Compute water bridges for a protein that has water molecules.

    :arg atoms: Atomic object from which atoms are considered
    :type atoms: :class:`.Atomic`

    :arg distDA: maximal distance between water/protein donor and acceptor
        default is 3.5
    :type distDA: int, float

    :arg anglePDWA: angle range where protein is donor and water is acceptor
        default is (100, 200)
    :type anglePDWA: (int, int)

    :arg anglePAWD: angle range where protein is acceptor and water is donor
        default is (100, 140)
    :type anglePDWA: (int, int)

    :arg angleWW: angle between water donor/acceptor
        default is (140, 180)
    :type angleWW: (int, int)

    :arg donors: which atoms to count as donors 
        default is ['N', 'O', 'S', 'F']
    :type donors: list

    :arg acceptors: which atoms to count as acceptors 
        default is ['N', 'O', 'S', 'F']
    :type acceptors: list
    """

    distDA = kwargs.pop('distDA', 3.5)
    anglePDWA = kwargs.pop('anglePDWA', (100, 200))
    anglePAWD = kwargs.pop('anglePAWD', (100, 140))
    angleWW = kwargs.pop('angleWW', (140, 180))
    donors = kwargs.pop('donors', ['nitrogen', 'oxygen', 'sulfur'])
    acceptors = kwargs.pop('acceptors', ['nitrogen', 'oxygen', 'sulfur'])

    waterHydrogens = atoms.select('water and hydrogen')
    waterOxygens = atoms.select('water and oxygen')
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

    waterBridges = []

    while len(hydrogenBonds):
        if not len(hydrogenBonds) % 500:
            LOGGER.info(len(hydrogenBonds))
        waterGroup = []
        firstWater = hydrogenBonds[0][0] if hydrogenBonds[0][0].type == ResType.WATER else hydrogenBonds[0][1]
        newWaters = [firstWater]

        while True:
            currentWater = newWaters[0]
            for bond in hydrogenBonds.copy():
                if currentWater.hydrophilic.getIndex() not in [bond[0].hydrophilic.getIndex(), bond[1].hydrophilic.getIndex()]:
                    continue

                hydrogenBonds.remove(bond)
                waterGroup.append(bond)
                otherAtom = bond[0] if currentWater == bond[1] else bond[1]
                if otherAtom.type == ResType.WATER:
                    newWaters.append(otherAtom)

            newWaters.pop(0)
            if not len(newWaters):
                if checkDoesGroupFormBridge(waterGroup):
                    atoms = getWaterGroupAtoms(waterGroup)
                    waterBridges.append(atoms)

                break
    LOGGER.info('Water bridges calculated.')

    return waterBridges
