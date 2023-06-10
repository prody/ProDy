# -*- coding: utf-8 -*-

"""This module defines functions for calculating different types of interactions 
in protein structure, between proteins or between protein and ligand.
The following interactions are available for protein interactions:
        (1) Hydrogen bonds
        (2) Salt Bridges
        (3) Repulsive Ionic Bonding 
        (4) Pi stacking interactions
        (5) Pi-cation interactions
        (6) Hydrophobic interactions
        (7) Disulfide Bonds

For protein-ligand interactions (3) is replaced by water bridges.
"""

__author__ = 'Karolina Mikulska-Ruminska'
__credits__ = ['James Krieger', 'Karolina Mikulska-Ruminska']
__email__ = ['karolamik@fizyka.umk.pl', 'jamesmkrieger@gmail.com']


import numpy as np
from numpy import *
from prody import LOGGER, SETTINGS
from prody.atomic import AtomGroup, Atom, Atomic, Selection, Select
from prody.atomic import flags
from prody.utilities import importLA, checkCoords, showFigure, getCoords
from prody.measure import calcDistance, calcAngle, calcCenter
from prody.measure.contacts import findNeighbors
from prody.proteins import writePDB, parsePDB
from collections import Counter

from prody.trajectory import TrajBase, Trajectory
from prody.ensemble import Ensemble

__all__ = ['calcHydrogenBonds', 'calcChHydrogenBonds', 'calcSaltBridges',
           'calcRepulsiveIonicBonding', 'calcPiStacking', 'calcPiCation',
           'calcHydrophobic', 'calcDisulfideBonds', 'calcMetalInteractions',
           'calcHydrogenBondsTrajectory', 'calcSaltBridgesTrajectory',
           'calcRepulsiveIonicBondingTrajectory', 'calcPiStackingTrajectory', 
           'calcPiCationTrajectory', 'calcHydrophobicTrajectory', 'calcDisulfideBondsTrajectory',
           'calcProteinInteractions', 'calcStatisticsInteractions', 'calcDistribution',
           'compareInteractions', 'showInteractionsGraph',
           'calcLigandInteractions', 'listLigandInteractions', 
           'showProteinInteractions_VMD', 'showLigandInteraction_VMD', 
           'calcHydrogenBondsTrajectory',
           'Interactions', 'InteractionsTrajectory']


def cleanNumbers(listContacts):
    """Provide short list with indices and value of distance."""
    
    shortList = [ [int(str(i[0]).split()[-1].strip(')')), 
                           int(str(i[1]).split()[-1].strip(')')), 
                           str(i[0]).split()[1], 
                           str(i[1]).split()[1], 
                           float(i[2])] for i in listContacts ]    
    
    return shortList


def calcPlane(atoms):
    """Function provide parameters of a plane for aromatic rings (based on 3 points).
    Used in calcPiStacking()"""
    
    coordinates = atoms.getCoords()
    p1, p2, p3 = coordinates[:3] # 3 points will be enough to obtain the plane
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    x3, y3, z3 = p3    
    vec1 = p3 - p1 # These two vectors are in the plane
    vec2 = p2 - p1
    cp = np.cross(vec1, vec2) # the cross product is a vector normal to the plane
    a, b, c = cp
    d = np.dot(cp, p3) # This evaluates a * x3 + b * y3 + c * z3 which equals d
    
    return a,b,c,d


def calcAngleBetweenPlanes(a1, b1, c1, a2, b2, c2):  
    """Find angle between two planes."""
    import math 
          
    d = ( a1 * a2 + b1 * b2 + c1 * c2 ) 
    eq1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1) 
    eq2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2) 
    d = d / (eq1 * eq2) 
    AngleBetweenPlanes = math.degrees(math.acos(d)) 
    
    return AngleBetweenPlanes
    
    
def removeDuplicates(list_of_interactions):
    ls=[]
    newList = []
    for no, i in enumerate(list_of_interactions):
       i = sorted(list(array(i).astype(str)))
       if i not in ls:
           ls.append(i)
           newList.append(list_of_interactions[no])
    return newList


def filterInteractions(list_of_interactions, atoms, **kwargs):
    """Return interactions based on selection."""
    
    if 'selection' in kwargs:
        if 'selection2' in kwargs:
            ch1 = kwargs['selection'].split()[-1] 
            ch2 = kwargs['selection2'].split()[-1] 
            final = [i for i in list_of_interactions if (i[2] == ch1 and i[5] == ch2) or (i[5] == ch1 and i[2] == ch2)]
        else:
            p = atoms.select('same residue as protein within 10 of ('+kwargs['selection']+')')
            x = p.select(kwargs['selection']).getResnames()
            y = p.select(kwargs['selection']).getResnums()
            listOfselection = np.unique(list(map(lambda x, y: x + str(y), x, y)))
            final = [i for i in list_of_interactions if i[0] in listOfselection or i[3] in listOfselection]
    else:
        final = list_of_interactions
    return final


def calcHydrogenBonds(atoms, **kwargs):
    """Compute hydrogen bonds for proteins and other molecules.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between donor and acceptor.
        default is 3.5
    :type distA: int, float
    
    :arg angle: non-zero value, maximal (180 - D-H-A angle) (donor, hydrogen, acceptor).
        default is 40.
    :type angle: int, float
    
    :arg seq_cutoff: non-zero value, interactions will be found between atoms with index differences
        that are higher than seq_cutoff.
        default is 25 atoms.
    :type seq_cutoff: int

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'

    Structure should contain hydrogens.
    If not they can be added using addMissingAtoms(pdb_name) function available in ProDy after Openbabel installation.
    `conda install -c conda-forge openbabel`
    
    Note that the angle which it is considering is 180-defined angle D-H-A (in a good agreement with VMD)
    Results can be displayed in VMD by using showVMDinteraction() """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    if atoms.hydrogen is None:
        raise ValueError('atoms should have hydrogens to calculate hydrogen bonds. '
                         'Use addMissingAtoms to add hydrogens')
    
    distA = kwargs.pop('distA', 3.5)
    angle = kwargs.pop('angle', 40)
    seq_cutoff = kwargs.pop('seq_cutoff', 25)
    
    donors = kwargs.get('donors', ['N', 'O', 'S', 'F'])
    acceptors = kwargs.get('acceptors', ['N', 'O', 'S', 'F'])
    
    if atoms.hydrogen == None or atoms.hydrogen.numAtoms() < 10:
        LOGGER.info("Provide structure with hydrogens or install Openbabel to add missing hydrogens using addMissingAtoms(pdb_name) first.")
    
    contacts = findNeighbors(atoms.heavy, distA)
    short_contacts = cleanNumbers(contacts)
    pairList = [] # list with Donor-Hydrogen-Acceptor(indices)-distance-Angle
    
    LOGGER.info('Calculating hydrogen bonds.')
    for nr_i,i in enumerate(short_contacts):
        # Removing those close contacts which are between neighbour atoms
        if i[1] - seq_cutoff < i[0] < i[1] + seq_cutoff:
            continue
        
        if (i[2][0] in donors and i[3][0] in acceptors) or (i[2][0] in acceptors and i[3][0] in donors): # First letter is checked
            listOfHydrogens1 = cleanNumbers(findNeighbors(atoms.hydrogen, 1.4, atoms.select('index '+str(i[0]))))
            listOfHydrogens2 = cleanNumbers(findNeighbors(atoms.hydrogen, 1.4, atoms.select('index '+str(i[1]))))
            AtomsForAngle = ['D','H','A', 'distance','angle']
            
            if not listOfHydrogens1:
                for j in listOfHydrogens2:
                    AtomsForAngle = [j[1], j[0], i[0], i[-1], calcAngle(atoms.select('index '+str(j[1])), 
                                                                    atoms.select('index '+str(j[0])), 
                                                                    atoms.select('index '+str(i[0])))[0]]                                                                                   
                    pairList.append(AtomsForAngle)            
            
            elif not listOfHydrogens2:
                for jj in listOfHydrogens1:
                    AtomsForAngle = [jj[1], jj[0], i[1], i[-1], calcAngle(atoms.select('index '+str(jj[1])), 
                                                                        atoms.select('index '+str(jj[0])), 
                                                                        atoms.select('index '+str(i[1])))[0]]
                    pairList.append(AtomsForAngle)            
    
            else:            
                for j in listOfHydrogens2:
                    AtomsForAngle = [j[1], j[0], i[0], i[-1], calcAngle(atoms.select('index '+str(j[1])), 
                                                                        atoms.select('index '+str(j[0])), 
                                                                        atoms.select('index '+str(i[0])))[0]]                                                                                   
                    pairList.append(AtomsForAngle)

                
                for jj in listOfHydrogens1:
                    AtomsForAngle = [jj[1], jj[0], i[1], i[-1], calcAngle(atoms.select('index '+str(jj[1])), 
                                                                            atoms.select('index '+str(jj[0])), 
                                                                            atoms.select('index '+str(i[1])))[0]]
                    pairList.append(AtomsForAngle)
    
    HBs_list = []
    for k in pairList:
        if 180-angle < float(k[-1]) < 180 and float(k[-2]) < distA:
            ag = atoms.getAtomGroup()
            aa_donor = ag.getResnames()[k[0]]+str(ag.getResnums()[k[0]])
            aa_donor_atom = ag.getNames()[k[0]]+'_'+str(k[0])
            aa_donor_chain = ag.getChids()[k[0]]
            aa_acceptor = ag.getResnames()[k[2]]+str(ag.getResnums()[k[2]])
            aa_acceptor_atom = ag.getNames()[k[2]]+'_'+str(k[2])
            aa_acceptor_chain = ag.getChids()[k[2]]
            
            HBs_list.append([str(aa_donor), str(aa_donor_atom), str(aa_donor_chain), str(aa_acceptor), str(aa_acceptor_atom), 
                             str(aa_acceptor_chain), np.round(float(k[-2]),2), np.round(180.0-float(k[-1]),2)])
    
    HBs_list = sorted(HBs_list, key=lambda x : x[-2])
    HBs_list_final = removeDuplicates(HBs_list)
    
    selection = kwargs.get('selection', None)
    selection2 = kwargs.get('selection2', None)    
    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
    HBs_list_final2 = filterInteractions(HBs_list_final, atoms, **sel_kwargs)
    
    LOGGER.info(("%26s   <---> %30s%12s%7s" % ('DONOR (res chid atom)','ACCEPTOR (res chid atom)','Distance','Angle')))
    for kk in HBs_list_final2:
        LOGGER.info("%10s%5s%14s  <---> %10s%5s%14s%8.1f%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6], kk[7]))
                                
    LOGGER.info("Number of detected hydrogen bonds: {0}.".format(len(HBs_list_final2)))
                
    return HBs_list_final2   
    
    
def calcChHydrogenBonds(atoms, **kwargs):
    """Finds hydrogen bonds between different chains.
    See more details in calcHydrogenBonds().
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between donor and acceptor.
        default is 3.0.
    :type distA: int, float

    :arg angle: non-zero value, D-H-A angle (donor, hydrogen, acceptor).
        default is 40.
    :type angle: int, float
    
    :arg seq_cutoff: non-zero value, interactions will be found between atoms with index differences
        that are higher than seq_cutoff.
        default is 25 atoms.
    :type seq_cutoff: int

    Structure should contain hydrogens.
    If not they can be added using addMissingAtoms(pdb_name) function available in ProDy after Openbabel installation.
    `conda install -c conda-forge openbabel`
    
    Note that the angle which it is considering is 180-defined angle D-H-A (in a good agreement with VMD)
    Results can be displayed in VMD. """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    distA = kwargs.pop('distA', 3.5)
    angle = kwargs.pop('angle', 40)
    seq_cutoff = kwargs.pop('seq_cutoff', 25)

    if len(np.unique(atoms.getChids())) > 1:
        HBS_calculations = calcHydrogenBonds(atoms, **kwargs)
    
        ChainsHBs = [ i for i in HBS_calculations if str(i[2]) != str(i[5]) ]
        if not ChainsHBs:
            ligand_name = list(set(atoms.select('all not protein and not ion').getResnames()))[0]
            ChainsHBs = [ ii for ii in HBS_calculations if ii[0][:3] == ligand_name or ii[3][:3] == ligand_name ]
        
        return ChainsHBs 
        

def calcSaltBridges(atoms, **kwargs):
    """Finds salt bridges in protein structure.
    Histidines are considered as charge residues.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between center of masses 
        of N and O atoms of negatively and positevely charged residues.
        default is 5.
    :type distA: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'
        
    Results can be displayed in VMD."""

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    
    distA = kwargs.pop('distA', 5.)
    atoms_KRED = atoms.select('protein and ((resname ASP GLU LYS ARG and not backbone and not name OXT NE "C.*" and noh) or (resname HIS HSE HSD HSP and name NE2))')
    charged_residues = list(set(zip(atoms_KRED.getResnums(), atoms_KRED.getChids())))
    
    LOGGER.info('Calculating salt bridges.')
    SaltBridges_list = []
    for i in charged_residues:
        sele1 = atoms_KRED.select('resid '+str(i[0])+' and chain '+i[1])
        try:
            sele1_center = calcCenter(sele1.getCoords())
            sele2 = atoms_KRED.select('same residue as exwithin '+str(distA)+' of center', center=sele1_center)
        except:
            sele1_center = sele1.getCoords()
            sele2 = atoms_KRED.select('same residue as exwithin '+str(distA)+' of center', center=sele1.getCoords())            
 
        if sele1 != None and sele2 != None:
            for ii in np.unique(sele2.getResnums()):                
                sele2_single = sele2.select('resid '+str(ii))
                try:
                    distance = calcDistance(sele1_center,calcCenter(sele2_single.getCoords()))
                except: 
                    distance = calcDistance(sele1_center,sele2_single.getCoords())
                
                if distance < distA and sele1.getNames()[0][0] != sele2_single.getNames()[0][0]:
                    SaltBridges_list.append([sele1.getResnames()[0]+str(sele1.getResnums()[0]), sele1.getNames()[0]+'_'+'_'.join(map(str,sele1.getIndices())), sele1.getChids()[0],
                                                  sele2_single.getResnames()[0]+str(sele2_single.getResnums()[0]), sele2_single.getNames()[0]+'_'+'_'.join(map(str,sele2_single.getIndices())), 
                                                  sele2_single.getChids()[0], round(distance,3)])
    
    SaltBridges_list = sorted(SaltBridges_list, key=lambda x : x[-1])
    [ SaltBridges_list.remove(j) for i in SaltBridges_list for j in SaltBridges_list if Counter(i) == Counter(j) ]
    SaltBridges_list_final = removeDuplicates(SaltBridges_list)
    
    selection = kwargs.get('selection', None)
    selection2 = kwargs.get('selection2', None)    
    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
    SaltBridges_list_final2 = filterInteractions(SaltBridges_list_final, atoms, **sel_kwargs)
    
    for kk in SaltBridges_list_final2:
        LOGGER.info("%10s%5s%16s  <---> %10s%5s%16s%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6]))
        
    LOGGER.info("Number of detected salt bridges: {0}.".format(len(SaltBridges_list_final2)))        

    return SaltBridges_list_final2
    

def calcRepulsiveIonicBonding(atoms, **kwargs):
    """Finds repulsive ionic bonding in protein structure
    i.e. between positive-positive or negative-negative residues.
    Histidine is not considered as a charged residue.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between center of masses 
            between N-N or O-O atoms of residues.
            default is 4.5.
    :type distA: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'

    Results can be displayed in VMD."""

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    
    distA = kwargs.pop('distA', 4.5)
    atoms_KRED = atoms.select('protein and resname ASP GLU LYS ARG and not backbone and not name OXT NE "C.*" and noh')
    charged_residues = list(set(zip(atoms_KRED.getResnums(), atoms_KRED.getChids())))
    
    LOGGER.info('Calculating repulsive ionic bonding.')
    RepulsiveIonicBonding_list = []
    for i in charged_residues:
        sele1 = atoms_KRED.select('resid '+str(i[0])+' and chain '+i[1])
        try:
            sele1_center = calcCenter(sele1.getCoords())
            sele2 = atoms_KRED.select('same residue as exwithin '+str(distA)+' of center', center=sele1_center)
        except:
            sele1_center = sele1.getCoords()
            sele2 = atoms_KRED.select('same residue as exwithin '+str(distA)+' of center', center=sele1_center)            
 
        if sele1 != None and sele2 != None:
            for ii in np.unique(sele2.getResnums()):                
                sele2_single = sele2.select('resid '+str(ii))
                try:
                    distance = calcDistance(sele1_center,calcCenter(sele2_single.getCoords()))
                except: 
                    distance = calcDistance(sele1_center,sele2_single.getCoords())
                
                if distance < distA and sele1.getNames()[0][0] == sele2_single.getNames()[0][0] and distance > 0:
                    RepulsiveIonicBonding_list.append([sele1.getResnames()[0]+str(sele1.getResnums()[0]), sele1.getNames()[0]+'_'+'_'.join(map(str,sele1.getIndices())), sele1.getChids()[0],
                                                  sele2_single.getResnames()[0]+str(sele2_single.getResnums()[0]), sele2_single.getNames()[0]+'_'+'_'.join(map(str,sele2_single.getIndices())), 
                                                  sele2_single.getChids()[0], round(distance,3)])
    
    [ RepulsiveIonicBonding_list.remove(j) for i in RepulsiveIonicBonding_list for j in RepulsiveIonicBonding_list if Counter(i) == Counter(j) ]
    RepulsiveIonicBonding_list = sorted(RepulsiveIonicBonding_list, key=lambda x : x[-1])
    RepulsiveIonicBonding_list_final = removeDuplicates(RepulsiveIonicBonding_list)

    selection = kwargs.get('selection', None)
    selection2 = kwargs.get('selection2', None)    
    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}    
    RepulsiveIonicBonding_list_final2 = filterInteractions(RepulsiveIonicBonding_list_final, atoms, **sel_kwargs)
    
    for kk in RepulsiveIonicBonding_list_final2:
        LOGGER.info("%10s%5s%16s  <---> %10s%5s%16s%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6]))
        
    LOGGER.info("Number of detected Repulsive Ionic Bonding interactions: {0}.".format(len(RepulsiveIonicBonding_list_final2)))
    
    return RepulsiveIonicBonding_list_final2


def calcPiStacking(atoms, **kwargs):
    """Finds π–π stacking interactions (between aromatic rings).
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between center of masses 
                of residues aromatic rings.
                default is 5.
    :type distA: int, float
    
    :arg angle_min: minimal angle between aromatic rings.
        default is 0.
    :type angle_min: int, float

    :arg angle_max: maximal angle between rings.
        default is 360.
    :type angle_max: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg non_standard: dictionary of non-standard residue in the protein structure
                        that need to be included in calculations
    :type non_standard: dictionary

    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'
    
    Results can be displayed in VMD.
    By default three residues are included TRP, PHE, TYR and HIS.
    Additional selection can be added: 
        >>> non_standard = {"HSE": "noh and not backbone and not name CB", 
                    "HSD": "noh and not backbone and not name CB"}
        >>> calcPiStacking(atoms, non_standard)
    
    Predictions for proteins only. 
    To compute protein-ligand interactions use calcLigandInteractions() or define **kwargs. """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    aromatic_dic = {'TRP':'noh and not backbone and not name CB NE1 CD1 CG',
                'PHE':'noh and not backbone and not name CB',
                'TYR':'noh and not backbone and not name CB and not name OH',
                'HIS':'noh and not backbone and not name CB'}
    
    distA = kwargs.pop('distA', 5.0)
    angle_min = kwargs.pop('angle_min', 0)
    angle_max = kwargs.pop('angle_max', 360)
    
    non_standard = kwargs.get('non_standard', {})
    for key, value in non_standard.items():
        aromatic_dic[key] = value
    
    atoms_cylic = atoms.select('resname TRP PHE TYR HIS')
    if atoms_cylic is None:
        return []
    
    aromatic_resids = list(set(zip(atoms_cylic.getResnums(), atoms_cylic.getChids())))

    LOGGER.info('Calculating Pi stacking interactions.')
    PiStack_calculations = []
    for i in aromatic_resids:
        for j in aromatic_resids:
            if i != j: 
                sele1_name = atoms.select('resid '+str(i[0])+' and chain '+i[1]+' and name CA').getResnames()
                sele1 = atoms.select('resid '+str(i[0])+' and chain '+i[1]+' and '+aromatic_dic[sele1_name[0]])
                
                sele2_name = atoms.select('resid '+str(j[0])+' and chain '+j[1]+' and name CA').getResnames()
                sele2 = atoms.select('resid '+str(j[0])+' and chain '+j[1]+' and '+aromatic_dic[sele2_name[0]])
                
                if sele1 != None and sele2 != None:
                    a1, b1, c1, a2, b2, c2 = calcPlane(sele1)[:3]+calcPlane(sele2)[:3]
                    RingRing_angle = calcAngleBetweenPlanes(a1, b1, c1, a2, b2, c2) # plane is computed based on 3 points of rings           
                    RingRing_distance = calcDistance(calcCenter(sele1.getCoords()),calcCenter(sele2.getCoords()))
                    if RingRing_distance < distA and angle_min < RingRing_angle < angle_max:
                        PiStack_calculations.append([str(sele1_name[0])+str(sele1.getResnums()[0]), '_'.join(map(str,sele1.getIndices())), str(sele1.getChids()[0]),
                                                     str(sele2_name[0])+str(sele2.getResnums()[0]), '_'.join(map(str,sele2.getIndices())), str(sele2.getChids()[0]),
                                                     round(RingRing_distance,3), round(RingRing_angle,3)])
    
    PiStack_calculations = sorted(PiStack_calculations, key=lambda x : x[-2])   
    PiStack_calculations_final = removeDuplicates(PiStack_calculations)
    
    selection = kwargs.get('selection', None)
    selection2 = kwargs.get('selection2', None)    
    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
    PiStack_calculations_final2 = filterInteractions(PiStack_calculations_final, atoms, **sel_kwargs)
    
    for kk in PiStack_calculations_final2:
        LOGGER.info("%10s%8s%32s  <---> %10s%8s%32s%8.1f%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6], kk[7]))
        
    LOGGER.info("Number of detected Pi stacking interactions: {0}.".format(len(PiStack_calculations_final2)))
    
    return PiStack_calculations_final2


def calcPiCation(atoms, **kwargs):
    """Finds cation-Pi interaction i.e. between aromatic ring and 
    positively charged residue (ARG and LYS).
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between center of masses 
                of aromatic ring and positively charge group.
                default is 5.
    :type distA: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg non_standard: dictionary of non-standard residue in the protein structure
                        that need to be included in calculations
    :type non_standard: dictionary

    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'

    By default three residues are included TRP, PHE, TYR and HIS.
    Additional selection can be added: 
        >>> calcPiCation(atoms, 'HSE'='noh and not backbone and not name CB')
        or
        >>> non_standard = {"HSE": "noh and not backbone and not name CB", 
                "HSD": "noh and not backbone and not name CB"}
        >>> calcPiCation(atoms, non_standard)
    
    Results can be displayed in VMD.
    Predictions for proteins only. To compute protein-ligand interactions use 
    calcLigandInteractions() or define **kwargs"""

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    
    aromatic_dic = {'TRP':'noh and not backbone and not name CB NE1 CD1 CG',
                'PHE':'noh and not backbone and not name CB',
                'TYR':'noh and not backbone and not name CB and not name OH',
                'HIS':'noh and not backbone and not name CB'}
    
    distA = kwargs.pop('distA', 5.0)
    
    non_standard = kwargs.get('non_standard', {})
    for key, value in non_standard.items():
        aromatic_dic[key] = value
        
    atoms_cylic = atoms.select('resname TRP PHE TYR HIS')
    if atoms_cylic is None:
        return []
    
    aromatic_resids = list(set(zip(atoms_cylic.getResnums(), atoms_cylic.getChids())))

    PiCation_calculations = []
    LOGGER.info('Calculating cation-Pi interactions.')
    
    for i in aromatic_resids:
        sele1_name = atoms.select('resid '+str(i[0])+' and chain '+i[1]+' and name CA').getResnames()
        
        try:
            sele1 = atoms.select('resid '+str(i[0])+' and chain '+i[1]+' and '+aromatic_dic[sele1_name[0]])
            sele2 = atoms.select('(same residue as exwithin '+str(distA)+' of center) and resname ARG LYS and noh and not backbone and not name NE "C.*"', 
                               center=calcCenter(sele1.getCoords()))
        except:
            raise ValueError("Missing atoms from the side chains of the structure. Use addMissingAtoms.")
        
        if sele1 != None and sele2 != None:
            for ii in np.unique(sele2.getResnums()):
                sele2_single = sele2.select('resid '+str(ii))
                try:
                    RingCation_distance = calcDistance(calcCenter(sele1.getCoords()),calcCenter(sele2_single.getCoords()))
                except: 
                    RingCation_distance = calcDistance(calcCenter(sele1.getCoords()),sele2_single.getCoords())
                
                if RingCation_distance < distA:
                    PiCation_calculations.append([str(sele1_name[0])+str(sele1.getResnums()[0]), '_'.join(map(str,sele1.getIndices())), str(sele1.getChids()[0]),
                                                  str(sele2_single.getResnames()[0])+str(sele2_single.getResnums()[0]), sele2_single.getNames()[0]+'_'+'_'.join(map(str,sele2_single.getIndices())), 
                                                  str(sele2_single.getChids()[0]), round(RingCation_distance,3)])
    
    PiCation_calculations = sorted(PiCation_calculations, key=lambda x : x[-1]) 
    PiCation_calculations_final = removeDuplicates(PiCation_calculations)
    
    selection = kwargs.get('selection', None)
    selection2 = kwargs.get('selection2', None)    
    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
    PiCation_calculations_final2 = filterInteractions(PiCation_calculations_final, atoms, **sel_kwargs)
    
    for kk in PiCation_calculations_final2:
        LOGGER.info("%10s%4s%32s  <---> %10s%4s%32s%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6]))
        
    LOGGER.info("Number of detected cation-pi interactions: {0}.".format(len(PiCation_calculations_final2)))
    
    return PiCation_calculations_final2


def calcHydrophobic(atoms, **kwargs): 
    """Prediction of hydrophobic interactions between hydrophobic residues 
    (ALA, ILE, LEU, MET, PHE, TRP, VAL).
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between atoms of hydrophobic residues.
        default is 4.5.
    :type distA: int, float
    
    :arg non_standard: dictionary of non-standard residue in the protein structure
                        that need to be included in calculations
    :type non_standard: dictionary

    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'
    
    Additional selection can be added as shown below (with selection that includes 
    only hydrophobic part): 
        >>> calcHydrophobic(atoms, non_standard={'XLE'='noh and not backbone', 
                                                'XLI'='noh and not backbone'})
    Predictions for proteins only. To compute protein-ligand interactions use 
    calcLigandInteractions().
    Results can be displayed in VMD by using showVMDinteraction() 
    
    Note that interactions between aromatic residues are omitted becasue they are 
    provided by calcPiStacking(). """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    
    distA = kwargs.pop('distA', 4.5)
    
    Hydrophobic_list = []  
    atoms_hydrophobic = atoms.select('resname ALA VAL ILE MET LEU PHE TYR TRP')
    hydrophobic_resids = list(set(zip(atoms_hydrophobic.getResnums(), atoms_hydrophobic.getChids())))

    if atoms.aromatic is None:
        return []
    
    aromatic_nr = list(set(zip(atoms.aromatic.getResnums(),atoms.aromatic.getChids())))   
    aromatic = list(set(zip(atoms.aromatic.getResnames())))
    
    hydrophobic_dic = {'ALA': 'noh and not backbone', 'VAL': 'noh and not (backbone or name CB)',
    'ILE': 'noh and not (backbone or name CB)', 'LEU': 'noh and not (backbone or name CB)',
    'MET': 'noh and not (backbone or name CB)', 'PHE': 'noh and not (backbone or name CB)',
    'TYR': 'noh and not (backbone or name CB)', 'TRP': 'noh and not (backbone or name CB)'}

    non_standard = kwargs.get('non_standard', {})
    for key, value in non_standard.items():
        aromatic_dic[key] = value
    
    LOGGER.info('Calculating hydrophobic interactions.')
    Hydrophobic_calculations = []
    for i in hydrophobic_resids:
        try:
            sele1_name = atoms.select('resid '+str(i[0])+' and chain '+i[1]+' and name CA').getResnames()
            sele1 = atoms.select('resid '+str(i[0])+' and '+' chain '+i[1]+' and '+hydrophobic_dic[sele1_name[0]]) 
            sele1_nr = sele1.getResnums()[0]  
            sele2 = atoms.select('(same residue as exwithin '+str(distA)+' of (resid '+str(sele1_nr)+' and chain '+i[1]+' and resname '+sele1_name[0]+
                               ')) and ('+' or '.join([ '(resname '+item[0]+' and '+item[1]+')' for item in hydrophobic_dic.items() ])+')')

        except:
            LOGGER.info("Missing atoms from the side chains of the structure. Use PDBFixer.")
            sele1 = None
            sele2 = None
        
        if sele2 != None:
            sele2_nr = list(set(zip(sele2.getResnums(), sele2.getChids())))

            if sele1_name[0] in aromatic:
                sele2_filter = sele2.select('all and not (resname TYR PHE TRP or resid '+str(i)+')')
                if sele2_filter != None:
                    listOfAtomToCompare = cleanNumbers(findNeighbors(sele1, distA, sele2_filter))
                
            elif sele1_name[0] not in aromatic and i in sele2_nr:
                sele2_filter = sele2.select(sele2.select('all and not (resid '+str(i[0])+' and chain '+i[1]+')'))
                if sele2_filter != None:
                    listOfAtomToCompare = cleanNumbers(findNeighbors(sele1, distA, sele2_filter))
            else:
                listOfAtomToCompare = cleanNumbers(findNeighbors(sele1, distA, sele2))
                                                           
            if listOfAtomToCompare != []:
                listOfAtomToCompare = sorted(listOfAtomToCompare, key=lambda x : x[-1])
                minDistancePair = listOfAtomToCompare[0]
                if minDistancePair[-1] < distA:
                    sele1_new = atoms.select('index '+str(minDistancePair[0])+' and name '+str(minDistancePair[2]))
                    sele2_new = atoms.select('index '+str(minDistancePair[1])+' and name '+str(minDistancePair[3]))
                    Hydrophobic_calculations.append([sele1_new.getResnames()[0]+str(sele1_new.getResnums()[0]), 
                                                             minDistancePair[2]+'_'+str(minDistancePair[0]), sele1_new.getChids()[0],
                                                             sele2_new.getResnames()[0]+str(sele2_new.getResnums()[0]), 
                                                             minDistancePair[3]+'_'+str(minDistancePair[1]), sele2_new.getChids()[0],
                                                             round(minDistancePair[-1],3)]) 
                    
    Hydrophobic_calculations = sorted(Hydrophobic_calculations, key=lambda x : x[-1])
    Hydrophobic_calculations_final = removeDuplicates(Hydrophobic_calculations)
    
    selection = kwargs.get('selection', None)
    selection2 = kwargs.get('selection2', None)    
    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
    Hydrophobic_calculations_final2 = filterInteractions(Hydrophobic_calculations_final, atoms, **sel_kwargs)
    
    for kk in Hydrophobic_calculations_final2:
        LOGGER.info("%10s%5s%14s  <---> %10s%5s%14s%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6]))
        
    LOGGER.info("Number of detected hydrophobic interactions: {0}.".format(len(Hydrophobic_calculations_final2)))
    
    return Hydrophobic_calculations_final2


def calcDisulfideBonds(atoms, **kwargs):
    """Prediction of disulfide bonds.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between atoms of hydrophobic residues.
        default is 3.
    :type distA: int, float
    """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    
    distA = kwargs.pop('distA', 3)
    
    try:
        atoms_SG = atoms.select('protein and resname CYS and name SG')
        atoms_SG_res = list(set(zip(atoms_SG.getResnums(), atoms_SG.getChids())))
    
        LOGGER.info('Calculating disulfide bonds.')
        DisulfideBonds_list = []
        for i in atoms_SG_res:
            CYS_pairs = atoms.select('(same residue as protein within '+str(distA)+' of ('+'resid '+str(i[0])+' and chain '+i[1]+' and name SG)) and (resname CYS and name SG)')
            if CYS_pairs.numAtoms() > 1:
                sele1 = CYS_pairs[0]
                sele2 = CYS_pairs[1]

                listOfAtomToCompare = cleanNumbers(findNeighbors(sele1, distA, sele2))
                if listOfAtomToCompare != []:
                    listOfAtomToCompare = sorted(listOfAtomToCompare, key=lambda x : x[-1])
                    minDistancePair = listOfAtomToCompare[0]
                    if minDistancePair[-1] < distA:
                        sele1_new = atoms.select('index '+str(minDistancePair[0])+' and name '+str(minDistancePair[2]))
                        sele2_new = atoms.select('index '+str(minDistancePair[1])+' and name '+str(minDistancePair[3]))
                        DisulfideBonds_list.append([sele1_new.getResnames()[0]+str(sele1_new.getResnums()[0]),
                                                                minDistancePair[2]+'_'+str(minDistancePair[0]), sele1_new.getChids()[0],
                                                                sele2_new.getResnames()[0]+str(sele2_new.getResnums()[0]),
                                                                minDistancePair[3]+'_'+str(minDistancePair[1]), sele2_new.getChids()[0],
                                                                round(minDistancePair[-1],3)])
    except:
        atoms_SG = atoms.select('protein and resname CYS')
        if atoms_SG is None:
            LOGGER.info('Lack of cysteines in the structure.')
            DisulfideBonds_list = []

    DisulfideBonds_list_final = removeDuplicates(DisulfideBonds_list)

    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
    DisulfideBonds_list_final2 = filterInteractions(DisulfideBonds_list_final, atoms, **sel_kwargs)

    for kk in DisulfideBonds_list_final2:
        LOGGER.info("%10s%5s%14s  <---> %10s%5s%14s%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6]))

    LOGGER.info("Number of detected disulfide bonds: {0}.".format(len(DisulfideBonds_list_final2)))

    return DisulfideBonds_list_final2


def calcMetalInteractions(atoms, distA=3.0, extraIons=['FE'], excluded_ions=['SOD', 'CLA']):
    """Interactions with metal ions (includes water, ligands and other ions).
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distA: non-zero value, maximal distance between ion and residue.
        default is 3.0
    :type distA: int, float
    
    :arg extraIons: ions to be included in the analysis.
    :type extraIons: list
    
    :arg excluded_ions: ions which should be excluded from the analysis.
    :type excluded_ions: list """
    
    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    
    try:
        atoms_ions = atoms.select('ion and not name '+' '.join(excluded_ions)+' or (name '+' '.join(map(str,extraIons))+')')
        MetalResList = []
        MetalRes_calculations = cleanNumbers(findNeighbors(atoms_ions, distA, atoms.select('all and noh')))
        for i in MetalRes_calculations:
            if i[-1] != 0:
                MetalResList.append([atoms.getResnames()[i[0]]+str(atoms.getResnums()[i[0]]), i[2], 
                                 atoms.getResnames()[i[1]]+str(atoms.getResnums()[i[1]]), i[3], i[-1]])

        return MetalResList
        
    except TypeError:
        raise TypeError('An object should contain ions')


def calcInteractionsMultipleFrames(atoms, interaction_type, trajectory, **kwargs):
    """Compute selected type interactions for DCD trajectory or multi-model PDB 
    using default parameters."""
    
    try:
        coords = getCoords(atoms)
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')    
    
    interactions_all = []
    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)

    interactions_dic = {
    "HBs": calcHydrogenBonds,
    "SBs": calcSaltBridges,
    "RIB": calcRepulsiveIonicBonding,
    "PiStack": calcPiStacking,
    "PiCat": calcPiCation,
    "HPh": calcHydrophobic,
    "DiB": calcDisulfideBonds
    }
    
    if trajectory is not None: 
        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)
        
        nfi = trajectory._nfi    
        trajectory.reset()
        numFrames = trajectory._n_csets
        
        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame+1]
        
        for j0, frame0 in enumerate(traj, start=start_frame):
            LOGGER.info('Frame: {0}'.format(j0))
            protein = atoms.select('protein')
            interactions = interactions_dic[interaction_type](protein, **kwargs)
            interactions_all.append(interactions)
        trajectory._nfi = nfi
    
    else:
        if atoms.numCoordsets() > 1:
            for i in range(len(atoms.getCoordsets()[start_frame:stop_frame])):
                LOGGER.info('Model: {0}'.format(i+start_frame))
                atoms.setACSIndex(i+start_frame)
                protein = atoms.select('protein')
                interactions = interactions_dic[interaction_type](protein, **kwargs)
                interactions_all.append(interactions)
        else:
            LOGGER.info('Include trajectory or use multi-model PDB file.')
    
    return interactions_all


def calcProteinInteractions(atoms, **kwargs):
    """Compute all protein interactions (shown below) using default parameters.
        (1) Hydrogen bonds
        (2) Salt Bridges
        (3) RepulsiveIonicBonding 
        (4) Pi stacking interactions
        (5) Pi-cation interactions
        (6) Hydrophobic interactions
        (7) Disulfide Bonds
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str 
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'  """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    LOGGER.info('Calculating all interations.') 
    HBs_calculations = calcHydrogenBonds(atoms.protein, **kwargs)               #1 in counting
    SBs_calculations = calcSaltBridges(atoms.protein, **kwargs)                 #2
    SameChargeResidues = calcRepulsiveIonicBonding(atoms.protein, **kwargs)     #3
    Pi_stacking = calcPiStacking(atoms.protein, **kwargs)                       #4
    Pi_cation = calcPiCation(atoms.protein, **kwargs)                           #5
    Hydroph_calculations = calcHydrophobic(atoms.protein, **kwargs)             #6
    Disulfide_Bonds = calcDisulfideBonds(atoms.protein, **kwargs)		#7
    AllInteractions = [HBs_calculations, SBs_calculations, SameChargeResidues, Pi_stacking, 
                        Pi_cation, Hydroph_calculations, Disulfide_Bonds]   
    
    return AllInteractions


def calcHydrogenBondsTrajectory(atoms, trajectory=None, **kwargs):   
    """Compute hydrogen bonds for DCD trajectory or multi-model PDB using default parameters.
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
        
    :arg trajectory: trajectory file
    :type trajectory: class:`.Trajectory`

    :arg distA: non-zero value, maximal distance between donor and acceptor.
        default is 3.5
    :type distA: int, float
    
    :arg angle: non-zero value, maximal (180 - D-H-A angle) (donor, hydrogen, acceptor).
        default is 40.
    :type angle: int, float
    
    :arg seq_cutoff: non-zero value, interactions will be found between atoms with index differences
        that are higher than seq_cutoff.
        default is 20 atoms.
    :type seq_cutoff: int

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg start_frame: index of first frame to read
    :type start_frame: int

    :arg stop_frame: index of last frame to read
    :type stop_frame: int
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B' """
    
    return calcInteractionsMultipleFrames(atoms, 'HBs', trajectory, **kwargs)


def calcSaltBridgesTrajectory(atoms, trajectory=None, **kwargs):
    """Compute salt bridges for DCD trajectory or multi-model PDB using default parameters.
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg trajectory: trajectory file
    :type trajectory: class:`.Trajectory`
    
    :arg distA: non-zero value, maximal distance between center of masses 
        of N and O atoms of negatively and positevely charged residues.
        default is 5.
    :type distA: int, float
    
    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg start_frame: index of first frame to read
    :type start_frame: int

    :arg stop_frame: index of last frame to read
    :type stop_frame: int
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'  """

    return calcInteractionsMultipleFrames(atoms, 'SBs', trajectory, **kwargs)
    

def calcRepulsiveIonicBondingTrajectory(atoms, trajectory=None, **kwargs):  
    """Compute repulsive ionic bonding for DCD trajectory or multi-model PDB 
    using default parameters.
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg trajectory: trajectory file
    :type trajectory: class:`.Trajectory`
    
    :arg distA: non-zero value, maximal distance between center of masses 
            between N-N or O-O atoms of residues.
            default is 4.5.
    :type distA: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg start_frame: index of first frame to read
    :type start_frame: int

    :arg stop_frame: index of last frame to read
    :type stop_frame: int
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'  """

    return calcInteractionsMultipleFrames(atoms, 'RIB', trajectory, **kwargs)


def calcPiStackingTrajectory(atoms, trajectory=None, **kwargs):   
    """Compute Pi-stacking interactions for DCD trajectory or multi-model PDB 
    using default parameters.
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
      
    :arg trajectory: trajectory file
    :type trajectory: class:`.Trajectory`
    
    :arg distA: non-zero value, maximal distance between center of masses 
                of residues aromatic rings.
                default is 5.
    :type distA: int, float
    
    :arg angle_min: minimal angle between aromatic rings.
        default is 0.
    :type angle_min: int

    :arg angle_max: maximal angle between rings.
        default is 360.
    :type angle_max: int, float
    
    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg start_frame: index of first frame to read
    :type start_frame: int

    :arg stop_frame: index of last frame to read
    :type stop_frame: int
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B' """            

    return calcInteractionsMultipleFrames(atoms, 'PiStack', trajectory, **kwargs)


def calcPiCationTrajectory(atoms, trajectory=None, **kwargs):  
    """Compute Pi-cation interactions for DCD trajectory or multi-model PDB 
    using default parameters.
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
        
    :arg trajectory: trajectory file
    :type trajectory: class:`.Trajectory`
    
    :arg distA: non-zero value, maximal distance between center of masses of aromatic ring 
                and positively charge group.
                default is 5.
    :type distA: int, float
    
    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg start_frame: index of first frame to read
    :type start_frame: int

    :arg stop_frame: index of last frame to read
    :type stop_frame: int

    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B' """

    return calcInteractionsMultipleFrames(atoms, 'PiCat', trajectory, **kwargs)


def calcHydrophobicTrajectory(atoms, trajectory=None, **kwargs):  
    """Compute hydrophobic interactions for DCD trajectory or multi-model PDB 
    using default parameters.
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
        
    :arg trajectory: trajectory file
    :type trajectory: class:`.Trajectory`
    
    :arg distA: non-zero value, maximal distance between atoms of hydrophobic residues.
        default is 4.5.
    :type distA: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg start_frame: index of first frame to read
    :type start_frame: int

    :arg stop_frame: index of last frame to read
    :type stop_frame: int
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B' """

    return calcInteractionsMultipleFrames(atoms, 'HPh', trajectory, **kwargs)


def calcDisulfideBondsTrajectory(atoms, trajectory=None, **kwargs):
    """Compute disulfide bonds for DCD trajectory or multi-model PDB using default parameters.
        
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
        
    :arg trajectory: trajectory file
    :type trajectory: class:`.Trajectory`
    
    :arg distA: non-zero value, maximal distance between atoms of hydrophobic residues.
        default is 2.5.
    :type distA: int, float

    :arg start_frame: index of first frame to read
    :type start_frame: int

    :arg stop_frame: index of last frame to read
    :type stop_frame: int """

    return calcInteractionsMultipleFrames(atoms, 'DiB', trajectory, **kwargs)


def compareInteractions(data1, data2, **kwargs):
    """Comparison of two outputs from interactions. 
    It will provide information about the disappearance and appearance of some interactions
    as well as the similarities in the interactions
    for the same system. Two conformations can be compared.
    
    :arg data1: list with interactions from calcHydrogenBonds() or other types
    :type data1: list
 
    :arg data2: list with interactions from calcHydrogenBonds() or other types
    :type data2: list
    
    :arg filename: name of text file in which the comparison between two sets of interactions 
                will be saved 
    type filename: str 
    
    Example of usage: 
    >>> atoms1 = parsePDB('PDBfile1.pdb').select('protein')
    >>> atoms2 = parsePDB('PDBfile2.pdb').select('protein')
    >>> compareInteractions(calcHydrogenBonds(atoms1), calcHydrogenBonds(atoms2))
    """
    
    if not isinstance(data1, list):
        raise TypeError('data1 must be a list of interactions.')

    if not isinstance(data2, list):
        raise TypeError('data2 must be a list of interactions.')        

    data1_tuple = [ tuple([i[0]+i[2], i[3]+i[5]]) for i in data1 ]
    data2_tuple = [ tuple([i[0]+i[2], i[3]+i[5]]) for i in data2 ]
    diff_21 = set(data2_tuple) - set(data1_tuple)
    diff_12 = set(data1_tuple) - set(data2_tuple)
    similar_12 = set(data1_tuple) & set(data2_tuple)
    
    LOGGER.info("Which interactions disappeared: {0}".format(len(diff_21)))
    for j in diff_21:
        LOGGER.info('{0} <---> {1}'.format(j[0],j[1]))
        
    LOGGER.info("\nWhich interactions appeared: {0}".format(len(diff_12)))  
    for j in diff_12:  
        LOGGER.info('{0} <---> {1}'.format(j[0],j[1]))
    
    LOGGER.info("Which interactions are the same: {0}".format(len(similar_12)))
    for j in similar_12:
        if len(similar_12) != 0:
            LOGGER.info('{0} <---> {1}'.format(j[0],j[1]))
        else: LOGGER.info("None")
    
    try:
        if 'filename' in kwargs:
            with open(kwargs['filename'], 'w') as f:  # what disapperaed from initial
                f.write("Which interactions disappeared:\n")
                for i in diff_21:
                    f.write(i[0]+'-'+i[1]+'\n')
                f.write("\nWhich interactions appeared:\n")
                for i in diff_12:
                    f.write(i[0]+'-'+i[1]+'\n')
                f.write("\nWhich interactions are the same:\n")
                for i in similar_12:
                    f.write(i[0]+'-'+i[1]+'\n')
            f.close()
    except: pass
    
    return diff_21, diff_12, similar_12


def showInteractionsGraph(statistics, **kwargs):
    """Return residue-residue interactions as graph/network.
    
    :arg statistics: Results obtained from calcStatisticsInteractions analysis
    :type statistics: list
    
    :arg cutoff: Minimal number of counts per residue in the trajectory
    :type cutoff: int, by default 3.

    :arg code: representation of the residues, 3-letter or 1-letter
    :type code: str, by default 3-letter.

    :arg edge_cmap: color of the residue connection
    :type edge_cmap: str, by default plt.cm.Blues (blue color).

    :arg node_size: size of the nodes which describes residues
    :type node_size: int, by default 300.
    
    :arg node_distance: value which will scale residue-residue interactions
    :type node_distance: int, by default 5.

    :arg font_size: size of the font
    :type font_size: int, by default 14.

    :arg seed: random number which affect the distribution of residues
    :type seed: int, by default 42.  """
    
    import networkx as nx
    import matplotlib.pyplot as plt
    
    amino_acid_colors_dic = {
        "ALA": "silver",     # non-polar
        "ILE": "silver",     # non-polar
        "LEU": "silver",     # non-polar
        "VAL": "silver",     # non-polar
        "PHE": "silver",     # non-polar
        "MET": "silver",     # non-polar
        "TRP": "silver",     # non-polar
        "GLY": "limegreen",     # polar
        "SER": "limegreen",     # polar
        "THR": "limegreen",     # polar
        "CYS": "limegreen",     # polar
        "TYR": "limegreen",     # polar
        "ASN": "limegreen",     # polar
        "GLN": "limegreen",     # polar
        "HIS": "deepskyblue",      # basic
        "HSE": "deepskyblue",      # basic
        "HSD": "deepskyblue",      # basic
        "LYS": "deepskyblue",      # basic
        "ARG": "deepskyblue",      # basic
        "ASP": "tomato",       # acidic
        "GLU": "tomato",       # acidic
        "PRO": "pink",
        "A": "silver",     # non-polar
        "I": "silver",     # non-polar
        "L": "silver",     # non-polar
        "V": "silver",     # non-polar
        "F": "silver",     # non-polar
        "M": "silver",     # non-polar
        "W": "silver",     # non-polar
        "G": "limegreen",     # polar
        "S": "limegreen",     # polar
        "T": "limegreen",     # polar
        "C": "limegreen",     # polar
        "Y": "limegreen",     # polar
        "N": "limegreen",     # polar
        "Q": "limegreen",     # polar
        "H": "deepskyblue",      # basic
        "K": "deepskyblue",      # basic
        "R": "deepskyblue",      # basic
        "D": "tomato",       # acidic
        "E": "tomato",       # acidic
        "P": "pink" }
        
    aa_dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
              'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 
              'HIS': 'H', 'HSD': 'H','HSE': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
              'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    if len(statistics[0]) != 4:
        raise TypeError('data must be a list obtained from calcStatisticsInteractions')
    else:
        if isinstance(statistics, int) or isinstance(statistics, str):
            raise TypeError('node_size must be a list')

    code = kwargs.pop('code', None)
    if code is None:
        code = '3-letter'
    elif isinstance(code, str):
        code = code
    elif isinstance(code, int):
        raise TypeError('code must be 3-letter or 1-letter')

    edge_cmap = kwargs.pop('edge_cmap', plt.cm.Blues)
    node_size = kwargs.pop('node_size', 300)
    node_distance = kwargs.pop('node_distance', 5)
    font_size = kwargs.pop('font_size', 14)
    seed = kwargs.pop('seed', 42)
    cutoff = kwargs.pop('cutoff', 3)
    
    X = [i for i in statistics if i[1] >= cutoff]   
    G = nx.Graph()
    
    for row in X:  
        if code == '1-letter':
            aa1 = aa_dic[row[0].split('-')[0][:3]] + row[0].split('-')[0][3:]
            aa2 = aa_dic[row[0].split('-')[1][:3]] + row[0].split('-')[1][3:]
        else:
            aa1 = row[0].split('-')[0]
            aa2 = row[0].split('-')[1]   
            
        G.add_node(aa1)
        G.add_node(aa2)
        weight = row[1]
        length = row[2]
        G.add_edge(aa1, aa2, weight=weight, length=length)
        
    try:
        node_colors = [ amino_acid_colors_dic[i[:3]] for i in G.nodes() ]
    except:
        node_colors = [ amino_acid_colors_dic[i[:1]] for i in G.nodes() ]

    if SETTINGS['auto_show']:
        plt.figure()

    pos = nx.spring_layout(G, k=node_distance, seed=seed)
    edges = G.edges()
    weights = [G[u][v]['weight'] for u,v in edges]
    lengths = [G[u][v]['length'] for u,v in edges]

    show = nx.draw(G, pos, edgelist=edges, edge_color=weights, width=lengths, edge_cmap=edge_cmap, 
            node_color=node_colors, node_size=node_size, font_size=font_size, with_labels=True)
    
    if SETTINGS['auto_show']:
        showFigure()
    return show
    
    
def calcStatisticsInteractions(data):
    """Return the statistics of interactions from DCD trajectory or multi-model PDB including
    the number of counts for each residue pair, 
    average distance of interactions for each pair [in Ang] and standard deviation.
        
    :arg data: list with interactions from calcHydrogenBondsTrajectory() or other types
    :type data: list
    
    Example of usage: 
    >>> atoms = parsePDB('PDBfile.pdb')
    >>> dcd = Trajectory('DCDfile.dcd')
    >>> dcd.link(atoms)
    >>> dcd.setCoords(atoms)
    
    >>> data = calcPiCationTrajectory(atoms, dcd, distA=5)
    or
    >>> interactionsTrajectory = InteractionsTrajectory()
    >>> data = interactionsTrajectory.getPiCation()
    """
    
    interactions_list = [ (jj[0]+jj[2]+'-'+jj[3]+jj[5], jj[6]) for ii in data for jj in ii]
    import numpy as np
    elements = [t[0] for t in interactions_list]
    stats = {}

    for element in elements:
        if element not in stats:
            values = [t[1] for t in interactions_list if t[0] == element]
            stats[element] = {
                "stddev": np.round(np.std(values),2),
                "mean": np.round(np.mean(values),2),
                "count": len(values)
            }

    statistic = []
    for key, value in stats.items():
        LOGGER.info("Statistics for {0}:".format(key))
        LOGGER.info("  Average [Ang.]: {}".format(value['mean']))
        LOGGER.info("  Standard deviation [Ang.]: {0}".format(value['stddev']))
        LOGGER.info("  Count: {0}".format(value['count']))
        statistic.append([key, value['count'], value['mean'], value['stddev']])
    
    statistic.sort(key=lambda x: x[1], reverse=True)
    return statistic


def calcDistribution(interactions, residue1, residue2=None, **kwargs):
    """Distributions/histograms of pairs of amino acids. 
    Histograms are normalized.

    :arg interactions: list of interactions
    :type interactions: list
    
    :arg residue1: residue name in 3-letter code and residue number
    :type residue1: str
    
    :arg residue2: residue name in 3-letter code and residue number
    :type residue2: str
    
    :arg metrics: name of the data type
        'distance' or 'angle' depends on the type of interaction
    :type metrics: str
    """
    import matplotlib
    import matplotlib.pyplot as plt
    metrics = kwargs.pop('metrics', 'distance')
    
    if residue2 == None:
        additional_residues = []
    
        for sublist in interactions:
            for line in sublist:
                if (residue1 in line[0] or residue1 in line[3]):
                    additional_residues.append(line[0])
                    additional_residues.append(line[3])
        
        additional_residues = list(set(additional_residues))
        additional_residues.remove(residue1)
        
        if (additional_residues) == []:
            pass
        else:
            LOGGER.info('Possible contacts for '+residue1+':')
            for i in additional_residues:
                LOGGER.info(i)
    
    else:
        values = []
        additional_residues = []

        for sublist in interactions:
            for line in sublist:
                if (residue1 in line[0] or residue1 in line[3]):
                    additional_residues.append(line[0])
                    additional_residues.append(line[3])
                    if residue2 in line:
                        if metrics == 'distance':
                            values.append(line[6])
                        elif metrics == 'angle':
                            values.append(line[7])

        plt.hist(values, rwidth=0.95, density=True)
        plt.ylabel('Probability')
        
        if metrics == 'distance':
            plt.xlabel('Distance [Å]')
        elif metrics == 'angle':
            plt.xlabel('Angle [deg.]')
        
        plt.show()

        additional_residues = list(set(additional_residues))
        additional_residues.remove(residue1)
        additional_residues.remove(residue2)

        if (additional_residues) == []:
            pass
        else:
            LOGGER.info('Additional contacts for '+residue1+':')
            for i in additional_residues:
                LOGGER.info(i)

    

def calcLigandInteractions(atoms, **kwargs):
    """Provide ligand interactions with other elements of the system including protein, 
    water and ions. Results are computed by PLIP [SS15]_ which should be installed.
    Note that PLIP will not recognize ligand unless it will be HETATM in the PDB file.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg select: a selection string for residues of interest
            default is 'all not (water or protein or ion)'
    :type select: str
    
    :arg ignore_ligs: List of ligands which will be excluded from the analysis.
    :type ignore_ligs: list
    
    To display results as a list of interactions use listLigandInteractions()
    and for visualization in VMD please use showLigandInteraction_VMD() 
    
    Requirements of usage:
    ## Instalation of Openbabel:
    >>> conda install -c conda-forge openbabel    
    ## https://anaconda.org/conda-forge/openbabel
    
    ## Instalation of PLIP
    >>> conda install -c conda-forge plip
    ## https://anaconda.org/conda-forge/plip
    # https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md

    .. [SS15] Salentin S., Schreiber S., Haupt V. J., Adasme M. F., Schroeder M.  
    PLIP: fully automated protein–ligand interaction profiler 
    *Nucl. Acids Res.* **2015** 43:W443-W447.  """
    
    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    try:
        from plip.structure.preparation import PDBComplex   
        
        pdb_name = atoms.getTitle()+'_sele.pdb'
        LOGGER.info("Writing PDB file with selection in the local directory.")
        writePDB(pdb_name, atoms)

        try:
            if atoms.hydrogen == None or atoms.hydrogen.numAtoms() < 30: # if there is no hydrogens in PDB structure
                addMissingAtoms(pdb_name)
                pdb_name = pdb_name[:-4]+'_addH.pdb'
                atoms = parsePDB(pdb_name)
                LOGGER.info("Lack of hydrogens in the structure. Hydrogens have been added.")
        except: 
            raise ValueError("Missing atoms from the side chains of the structure. Use addMissingAtoms.")
    
        Ligands = [] # Ligands can be more than one
        my_mol = PDBComplex()
        my_mol.load_pdb(pdb_name) # Load the PDB file into PLIP class
        
        if 'select' in kwargs:
            select = kwargs['select']
            LOGGER.info('Selection will be replaced.')
        else:
            select='all not (water or protein or ion)'
            LOGGER.info('Default selection will be used.')

        if 'ignore_ligs' in kwargs:
            ignore_ligs = kwargs['ignore_ligs']
            LOGGER.info('Ignoring list of ligands is uploaded.')
        else:
            ignore_ligs=['NAG','BMA','MAN']
            LOGGER.info('Three molecules will be ignored from analysis: NAG, BMA and MAN.')
        
        select = select+' and not (resname '+' '.join(ignore_ligs)+')'
        ligand_select = atoms.select(select)
        analyzedLigand = []
        LOGGER.info("Detected ligands: ")
        for i in range(len(ligand_select.getResnums())): # It has to be done by each atom
            try:
                ResID = ligand_select.getResnames()[i]
                ChainID = ligand_select.getChids()[i]
                ResNames = ligand_select.getResnums()[i]
                my_bsid = str(ResID)+':'+str(ChainID)+':'+str(ResNames)
                if my_bsid not in analyzedLigand: 
                    LOGGER.info(my_bsid)
                    analyzedLigand.append(my_bsid)
                    my_mol.analyze()
                    my_interactions = my_mol.interaction_sets[my_bsid] # Contains all interaction data      
                    Ligands.append(my_interactions)
            except: 
                LOGGER.info(my_bsid+" not analyzed")

        return Ligands, analyzedLigand

    except:
        LOGGER.info("Install Openbabel and PLIP.")


def listLigandInteractions(PLIP_output):
    """Create a list of interactions from PLIP output created using calcLigandInteractions().
    Results can be displayed in VMD. 
    
    :arg PLIP_output: Results from PLIP for protein-ligand interactions.
    :type PLIP_output: PLIP object obtained from calcLigandInteractions() 
    
    Note that five types of interactions are considered: hydrogen bonds, salt bridges, 
    pi-stacking, cation-pi, hydrophobic and water bridges."""
    
    Inter_list_all = []
    for i in PLIP_output.all_itypes:
        param_inter = [method for method in dir(i) if method.startswith('_') is False]
        
        if str(type(i)).split('.')[-1].strip("'>") == 'hbond':
            Inter_list = ['hbond',i.restype+str(i.resnr), i[0].type+'_'+str(i.d_orig_idx), i.reschain,
                          i.restype+str(i.resnr_l), i[2].type+'_'+str(i.a_orig_idx), i.reschain_l, 
                          i.distance_ad, i.angle, i[0].coords, i[2].coords]
     
        if str(type(i)).split('.')[-1].strip("'>") == 'saltbridge':
            Inter_list = ['saltbridge',i.restype+str(i.resnr), '_'.join(map(str,i.negative.atoms_orig_idx)), i.reschain,
                          i.restype+str(i.resnr_l), '_'.join(map(str,i.positive.atoms_orig_idx)), i.reschain_l, 
                          i.distance, None, i.negative.center, i.positive.center]
                 
        if str(type(i)).split('.')[-1].strip("'>") == 'pistack':
             Inter_list = ['pistack',i.restype+str(i.resnr), '_'.join(map(str,i[0].atoms_orig_idx)), i.reschain,
                          i.restype+str(i.resnr_l), '_'.join(map(str,i[1].atoms_orig_idx)), i.reschain_l, 
                          i.distance, i.angle, i[0].center, i[1].center]           
        
        if str(type(i)).split('.')[-1].strip("'>") == 'pication':
             Inter_list = ['pication',i.restype+str(i.resnr), '_'.join(map(str,i[0].atoms_orig_idx)), i.reschain,
                          i.restype+str(i.resnr_l), '_'.join(map(str,i[1].atoms_orig_idx)), i.reschain_l, 
                          i.distance, None, i[0].center, i[1].center]                       
        
        if str(type(i)).split('.')[-1].strip("'>") == 'hydroph_interaction':
            Inter_list = ['hydroph_interaction',i.restype+str(i.resnr), i[0].type+'_'+str(i[0].idx), i.reschain,
                          i.restype+str(i.resnr_l), i[2].type+'_'+str(i[2].idx), i.reschain_l, 
                          i.distance, None, i[0].coords, i[2].coords]           
             
        if str(type(i)).split('.')[-1].strip("'>") == 'waterbridge':
            water = i.water
            Inter_list = ['waterbridge',i.restype+str(i.resnr), i[0].type+'_'+str(i[0].idx), i.reschain,
                          i.restype+str(i.resnr_l), i[3].type+'_'+str(i[3].idx), i.reschain_l, 
                          [i.distance_aw, i.distance_dw], [i.d_angle, i.w_angle], i[0].coords, i[3].coords, 
                          i.water.coords, i[7].residue.name+'_'+str(i[7].residue.idx)]
                      
        Inter_list_all.append(Inter_list)               
    
    for nr_k,k in enumerate(Inter_list_all):
        LOGGER.info("%3i%22s%10s%26s%4s  <---> %8s%12s%4s%6.1f" % (nr_k,k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]))
    
    return Inter_list_all


def showProteinInteractions_VMD(atoms, interactions, color='red',**kwargs):
    """Save information about protein interactions to a TCL file (filename)
    which can be further use in VMD to display all intercations in a graphical interface
    (in TKConsole: play script_name.tcl).
    Different types of interactions can be saved separately (color can be selected) 
    or all at once for all types of interactions (hydrogen bonds - blue, salt bridges - yellow,
    pi stacking - green, cation-pi - orangem, hydrophobic - silver, and disulfide bonds - black).
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg interactions: List of interactions for protein interactions.
    :type interactions: List of lists
    
    :arg color: color to draw interactions in VMD,
                not used only for single interaction type.
    :type color: str or **None**, by default `red`.
    
    :arg filename: name of TCL file where interactions will be saved.
    :type filename: str
        
    Example (hydrogen bonds for protein only): 
    >>> interactions = calcHydrogenBonds(atoms.protein, distA=3.2, angle=30)
    or all interactions at once:
    >>> interactions = calcProteinInteractions(atoms.protein) """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    if not isinstance(interactions, list):
        raise TypeError('interactions must be a list of interactions.')
    
    try:
        filename = kwargs['filename']
    except:
        filename = atoms.getTitle()+'_interaction.tcl'
    
    tcl_file = open(filename, 'w') 
    
    
    def TCLforSingleInteraction(interaction, color='blue', tcl_file=tcl_file):
        """Creates TCL file for the VMD program based on the interactions
        computed by any function which returns interactions.
        
        :arg interactions: List of interactions for protein interactions.
        :type interactions: List of lists
        
        :arg color: Name of the color which will be used for the visualization of 
                    interactions in VMD
        :type color: str
        
        :arg tcl_file: name of the TCL file which will be saved for visualization                
        :type tcl_file: str """
        
        tcl_file.write('draw color '+color+'\n')
        for nr_i,i in enumerate(interaction):
            try:
                at1 = atoms.select('index '+' '.join([k for k in i[1].split('_') if k.isdigit() ] ))
                at1_atoms = ' '.join(map(str,list(calcCenter(at1.getCoords()))))
                at2 = atoms.select('index '+' '.join([kk for kk in i[4].split('_') if kk.isdigit() ] ))
                at2_atoms = ' '.join(map(str,list(calcCenter(at2.getCoords()))))
                            
                tcl_file.write('draw line {'+at1_atoms+'} {'+at2_atoms+'} style dashed width 4\n')
                
                tcl_file.write('mol color Name\n')
                tcl_file.write('mol representation Licorice 0.100000 12.000000 12.000000\n')
                tcl_file.write('mol selection (resname '+at1.getResnames()[0]+' and resid '+str(at1.getResnums()[0])
                               +' and chain '+at1.getChids()[0]+' and noh) or (resname '+at2.getResnames()[0]+' and resid '
                               +str(at2.getResnums()[0])+' and chain '+at2.getChids()[0]+' and noh)\n')
                tcl_file.write('mol material Opaque\n')
                tcl_file.write('mol addrep 0 \n')
            except: LOGGER.info("There was a problem.")
     
    if len(interactions) == 7:   
        # For all seven types of interactions at once
        # HBs_calculations, SBs_calculations, SameChargeResidues, Pi_stacking, Pi_cation, Hydroph_calculations, Disulfide Bonds
        colors = ['blue', 'yellow', 'red', 'green', 'orange', 'silver', 'black']
        
        for nr_inter,inter in enumerate(interactions):
            TCLforSingleInteraction(inter, color=colors[nr_inter], tcl_file=tcl_file)

    elif len(interactions[0]) == 0 or interactions == []:
        LOGGER.info("Lack of results")
        
    else:
        TCLforSingleInteraction(interactions,color)

    tcl_file.write('draw materials off')
    tcl_file.close()   
    LOGGER.info("TCL file saved")


def showLigandInteraction_VMD(atoms, interactions, **kwargs):
    """Save information from PLIP for ligand-protein interactions in a TCL file
    which can be further used in VMD to display all intercations in a graphical 
    interface (hydrogen bonds - `blue`, salt bridges - `yellow`,
    pi stacking - `green`, cation-pi - `orange`, hydrophobic - `silver` 
    and water bridges - `cyan`).
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg interactions: List of interactions for protein-ligand interactions.
    :type interactions: List of lists
    
    :arg filename: name of TCL file where interactions will be saved.
    :type filename: str

    To obtain protein-ligand interactions:
    >>> calculations = calcLigandInteractions(atoms)
    >>> interactions = listLigandInteractions(calculations) """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    if not isinstance(interactions, list):
        raise TypeError('interactions must be a list of interactions.')
    
    try:
        filename = kwargs['filename']
    except:
        filename = atoms.getTitle()+'_interaction.tcl'
    
    tcl_file = open(filename, 'w') 
    
    if len(interactions[0]) >= 10: 
        dic_color = {'hbond':'blue','pistack':'green','saltbridge':'yellow','pication':'orange',
                     'hydroph_interaction':'silver','waterbridge':'cyan'}
        
        for i in interactions:
            tcl_file.write('draw color '+dic_color[i[0]]+'\n')
            
            if i[0] == 'waterbridge':
                hoh_id = atoms.select('x `'+str(i[11][0])+'` and y `'+str(i[11][1])+'` and z `'+str(i[11][2])+'`').getResnums()[0]
                tcl_file.write('draw line {'+str(' '.join(map(str,i[9])))+'} {'+str(' '.join(map(str,i[11])))+'} style dashed width 4\n')
                tcl_file.write('draw line {'+str(' '.join(map(str,i[10])))+'} {'+str(' '.join(map(str,i[11])))+'} style dashed width 4\n')
                tcl_file.write('mol color Name\n')
                tcl_file.write('mol representation Licorice 0.100000 12.000000 12.000000\n')
                tcl_file.write('mol selection (resname '+i[1][:3]+' and resid '+str(i[1][3:])
                               +' and chain '+i[3]+' and noh) or (resname '+i[4][:3]+' and resid '
                               +str(i[4][3:])+' and chain '+i[6]+' and noh) or (water and resid '+str(hoh_id)+')\n')
                
            else:
                tcl_file.write('draw line {'+str(' '.join(map(str,i[9])))+'} {'+str(' '.join(map(str,i[10])))+'} style dashed width 4\n')
                tcl_file.write('mol color Name\n')
                tcl_file.write('mol representation Licorice 0.100000 12.000000 12.000000\n')
                tcl_file.write('mol selection (resname '+i[1][:3]+' and resid '+str(i[1][3:])
                               +' and chain '+i[3]+' and noh) or (resname '+i[4][:3]+' and resid '
                               +str(i[4][3:])+' and chain '+i[6]+' and noh)\n')
            tcl_file.write('mol material Opaque\n')
            tcl_file.write('mol addrep 0 \n')            

    tcl_file.write('draw materials off')
    tcl_file.close()   
    LOGGER.info("TCL file saved")


class Interactions(object):

    """Class for Interaction analysis of proteins."""

    def __init__(self, title='Unknown'):
        self._title = str(title).strip()
        self._atoms = None
        self._interactions = None
        self._interactions_matrix = None
        self._hbs = None
        self._sbs = None
        self._rib = None
        self._piStack = None
        self._piCat = None
        self._hps = None
        self._dibs = None
        #super(Interactions, self).__init__(name)


    def setTitle(self, title):
        """Set title of the model."""

        self._title = str(title)
        
           
    def calcProteinInteractions(self, atoms, **kwargs):
        """Compute all protein interactions (shown below) using default parameters.
            (1) Hydrogen bonds
            (2) Salt Bridges
            (3) RepulsiveIonicBonding 
            (4) Pi stacking interactions
            (5) Pi-cation interactions
            (6) Hydrophobic interactions
            (7) Disulfide Bonds
        
        :arg atoms: an Atomic object from which residues are selected
        :type atoms: :class:`.Atomic` """

        try:
            coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                        atoms.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be an object '
                                'with `getCoords` method')

        LOGGER.info('Calculating all interations.') 
        HBs_calculations = calcHydrogenBonds(atoms.protein, **kwargs)               #1 in scoring
        SBs_calculations = calcSaltBridges(atoms.protein, **kwargs)                 #2
        SameChargeResidues = calcRepulsiveIonicBonding(atoms.protein, **kwargs)     #3
        Pi_stacking = calcPiStacking(atoms.protein, **kwargs)                       #4
        Pi_cation = calcPiCation(atoms.protein, **kwargs)                           #5
        Hydroph_calculations = calcHydrophobic(atoms.protein, **kwargs)             #6
        Disulfide_Bonds = calcDisulfideBonds(atoms.protein, **kwargs)               #7
        AllInteractions = [HBs_calculations, SBs_calculations, SameChargeResidues, Pi_stacking, 
                            Pi_cation, Hydroph_calculations, Disulfide_Bonds]   
        
        self._atoms = atoms
        self._interactions = AllInteractions
        
        self._hbs = HBs_calculations
        self._sbs = SBs_calculations
        self._rib = SameChargeResidues
        self._piStack = Pi_stacking
        self._piCat = Pi_cation
        self._hps = Hydroph_calculations
        self._dibs = Disulfide_Bonds
        
        return self._interactions

    
    def getAtoms(self):
        """Returns associated atoms. """

        return self._atoms


    def getInteractions(self, **kwargs):
        """Returns the list of all interactions.
        
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            results = [filterInteractions(j, self._atoms, **kwargs) for j in self._interactions]
        else: 
            results = self._interactions
        
        return results


    def getHydrogenBonds(self, **kwargs):
        """Returns the list of hydrogen bonds.

        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """
        
        if len(kwargs) != 0:
            results = filterInteractions(self._hbs, self._atoms, **kwargs)
        else: 
            results = self._hbs

        return results


    def getSaltBridges(self, **kwargs):
        """Returns the list of salt bridges.

        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """ 
        
        if len(kwargs) != 0:
            results = filterInteractions(self._sbs, self._atoms, **kwargs)
        else: 
            results = self._sbs

        return results
        

    def getRepulsiveIonicBonding(self, **kwargs):
        """Returns the list of repulsive ionic bonding.

        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """


        if len(kwargs) != 0:
            results = filterInteractions(self._rib, self._atoms, **kwargs)
        else: 
            results = self._rib

        return results
        

    def getPiStacking(self, **kwargs):
        """Returns the list of Pi-stacking interactions.

        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """
        
        if len(kwargs) != 0:
            results = filterInteractions(self._piStack, self._atoms, **kwargs)
        else: 
            results = self._piStack

        return results


    def getPiCation(self, **kwargs):
        """Returns the list of Pi-cation interactions.
        
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            results = filterInteractions(self._piCat, self._atoms, **kwargs)
        else: 
            results = self._piCat

        return results

        
    def getHydrophobic(self, **kwargs):
        """Returns the list of hydrophobic interactions.

        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """
        
        if len(kwargs) != 0:
            results = filterInteractions(self._hps, self._atoms, **kwargs)
        else: 
            results = self._hps

        return results


    def getDisulfideBonds(self, **kwargs):
        """Returns the list of disulfide bonds.

        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """
        
        if len(kwargs) != 0:
            results = filterInteractions(self._dibs, self._atoms, **kwargs)
        else: 
            results = self._dibs

        return results


    def setNewHydrogenBonds(self, interaction):
        """Replace default calculation of hydrogen bonds by the one provided by user."""

        self._interactions[0] = interaction
        self._hbs = self._interactions[0]    
        LOGGER.info('Hydrogen Bonds are replaced')


    def setNewSaltBridges(self, interaction):
        """Replace default calculation of salt bridges by the one provided by user."""

        self._interactions[1] = interaction
        self._sbs = self._interactions[1]  
        LOGGER.info('Salt Bridges are replaced')


    def setNewRepulsiveIonicBonding(self, interaction):
        """Replace default calculation of repulsive ionic bonding by the one provided by user."""

        self._interactions[2] = interaction
        self._rib = self._interactions[2]   
        LOGGER.info('Repulsive Ionic Bonding are replaced')

        
    def setNewPiStacking(self, interaction):
        """Replace default calculation of pi-stacking interactions by the one provided by user."""

        self._interactions[3] = interaction
        self._piStack = self._interactions[3]   
        LOGGER.info('Pi-Stacking interactions are replaced')


    def setNewPiCation(self, interaction):
        """Replace default calculation of pi-cation interactions by the one provided by user."""

        self._interactions[4] = interaction
        self._piCat = self._interactions[4]   
        LOGGER.info('Pi-Cation interactions are replaced')


    def setNewHydrophobic(self, interaction):
        """Replace default calculation of hydrophobic interactions by the one provided by user."""

        self._interactions[5] = interaction
        self._hps = self._interactions[5]  
        LOGGER.info('Hydrophobic interactions are replaced')


    def setNewDisulfideBonds(self, interaction):
        """Replace default calculation of hydrophobic interactions by the one provided by user."""

        self._interactions[6] = interaction
        self._dibs = self._interactions[6]  
        LOGGER.info('Disulfide bonds are replaced')


    def buildInteractionMatrix(self, **kwargs):
        """Build matrix with protein interactions. Interactions are counted as follows:
            (1) Hydrogen bonds (HBs) +1
            (2) Salt Bridges (SBs) +1 (Salt bridges might be included in hydrogen bonds)
            (3) Repulsive Ionic Bonding (RIB) +1 
            (4) Pi stacking interactions (PiStack) +1
            (5) Pi-cation interactions (PiCat) +1
            (6) Hydrophobic interactions (HPh) +1
            (7) Disulfide bonds (DiBs) +1
                 
        Some type of interactions could be removed from the analysis by replacing value 1 by 0.
        Example: 
        >>> interactions = Interactions()
        >>> atoms = parsePDB(PDBfile).select('protein')
        >>> interactions.calcProteinInteractions(atoms)
        >>> interactions.buildInteractionMatrix(RIB=0, HBs=0, HPh=0)
                 
        :arg HBs: score per single hydrogen bond
        :type HBs: int, float

        :arg SBs: score per single salt bridge
        :type SBs: int, float

        :arg RIB: score per single repulsive ionic bonding
        :type RIB: int, float

        :arg PiStack: score per pi-stacking interaction
        :type PiStack: int, float

        :arg PiCat: score per pi-cation interaction
        :type PiCat: int, float

        :arg HPh: score per hydrophobic interaction
        :type HPh: int, float

        :arg DiBs: score per disulfide bond
        :type DiBs: int, float
        """
        
        atoms = self._atoms   
        interactions = self._interactions
        
        LOGGER.info('Calculating all interactions')
        InteractionsMap = np.zeros([atoms.select('name CA').numAtoms(),atoms.select('name CA').numAtoms()])
        resIDs = list(atoms.select('name CA').getResnums())
        resChIDs = list(atoms.select('name CA').getChids())
        resIDs_with_resChIDs = list(zip(resIDs, resChIDs))
        
        dic_interactions = {'HBs':'Hydrogen Bonds', 'SBs':'Salt Bridges', 'RIB':'Repulsive Ionic Bonding', 
        'PiStack':'Pi-stacking interactions', 'PiCat':'Pi-cation interactions', 'HPh':'Hydrophobic interactions',
        'DiBs':'Disulfide Bonds'}

        if not 'HBs' in kwargs:
            kwargs['HBs'] = 1
        if not 'SBs' in kwargs:
            kwargs['SBs'] = 1
        if not 'RIB' in kwargs:
            kwargs['RIB'] = 1
        if not 'PiStack' in kwargs:
            kwargs['PiStack'] = 1
        if not 'PiCat' in kwargs:
            kwargs['PiCat'] = 1
        if not 'HPh' in kwargs:
            kwargs['HPh'] = 1            
        if not 'DiBs' in kwargs:
            kwargs['DiBs'] = 1            
        
        scoring = [kwargs['HBs'], kwargs['SBs'], kwargs['RIB'], kwargs['PiStack'], kwargs['PiCat'], kwargs['HPh'], kwargs['DiBs']]        

        LOGGER.info('Following scores will be used:')        
        for key,value in kwargs.items(): 
            LOGGER.info('{0} = {1}'.format(dic_interactions[key], value))
        
        for nr_i,i in enumerate(interactions):
            if i != []:
                for ii in i: 
                    m1 = resIDs_with_resChIDs.index((int(ii[0][3:]),ii[2]))
                    m2 = resIDs_with_resChIDs.index((int(ii[3][3:]),ii[5]))
                    InteractionsMap[m1][m2] = InteractionsMap[m2][m1] = InteractionsMap[m1][m2] + scoring[nr_i]
        
        self._interactions_matrix = InteractionsMap
        
        return InteractionsMap


    def showInteractors(self, **kwargs):
        """Display protein residues and their number of potential interactions
        with other residues from protein structure. """
        
        import numpy as np
        import matplotlib
        import matplotlib.pylab as plt        
        from prody.dynamics.plotting import showAtomicLines
        
        if not hasattr(self, '_interactions_matrix') or self._interactions_matrix is None:
            raise ValueError('Please calculate interactions matrix first')

        interaction_matrix = self._interactions_matrix
        atoms = self._atoms 
        
        freq_contacts_residues = np.sum(interaction_matrix, axis=0)
        ResNumb = atoms.select('protein and name CA').getResnums()
        ResName = atoms.select('protein and name CA').getResnames()
        ResChid = atoms.select('protein and name CA').getChids()

        ResList = [ i[0]+str(i[1])+i[2] for i in list(zip(ResName, ResNumb, ResChid)) ]
        
        if SETTINGS['auto_show']:
            matplotlib.rcParams['font.size'] = '20' 
            fig = plt.figure(num=None, figsize=(12,6), facecolor='w')
        show = showAtomicLines(freq_contacts_residues, atoms=atoms.select('name CA'), **kwargs)
        plt.ylabel('Score of interactions')
        plt.xlabel('Residue')
        plt.tight_layout()
        
        if SETTINGS['auto_show']:
            showFigure()
        return show
        
    
    def saveInteractionsPDB(self, **kwargs):
        """Save the number of potential interactions to PDB file in occupancy column.
        
        :arg filename: name of the PDB file which will be saved for visualization,
                     it will contain the results in occupancy column.
        :type filename: str  """
        
        if not hasattr(self, '_interactions_matrix') or self._interactions_matrix is None:
            raise ValueError('Please calculate interactions matrix first.')

        import numpy as np
        interaction_matrix = self._interactions_matrix
        atoms = self._atoms     
        freq_contacts_residues = np.sum(interaction_matrix, axis=0)
        
        try:
            from collections import Counter
            lista_ext = []
            atoms = atoms.select("all and noh")
            aa_counter = Counter(atoms.getResindices())
            calphas = atoms.select('name CA')
            for i in range(calphas.numAtoms()):
                lista_ext.extend(list(aa_counter.values())[i]*[round(freq_contacts_residues[i], 8)])
            
            kw = {'occupancy': lista_ext}
            if 'filename' in kwargs:
                writePDB(kwargs['filename'], atoms, **kw)  
                LOGGER.info('PDB file saved.')
            else:
                writePDB('filename', atoms, **kw)
                LOGGER.info('PDB file saved.')
        except: LOGGER.info('There is a problem.')
        

    def getFrequentInteractors(self, contacts_min=3):
        """Provide a list of residues with the most frequent interactions based 
        on the following interactions:
            (1) Hydrogen bonds (hb)
            (2) Salt Bridges (sb)
            (3) Repulsive Ionic Bonding (rb) 
            (4) Pi stacking interactions (ps)
            (5) Pi-cation interactions (pc)
            (6) Hydrophobic interactions (hp)
            (7) Disulfide bonds (disb)
        
        :arg contacts_min: Minimal number of contacts which residue may form with other residues. 
        :type contacts_min: int, be default 3.  """

        atoms = self._atoms   
        interactions = self._interactions
        
        InteractionsMap = np.empty([atoms.select('name CA').numAtoms(),atoms.select('name CA').numAtoms()], dtype='S256')
        resIDs = list(atoms.select('name CA').getResnums())
        resChIDs = list(atoms.select('name CA').getChids())
        resIDs_with_resChIDs = list(zip(resIDs, resChIDs))
        interaction_type = ['hb','sb','rb','ps','pc','hp','dibs']

        for nr,i in enumerate(interactions):
            if i != []:
                for ii in i: 
                    m1 = resIDs_with_resChIDs.index((int(ii[0][3:]),ii[2]))
                    m2 = resIDs_with_resChIDs.index((int(ii[3][3:]),ii[5]))
                    InteractionsMap[m1][m2] = interaction_type[nr]+':'+ii[0]+ii[2]+'-'+ii[3]+ii[5]
            
        InteractionsMap = InteractionsMap.astype(str)

        ListOfInteractions = [ list(filter(None, InteractionsMap[:,j])) for j in range(len(interactions[0])) ]
        ListOfInteractions = list(filter(lambda x : x != [], ListOfInteractions))
        ListOfInteractions = [k for k in ListOfInteractions if len(k) >= contacts_min ]
        ListOfInteractions_list = [ (i[0].split('-')[-1], [ j.split('-')[0] for j in i]) for i in ListOfInteractions ]
        LOGGER.info('The most frequent interactions between:')
        for res in ListOfInteractions_list:
            LOGGER.info('{0}  <--->  {1}'.format(res[0], '  '.join(res[1])))

        LOGGER.info('Legend: hb-hydrogen bond, sb-salt bridge, rb-repulsive ionic bond, ps-Pi stacking interaction,'
                             'pc-Cation-Pi interaction, hp-hydrophobic interaction, dibs - disulfide bonds')
        
        try:
            from toolz.curried import count
        except ImportError:
            LOGGER.warn('This function requires the module toolz')
            return
        
        LOGGER.info('The biggest number of interactions: {}'.format(max(map(count, ListOfInteractions))))
        
        return ListOfInteractions_list
        

    def showFrequentInteractors(self, cutoff=5, **kwargs):
        """Plots regions with the most frequent interactions.
        
        :arg cutoff: minimal score per residue which will be displayed.
                     If cutoff value is to big, top 30% with the higest values will be returned.
                     Default is 5.
        :type cutoff: int, float

        Nonstandard resiudes can be updated in a following way:
        d = {'CYX': 'X', 'CEA': 'Z'}
        >>> name.showFrequentInteractors(d)  """
        
        if not hasattr(self, '_interactions_matrix') or self._interactions_matrix is None:
            raise ValueError('Please calculate interactions matrix first')

        import numpy as np
        import matplotlib
        import matplotlib.pyplot as plt
        
        atoms = self._atoms
        interaction_matrix = self._interactions_matrix        
        
        aa_dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                       'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HSE': 'H', 'HSD': 'H'}

        for key, value in kwargs.items():
            aa_dict[key] = value

        freq_contacts_residues = np.sum(interaction_matrix, axis=0)
        ResNumb = atoms.select('protein and name CA').getResnums()
        ResName = atoms.select('protein and name CA').getResnames()
        ResChid = atoms.select('protein and name CA').getChids()
        ResList = [ i[0]+str(i[1])+i[2] for i in list(zip(ResName, ResNumb, ResChid)) ]

        all_y = [ aa_dic[i[:3]]+i[3:] for i in  ResList]

        if cutoff > np.max(freq_contacts_residues):
            cutoff = round(np.max(freq_contacts_residues)*0.7)

        y = []
        x = []
        for nr_ii, ii in enumerate(freq_contacts_residues):
            if ii >= cutoff:
                x.append(ii)
                y.append(all_y[nr_ii])

        if SETTINGS['auto_show']:
            matplotlib.rcParams['font.size'] = '20' 
            fig = plt.figure(num=None, figsize=(12,6), facecolor='w')
        
        y_pos = np.arange(len(y))
        show = plt.bar(y_pos, x, align='center', alpha=0.5, color='blue')
        plt.xticks(y_pos, y, rotation=45, fontsize=20)
        plt.ylabel('Number of interactions')
        plt.tight_layout()

        if SETTINGS['auto_show']:
            showFigure()
        return show        


    def showCumulativeInteractionTypes(self, **kwargs):
        
        """Bar plot with the number of potential inetractions per residue."""

        import numpy as np
        import matplotlib
        import matplotlib.pyplot as plt
        from prody.dynamics.plotting import pplot
        
        aa_dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
               'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HSE': 'H', 'HSD': 'H'}

        atoms = self._atoms

        ResNumb = atoms.select('protein and name CA').getResnums()
        ResName = atoms.select('protein and name CA').getResnames()
        ResChid = atoms.select('protein and name CA').getChids()
        ResList = [ i[0]+str(i[1])+i[2] for i in list(zip([ aa_dic[i] for i in ResName ], ResNumb, ResChid)) ]

        matrix_hbs = self.buildInteractionMatrix(HBs=1, SBs=0, RIB=0,PiStack=0,PiCat=0,HPh=0,DiBs=0)
        matrix_sbs = self.buildInteractionMatrix(HBs=0, SBs=1, RIB=0,PiStack=0,PiCat=0,HPh=0,DiBs=0)
        matrix_rib = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=1,PiStack=0,PiCat=0,HPh=0,DiBs=0)
        matrix_pistack = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=1,PiCat=0,HPh=0,DiBs=0)
        matrix_picat = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=0,PiCat=1,HPh=0,DiBs=0)
        matrix_hph = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=0,PiCat=0,HPh=1,DiBs=0)
        matrix_dibs = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=0,PiCat=0,HPh=0,DiBs=1)

        matrix_hbs_sum = np.sum(matrix_hbs, axis=0)
        matrix_sbs_sum = np.sum(matrix_sbs, axis=0)
        matrix_rib_sum = np.sum(matrix_rib, axis=0)
        matrix_pistack_sum = np.sum(matrix_pistack, axis=0)
        matrix_picat_sum = np.sum(matrix_picat, axis=0)
        matrix_hph_sum = np.sum(matrix_hph, axis=0)
        matrix_dibs_sum = np.sum(matrix_dibs, axis=0)

        width = 0.8
        fig, ax = plt.subplots(num=None, figsize=(20,6), facecolor='w')
        matplotlib.rcParams['font.size'] = '24'

        sum_matrix = np.zeros(matrix_hbs_sum.shape)
        pplot(sum_matrix, atoms=atoms.ca)

        ax.bar(ResList, matrix_hbs_sum, width, color = 'blue', bottom = 0, label='HBs')
        sum_matrix += matrix_hbs_sum

        ax.bar(ResList, matrix_sbs_sum, width, color = 'yellow', bottom = sum_matrix, label='SBs')
        sum_matrix += matrix_sbs_sum

        ax.bar(ResList, matrix_hph_sum, width, color = 'silver', bottom = sum_matrix, label='HPh')
        sum_matrix += matrix_hph_sum

        ax.bar(ResList, matrix_rib_sum, width, color = 'red', bottom = sum_matrix, label='RIB')
        sum_matrix += matrix_rib_sum

        ax.bar(ResList, matrix_pistack_sum, width, color = 'green', bottom = sum_matrix, label='PiStack')
        sum_matrix += matrix_pistack_sum

        ax.bar(ResList, matrix_picat_sum, width, color = 'orange', bottom = sum_matrix, label='PiCat')
        sum_matrix += matrix_picat_sum

        ax.bar(ResList, matrix_dibs_sum, width, color = 'black', bottom = sum_matrix, label='DiBs')
        sum_matrix += matrix_dibs_sum

        ax.legend(ncol=7, loc='upper center')
        plt.ylim([0,max(sum_matrix)+3])
        plt.tight_layout()    
        plt.xlabel('Residue')
        plt.ylabel('Number of counts')
       
        return matrix_hbs_sum, matrix_sbs_sum, matrix_rib_sum, matrix_pistack_sum, matrix_picat_sum, matrix_hph_sum, matrix_dibs_sum 
        
class InteractionsTrajectory(object):

    """Class for Interaction analysis of DCD trajectory or multi-model PDB (Ensemble PDB)."""

    def __init__(self, name='Unknown'):
        
        self._atoms = None
        self._traj = None
        self._interactions_traj = None
        self._interactions_nb_traj = None
        self._interactions_matrix_traj = None
        self._hbs_traj = None
        self._sbs_traj = None
        self._rib_traj = None
        self._piStack_traj = None
        self._piCat_traj = None
        self._hps_traj = None
        self._dibs_traj = None


    def calcProteinInteractionsTrajectory(self, atoms, trajectory=None, filename=None, **kwargs):
        """Compute all protein interactions (shown below) for DCD trajectory or multi-model PDB 
            using default parameters.
            (1) Hydrogen bonds
            (2) Salt Bridges
            (3) RepulsiveIonicBonding 
            (4) Pi stacking interactions
            (5) Pi-cation interactions
            (6) Hydrophobic interactions
            (7) Disulfide Bonds
        
        :arg atoms: an Atomic object from which residues are selected
        :type atoms: :class:`.Atomic`
        
        :arg trajectory: trajectory file
        :type trajectory: class:`.Trajectory`

        :arg filename: Name of pkl filename in which interactions will be storage
        :type filename: pkl 
        
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str

        :arg start_frame: index of first frame to read
        :type start_frame: int

        :arg stop_frame: index of last frame to read
        :type stop_frame: int
    
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B' """

        try:
            coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                        atoms.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be an object '
                                'with `getCoords` method')
        
        HBs_all = []
        SBs_all = []
        RIB_all = []
        PiStack_all = []
        PiCat_all = []
        HPh_all = []
        DiBs_all = []

        HBs_nb = []
        SBs_nb = []
        RIB_nb = []
        PiStack_nb = []
        PiCat_nb = []
        HPh_nb = []
        DiBs_nb = []

        start_frame = kwargs.pop('start_frame', 0)
        stop_frame = kwargs.pop('stop_frame', -1)

        if trajectory is not None:
            if isinstance(trajectory, Atomic):
                trajectory = Ensemble(trajectory)
        
            nfi = trajectory._nfi
            trajectory.reset()
            numFrames = trajectory._n_csets

            if stop_frame == -1:
                traj = trajectory[start_frame:]
            else:
                traj = trajectory[start_frame:stop_frame+1]

            for j0, frame0 in enumerate(traj, start=start_frame):
                LOGGER.info('Frame: {0}'.format(j0))
                hydrogen_bonds = calcHydrogenBonds(atoms.protein, **kwargs)
                salt_bridges = calcSaltBridges(atoms.protein, **kwargs)
                RepulsiveIonicBonding = calcRepulsiveIonicBonding(atoms.protein, **kwargs)
                Pi_stacking = calcPiStacking(atoms.protein, **kwargs)
                Pi_cation = calcPiCation(atoms.protein, **kwargs)
                hydrophobic = calcHydrophobic(atoms.protein, **kwargs)
                Disulfide_Bonds = calcDisulfideBonds(atoms.protein, **kwargs)

                HBs_all.append(hydrogen_bonds)
                SBs_all.append(salt_bridges)
                RIB_all.append(RepulsiveIonicBonding)
                PiStack_all.append(Pi_stacking)
                PiCat_all.append(Pi_cation)
                HPh_all.append(hydrophobic)
                DiBs_all.append(Disulfide_Bonds)
            
                HBs_nb.append(len(hydrogen_bonds))
                SBs_nb.append(len(salt_bridges))
                RIB_nb.append(len(RepulsiveIonicBonding))
                PiStack_nb.append(len(Pi_stacking))
                PiCat_nb.append(len(Pi_cation))
                HPh_nb.append(len(hydrophobic))
                DiBs_nb.append(len(Disulfide_Bonds))
      
        else:
            if atoms.numCoordsets() > 1:
                for i in range(len(atoms.getCoordsets()[start_frame:stop_frame])):
                    LOGGER.info('Model: {0}'.format(i+start_frame))
                    atoms.setACSIndex(i+start_frame) 
                    protein = atoms.select('protein')
                    
                    hydrogen_bonds = calcHydrogenBonds(atoms.protein, **kwargs)
                    salt_bridges = calcSaltBridges(atoms.protein, **kwargs)
                    RepulsiveIonicBonding = calcRepulsiveIonicBonding(atoms.protein, **kwargs)
                    Pi_stacking = calcPiStacking(atoms.protein, **kwargs)
                    Pi_cation = calcPiCation(atoms.protein, **kwargs)
                    hydrophobic = calcHydrophobic(atoms.protein, **kwargs)
                    Disulfide_Bonds = calcDisulfideBonds(atoms.protein, **kwargs)

                    HBs_all.append(hydrogen_bonds)
                    SBs_all.append(salt_bridges)
                    RIB_all.append(RepulsiveIonicBonding)
                    PiStack_all.append(Pi_stacking)
                    PiCat_all.append(Pi_cation)
                    HPh_all.append(hydrophobic)
                    DiBs_all.append(Disulfide_Bonds)
            
                    HBs_nb.append(len(hydrogen_bonds))
                    SBs_nb.append(len(salt_bridges))
                    RIB_nb.append(len(RepulsiveIonicBonding))
                    PiStack_nb.append(len(Pi_stacking))
                    PiCat_nb.append(len(Pi_cation))
                    HPh_nb.append(len(hydrophobic))
                    DiBs_nb.append(len(Disulfide_Bonds))
            else:
                LOGGER.info('Include trajectory or use multi-model PDB file.') 
        
        self._atoms = atoms
        self._traj = trajectory
        self._interactions_traj = [HBs_all, SBs_all, RIB_all, PiStack_all, PiCat_all, HPh_all, DiBs_all]
        self._interactions_nb_traj = [HBs_nb, SBs_nb, RIB_nb, PiStack_nb, PiCat_nb, HPh_nb, DiBs_nb]
        self._hbs_traj = HBs_all
        self._sbs_traj = SBs_all  
        self._rib_traj = RIB_all
        self._piStack_traj = PiStack_all
        self._piCat_traj = PiCat_all
        self._hps_traj = HPh_all
        self._dibs_traj = DiBs_all
        
        if filename is not None:
            import pickle
            with open(str(filename)+'.pkl', 'wb') as f:
                pickle.dump(self._interactions_traj, f)  
            LOGGER.info('File with interactions saved.')
            
        return HBs_nb, SBs_nb, RIB_nb, PiStack_nb, PiCat_nb, HPh_nb, DiBs_nb


    def getInteractions(self, **kwargs):
        """Return the list of all interactions.
        
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
            
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            sele_inter = []
            for i in self._interactions_traj:
                for nr_j,j in enumerate(i):
                    sele_inter.append(filterInteractions(i[nr_j], self._atoms, **kwargs))
            results = sele_inter
        else: 
            results = self._interactions_traj
        
        return results


    def getAtoms(self):
        """Returns associated atoms."""

        return self._atoms

    
    def getInteractionsNumber(self):
        """Return the number of interactions in each frame."""
        
        return self._interactions_nb_traj 
    
    
    def getHydrogenBonds(self, **kwargs):
        """Return the list of hydrogen bonds computed from DCD trajectory.
              
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B' """

        if len(kwargs) != 0:
            sele_inter = []
            for nr_i,i in enumerate(self._hbs_traj):
                sele_inter.append(filterInteractions(self._hbs_traj[nr_i], self._atoms, **kwargs))
            results = sele_inter
        else: 
            results = self._hbs_traj
        
        return results

        
    def getSaltBridges(self, **kwargs):
        """Return the list of salt bridges computed from DCD trajectory.
              
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B' """
        
        if len(kwargs) != 0:
            sele_inter = []
            for nr_i,i in enumerate(self._sbs_traj):
                sele_inter.append(filterInteractions(self._sbs_traj[nr_i], self._atoms, **kwargs))
            results = sele_inter
        else:
            results = self._sbs_traj
            
        return results
            

    def getRepulsiveIonicBonding(self, **kwargs):
        """Return the list of repulsive ionic bonding computed from DCD trajectory.
              
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            sele_inter = []
            for nr_i,i in enumerate(self._rib_traj):
                sele_inter.append(filterInteractions(self._rib_traj[nr_i], self._atoms, **kwargs))
            results = sele_inter
        else:
            results = self._rib_traj
        
        return results
        

    def getPiStacking(self, **kwargs):
        """Return the list of Pi-stacking interactions computed from DCD trajectory.
              
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            sele_inter = []
            for nr_i,i in enumerate(self._piStack_traj):
                sele_inter.append(filterInteractions(self._piStack_traj[nr_i], self._atoms, **kwargs))
            results = sele_inter
        else:
            results =  self._piStack_traj
            
        return results


    def getPiCation(self, **kwargs):
        """Return the list of Pi-cation interactions computed from DCD trajectory.
              
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            sele_inter = []
            for nr_i,i in enumerate(self._piCat_traj):
                sele_inter.append(filterInteractions(self._piCat_traj[nr_i], self._atoms, **kwargs))
            results = sele_inter
        else: 
            results = self._piCat_traj

        return results
        

    def getHydrophobic(self, **kwargs):
        """Return the list of hydrophobic interactions computed from DCD trajectory.
              
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            sele_inter = []
            for nr_i,i in enumerate(self._hps_traj):
                sele_inter.append(filterInteractions(self._hps_traj[nr_i], self._atoms, **kwargs))
            results = sele_inter
        else:
            results = self._hps_traj
            
        return results


    def getDisulfideBonds(self, **kwargs):
        """Return the list of disulfide bonds computed from DCD trajectory.
              
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str
        
        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """

        if len(kwargs) != 0:
            sele_inter = []
            for nr_i,i in enumerate(self._dibs_traj):
                sele_inter.append(filterInteractions(self._dibs_traj[nr_i], self._atoms, **kwargs))
            results = sele_inter
        else:
            results = self._dibs_traj
            
        return results


    def setNewHydrogenBondsTrajectory(self, interaction):
        """Replace default calculation of hydrogen bonds by the one provided by user."""

        self._interactions_traj[0] = interaction
        self._hbs_traj = self._interactions_traj[0]    
        self._interactions_nb_traj[0] = [ len(i) for i in interaction ]
        LOGGER.info('Hydrogen Bonds are replaced')


    def setNewSaltBridgesTrajectory(self, interaction):
        """Replace default calculation of salt bridges by the one provided by user."""

        self._interactions_traj[1] = interaction
        self._sbs_traj = self._interactions_traj[1]  
        self._interactions_nb_traj[1] = [ len(i) for i in interaction ]
        LOGGER.info('Salt Bridges are replaced')


    def setNewRepulsiveIonicBondingTrajectory(self, interaction):
        """Replace default calculation of repulsive ionic bonding by the one provided by user. """

        self._interactions_traj[2] = interaction
        self._rib_traj = self._interactions_traj[2]   
        self._interactions_nb_traj[2] = [ len(i) for i in interaction ]
        LOGGER.info('Repulsive Ionic Bonding are replaced')

        
    def setNewPiStackingTrajectory(self, interaction):
        """Replace default calculation of pi-stacking interactions by the one provided by user."""

        self._interactions_traj[3] = interaction
        self._piStack_traj = self._interactions_traj[3]   
        self._interactions_nb_traj[3] = [ len(i) for i in interaction ]
        LOGGER.info('Pi-Stacking interactions are replaced')


    def setNewPiCationTrajectory(self, interaction):
        """Replace default calculation of pi-cation interactions by the one provided by user."""

        self._interactions_traj[4] = interaction
        self._piCat_traj = self._interactions_traj[4]   
        self._interactions_nb_traj[4] = [ len(i) for i in interaction ]
        LOGGER.info('Pi-Cation interactions are replaced')


    def setNewHydrophobicTrajectory(self, interaction):
        """Replace default calculation of hydrophobic interactions by the one provided by user."""

        self._interactions_traj[5] = interaction
        self._hps_traj = self._interactions_traj[5]  
        self._interactions_nb_traj[5] = [ len(i) for i in interaction ]
        LOGGER.info('Hydrophobic interactions are replaced')


    def setNewDisulfideBondsTrajectory(self, interaction):
        """Replace default calculation of disulfide bonds by the one provided by user."""

        self._interactions_traj[6] = interaction
        self._dibs_traj = self._interactions_traj[6]  
        self._interactions_nb_traj[6] = [ len(i) for i in interaction ]
        LOGGER.info('Disulfide bonds are replaced')

    
    def parseInteractions(self, filename):
        """Import interactions from analysis of trajectory which was saved via
        calcProteinInteractionsTrajectory().
        
        :arg filename: Name of pkl file in which interactions will be storage
        :type filename: pkl """
        
        import pickle
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        
        self._interactions_traj = data
        self._interactions_nb_traj = [[len(sublist) if sublist else 0 for sublist in sublist] for sublist in data]
        self._hbs_traj = data[0]
        self._sbs_traj = data[1]
        self._rib_traj = data[2]
        self._piStack_traj = data[3]
        self._piCat_traj = data[4]
        self._hps_traj = data[5]
        self._dibs_traj = data[6]
        
        return data

    
    def getTimeInteractions(self, filename=None, **kwargs):    
        """Return a bar plots with the number of interactions per each frame.
        
        :arg filename: PNG file name
        :type filename: str """
        
        HBs = self._interactions_nb_traj[0]
        SBs = self._interactions_nb_traj[1]
        RIB = self._interactions_nb_traj[2]
        PiStack = self._interactions_nb_traj[3]
        PiCat = self._interactions_nb_traj[4]
        HPh = self._interactions_nb_traj[5]
        DiBs = self._interactions_nb_traj[6]
        
        import numpy as np
        import matplotlib
        import matplotlib.pyplot as plt
        matplotlib.rcParams['font.size'] = '20' 

        fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, num=None, figsize=(12,8), facecolor='w', sharex='all', **kwargs)
        hspace = 0.1
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.35)
        ax1.bar(np.arange(len(HBs)), HBs, color='deepskyblue')
        ax2.bar(np.arange(len(SBs)),SBs, color='yellow')
        ax3.bar(np.arange(len(HPh)), HPh, color='silver')
        ax4.bar(np.arange(len(PiStack)), PiStack, color='lightgreen')
        ax5.bar(np.arange(len(PiCat)), PiCat, color='orange')
        ax6.bar(np.arange(len(RIB)), RIB, color='red')
        ax7.bar(np.arange(len(DiBs)), DiBs, color='black')

        ax1.plot(HBs, 'k:')
        ax2.plot(SBs, 'k:')
        ax3.plot(HPh, 'k:')
        ax4.plot(PiStack, 'k:')
        ax5.plot(PiCat, 'k:')
        ax6.plot(RIB, 'k:')
        ax7.plot(DiBs, 'k:')

        plt.xlabel('Frame')
        
        if filename is not None:
            plt.savefig(filename+'.png', dpi=300)
        plt.show()
        
        return HBs, SBs, HPh, PiStack, PiCat, HPh, DiBs

