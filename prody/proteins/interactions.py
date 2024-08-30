# -*- coding: utf-8 -*-

"""This module called InSty defines functions for calculating different types of interactions 
in protein structure, between proteins or between protein and ligand.
The following interactions are available for protein interactions:
        (1) Hydrogen bonds
        (2) Salt Bridges
        (3) Repulsive Ionic Bonding 
        (4) Pi stacking interactions
        (5) Pi-cation interactions
        (6) Hydrophobic interactions
        (7) Disulfide Bonds
"""

__author__ = 'Karolina Mikulska-Ruminska'
__credits__ = ['James Krieger', 'Karolina Mikulska-Ruminska']
__email__ = ['karolamik@fizyka.umk.pl', 'jamesmkrieger@gmail.com']


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

import multiprocessing

__all__ = ['calcHydrogenBonds', 'calcChHydrogenBonds', 'calcSaltBridges',
           'calcRepulsiveIonicBonding', 'calcPiStacking', 'calcPiCation',
           'calcHydrophobic', 'calcDisulfideBonds', 'calcMetalInteractions',
           'calcHydrogenBondsTrajectory', 'calcSaltBridgesTrajectory',
           'calcRepulsiveIonicBondingTrajectory', 'calcPiStackingTrajectory', 
           'calcPiCationTrajectory', 'calcHydrophobicTrajectory', 'calcDisulfideBondsTrajectory',
           'calcProteinInteractions', 'calcStatisticsInteractions', 'calcDistribution',
           'calcSASA', 'calcVolume','compareInteractions', 'showInteractionsGraph',
           'calcLigandInteractions', 'listLigandInteractions', 
           'showProteinInteractions_VMD', 'showLigandInteraction_VMD', 
           'calcHydrogenBondsTrajectory', 'calcHydrophobicOverlapingAreas',
           'Interactions', 'InteractionsTrajectory', 'LigandInteractionsTrajectory',
           'calcSminaBindingAffinity', 'calcSminaPerAtomInteractions', 'calcSminaTermValues',
           'showSminaTermValues', 'showPairEnergy']


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
    
    coordinates = getCoords(atoms)
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
    """Remove duplicates from interactions."""
    ls=[]
    newList = []
    for no, i in enumerate(list_of_interactions):
       i = sorted(list(array(i).astype(str)))
       if i not in ls:
           ls.append(i)
           newList.append(list_of_interactions[no])
    return newList


def get_permutation_from_dic(dictionary, key):
    """Check permutations of residue pairs in a dictionary
    format: key=('VAL8', 'LEU89')
    ('VAL8', 'LEU89') and ('LEU89', 'VAL8') will be check and corresponding
    value will be returned."""
    
    if key in dictionary:
        return dictionary[key]

    reversed_key = tuple(reversed(key))
    if reversed_key in dictionary:
        return dictionary[reversed_key]

    return 0


def filterInteractions(list_of_interactions, atoms, **kwargs):
    """Return interactions based on selection."""
    
    if 'selection' in kwargs:
        if 'selection2' in kwargs:
            if not 'chid' in kwargs['selection'] and not 'chain' in kwargs['selection']:
                LOGGER.warn('selection does not include chid or chain, so no filtering is performed')
                return list_of_interactions

            if not 'chid' in kwargs['selection2'] and not 'chain' in kwargs['selection2']:
                LOGGER.warn('selection2 does not include chid or chain, so no filtering is performed')
                return list_of_interactions

            ch1 = kwargs['selection'].split()[-1] 
            ch2 = kwargs['selection2'].split()[-1] 
            final = [i for i in list_of_interactions if (i[2] == ch1 and i[5] == ch2) or (i[5] == ch1 and i[2] == ch2)]
        else:
            p = atoms.select('same residue as protein within 10 of ('+kwargs['selection']+')')
            if p is None:
                LOGGER.warn('selection did not work, so no filtering is performed')
                return list_of_interactions

            x = p.select(kwargs['selection']).getResnames()
            y = p.select(kwargs['selection']).getResnums()
            listOfselection = np.unique(list(map(lambda x, y: x + str(y), x, y)))
            final = [i for i in list_of_interactions if i[0] in listOfselection or i[3] in listOfselection]
    elif 'selection2' in kwargs:
        LOGGER.warn('selection2 by itself is ignored')
    else:
        final = list_of_interactions
    return final


def get_energy(pair, source):
    """Return energies based on the pairs of interacting residues (without distance criteria)
    Taking information from tabulated_energies.txt file"""

    import numpy as np
    import importlib.resources as pkg_resources    
    
    try:
        # Python 3
        with pkg_resources.path('prody.proteins', 'tabulated_energies.txt') as file_path:
            data = np.loadtxt(file_path, skiprows=1, dtype=str)
    except: 
        # Python 2.7
        import pkg_resources
        file_path = pkg_resources.resource_filename('prody.proteins', 'tabulated_energies.txt')
        with open(file_path) as f:
            data = np.loadtxt(f, skiprows=1, dtype=str)

    
    sources = ["IB_nosolv", "IB_solv", "CS"]
    aa_pairs = []
    
    for row in data:
        aa_pairs.append(row[0]+row[1])
    
    lookup = pair[0]+pair[1]
    
    return data[np.where(np.array(aa_pairs)==lookup)[0]][0][2:][np.where(np.array(sources)==source)][0]


def showPairEnergy(data, **kwargs):
    """Return energies when a list of interactions is given. Energies will be added to each pair of residues 
    at the last position in the list. Energy is based on the residue types and not on the distances.
    The unit of energy is kcal/mol. The energies defined as 'IB_nosolv', 'IB_solv' are taken from XX and 
    'CS' from YY.
    
    :arg data: list with interactions from calcHydrogenBonds() or other types
    :type data: list
    
    :arg energy_list_type: name of the list with energies 
                            default is 'IB_solv'
    :type energy_list_type: 'IB_nosolv', 'IB_solv', 'CS'
    """
    
    if not isinstance(data, list):
        raise TypeError('list_of_interactions must be a list of interactions.')

    energy_list_type = kwargs.pop('energy_list_type', 'IB_solv')
    
    for i in data:
        energy = get_energy([i[0][:3], i[3][:3]], energy_list_type)
        i.append(float(energy))
        
    return data



def calcHydrophobicOverlapingAreas(atoms, **kwargs):
    """Provide information about hydrophobic contacts between pairs of residues based on 
    the regsurf program. To use this function compiled hpb.so is needed.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg selection: selection string of hydrophobic residues
    :type selection: str
    
    :arg hpb_cutoff: cutoff for hydrophobic overlapping area values
        default is 0.0
    :type hpb_cutoff: float, int
    
    :arg cumulative_values: sum of results for pairs of residues or
        for single residues without distinguishing pairs,
        default is None
    :type cumulative_values: 'pairs' or 'single_residues'
    """
    
    imported_hpb = False

    try:
        import prody.proteins.hpb as hpb
        imported_hpb = True
    except ImportError:
        try:
            import hpb
            imported_hpb = True
        except ImportError:
            raise ImportError('Please provide hpb.so file.')

    if imported_hpb:
        selection = kwargs.pop('selection', 'protein and noh')
        hpb_cutoff = kwargs.pop('hpb_cutoff', 0.0)
        cumulative_values = kwargs.pop('cumulative_values', None)
        sele = atoms.select(selection)
        lB = sele.getCoords().tolist()
        
        x = sele.getResnames()
        y = sele.getResnums()
        z = sele.getNames()
        w = sele.getIndices()
        ch = sele.getChids()
        lA = [ [x[i] + str(y[i]), z[i] +'_'+ str(w[i]), ch[i]] for i in range(len(x))]
            
        output = hpb.hpb((lB,lA))
        LOGGER.info("Hydrophobic Overlapping Areas are computed.")
        output_final = [i for i in output if i[-1] >= hpb_cutoff]
        
        if cumulative_values == None:            
            return output_final
            
        if cumulative_values == 'pairs':
            sums = {}
            for item in output_final:
                k = (item[0]+item[2], item[3]+item[5])
                v = item[-1]
                if k in sums:
                    sums[k] += v
                else:
                    sums[k] = v
            return sums 
            
        if cumulative_values == 'single_residues':
            sums = {}
            for j in output_final:
                k = j[0]+j[2]
                w = j[-1]
                sums[k] = sums.get(k, 0) + w
            return sums
    

def calcSASA(atoms, **kwargs):
    """Provide information about solvent accessible surface area (SASA) based on 
    the sasa program. To use this function compiled hpb.so is needed.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg selection: selection string
        by default all the protein structure is used
    :type selection: str
    
    :arg sasa_cutoff: cutoff for SASA values
        default is 0.0
    :type sasa_cutoff: float, int

    :arg resnames: residues name included
        default is False
    :type resnames: bool    
    """
    
    imported_hpb = False

    try:
        import prody.proteins.hpb as hpb
        imported_hpb = True
    except ImportError:
        try:
            import hpb
            imported_hpb = True
        except ImportError:
            raise ImportError('Please provide hpb.so file.')

    if imported_hpb:
        selection = kwargs.pop('selection', 'protein and noh')
        resnames = kwargs.pop('resnames', False)
        sasa_cutoff = kwargs.pop('sasa_cutoff', 0.0)
        
        sele = atoms.select(selection)
        lB = sele.getCoords().tolist()
        
        x = sele.getResnames()
        y = sele.getResnums()
        z = sele.getNames()
        w = sele.getIndices()
        
        ch = sele.getChids()
        lA = [ [x[i] + str(y[i]), z[i] +'_'+ str(w[i]), ch[i]] for i in range(len(x))]
            
        output = hpb.sasa((lB,lA))
        LOGGER.info("Solvent Accessible Surface Area (SASA) is computed.")
        output_final = [i for i in output if i[-1] >= sasa_cutoff]
        
        if resnames == True:            
            return output_final
        else:
            return [ float(i[-1]) for i in output_final ]


def calcVolume(atoms, **kwargs):
    """Provide information about volume for each residue/molecule/chain
    or other selection". To use this function compiled hpb.so is needed.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg selection: selection string
        by default all the protein structure is used
    :type selection: str
    
    :arg volume_cutoff: cutoff for volume
        default is 0.0 to include all the results
    :type sasa_volume: float, int

    :arg split_residues: it will provide values for each residue
        default is False
    :type split_residues: bool

    :arg resnames: residues name included
        default is False
        True - will give residue names and values for each residue
        False - will give only the values for each residues
    :type resnames: bool
    """
    
    imported_hpb = False

    try:
        import prody.proteins.hpb as hpb
        imported_hpb = True
    except ImportError:
        try:
            import hpb
            imported_hpb = True
        except ImportError:
            raise ImportError('Please provide hpb.so file.')

    if imported_hpb:
        selection = kwargs.pop('selection', 'protein and noh')
        resnames = kwargs.pop('resnames', False)
        volume_cutoff = kwargs.pop('volume_cutoff', 0.0)
        split_residues = kwargs.pop('split_residues', False)
        
        sele = atoms.select(selection)
        lB = sele.getCoords().tolist()
        
        x = sele.getResnames()
        y = sele.getResnums()
        z = sele.getNames()
        w = sele.getIndices()
        
        ch = sele.getChids()
        lA = [ [x[i] + str(y[i]), z[i] +'_'+ str(w[i]), ch[i]] for i in range(len(x))]
            
        output = hpb.volume((lB,lA))
        LOGGER.info("Volume is computed.")
        output_final = [i for i in output if i[-1] >= volume_cutoff]
        
        if resnames == True:            
            return output_final

        if split_residues == True and resnames == False:
            return [ float(i[-1]) for i in output_final ]
            
        else:
            return sum( [float(i[-1]) for i in output_final] )


def calcHydrogenBonds(atoms, **kwargs):
    """Compute hydrogen bonds for proteins and other molecules.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distDA: non-zero value, maximal distance between donor and acceptor.
        default is 3.5
        distA also works
    :type distDA: int, float
    
    :arg angleDHA: non-zero value, maximal (180 - D-H-A angle) (donor, hydrogen, acceptor).
        default is 40.
        angle also works
    :type angleDHA: int, float
    
    :arg seq_cutoff_HB: non-zero value, interactions will be found between atoms with index differences
        that are higher than seq_cutoff_HB. 
        default is 25 atoms.
        seq_cutoff also works
    :type seq_cutoff_HB: int

    :arg donors: which atoms to count as donors 
        default is ['N', 'O', 'S', 'F']
    :type donors: list

    :arg acceptors: which atoms to count as acceptors 
        default is ['N', 'O', 'S', 'F']
    :type acceptors: list 

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
    
    distDA = kwargs.pop('distDA', 3.5)
    distA = kwargs.pop('distA', distDA)

    angleDHA = kwargs.pop('angleDA', 40)
    angle = kwargs.pop('angle', angleDHA)
    
    seq_cutoff_HB = kwargs.pop('seq_cutoff_HB', 25)
    seq_cutoff = kwargs.pop('seq_cutoff', seq_cutoff_HB)
    
    donors = kwargs.get('donors', ['N', 'O', 'S', 'F'])
    acceptors = kwargs.get('acceptors', ['N', 'O', 'S', 'F'])
    
    if atoms.hydrogen == None or atoms.hydrogen.numAtoms() < 10:
        LOGGER.info("Provide structure with hydrogens or install Openbabel to add missing hydrogens using addMissingAtoms(pdb_name) first.")
    
    contacts = findNeighbors(atoms.heavy, distA)
    short_contacts = cleanNumbers(contacts)
    pairList = [] # list with Donor-Hydrogen-Acceptor(indices)-distance-Angle
    
    LOGGER.info('Calculating hydrogen bonds.')
    hydrogens = atoms.hydrogen
    for nr_i,i in enumerate(short_contacts):
        # Removing those close contacts which are between neighbour atoms
        if i[1] - seq_cutoff < i[0] < i[1] + seq_cutoff:
            continue
        
        if (i[2][0] in donors and i[3][0] in acceptors) or (i[2][0] in acceptors and i[3][0] in donors): # First letter is checked
            listOfHydrogens1 = cleanNumbers(findNeighbors(hydrogens, 1.4, atoms.select('index '+str(i[0]))))
            listOfHydrogens2 = cleanNumbers(findNeighbors(hydrogens, 1.4, atoms.select('index '+str(i[1]))))
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
    ag = atoms.getAtomGroup()
    resnames = ag.getResnames()
    resnums = ag.getResnums()
    names = ag.getNames()
    chids = ag.getChids()
    for k in pairList:
        if 180-angle < float(k[-1]) < 180 and float(k[-2]) < distA:
            aa_donor = resnames[k[0]]+str(resnums[k[0]])
            aa_donor_atom = names[k[0]]+'_'+str(k[0])
            aa_donor_chain = chids[k[0]]
            aa_acceptor = resnames[k[2]]+str(resnums[k[2]])
            aa_acceptor_atom = names[k[2]]+'_'+str(k[2])
            aa_acceptor_chain = chids[k[2]]
            
            HBs_list.append([str(aa_donor), str(aa_donor_atom), str(aa_donor_chain), str(aa_acceptor), str(aa_acceptor_atom), 
                             str(aa_acceptor_chain), np.round(float(k[-2]),4), np.round(180.0-float(k[-1]),4)])
    
    HBs_list = sorted(HBs_list, key=lambda x : x[-2])
    HBs_list_final = removeDuplicates(HBs_list)
    
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
        default is 3.5
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
        HBS_calculations = calcHydrogenBonds(atoms, distA=distA, angle=angle, seq_cutoff=seq_cutoff)
    
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
    
    :arg distSB: non-zero value, maximal distance between center of masses 
        of N and O atoms of negatively and positevely charged residues.
        default is 5.
        distA also works
    :type distSB: int, float

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
    
    distSB = kwargs.pop('distSB', 5.)
    distA = kwargs.pop('distA', distSB)

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
                                                  sele2_single.getChids()[0], round(distance,4)])
    
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
    
    :arg distRB: non-zero value, maximal distance between center of masses 
            between N-N or O-O atoms of residues.
            default is 4.5.
            distA works too
    :type distRB: int, float

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
    
    distRB = kwargs.pop('distRB', 4.5)
    distA = kwargs.pop('distA', distRB)

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
                                                  sele2_single.getChids()[0], round(distance,4)])
    
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


def calcPiStacking_once(sele1, sele2, distA, angle_min, angle_max, data):
    """Helper function to be used by calcPiStacking"""
    a1, b1, c1, a2, b2, c2 = calcPlane(sele1)[:3]+calcPlane(sele2)[:3]
    RingRing_angle = calcAngleBetweenPlanes(a1, b1, c1, a2, b2, c2) # plane is computed based on 3 points of rings
    RingRing_distance = calcDistance(calcCenter(sele1), calcCenter(sele2))
    if RingRing_distance < distA and angle_min < RingRing_angle < angle_max:
        return data+[round(RingRing_distance,4), round(RingRing_angle,4)]


def calcPiStacking(atoms, **kwargs):
    """Finds π–π stacking interactions (between aromatic rings).
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distPS: non-zero value, maximal distance between center of masses 
                of residues aromatic rings.
                default is 5.
                distA works too
    :type distPS: int, float
    
    :arg angle_min_PS: minimal angle between aromatic rings.
        default is 0.
        angle_min works too
    :type angle_min_PS: int, float

    :arg angle_max_PS: maximal angle between rings.
        default is 360.
        angle_max works too
    :type angle_max_PS: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg non_standard_PS: dictionary of non-standard residue in the protein structure
                        that need to be included in calculations
                        non_standard works too
    :type non_standard_PS: dict

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
                'HIS':'noh and not backbone and not name CB',
                'HSE':'noh and not backbone and not name CB',
                'HSD':'noh and not backbone and not name CB'}
    
    distPS = kwargs.pop('distPS', 5.0)
    distA = kwargs.pop('distA', distPS)

    angle_min_RB = kwargs.pop('angle_min_RB', 0)
    angle_min = kwargs.pop('angle_min', angle_min_RB)

    angle_max_RB = kwargs.pop('angle_max_RB', 360)
    angle_max = kwargs.pop('angle_max', angle_max_RB)
    
    non_standard_PS = kwargs.get('non_standard_PS', {})
    non_standard = kwargs.get('non_standard', non_standard_PS)

    for key, value in non_standard.items():
        aromatic_dic[key] = value
    
    atoms_cylic = atoms.select('resname TRP PHE TYR HIS HSE HSD')
    if atoms_cylic is None:
        return []
    
    aromatic_resids = list(set(zip(atoms_cylic.getResnums(), atoms_cylic.getChids())))

    LOGGER.info('Calculating Pi stacking interactions.')
    PiStack_calculations = []
    items = []
    for i in aromatic_resids:
        for j in aromatic_resids:
            if i != j: 
                sele1_full = atoms.select('resid '+str(i[0])+' and chain '+i[1])
                sele1_name = sele1_full.getResnames()
                sele1 = sele1_full.select(aromatic_dic[sele1_name[0]])
                
                sele2_full = atoms.select('resid '+str(j[0])+' and chain '+j[1])
                sele2_name = sele2_full.getResnames()
                sele2 = sele2_full.select(aromatic_dic[sele2_name[0]])

                if sele1 != None and sele2 != None:
                    items.append([sele1.getCoords(), sele2.getCoords(), distA, angle_min, angle_max,
                                  [str(sele1.getResnames()[0])+str(sele1.getResnums()[0]), '_'.join(map(str,sele1.getIndices())), str(sele1.getChids()[0]),
                                   str(sele2.getResnames()[0])+str(sele2.getResnums()[0]), '_'.join(map(str,sele2.getIndices())), str(sele2.getChids()[0])]])

    # create a process pool that uses all cpus
    with multiprocessing.Pool() as pool:
        # call the function for each item in parallel with multiple arguments
        for result in pool.starmap(calcPiStacking_once, items):
            if result is not None:
                PiStack_calculations.append(result)
    
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
    
    :arg distPC: non-zero value, maximal distance between center of masses 
                of aromatic ring and positively charge group.
                default is 5.
                distA works too
    :type distPC: int, float

    :arg selection: selection string
    :type selection: str
    
    :arg selection2: selection string
    :type selection2: str

    :arg non_standard_PC: dictionary of non-standard residue in the protein structure
                        that need to be included in calculations
                        non_standard also works
    :type non_standard_PC: dict

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
                'HIS':'noh and not backbone and not name CB',
                'HSE':'noh and not backbone and not name CB',
                'HSD':'noh and not backbone and not name CB'}
    
    distPC = kwargs.pop('distPC', 5.0)
    distA = kwargs.pop('distA', distPC)
    
    non_standard_PC = kwargs.get('non_standard_PC', {})
    non_standard = kwargs.get('non_standard', non_standard_PC)
    
    for key, value in non_standard.items():
        aromatic_dic[key] = value
        
    atoms_cylic = atoms.select('resname TRP PHE TYR HIS HSE HSD')
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
                                                  str(sele2_single.getChids()[0]), round(RingCation_distance,4)])
    
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
    (ALA, ILE, LEU, MET, PHE, TRP, VAL, CG of ARG, and CG and CD of LYS).
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distHPh: non-zero value, maximal distance between atoms of hydrophobic residues.
        default is 4.5.
        distA works too
    :type distHPh: int, float
    
    :arg non_standard_Hph: dictionary of non-standard residue in the protein structure
                        that need to be included in calculations
                        non_standard works too
    :type non_standard_Hph: dict

    :arg zerosHPh: zero values of hydrophobic overlapping areas included
        default is False
    :type zerosHPh: bool

    Last value in the output corresponds to the total hydrophobic overlapping area for two residues
    not only for the atoms that are included in the list. Atoms that which are listed are the closest
    between two residues and they will be inluded to draw the line in VMD_.
    
    Selection:
    If we want to select interactions for the particular residue or group of residues: 
        selection='chain A and resid 1 to 50'
    If we want to study chain-chain interactions:
        selection='chain A', selection2='chain B'
    
    Additional selection can be added as shown below (with selection that includes 
    only hydrophobic part): 
        >>> calcHydrophobic(atoms, non_standard_Hph={'XLE'='noh and not backbone', 
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
    
    distHph = kwargs.pop('distHph', 4.5)
    distA = kwargs.pop('distA', distHph)
    zerosHPh = kwargs.pop('zerosHPh', False)

    # Nonstandard residues:    
    hydrophobic_dic = {'ALA': 'noh and not backbone', 'VAL': 'noh and not (backbone or name CB)',
    'ILE': 'noh and not (backbone or name CB)', 'LEU': 'noh and not (backbone or name CB)',
    'MET': 'noh and not (backbone or name CB)', 'PHE': 'noh and not (backbone or name CB)',
    'TYR': 'noh and not (backbone or name CB)', 'TRP': 'noh and not (backbone or name CB)',
    'ARG': 'name CG', 'LYS': 'name CG CD'}

    non_standard_Hph = kwargs.get('non_standard_Hph', {})
    non_standard = kwargs.get('non_standard', non_standard_Hph)

    for key, value in non_standard.items():
        hydrophobic_dic[key] = value

    
    Hydrophobic_list = []  
    # All residues, also non-standard will be included in the selection:
    residue_list = list(hydrophobic_dic.keys())
    atoms_hydrophobic = atoms.select('resname '+' '.join(residue_list))
    hydrophobic_resids = list(set(zip(atoms_hydrophobic.getResnums(), atoms_hydrophobic.getChids())))

    if atoms.aromatic is None:
        return []
    
    aromatic_nr = list(set(zip(atoms.aromatic.getResnums(),atoms.aromatic.getChids())))   
    aromatic = list(set(atoms.aromatic.getResnames()))
    
    # Computing hydrophobic overlapping areas for pairs of residues:
    try:
        hpb_overlaping_results = calcHydrophobicOverlapingAreas(atoms_hydrophobic, cumulative_values='pairs')
    except: 
        LOGGER.info('Please provide hpb.so file to obtain additional data.')
    
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
                # avoid double counting pi stacking and don't include same residue interactions
                sele2_filter = sele2.select('all and not (resname TYR PHE TRP or resid '+str(i[0])+' and chain '+i[1]+')')
                if sele2_filter != None:
                    listOfAtomToCompare = cleanNumbers(findNeighbors(sele1, distA, sele2_filter))
                
            elif sele1_name[0] not in aromatic and i in sele2_nr:
                # don't include same residue interactions but don't worry about double counting pi stacking
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
                    residue1 = sele1_new.getResnames()[0]+str(sele1_new.getResnums()[0]) 
                    residue2 = sele2_new.getResnames()[0]+str(sele2_new.getResnums()[0])
                    try:
                        Hydrophobic_calculations.append([residue1, 
                                                    minDistancePair[2]+'_'+str(minDistancePair[0]), sele1_new.getChids()[0],
                                                    residue2, 
                                                    minDistancePair[3]+'_'+str(minDistancePair[1]), sele2_new.getChids()[0],
                                                    round(minDistancePair[-1],4),
                                                    round(get_permutation_from_dic(hpb_overlaping_results,(residue1+sele1_new.getChids()[0],
                                                    residue2+sele2_new.getChids()[0])),4)])
                    except:
                        Hydrophobic_calculations.append([residue1, 
                                                    minDistancePair[2]+'_'+str(minDistancePair[0]), sele1_new.getChids()[0],
                                                    residue2, 
                                                    minDistancePair[3]+'_'+str(minDistancePair[1]), sele2_new.getChids()[0],
                                                    round(minDistancePair[-1],4)])                         
    
    selection = kwargs.get('selection', None)
    selection2 = kwargs.get('selection2', None) 
    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
        
    imported_hpb = False

    try:
        import prody.proteins.hpb as hpb
        imported_hpb = True
    except ImportError:
        try:
            import hpb
            imported_hpb = True
        except ImportError:
            LOGGER.info('Please provide hpb.so file.')

    if imported_hpb:
        Hydrophobic_calculations = sorted(Hydrophobic_calculations, key=lambda x : x[-2])
        Hydrophobic_calculations_final = removeDuplicates(Hydrophobic_calculations)
        Hydrophobic_calculations_final2 = filterInteractions(Hydrophobic_calculations_final, atoms, **sel_kwargs)
        
        if zerosHPh == False:
            Hydrophobic_calculations_final3 = [ i for i in Hydrophobic_calculations_final2 if i[-1] != 0 ]
        else:
            Hydrophobic_calculations_final3 = Hydrophobic_calculations_final2
        
        for kk in Hydrophobic_calculations_final3:
            LOGGER.info("%10s%5s%14s14s  <---> %10s%5s%14s%8.1f%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6], kk[7]))
            
        LOGGER.info("Number of detected hydrophobic interactions: {0}.".format(len(Hydrophobic_calculations_final3)))
        
        return Hydrophobic_calculations_final3

    else:
        Hydrophobic_calculations = sorted(Hydrophobic_calculations, key=lambda x : x[-1])
        Hydrophobic_calculations_final = removeDuplicates(Hydrophobic_calculations)
        Hydrophobic_calculations_final2 = filterInteractions(Hydrophobic_calculations_final, atoms, **sel_kwargs)
        
        for kk in Hydrophobic_calculations_final2:
            LOGGER.info("%10s%5s%14s  <---> %10s%5s%14s%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6]))
            
        LOGGER.info("Number of detected hydrophobic interactions: {0}.".format(len(Hydrophobic_calculations_final2)))
        
        return Hydrophobic_calculations_final2
    

def calcDisulfideBonds(atoms, **kwargs):
    """Prediction of disulfide bonds.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg distDB: non-zero value, maximal distance between atoms of hydrophobic residues.
        default is 3.
        distA works too
    :type distDB: int, float
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
    
    distDB = kwargs.pop('distDB', 3)
    distA = kwargs.pop('distA', distDB)
    
    from prody.measure import calcDihedral
    
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
                        sele1_CB = atoms.select('resname CYS and name CB and resid '+str(sele1_new.getResnums()[0])+
                            ' and chain '+str(sele1_new.getChids()[0]))
                        sele2_CB = atoms.select('resname CYS and name CB and resid '+str(sele2_new.getResnums()[0])+
                            ' and chain '+str(sele2_new.getChids()[0]))
                        diheAng = calcDihedral(sele1_CB, sele1_new, sele2_new, sele2_CB)
                        DisulfideBonds_list.append([sele1_new.getResnames()[0]+str(sele1_new.getResnums()[0]),
                                                                minDistancePair[2]+'_'+str(minDistancePair[0]), sele1_new.getChids()[0],
                                                                sele2_new.getResnames()[0]+str(sele2_new.getResnums()[0]),
                                                                minDistancePair[3]+'_'+str(minDistancePair[1]), sele2_new.getChids()[0],
                                                                round(minDistancePair[-1],4), round(float(diheAng),4)])
    except:
        atoms_SG = atoms.select('protein and resname CYS')
        if atoms_SG is None:
            LOGGER.info('Lack of cysteines in the structure.')
            DisulfideBonds_list = []

    DisulfideBonds_list_final = removeDuplicates(DisulfideBonds_list)

    sel_kwargs = {k: v for k, v in kwargs.items() if k.startswith('selection')}
    DisulfideBonds_list_final2 = filterInteractions(DisulfideBonds_list_final, atoms, **sel_kwargs)

    for kk in DisulfideBonds_list_final2:
        LOGGER.info("%10s%5s%14s  <---> %10s%5s%14s%8.1f%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6], kk[7]))

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
    :type excluded_ions: list
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
    
    try:
        atoms_ions = atoms.select('ion and not name '+' '.join(excluded_ions)+' or (name '+' '.join(map(str,extraIons))+')')
        MetalResList = []
        MetalRes_calculations = cleanNumbers(findNeighbors(atoms_ions, distA, atoms.select('all and noh')))
        for i in MetalRes_calculations:
            if i[-1] != 0:
                MetalResList.append([atoms.getResnames()[i[0]]+str(atoms.getResnums()[i[0]]), i[2], atoms.getChids()[i[0]],
                                     atoms.getResnames()[i[1]]+str(atoms.getResnums()[i[1]]), i[3], atoms.getChids()[i[1]],
                                     i[-1]])
                
        for kk in MetalResList:
            LOGGER.info("%10s%5s%14s  <---> %10s%5s%14s%8.1f" % (kk[0], kk[2], kk[1], kk[3], kk[5], kk[4], kk[6]))

        LOGGER.info("Number of detected metal bonds: {0}.".format(len(MetalResList)))

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
        
        atoms_copy = atoms.copy()
        for j0, frame0 in enumerate(traj, start=start_frame):
            LOGGER.info('Frame: {0}'.format(j0))
            atoms_copy.setCoords(frame0.getCoords())
            protein = atoms_copy.select('protein')
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

    LOGGER.info('Calculating interations.') 
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
    :type stop_frame: int
    """

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
        by default 0.1.
    :type cutoff: int, float

    :arg code: representation of the residues, 3-letter or 1-letter
        by default 3-letter
    :type code: str

    :arg multiple_chains: display chain name for structure with many chains
        by default False
    :type multiple_chains: bool
    
    :arg edge_cmap: color of the residue connection
        by default plt.cm.Blues (blue color)
    :type edge_cmap: str

    :arg node_size: size of the nodes which describes residues
        by default 300
    :type node_size: int
    
    :arg node_distance: value which will scale residue-residue interactions
        by default 5
    :type node_distance: int

    :arg font_size: size of the font
        by default 14
    :type font_size: int

    :arg seed: random number which affect the distribution of residues
        by default 42
    :type seed: int
    """
    
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
    
    if len(statistics[0]) != 5:
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
    cutoff = kwargs.pop('cutoff', 0.1)
    multiple_chains = kwargs.pop('multiple_chains', False)
    
    X = [i for i in statistics if i[1] >= cutoff]   
    G = nx.Graph()
    
    for row in X:  
        if multiple_chains == False:
            if code == '1-letter':
                aa1 = aa_dic[row[0].split('-')[0][:3]] + row[0].split('-')[0][3:-1]
                aa2 = aa_dic[row[0].split('-')[1][:3]] + row[0].split('-')[1][3:-1]
            else:
                aa1 = row[0].split('-')[0][:-1]
                aa2 = row[0].split('-')[1][:-1]   
        else:
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
    
    
def calcStatisticsInteractions(data, **kwargs):
    """Return the statistics of interactions from PDB Ensemble or trajectory including:
    (1) the weight for each residue pair: corresponds to the number of counts divided by the 
    number of frames (values >1 are obtained when residue pair creates multiple contacts); 
    (2) average distance of interactions for each pair [in Ang] and (3) standard deviation [Ang.].
        
    :arg data: list with interactions from calcHydrogenBondsTrajectory() or other types
    :type data: list
    
    :arg weight_cutoff: value above which results will be displayed
        1 or more means that residue contact is present in all conformations/frames
        default value is 0.2 (in 20% of conformations contact appeared)
    :type weight_cutoff: int, float 
    
    :arg energy_list_type: name of the list with energies 
                            default is 'IB_solv'
    :type energy_list_type: 'IB_nosolv', 'IB_solv', 'CS'

    
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
    weight_cutoff = kwargs.pop('weight_cutoff', 0.2)
    energy_list_type = kwargs.pop('energy_list_type', 'IB_solv')
    
    import numpy as np
    elements = [t[0] for t in interactions_list]
    stats = {}

    for element in elements:
        if element not in stats:
            values = [t[1] for t in interactions_list if t[0] == element]
            stats[element] = {
                "stddev": np.round(np.std(values),6),
                "mean": np.round(np.mean(values),6),
                "weight": np.round(float(len(values))/len(data), 6),
                "energy": get_energy([element.split('-')[0][:3], element.split('-')[1][:3]], energy_list_type)
            }

    statistic = []
    for key, value in stats.items():
        if float(value['weight']) > weight_cutoff:
            LOGGER.info("Statistics for {0}:".format(key))
            LOGGER.info("  Average [Ang.]: {}".format(value['mean']))
            LOGGER.info("  Standard deviation [Ang.]: {0}".format(value['stddev']))
            LOGGER.info("  Weight: {0}".format(value['weight']))
            LOGGER.info("  Energy [kcal/mol]: {0}".format(value['energy']))
            statistic.append([key, value['weight'], value['mean'], value['stddev'], value['energy']])
        else: pass
    
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


def listLigandInteractions(PLIP_output, **kwargs):
    """Create a list of interactions from PLIP output created using calcLigandInteractions().
    Results can be displayed in VMD. 
    
    :arg PLIP_output: Results from PLIP for protein-ligand interactions.
    :type PLIP_output: PLIP object obtained from calcLigandInteractions()
    
    :arg output: parameter to print the interactions on the screen
                 while analyzing the structure (True | False)
                 by default is False
    :type output: bool     
    
    Note that five types of interactions are considered: hydrogen bonds, salt bridges, 
    pi-stacking, cation-pi, hydrophobic and water bridges."""
    
    Inter_list_all = []
    output = kwargs.pop('output', False)
    
    for i in PLIP_output.all_itypes:
        param_inter = [method for method in dir(i) if method.startswith('_') is False]
        
        if str(type(i)).split('.')[-1].strip("'>") == 'hbond':
            Inter_list = ['HBs',i.restype+str(i.resnr), i[0].type+'_'+str(i.d_orig_idx), i.reschain,
                          i.restype_l+str(i.resnr_l), i[2].type+'_'+str(i.a_orig_idx), i.reschain_l, 
                          i.distance_ad, i.angle, i[0].coords, i[2].coords]
     
        if str(type(i)).split('.')[-1].strip("'>") == 'saltbridge':
            Inter_list = ['SBs',i.restype+str(i.resnr), '_'.join(map(str,i.negative.atoms_orig_idx)), i.reschain,
                          i.restype_l+str(i.resnr_l), '_'.join(map(str,i.positive.atoms_orig_idx)), i.reschain_l, 
                          i.distance, None, i.negative.center, i.positive.center]
                 
        if str(type(i)).split('.')[-1].strip("'>") == 'pistack':
             Inter_list = ['PiStack',i.restype+str(i.resnr), '_'.join(map(str,i[0].atoms_orig_idx)), i.reschain,
                          i.restype_l+str(i.resnr_l), '_'.join(map(str,i[1].atoms_orig_idx)), i.reschain_l, 
                          i.distance, i.angle, i[0].center, i[1].center]           
        
        if str(type(i)).split('.')[-1].strip("'>") == 'pication':
             Inter_list = ['PiCat',i.restype+str(i.resnr), '_'.join(map(str,i[0].atoms_orig_idx)), i.reschain,
                          i.restype_l+str(i.resnr_l), '_'.join(map(str,i[1].atoms_orig_idx)), i.reschain_l, 
                          i.distance, None, i[0].center, i[1].center]                       
        
        if str(type(i)).split('.')[-1].strip("'>") == 'hydroph_interaction':
            Inter_list = ['HPh',i.restype+str(i.resnr), i[0].type+'_'+str(i[0].idx), i.reschain,
                          i.restype_l+str(i.resnr_l), i[2].type+'_'+str(i[2].idx), i.reschain_l, 
                          i.distance, None, i[0].coords, i[2].coords]           
             
        if str(type(i)).split('.')[-1].strip("'>") == 'waterbridge':
            water = i.water
            Inter_list = ['watBridge',i.restype+str(i.resnr), i[0].type+'_'+str(i[0].idx), i.reschain,
                          i.restype_l+str(i.resnr_l), i[3].type+'_'+str(i[3].idx), i.reschain_l, 
                          [i.distance_aw, i.distance_dw], [i.d_angle, i.w_angle], i[0].coords, i[3].coords, 
                          i.water.coords, i[7].residue.name+'_'+str(i[7].residue.idx)]
                      
        Inter_list_all.append(Inter_list)               
    
    if output == True:
        LOGGER.info("%3s%12s%10s%20s%8s  <---> %6s%10s%6s%10s%16s" % ('#','Type','Residue','Atoms','Chain','','Ligand','Atoms','Chain','Distance/Angle'))
        for nr_k,k in enumerate(Inter_list_all):
            if k[0] == 'watBridge':
                LOGGER.info("%3i%12s%10s%26s%4s  <---> %8s%12s%4s%12s%14s" % (nr_k+1,k[0],k[1],k[2],k[3],k[4],k[5],k[6], 
                                    ' '.join(str(np.round(x, 2)) for x in k[7]), ' '.join(str(np.round(x, 2)) for x in k[8])))
            else:
                LOGGER.info("%3i%12s%10s%26s%4s  <---> %8s%12s%4s%6.1f" % (nr_k+1,k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]))
    
    return Inter_list_all


def calcLigandInteractions(atoms, **kwargs):
    """Provide ligand interactions with other elements of the system including protein, 
    water and ions. Results are computed by PLIP [SS15]_ which should be installed.
    Note that PLIP will not recognize ligand unless it will be HETATM in the PDB file.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg sele: a selection string for residues of interest
            default is 'all not (water or protein or ion)'
    :type sele: str
    
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

    except ImportError:
        raise ImportError("Install Openbabel and PLIP.")
           
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb_file:
        writePDB(temp_pdb_file.name, atoms, csets=atoms.getACSIndex())
        temp_pdb_file_name = temp_pdb_file.name

    try:
        if atoms.hydrogen == None or atoms.hydrogen.numAtoms() < 30: # if there is no hydrogens in PDB structure
            LOGGER.info("Lack of hydrogens in the structure. Use addMissingAtoms to add them.")
    except: 
        pass

    Ligands = [] # Ligands can be more than one
    my_mol = PDBComplex()
    my_mol.load_pdb(temp_pdb_file_name) # Load the PDB file into PLIP class
    
    if 'sele' in kwargs:
        sele = kwargs['sele']
    else:
        sele='all not (water or protein or ion)'

    if 'ignore_ligs' in kwargs:
        ignore_ligs = kwargs['ignore_ligs']
    else:
        ignore_ligs=['MAN', 'SOD', 'CLA']
    
    sele = sele+' and not (resname '+' '.join(ignore_ligs)+')'
    ligand_select = atoms.select(sele)
    analyzedLigand = []
    
    try:
        for i in range(len(ligand_select.getResnums())):
            ResID = ligand_select.getResnames()[i]
            ChainID = ligand_select.getChids()[i]
            ResNames = ligand_select.getResnums()[i]
            my_bsid = str(ResID)+':'+str(ChainID)+':'+str(ResNames)
            if my_bsid not in analyzedLigand:
                LOGGER.info("LIGAND:  {0}".format(my_bsid))
                analyzedLigand.append(my_bsid)
                my_mol.analyze()
                my_interactions = my_mol.interaction_sets[my_bsid] # Contains all interaction data      
                Ligands.append(my_interactions)
                listLigandInteractions(my_interactions)
    
        return Ligands, analyzedLigand
    
    except:
        LOGGER.info("Ligand not found.")


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
                default **"red"**
    :type color: str
    
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
        
        :arg interactions: List of interaction lists for protein interactions.
        :type interactions: list
        
        :arg color: Name of the color which will be used for the visualization of 
                    interactions in VMD
        :type color: str
        
        :arg tcl_file: name of the TCL file which will be saved for visualization                
        :type tcl_file: str
        """
        
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

    elif interactions == []:
        LOGGER.info("Lack of results")

    elif len(interactions[0]) == 0:
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
    
    :arg interactions: List of interactions lists for protein-ligand interactions.
    :type interactions: list
    
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
    
    pdb_name = atoms.getTitle()+'_sele.pdb'
    writePDB(pdb_name, atoms)
    
    tcl_file = open(filename, 'w') 
    
    if len(interactions[0]) >= 10: 
        dic_color = {'HBs':'blue','PiStack':'green','SBs':'yellow','PiCat':'orange',
                     'HPh':'silver','watBridge':'cyan'}
        
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


def SMINA_extract_data(result):
    """Supporting function for analysis of the SMINA results."""
    import re
    
    data = {}

    # Extracting Affinity from SMINA output
    affinity_match = re.search(r'Affinity: ([0-9.-]+) \(kcal/mol\)', result)
    if affinity_match:
        data['Affinity'] = float(affinity_match.group(1))

    # Extracting Intramolecular energy from SMINA output
    intramolecular_energy_match = re.search(r'Intramolecular energy: ([0-9.-]+)', result)
    if intramolecular_energy_match:
        data['Intramolecular energy'] = float(intramolecular_energy_match.group(1))

    # Extracting Weights and Terms from SMINA output
    weights_terms_match = re.search(r'Weights\s+Terms\s*([\s\S]*?)## Name', result, re.DOTALL)
    if weights_terms_match:
        weights_terms_text = weights_terms_match.group(1)
        term_lines = weights_terms_text.strip().split('\n')
        term_dict = {}
        for line in term_lines:
            parts = line.split()
            if len(parts) >= 2: 
                weight = float(parts[0])
                term = ' '.join(parts[1:])
                term_dict[term] = weight
        data.update(term_dict)
        
    term_values_match = re.search(r'Term values, before weighting:\n##\s+(.*?)\n', result, re.DOTALL)
    if term_values_match:
        term_values_text = term_values_match.group(1)
        term_values_array = np.array([float(value) for value in term_values_text.split()])
        data['Term values, before weighting'] = term_values_array.tolist()           
        
    return data


def calcSminaBindingAffinity(atoms, trajectory=None, **kwargs):
    """Computing binding affinity of ligand toward protein structure
    using SMINA package [DRK13]_.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`
    
    :arg protein_selection: selection string for the protein and other compoment
                            of the system that should be included,
                            e.g. "protein and chain A",
                            by default "protein" 
    :type protein_selection: str
    
    :arg ligand_selection: selection string for ligand,
                           e.g. "resname ADP",
                           by default "all not protein"
    :type ligand_selection: str
    
    :arg ligand_selection: scoring function (vina or vinardo)
                           by default is "vina"
    
    :type ligand_selection: str
    
    :arg atom_terms: write per-atom interaction term values.
                     It will provide more information as dictionary for each frame/model,
                     and details for atom terms will be saved in terms_*frame_number*.txt,    
                     by default is False

    :type atom_terms: bool
     

    SMINA installation is required to compute ligand binding affinity:
    >> conda install -c conda-forge smina       (for Anaconda)
    
    For more information on SMINA see https://sourceforge.net/projects/smina/.
    If you benefited from SMINA, please consider citing [DRK13]_.

    .. [DRK13] Koes D. R., Baumgartner M. P., Camacho C. J., Lessons Learned in 
    Empirical Scoring with smina from the CSAR 2011 Benchmarking Exercise,
    *J. Chem. Inf. Model.* **2013** 53: 1893–1904. """

    import tempfile
    import subprocess
    import re

    if isinstance(atoms, LigandInteractionsTrajectory):
        atoms = atoms._atoms
        trajectory = atoms._traj
    
    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)
    protein_selection = kwargs.pop('protein_selection', "protein")
    ligand_selection = kwargs.pop('ligand_selection', "all not protein")
    scoring_function = kwargs.pop('scoring_function', 'vina')
    atom_terms = kwargs.pop('atom_terms', False)
    bindingAffinity = []

    if trajectory is not None:
        # Trajectory
        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)
    
        nfi = trajectory._nfi
        trajectory.reset()
        numFrames = trajectory._n_csets

        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame+1]

        atoms_copy = atoms.copy()
        data_final = []
        
        for j0, frame0 in enumerate(traj, start=start_frame):
            atoms_copy.setCoords(frame0.getCoords())        
            protein = atoms_copy.select(protein_selection)
            ligand = atoms_copy.select(ligand_selection)        

            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb_file, \
                 tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb_file_lig:
                
                # Files are starage in the memory
                writePDB(temp_pdb_file.name, protein)
                writePDB(temp_pdb_file_lig.name, ligand)

                if atom_terms == False:
                    command = "smina -r {} -l {} --score_only --scoring {}".format(temp_pdb_file.name, 
                                                        temp_pdb_file_lig.name, scoring_function)
                else:
                    command = "smina -r {} -l {} --score_only --scoring {} --atom_terms terms_{}.txt".format(temp_pdb_file.name, 
                                                        temp_pdb_file_lig.name, scoring_function, j0)
                
                result = subprocess.check_output(command, shell=True, text=True)
                data = SMINA_extract_data(result)
                LOGGER.info('Frame {0}: {1} kcal/mol'.format(j0, data['Affinity']))
                bindingAffinity.append(data['Affinity'])
                data_final.append(data)
        
        trajectory._nfi = nfi
                
    else:
        if atoms.numCoordsets() == 1:
            # Single PDB
            protein = atoms.select(protein_selection)
            ligand = atoms.select(ligand_selection)        

            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb_file, \
                 tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb_file_lig:
                
                writePDB(temp_pdb_file.name, protein)
                writePDB(temp_pdb_file_lig.name, ligand)

                if atom_terms == False:
                    command = "smina -r {} -l {} --score_only --scoring {}".format(temp_pdb_file.name, 
                                                        temp_pdb_file_lig.name, scoring_function)
                else:
                    command = "smina -r {} -l {} --score_only --scoring {} --atom_terms terms.txt".format(temp_pdb_file.name, 
                                                        temp_pdb_file_lig.name, scoring_function)
               
                result = subprocess.check_output(command, shell=True, text=True)
                data = SMINA_extract_data(result)
                data_final = data
                bindingAffinity.append(data['Affinity'])

        if atoms.numCoordsets() > 1:
            # Multi-model PDB
            data_final = []
            for i in range(len(atoms.getCoordsets()[start_frame:stop_frame+1])):
                atoms.setACSIndex(i+start_frame) 
                protein = atoms.select(protein_selection)
                ligand = atoms.select(ligand_selection)        

                with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb_file, \
                     tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb_file_lig:
                    
                    writePDB(temp_pdb_file.name, protein, csets=atoms.getACSIndex())
                    writePDB(temp_pdb_file_lig.name, ligand, csets=atoms.getACSIndex())

                    if atom_terms == False:
                        command = "smina -r {} -l {} --score_only --scoring {}".format(temp_pdb_file.name, 
                                                            temp_pdb_file_lig.name, scoring_function)
                    else:
                        command = "smina -r {} -l {} --score_only --scoring {} --atom_terms terms_{}.txt".format(temp_pdb_file.name, 
                                                        temp_pdb_file_lig.name, scoring_function, i+start_frame)
               
                    result = subprocess.check_output(command, shell=True, text=True)
                    data = SMINA_extract_data(result)
                    LOGGER.info('Model {0}: {1} kcal/mol'.format(i+start_frame, data['Affinity']))
                    bindingAffinity.append(data['Affinity'])
                    data_final.append(data)

        else:
            LOGGER.info('Include trajectory or use multi-model PDB file.') 

    if atom_terms == False:
        return bindingAffinity
    else:
        return data_final


def calcSminaPerAtomInteractions(atoms, list_terms):
    """Computing the summary of per-atom interaction term values
    using SMINA package [DRK13]_. It will return dictionary with per-atom interaction term values.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`
    
    :arg list_terms: List of *terms.txt* files obtained from meth:`.calcSminaBindingAffinity`
                     using *atom_terms = True*
    :type list_terms: list 
    
    Important: First text file in the list should be reference structure which correspond to the 
    provided coordinates as atoms.
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
    
    if not isinstance(list_terms, list):
        raise TypeError('list_terms must be a list of text files with per-atom interaction term values.')
    
    LOGGER.info('Reference file: {}'.format(list_terms[0]))
    ref_file = open(list_terms[0], 'r').readlines()
    dict_terms = {}
    for nr_j, j in enumerate(ref_file[1:-2]):
            inter_type = j.split()[0]
            xyz = j.split()[1].strip('<').strip('>').split(',')
            try:
                xyz_atom = (atoms.select('x `'+str(xyz[0])+'` y `'+str(xyz[1])+'` z `'+str(xyz[2])+'`')).getNames()[0]
            except:
                xyz_atom = ' '

            sum_of_energy = np.sum([float(i) for i in j.split()[2:]])
            atom_name_with_type = inter_type+ ' '+xyz_atom
            dict_terms[atom_name_with_type] = [sum_of_energy]            

            # Checking by specific line each file
            for i in list_terms[1:]:
                an_line = open(i, 'r').readlines()[nr_j+1]
                sum_of_energy_line = np.sum([float(i) for i in an_line.split()[2:]])
                dict_terms[atom_name_with_type].append(sum_of_energy_line)
    
    return dict_terms


def calcSminaTermValues(data):
    """Computing weights multiplied by term values, before weighting for each Term.
    As a results will are obtaining a dictionary.
    
    :arg data: List of results provided by Smina using meth:`.calcSminaBindingAffinity`
                     with *atom_terms = True*
    :type data: list 
    """
    
    if not isinstance(data, list):
        raise TypeError('data must be a list')

    result_dict = {key: [] for key in list(data[0].keys())[2:-1]}

    for i in data:
        weights = i['Term values, before weighting']
        for idx, key in enumerate(result_dict.keys()):
            result_dict[key].append(i[key] * weights[idx] if key in i else None)

    return result_dict


def showSminaTermValues(data):
    """Display a histogram of weights multiplied by term values, before weighting for each Term.
    As a results will are obtaining a dictionary.
    
    :arg data: List of results provided by Smina using meth:`.calcSminaBindingAffinity`
                     with *atom_terms = True*
    :type data: list 
    """

    import matplotlib.pyplot as plt
    import numpy as np
    
    if not isinstance(data, list):
        raise TypeError('data must be a list')

    term_values = calcSminaTermValues(data)
    non_zero_values = {key: [v for v in value if v != 0] for key, value in term_values.items()}

    fig, ax = plt.subplots()
    colors = ['blue', 'orange', 'red', 'green', 'purple', 'silver', 'cyan', 'magenta', 'yellow']
    alpha = 0.5

    for i, (key, values) in enumerate(non_zero_values.items()):
        if values:
            show = ax.hist(values, bins=10, alpha=alpha, label=key, color=colors[i % len(colors)])

    ax.legend()
    ax.set_xlabel('Energy [kcal/mol]')
    ax.set_ylabel('# counts')

    if SETTINGS['auto_show']:
        showFigure()
    return show



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

        LOGGER.info('Calculating interations.') 
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
        
        LOGGER.info('Calculating interactions')
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

        if all(x in [0, 1] for x in scoring):
            pass
        else:
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
        
        from collections import Counter
        lista_ext = []
        atoms = atoms.select("protein and noh")
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
        
        :arg contacts_min: Minimal number of contacts which residue may form with other residues, 
                           by default 3.
        :type contacts_min: int  """

        atoms = self._atoms   
        interactions = self._interactions
        
        InteractionsMap = np.empty([atoms.select('name CA').numAtoms(),atoms.select('name CA').numAtoms()], dtype=object)
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
            
        ListOfInteractions = [ list(filter(None, InteractionsMap[:,j])) for j in range(len(interactions[0])) ]
        ListOfInteractions = list(filter(lambda x : x != [], ListOfInteractions))
        ListOfInteractions = [k for k in ListOfInteractions if len(k) >= contacts_min ]
        ListOfInteractions_list = [ (i[0].split('-')[-1], [ j.split('-')[0] for j in i]) for i in ListOfInteractions ]
        LOGGER.info('The most frequent interactions between:')
        for res in ListOfInteractions_list:
            LOGGER.info('{0}  <--->  {1}'.format(res[0], '  '.join(res[1])))

        LOGGER.info('Legend: hb-hydrogen bond, sb-salt bridge, rb-repulsive ionic bond, ps-Pi stacking interaction,'
                             'pc-Cation-Pi interaction, hp-hydrophobic interaction, dibs-disulfide bonds')
        
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
        """Bar plot with the number of potential inetractions per residue.
        Particular type of interactions can be excluded by using keywords HBs=0, RIB=0, etc.
        
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
        :type DiBs: int, float """

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
        
        replace_matrix = kwargs.get('replace_matrix', False)
        matrix_all = self._interactions_matrix

        HBs = kwargs.get('HBs', 1)
        SBs = kwargs.get('SBs', 1)
        RIB = kwargs.get('RIB', 1)
        PiStack = kwargs.get('PiStack', 1)
        PiCat = kwargs.get('PiCat', 1)
        HPh = kwargs.get('HPh', 1)
        DiBs = kwargs.get('DiBs', 1)
    
        matrix_hbs = self.buildInteractionMatrix(HBs=HBs, SBs=0, RIB=0,PiStack=0,PiCat=0,HPh=0,DiBs=0)
        matrix_sbs = self.buildInteractionMatrix(HBs=0, SBs=SBs, RIB=0,PiStack=0,PiCat=0,HPh=0,DiBs=0)
        matrix_rib = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=RIB,PiStack=0,PiCat=0,HPh=0,DiBs=0)
        matrix_pistack = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=PiStack,PiCat=0,HPh=0,DiBs=0)
        matrix_picat = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=0,PiCat=PiCat,HPh=0,DiBs=0)
        matrix_hph = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=0,PiCat=0,HPh=HPh,DiBs=0)
        matrix_dibs = self.buildInteractionMatrix(HBs=0, SBs=0, RIB=0,PiStack=0,PiCat=0,HPh=0,DiBs=DiBs)

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

        if HBs != 0:
            ax.bar(ResList, matrix_hbs_sum, width, color = 'blue', bottom = 0, label='HBs')
        sum_matrix += matrix_hbs_sum

        if SBs != 0:
            ax.bar(ResList, matrix_sbs_sum, width, color = 'yellow', bottom = sum_matrix, label='SBs')
        sum_matrix += matrix_sbs_sum

        if HPh != 0:
            ax.bar(ResList, matrix_hph_sum, width, color = 'silver', bottom = sum_matrix, label='HPh')
        sum_matrix += matrix_hph_sum
        
        if RIB != 0:
            ax.bar(ResList, matrix_rib_sum, width, color = 'red', bottom = sum_matrix, label='RIB')
        sum_matrix += matrix_rib_sum

        if PiStack != 0:
            ax.bar(ResList, matrix_pistack_sum, width, color = 'green', bottom = sum_matrix, label='PiStack')
        sum_matrix += matrix_pistack_sum

        if PiCat != 0:
            ax.bar(ResList, matrix_picat_sum, width, color = 'orange', bottom = sum_matrix, label='PiCat')
        sum_matrix += matrix_picat_sum
        
        if DiBs != 0:
            ax.bar(ResList, matrix_dibs_sum, width, color = 'black', bottom = sum_matrix, label='DiBs')
        sum_matrix += matrix_dibs_sum

        if replace_matrix:
            self._interactions_matrix = np.sum([matrix_hbs, matrix_sbs, matrix_rib, matrix_pistack,
                                                matrix_picat, matrix_hph, matrix_dibs], axis=0)
        else:
            self._interactions_matrix = matrix_all

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

            atoms_copy = atoms.copy()
            protein = atoms_copy.protein

            for j0, frame0 in enumerate(traj, start=start_frame):
                LOGGER.info('Frame: {0}'.format(j0))
                atoms_copy.setCoords(frame0.getCoords())
                
                hydrogen_bonds = calcHydrogenBonds(protein, **kwargs)
                salt_bridges = calcSaltBridges(protein, **kwargs)
                RepulsiveIonicBonding = calcRepulsiveIonicBonding(protein, **kwargs)
                Pi_stacking = calcPiStacking(protein, **kwargs)
                Pi_cation = calcPiCation(protein, **kwargs)
                hydrophobic = calcHydrophobic(protein, **kwargs)
                Disulfide_Bonds = calcDisulfideBonds(protein, **kwargs)

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

        data = [HBs, SBs, HPh, PiStack, PiCat, RIB, DiBs]
        colors = ['deepskyblue', 'yellow', 'silver', 'lightgreen', 'orange', 'red', 'black']
        labels = ['HBs', 'SBs', 'HPh', 'PiStack', 'PiCat', 'RIB', 'DiBs']
        
        # Removing empty interaction plots
        data_filtered = []
        colors_filtered = []
        labels_filtered = []

        for i, lst in enumerate(data):
            if any(lst):
                data_filtered.append(lst)
                colors_filtered.append(colors[i])
                labels_filtered.append(labels[i])

        data = data_filtered
        colors = colors_filtered
        labels = labels_filtered

        fig, axes = plt.subplots(len(data), num=None, figsize=(12, 8),
                                facecolor='w', sharex='all', **kwargs)
        hspace = 0.1
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.35)

        last_nonempty = None
        for i, ax in enumerate(axes):
            if any(data[i]):
                last_nonempty = i
                rects = ax.bar(np.arange(len(data[i])), data[i], color=colors[i])
                ax.plot(data[i], 'k:')
                
                data_min = min(data[i])
                data_max = max(data[i])
                if data_max > 15:
                    ax.set_ylim(data_min - 1, data_max + 1)
                else: pass
            
                if i != len(data) - 1:
                    ax.set_xlabel('')
            else:
                ax.axis('off')

        axes[last_nonempty].set_xlabel('Conformation')

        handles = []
        for i in range(len(data)):
            if any(data[i]):
                handles.append(matplotlib.patches.Rectangle((0, 0), 0.6, 0.6, color=colors[i]))

        if handles:
            legend = fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=len(handles), prop={'size': 14}, handlelength=0.8)
            legend.set_title("Interactions")

        if filename is not None:
            plt.savefig(filename + '.png', dpi=300)
        plt.show()

        return HBs, SBs, HPh, PiStack, PiCat, HPh, DiBs


 
class LigandInteractionsTrajectory(object):
    """Class for protein-ligand interaction analysis of DCD trajectory or multi-model PDB (Ensemble PDB).
    This class is using PLIP to provide the interactions. Install PLIP before using it.

    ## Instalation of PLIP using conda:
    >>> conda install -c conda-forge plip
    ## https://anaconda.org/conda-forge/plip
    # https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md
    
    ## Instalation using PIP:
    >>> pip install plip

    .. [SS15] Salentin S., Schreiber S., Haupt V. J., Adasme M. F., Schroeder M.  
    PLIP: fully automated protein–ligand interaction profiler 
    *Nucl. Acids Res.* **2015** 43:W443-W447.     """

    def __init__(self, name='Unknown'):
        
        self._atoms = None
        self._traj = None
        self._interactions_traj = None
        self._freq_interactors = None

    
    def calcLigandInteractionsTrajectory(self, atoms, trajectory=None, **kwargs):
        """Compute protein-ligand interactions for DCD trajectory or multi-model PDB 
            using PLIP library.
        
        :arg atoms: an Atomic object from which residues are selected
        :type atoms: :class:`.Atomic`
        
        :arg trajectory: trajectory file
        :type trajectory: class:`.Trajectory`

        :arg filename: Name of pkl filename in which interactions will be storage
        :type filename: pkl 
        
        :arg start_frame: index of first frame to read
        :type start_frame: int

        :arg stop_frame: index of last frame to read
        :type stop_frame: int 
        
        :arg output: parameter to print the interactions on the screen
                    while analyzing the structure,
                    use 'info'.
        :type output: bool
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
        
        interactions_all = []
        interactions_all_nb = []

        start_frame = kwargs.pop('start_frame', 0)
        stop_frame = kwargs.pop('stop_frame', -1)
        output = kwargs.pop('output', False)
        filename = kwargs.pop('filename', None)

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

            atoms_copy = atoms.copy()
            
            for j0, frame0 in enumerate(traj, start=start_frame):
                LOGGER.info(' ')
                LOGGER.info('Frame: {0}'.format(j0))
                atoms_copy.setCoords(frame0.getCoords())

                ligand_interactions, ligand = calcLigandInteractions(atoms_copy)
                
                ligs_per_frame_interactions = []
                for ligs in ligand_interactions:
                    LP_interactions = listLigandInteractions(ligs, output=output) 
                    ligs_per_frame_interactions.extend(LP_interactions)
                
                interactions_all.append(ligs_per_frame_interactions)
                interactions_all_nb.append(len(ligs_per_frame_interactions))

        else:
            if atoms.numCoordsets() > 1:
                for i in range(len(atoms.getCoordsets()[start_frame:stop_frame+1])):
                    LOGGER.info(' ')
                    LOGGER.info('Model: {0}'.format(i+start_frame))
                    atoms.setACSIndex(i+start_frame) 

                    ligand_interactions, ligand = calcLigandInteractions(atoms)
                    
                    ligs_per_frame_interactions = []
                    for ligs in ligand_interactions:
                        LP_interactions = listLigandInteractions(ligs, output=output) 
                        ligs_per_frame_interactions.extend(LP_interactions)
                    
                    interactions_all.append(ligs_per_frame_interactions)
                    interactions_all_nb.append(len(ligs_per_frame_interactions))

            else:
                LOGGER.info('Include trajectory or use multi-model PDB file.') 
        
        self._atoms = atoms
        self._traj = trajectory
        self._interactions_traj = interactions_all
        self._interactions_nb_traj = interactions_all_nb
        
        if filename is not None:
            import pickle
            with open(str(filename)+'.pkl', 'wb') as f:
                pickle.dump(self._interactions_traj, f)  
            LOGGER.info('File with interactions saved.')
            
        return interactions_all

    
    def getLigandInteractions(self, **kwargs):
        """Return the list of protein-ligand interactions.
                
        :arg filters: selection string of ligand with chain ID or interaction type
                     e.g. 'SBs' (HBs, SBs, HPh, PiStack, PiCat, HPh, watBridge)
        :type filters: str    
        
        :arg include_frames: used with filters, it will leave selected keyword in orginal 
                    lists, if False it will collect selected interactions in one list,
                    Use True to assign new selection using setLigandInteractions.
                    by default True
        :type include_frames: bool            
        """
        
        filters = kwargs.pop('filters', None)
        include_frames = kwargs.pop('include_frames', True)
        filtered_lists = []
                                 
        if filters != None:
            if include_frames == False:
                filtered_lists = [element for group in self._interactions_traj for element in group 
                                if filters in element]
            if include_frames == True:
                filtered_lists = []
                for i in self._interactions_traj:
                    filtered_lists.append([item for item in i if filters in item])
        
            return filtered_lists
            
        else:
            return self._interactions_traj
            

    def getAtoms(self):
        """Returns associated atoms."""

        return self._atoms


    def setLigandInteractions(self, atoms, interaction):
        """Replace protein-ligand interactions
        for example byb using getLigandInteractions() with filters to select particular ligand. 
        
        :arg atoms: an Atomic object from which residues are selected
        :type atoms: :class:`.Atomic`
        
        :arg interactions: list of interactions
        :type interactions: list
        """

        self._interactions_traj = interaction
        self._atoms = atoms
        LOGGER.info('Protein-ligand interactions are replaced.')
    
    
    def getLigandInteractionsNumber(self, **kwargs):
        """Return the number of interactions per each frame. Number of interactions can
        be a total number of interactions or it can be divided into interaction types.
        
        :arg types: Interaction types can be included (True) or not (False).
                    by default is True. 
        :type types: bool
        """

        types = kwargs.pop('types', True)
        
        if types == True:
            interactions = self._interactions_traj 
            unique_keywords = set()

            for sublist in interactions:
                for sublist_item in sublist:
                    keyword = sublist_item[0]
                    unique_keywords.add(keyword)
            unique_keywords_list = list(unique_keywords)
        
            keyword_counts = {keyword: [0] * len(interactions) for keyword in unique_keywords_list}

            for i, sublist in enumerate(interactions):
                for sublist_item in sublist:
                    keyword = sublist_item[0]  
                    keyword_counts[keyword][i] += 1

            return keyword_counts
        
        else:         
            return self._interactions_nb_traj 

    
    def parseLigandInteractions(self, filename):
        """Import interactions from analysis of trajectory which was saved via
        calcLigandInteractionsTrajectory().
        
        :arg filename: Name of pkl file in which interactions will be storage
        :type filename: pkl """
        
        import pickle
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        
        self._interactions_traj = data
        self._interactions_nb_traj = [[len(sublist) if sublist else 0 for sublist in sublist] for sublist in data]
        
        return data


    def getInteractionTypes(self):
        """Show which interaction types were detected for ligands."""
        
        interactions = self._interactions_traj 
        unique_keywords = set()

        for sublist in interactions:
            for sublist_item in sublist:
                keyword = sublist_item[0]
                unique_keywords.add(keyword)

        unique_keywords_list = list(unique_keywords)
        LOGGER.info("Interaction types: {0}".format(unique_keywords_list))
        
        return unique_keywords_list


    def getLigandsNames(self):
        """Show which ligands are in a system."""
        
        interactions = self._interactions_traj 
        ligands = set()

        for sublist in interactions:
            for sublist_item in sublist:
                keyword = sublist_item[4]
                ligands.add(keyword)

        ligands_list = list(ligands)
        
        return ligands_list


    def calcFrequentInteractors(self, **kwargs):
        """Returns a dictonary with residues involved in the interaction with ligand
        and their number of counts. 

        :arg selection: selection string of ligand with chain ID
                        e.g. "MESA" where MES is ligand resname and A is chain ID.
                        Selection pointed as None will return all interactions together
                        without ligands separation.
        :type selection: str
        """

        atoms = self._atoms   
        interactions = self._interactions_traj
        selection = kwargs.pop('selection', None)
        
        from collections import Counter

        if selection == None:  # Compute all interactions without distinguishing ligands
            all_residues = [ j[1]+j[3] for i in interactions for j in i ]
            dictOfInteractions = Counter(all_residues)
        
        else:
            interactions2 = [element for group in interactions for element in group]
            ligs = {}
            dictOfInteractions = []
            for i in interactions2:
                ligs_names = i[4]+i[6]
                if ligs_names not in ligs:
                    ligs[ligs_names] = []
                    
                res_name = i[1]+i[3]
                ligs[ligs_names].append(res_name)

            for i in ligs.keys():

                if selection == None:
                    LOGGER.info('LIGAND: {0}'.format(i))
                    aa_counter = Counter(ligs[i])
                    dictOfInteractions.append(aa_counter)
                    
                else:
                    if selection not in ligs.keys():
                        LOGGER.info('Wrong selection. Please provide ligand name with chain ID.')
                    if i == selection:
                        LOGGER.info('LIGAND: {0}'.format(selection))
                        aa_counter = Counter(ligs[selection])
                        dictOfInteractions.append(aa_counter)                    
                
        self._freq_interactors = dictOfInteractions
        
        return dictOfInteractions


    def saveInteractionsPDB(self, **kwargs):
        """Save the number of interactions with ligand to PDB file in occupancy column
        It will recognize the chains. If the system will contain one chain and many segments
        the PDB file will not be created in a correct way.
        
        :arg filename: name of the PDB file which will be saved for visualization,
                     it will contain the results in occupancy column.
        :type filename: str  
        
        :arg ligand_sele: ligand selection,
                          by default is 'all not (protein or water or ion)'.
        :type ligand_sele: str          
        """
        
        if self._freq_interactors is None:
            raise ValueError('Please calculate frequent interactors using getFrequentInteractors.')

        atoms = self._atoms     
        dictOfInteractions = self._freq_interactors
        ligand_sele = kwargs.pop('ligand_sele', 'all not (protein or water or ion)')

        chids_list = np.unique(atoms.getChids())
        freq_contacts_list = []
        
        for nr_chid, chid in enumerate(chids_list):
            atoms_chid = atoms.ca.select('protein and chain '+chid)
            freq_contacts_list_chids = np.zeros(atoms_chid.ca.numAtoms(), dtype=int)
            firstElement = atoms_chid.ca.getResnums()[0]
            dictOfInteractions_chids = dictOfInteractions[nr_chid]
            
            for k, v in dictOfInteractions_chids.items():
                res_index = int(k[3:-1])
                freq_contacts_list_chids[res_index - firstElement] = v

            freq_contacts_list.extend(freq_contacts_list_chids)

        freq_contacts_list = np.array(freq_contacts_list)

        from collections import Counter
        lista_ext = []
        ligands = atoms.select(ligand_sele)
        atoms = atoms.select("protein and noh")
        ligand_occupancy = np.zeros(len(ligands.getResnums()))
        
        aa_counter = Counter(atoms.getResindices())
        calphas = atoms.select('name CA')
        for i in range(calphas.numAtoms()):
            # in PDB values are normalized to 100 (max value)
            lista_ext.extend(list(aa_counter.values())[i]*[round((freq_contacts_list[i]/np.max(freq_contacts_list)*100), 8)])
        
        lista_ext.extend(ligand_occupancy)
        
        kw = {'occupancy': lista_ext}
        if 'filename' in kwargs:
            writePDB(kwargs['filename'], atoms+ligands, **kw)
            LOGGER.info('PDB file saved.')
        else:
            writePDB('filename', atoms+ligands, **kw)
            LOGGER.info('PDB file saved.')

        return freq_contacts_list

