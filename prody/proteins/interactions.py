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
__credits__ = ['James Krieger', 'Karolina Mikulska-Ruminska', 'Anupam Banerjee']
__email__ = ['karolamik@fizyka.umk.pl', 'jamesmkrieger@gmail.com', 'anupam.banerjee@stonybrook.edu']


import numpy as np
from numpy import *
from prody import LOGGER, SETTINGS, PY3K
from prody.atomic import AtomGroup, Atom, Atomic, Selection, Select
from prody.atomic import flags, sliceAtomicData
from prody.utilities import importLA, checkCoords, showFigure, getCoords
from prody.measure import calcDistance, calcAngle, calcCenter
from prody.measure.contacts import findNeighbors
from prody.proteins import writePDB, parsePDB, showProtein
from collections import Counter

from prody.trajectory import TrajBase, Trajectory, Frame
from prody.ensemble import Ensemble

import multiprocessing as mp
from .fixer import *
from .compare import *
from prody.measure import calcTransformation, calcDistance, calcRMSD, superpose


__all__ = ['calcHydrogenBonds', 'calcChHydrogenBonds', 'calcSaltBridges',
           'calcRepulsiveIonicBonding', 'calcPiStacking', 'calcPiCation',
           'calcHydrophobic', 'calcDisulfideBonds', 'calcMetalInteractions',
           'calcHydrogenBondsTrajectory', 'calcSaltBridgesTrajectory',
           'calcRepulsiveIonicBondingTrajectory', 'calcPiStackingTrajectory', 
           'calcPiCationTrajectory', 'calcHydrophobicTrajectory', 'calcDisulfideBondsTrajectory',
           'calcProteinInteractions', 'calcStatisticsInteractions', 'calcDistribution',
           'calcSASA', 'calcBSA', 'calcVolume','compareInteractions', 'showInteractionsGraph',
           'showInteractionsHist', 'calcLigandInteractions', 'listLigandInteractions', 
           'showProteinInteractions', 'showLigandInteraction',
           'showProteinInteractions_VMD', 'showLigandInteraction_VMD', 
           'calcHydrophobicOverlapingAreas',
           'Interactions', 'InteractionsTrajectory', 'LigandInteractionsTrajectory',
           'calcSminaBindingAffinity', 'calcSminaPerAtomInteractions', 'calcSminaTermValues',
           'showSminaPerAtomInteractions', 'calcSminaResidueContributions',
           'listSminaResidueContributions', 'showSminaResidueContributions',
           'showSminaLigandStructure',
           'showSminaTermValues', 'showPairEnergy', 'checkNonstandardResidues',
           'saveInteractionsAsDummyAtoms', 'createFoldseekAlignment', 'runFoldseek', 'runDali', 
           'runBLAST', 'extractMultiModelPDB', 'calcSignatureInteractions',
           'calcLigandResidueDistances', 'showLigandResidueDistances']


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
    """Return interactions based on *selection* and *selection2*."""
    
    if 'selection1' in kwargs:
        kwargs['selection'] = kwargs['selection1']

    if 'selection' in kwargs:
        selection = atoms.select(kwargs['selection'])
        if selection is None:
            LOGGER.warn('selection did not work, so no filtering is performed')
            return list_of_interactions

        ch1 = selection.getChids()
        x1 = selection.getResnames()
        y1 = selection.getResnums()
        listOfselection = np.unique(list(map(lambda x1, y1, ch1: (ch1, x1 + str(y1)),
                                             x1, y1, ch1)),
                                    axis=0)
        listOfselection = [list(i) for i in listOfselection] # needed for in check to work

        if 'selection2' in kwargs:
            selection2 = atoms.select(kwargs['selection2'])
            if selection2 is None:
                LOGGER.warn('selection2 did not work, so no filtering is performed')
                return list_of_interactions
            
            ch2 = selection2.getChids()
            x2 = selection2.getResnames()
            y2 = selection2.getResnums()
            listOfselection2 = np.unique(list(map(lambda x2, y2, ch2: (ch2, x2 + str(y2)),
                                                  x2, y2, ch2)),
                                         axis=0)
            listOfselection2 = [list(i) for i in listOfselection2] # needed for in check to work

            final = [i for i in list_of_interactions if (([i[2], i[0]] in listOfselection)
                                                         and ([i[5], i[3]] in listOfselection2)
                                                         or ([i[2], i[0]] in listOfselection2)
                                                         and ([i[5], i[3]] in listOfselection))]
        else:
            final = [i for i in list_of_interactions
                     if (([i[2], i[0]] in listOfselection)
                         or ([i[5], i[3]] in listOfselection))]

    elif 'selection2' in kwargs:
        LOGGER.warn('selection2 by itself is ignored')
        final = list_of_interactions
    else:
        final = list_of_interactions
    return final


def get_energy(pair, source):
    """Return energies based on the pairs of interacting residues (without distance criteria)
    Taking information from tabulated_energies.txt file"""

    import numpy as np
    import importlib.resources as pkg_resources    
    
    aa_correction = {
        # Histidine (His)
        'HSD': 'HIS',   # NAMD, protonated at ND1 (HID in AMBER)
        'HSE': 'HIS',   # NAMD, protonated at NE2 (HIE in AMBER)
        'HSP': 'HIS',   # NAMD, doubly protonated (HIP in AMBER)
        'HID': 'HIS',   # AMBER name, protonated at ND1
        'HIE': 'HIS',   # AMBER name, protonated at NE2
        'HIP': 'HIS',   # AMBER name, doubly protonated
        'HISD': 'HIS',  # GROMACS: protonated at ND1
        'HISE': 'HIS',  # GROMACS: protonated at NE2
        'HISP': 'HIS',  # GROMACS: doubly protonated

        # Cysteine (Cys)
        'CYX': 'CYS',   # Cystine (disulfide bridge)
        'CYM': 'CYS',   # Deprotonated cysteine, anion

        # Aspartic acid (Asp)
        'ASH': 'ASP',   # Protonated Asp
        'ASPP': 'ASP',

        # Glutamic acid (Glu)
        'GLH': 'GLU',   # Protonated Glu
        'GLUP': 'GLU',  # Protonated Glu

        # Lysine (Lys)
        'LYN': 'LYS',   # Deprotonated lysine (neutral)

        # Arginine (Arg)
        'ARN': 'ARG',   # Deprotonated arginine (rare, GROMACS)

        # Tyrosine (Tyr)
        'TYM': 'TYR',   # Deprotonated tyrosine (GROMACS)

        # Serine (Ser)
        'SEP': 'SER',   # Phosphorylated serine (GROMACS/AMBER)

        # Threonine (Thr)
        'TPO': 'THR',   # Phosphorylated threonine (GROMACS/AMBER)

        # Tyrosine (Tyr)
        'PTR': 'TYR',   # Phosphorylated tyrosine (GROMACS/AMBER)

        # Non-standard names for aspartic and glutamic acids in low pH environments
        'ASH': 'ASP',   # Protonated Asp
        'GLH': 'GLU',   # Protonated Glu
    }
    
    pair = [aa_correction.get(aa, aa) for aa in pair]    
    
    if PY3K:
        with pkg_resources.path('prody.proteins', 'tabulated_energies.txt') as file_path:
            with open(file_path) as f:
                data = np.loadtxt(f, dtype=str)
    else:
        file_path = pkg_resources.resource_filename('prody.proteins', 'tabulated_energies.txt')
        with open(file_path) as f:
            data = np.loadtxt(f, dtype=str)
    
    sources = ["IB_nosolv", "IB_solv", "CS"]
    aa_pairs = []
    
    for row in data:
        aa_pairs.append(row[0]+row[1])
    
    lookup = pair[0]+pair[1]
    
    try:
        data_results = data[np.nonzero(np.array(aa_pairs)==lookup)[0]][0][2:][sources.index(source)]
    except TypeError:
        raise TypeError('Please replace non-standard names of residues with standard names.')

    return data_results


def checkNonstandardResidues(atoms):
    """Check whether the atomic structure contains non-standard residues and inform to replace the name
    to the standard one so that non-standard residues are treated in a correct way while computing
    interactions.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
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
    
    amino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", 
                   "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    aa_list = atoms.select('name CA').getResnames()
    aa_list_nr = atoms.select('name CA').getResnums()
    nonstandard = []
    
    for nr_i,i in enumerate(aa_list):
        if i not in amino_acids:
            nonstandard.append(aa_list[nr_i] + str(aa_list_nr[nr_i]))
    
    if len(nonstandard) > 0:
        LOGGER.info('There are several non-standard residues in the structure.')
        LOGGER.info('Replace the non-standard name in the PDB file with the equivalent name from the standard one if you want to include them in the interactions.')
        LOGGER.info("Residues: {0}.".format(' '.join(nonstandard)))
        return True

    return False


def showPairEnergy(data, **kwargs):
    """Return energies when a list of interactions is given. Energies will be added to each pair of residues 
    at the last position in the list. Energy is based on the residue types and not on the distances.
    The unit of energy is kcal/mol. The energies defined as 'IB_nosolv' (non-solvent-mediated) and
    'IB_solv' (solvent-mediated) are taken from [OK98]_ and 'CS' from InSty paper (under preparation).
    Protonation of residues is not distinguished. The protonation of residues is not distinguished. 
    'IB_solv' and 'IB_nosolv' have RT units and 'CS' has units of kcal/mol.

    Known residues such as HSD, HSE, HIE, and HID (used in MD simulations) are treated as HIS.
    
    :arg data: list with interactions from calcHydrogenBonds() or other types
    :type data: list
    
    :arg energy_list_type: name of the list with energies 
                            default is 'IB_solv'
    :type energy_list_type: 'IB_nosolv', 'IB_solv', 'CS'
    
    
    .. [OK98] Keskin O., Bahar I., Badretdinov A.Y., Ptitsyn O.B., Jernigan R.L., 
    Empirical solvet-mediated potentials hold for both intra-molecular and 
    inter-molecular inter-residues interactions,
    *Protein Science* **1998** 7: 2578–2586.
    """
    
    if not isinstance(data, list):
        raise TypeError('list_of_interactions must be a list of interactions.')

    energy_list_type = kwargs.pop('energy_list_type', 'IB_solv')
    
    for i in data:
        energy = get_energy([i[0][:3], i[3][:3]], energy_list_type)
        i.append(float(energy))
        
    return data


# SignatureInteractions supporting functions
def remove_empty_strings(row):
    """Remove empty strings from a list."""
    return [elem for elem in row if elem != '']

def log_message(message, level="INFO"):
    """Log a message with a specified log level."""
    LOGGER.info("[{}] {}".format(level, message))

def is_module_installed(module_name):
    """Check if a Python module is installed."""
    import importlib.util
    spec = importlib.util.find_spec(module_name)
    return spec is not None

def is_command_installed(command):
    """Check if a command-line tool is installed."""
    import shutil
    return shutil.which(command) is not None

def load_residues_from_pdb(pdb_file):
    """Extract residue numbers and their corresponding one-letter amino acid codes from a PDB file."""
    from prody.atomic.atomic import AAMAP
    structure = parsePDB(pdb_file)
    residues = structure.iterResidues()
    residue_dict = {}

    for res in residues:
        resnum = res.getResnum()
        resname = res.getResname()  # Three-letter amino acid code
        try:
            one_letter_code = AAMAP[resname]
            residue_dict[resnum] = one_letter_code
        except KeyError:
            log_message("Unknown residue: {} at position {}".format(resname, resnum), "WARNING")
            
    return residue_dict

def append_residue_code(residue_num, residue_dict):
    """Return a string with one-letter amino acid code and residue number."""
    aa_code = residue_dict.get(residue_num, "X")  # Use "X" for unknown residues
    return "{}{}".format(aa_code, residue_num)

def process_data(mapping_file, pdb_folder, interaction_func, bond_type, files_for_analysis):
    """Process the mapping file and the PDB folder to compute interaction counts and percentages."""
    log_message("Loading mapping file: {}".format(mapping_file))

    # Load and clean the mapping file
    try:
        mapping = np.loadtxt(mapping_file, delimiter=' ', dtype=str)
    except Exception as e:
        log_message("Error loading mapping file: {}".format(e), "ERROR")
        return None

    mapping = np.where(mapping == '-', np.nan, mapping)
    filtered_mapping = np.array([remove_empty_strings(row) for row in mapping])
    mapping_num = filtered_mapping.astype(float)

    # Load the one-letter amino acid codes from model1.pdb
    import os 
    pdb_model_path = os.path.join(pdb_folder, 'model1.pdb')
    residue_dict = load_residues_from_pdb(pdb_model_path)

    log_message("Processing PDB files in folder: {}".format(pdb_folder))

    tar_bond_ind = []
    processed_files = 0  # To track the number of files successfully processed
    fixed_files = []  # Store paths of fixed files to remove at the end

    for i, files in enumerate(files_for_analysis):
        log_message("Processing file {}: {}".format(i + 1, files))
        
        try:
            coords = parsePDB(files)
            atoms = coords.select('protein')
            bonds = interaction_func(atoms)
            saveInteractionsAsDummyAtoms(atoms, bonds, filename=files.rstrip(files.split('/')[-1])+'INT_'+bond_type+'_'+files.split('/')[-1])
        except Exception as e:
            log_message("Error processing PDB file {}: {}".format(files, e), "ERROR")
            continue

        # If no bonds were found, skip this file
        if len(bonds) == 0:
            log_message("No {} found in file {}, skipping.".format(bond_type, files), "WARNING")
            continue

        processed_files += 1  # Increment successfully processed files
        fixed_files.append(files)

        if processed_files == 1:  # First valid file with bonds determines target indices
            for entries in bonds:
                ent = list(np.sort((int(entries[0][3:]), int(entries[3][3:]))))  # Ensure integers
                tar_bond_ind.append(ent)
            tar_bond_ind = np.unique(np.array(tar_bond_ind), axis=0).astype(int)
            count = np.zeros(tar_bond_ind.shape[0], dtype=int)

        bond_ind = []
        for entries in bonds:
            ent = list(np.sort((int(entries[0][3:]), int(entries[3][3:]))))  # Ensure integers
            bond_ind.append(ent)
        bond_ind = np.unique(np.array(bond_ind), axis=0)

        # Ensure bond_ind is a 2D array
        if bond_ind.ndim == 1:
            bond_ind = bond_ind.reshape(-1, 2)  # Reshape to (n, 2) if it's 1D

        for j, pairs in enumerate(tar_bond_ind):
            ind1_matches = np.where(mapping_num[0] == pairs[0])[0]
            ind2_matches = np.where(mapping_num[0] == pairs[1])[0]

            if ind1_matches.size > 0 and ind2_matches.size > 0:
                if processed_files - 1 < mapping_num.shape[0]:
                    ind1 = mapping_num[processed_files - 1, ind1_matches[0]]
                    ind2 = mapping_num[processed_files - 1, ind2_matches[0]]

                    if not (np.isnan(ind1) or np.isnan(ind2)):
                        index = np.where(np.logical_and(bond_ind[:, 0] == int(ind1), bond_ind[:, 1] == int(ind2)))[0]
                        if index.size != 0:
                            count[j] += 1
                else:
                    log_message("Skipping file {} due to index out of bounds error".format(files), "WARNING")
            else:
                log_message("No matching indices found for {} in {}".format(pairs, files), "WARNING")

    # If no files were successfully processed or no bonds were found
    if processed_files == 0 or len(tar_bond_ind) == 0:
        log_message("No valid {} entries found across all PDB files.".format(bond_type), "ERROR")
        return None

    count_reshaped = count.reshape(-1, 1)
    count_normalized = (count / processed_files * 100).reshape(-1, 1)

    # Modify tar_bond_ind to append the amino acid code before residue index
    tar_bond_with_aa = []
    for bond in tar_bond_ind:
        res1 = append_residue_code(bond[0], residue_dict)  # Append AA code for Res1
        res2 = append_residue_code(bond[1], residue_dict)  # Append AA code for Res2
        tar_bond_with_aa.append([res1, res2])

    tar_bond_with_aa = np.array(tar_bond_with_aa)

    # Combine tar_bond_with_aa with count and percentage
    output_data = np.hstack((tar_bond_with_aa, count_reshaped, count_normalized))

    log_message("Finished processing {} PDB files.".format(processed_files))
    output_filename = '{}_consensus.txt'.format(bond_type)

    # Save the result with amino acid codes and numeric values
    np.savetxt(output_filename, output_data, fmt='%s %s %s %s', delimiter=' ', 
           header='Res1 Res2 Count Percentage', comments='')

    return output_data, fixed_files  # Return fixed_files to remove later


def plot_barh(result, bond_type, **kwargs):
    """Plot horizontal bar plots of percentages of interactions, splitting the data into fixed-sized plots.

    :arg n_per_plot: The number of results per one plot
        Default is 20 if set to None 
    :type n_per_plot: int
    
    :arg min_height: Size of the bar plot
        Default is 8
    :type min_height: int
    """

    import matplotlib.pylab as plt
    plt.rcParams.update({'font.size': 20})

    n_per_plot = kwargs.pop('n_per_plot', None)
    min_height = kwargs.pop('min_height', 8)
    
    if n_per_plot is None:
        n_per_plot = 20

    # Error handling if result is None or empty
    if result is None or len(result) == 0:
        log_message("Skipping plot for {} due to insufficient data.".format(bond_type), "ERROR")
        return

    num_entries = result.shape[0]
    num_plots = (num_entries + n_per_plot - 1) // n_per_plot  # Number of plots required

    for plot_idx in range(num_plots):
        # Slice the data for the current plot (take the next 'n_per_plot' entries)
        start_idx = plot_idx * n_per_plot
        end_idx = min((plot_idx + 1) * n_per_plot, num_entries)
        result_chunk = result[start_idx:end_idx]

        log_message("Plotting entries {} to {} for {}.".format(start_idx + 1, end_idx, bond_type))
        # Use residue numbers for y-axis labels
        y_labels = ["{}-{}".format(str(row[0]), str(row[1])) for row in result_chunk]
        percentage_values = result_chunk[:, 3].astype('float')

        norm = plt.Normalize(vmin=0, vmax=100)
        cmap = plt.cm.get_cmap('coolwarm')

        # Set the figure height with a minimum height of `min_height` for small datasets
        fig_height = max(min_height, len(result_chunk) * 0.4)
        plt.figure(figsize=(18, fig_height))
        bars = plt.barh(y_labels, percentage_values, color=cmap(norm(percentage_values)))

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        # Adjust color bar height to match the number of entries
        plt.colorbar(sm, label='Percentage', fraction=0.02, pad=0.04)

        plt.ylim(-1, result_chunk.shape[0])
        plt.ylabel('{} Pairs of Residue Numbers'.format(bond_type))
        plt.xlabel('Percentage')
        plt.title('Persistence of {} (entries {}-{})'.format(bond_type, start_idx + 1, end_idx))
        # Save each plot with an incremented filename for multiple plots
        output_plot_file = '{}_plot_part{}.png'.format(bond_type, plot_idx + 1)
        log_message("Saving plot to: {}".format(output_plot_file))
        plt.savefig(output_plot_file)
        plt.close()
        log_message("Plot saved successfully.")


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
        Default all the protein structure is used
    :type selection: str
    
    :arg sasa_cutoff: cutoff for SASA values
        Default is 0.0
    :type sasa_cutoff: float, int

    :arg resnames: residues name included
        Default is False
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


def calcBSA(atoms, separate=False, selection='protein and noh', groups=None):
    """Calculate buried solvent-accessible surface area (BSA) of protein
    complexes using :func:`.calcSASA`.

    BSA is calculated as half of the difference between the sum of solvent-
    accessible surface areas (SASA) of isolated components and the SASA of
    the complete complex.
    The global BSA for N groups is (sum(SASA_i) - SASA_complex)/2.
    This represents half the sum of lost solvent accessibility on all partners.

    By default, individual components are identified using chain IDs.
    Alternatively, molecular components can be defined using the *groups*
    argument. At least two non-overlapping groups are required.

    If separate=False, return the global BSA (float, Å²).
    If separate=True, return a dict containing bsa_total, delta_sasa,
    sasa_complex, sasa_isolated, pairwise_bsa (isolated pairs) and
    bsa_vs_rest (each group versus all others considered together). Pairwise
    BSA must not be summed to reproduce the global BSA.

    By default solvent, ions, and hydrogens are excluded, just as for a
    protein-only calcSASA calculation.
    
    For complexes containing more than two components, the sum of
    pairwise BSA values does not necessarily equal the total BSA.
    Pairwise BSA is calculated for each isolated pair, whereas
    total BSA considers the complete molecular complex.

    SASA calculations require the compiled ProDy hpb.so extension.
    All surface areas are reported in Å².
    
    :arg atoms: Atomic object containing the molecular complex.
    :type atoms: :class:`.Atomic`

    :arg separate: If True, return detailed SASA and BSA values for
        individual components and interfaces. If False, return only
        the total BSA. Default is False.
    :type separate: bool

    :arg selection: Atom selection used for SASA calculations.
        Default is 'protein and noh'.
    :type selection: str

    :arg groups: Dictionary mapping component names to ProDy selection
        strings. If None, components are identified automatically using
        chain IDs. Groups must be non-overlapping and cover all atoms
        included in the selection. Default is None.
    :type groups: dict or None

    :returns: Total BSA (float) if *separate* is False. Otherwise,
        a dictionary containing:

        * **bsa_total** -- total buried surface area of the complex.
        * **delta_sasa** -- difference between the sum of isolated
          SASAs and the SASA of the complex (twice the total BSA).
        * **sasa_complex** -- SASA of the complete complex.
        * **sasa_isolated** -- SASA of each isolated component.
        * **pairwise_bsa** -- BSA between each pair of isolated components.
        * **bsa_vs_rest** -- BSA between each component and the remaining
          components considered together.
        * **groups** -- component names and their selection strings.
        * **unit** -- surface area unit ('A^2').

    Examples:

    Calculate the total BSA of a protein complex:
    >>> from prody import *
    >>> atoms = parsePDB('1a2k')
    >>> bsa = calcBSA(atoms)

    Calculate BSA for individual interfaces:
    >>> results = calcBSA(atoms, separate=True)
    >>> results['bsa_total']
    >>> results['pairwise_bsa']
    >>> results['bsa_vs_rest']

    Calculate BSA between two user-defined molecular components:
    >>> groups = {'component1': 'chain A B',
                  'component2': 'chain C D E'}
    >>> bsa = calcBSA(atoms, groups=groups) """
    
    from math import fsum, isfinite
    from itertools import combinations
    from prody.proteins.interactions import calcSASA
    from prody.atomic import Selection

    base = atoms.select(selection)
    if base is None or base.numAtoms() == 0:
        raise ValueError('The requested selection contains no atoms.')

    if groups is None:
        chids = [str(ch).strip() for ch in base.getChids()]
        if not chids or '' in chids:
            raise ValueError('Blank chain IDs: provide explicit groups.')
        groups = {ch: 'chain ' + ch for ch in sorted(set(chids))}
    elif not isinstance(groups, dict) or not groups:
        raise TypeError('groups must be a nonempty dict: label -> selection string.')

    if len(groups) < 2:
        raise ValueError('At least two groups are needed to calculate BSA.')

    chosen = {}
    chosen_ids = {}
    used = set()
    for label, expr in groups.items():
        if not isinstance(label, str) or not label.strip():
            raise ValueError('Group labels must be nonempty strings.')
        if not isinstance(expr, str) or not expr.strip():
            raise ValueError('Each group must have a nonempty selection expression.')
        sel = base.select('(' + expr + ')')
        if sel is None or sel.numAtoms() == 0:
            raise ValueError('Group %r is empty under selection %r.' % (label, selection))
        ids = set(int(i) for i in sel.getIndices())
        if used.intersection(ids):
            raise ValueError('Group %r overlaps another group.' % label)
        used.update(ids)
        chosen[label] = expr
        chosen_ids[label] = ids

    all_ids = set(int(i) for i in base.getIndices())
    if used != all_ids:
        raise ValueError('%d atoms in selection are not assigned to any group.' %
                         len(all_ids - used))

    parent = base.getAtomGroup()
    acsi = base.getACSIndex()

    def calculate(indices):
        if not indices:
            raise ValueError('Cannot calculate SASA for an empty atom subset.')
        selected = Selection(parent, sorted(indices), 'all', acsi=acsi)
        subset = selected.copy()
        values = calcSASA(subset, selection='all', sasa_cutoff=0.0)
        if values is None or len(values) == 0:
            raise ValueError('calcSASA returned no values for %d atoms' % len(indices))
        result = float(fsum(float(x) for x in values))
        if not isfinite(result) or result < 0:
            raise ValueError('Invalid SASA for %d atoms' % len(indices))
        return result

    def ids_of(names):
        result = set()
        for name in names:
            result.update(chosen_ids[name])
        return result

    names = list(chosen)
    isolated = {name: calculate(chosen_ids[name]) for name in names}
    complex_sasa = calculate(all_ids)
    delta = float(fsum(isolated.values()) - complex_sasa)

    if delta < -0.01:
        raise ValueError('Negative buried area (%.3f Å²); check SASA and groups.' % delta)
    if -0.01 <= delta < 0:
        delta = 0.0 
    total = delta / 2.0
    if not separate:
        return total

    pairwise = {}
    for a, b in combinations(names, 2):
        sasab = calculate(ids_of((a, b)))
        pairwise[(a, b)] = (isolated[a] + isolated[b] - sasab) / 2.0

    vs_rest = {}
    for name in names:
        rest = [n for n in names if n != name]
        sasa_rest = calculate(ids_of(rest)) if len(rest) > 1 else isolated[rest[0]]
        vs_rest[name] = (isolated[name] + sasa_rest - complex_sasa) / 2.0

    return {'bsa_total': total,
            'delta_sasa': delta,
            'sasa_complex': complex_sasa,
            'sasa_isolated': isolated,
            'pairwise_bsa': pairwise,
            'bsa_vs_rest': vs_rest,
            'groups': chosen,
            'unit': 'A^2'}


def calcVolume(atoms, **kwargs):
    """Provide information about volume for each residue/molecule/chain
    or other selection". To use this function compiled hpb.so is needed.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg selection: selection string
        Default all the protein structure is used
    :type selection: str
    
    :arg volume_cutoff: cutoff for volume
        Default is 0.0 to include all the results
    :type sasa_volume: float, int

    :arg split_residues: it will provide values for each residue
        Default is False
    :type split_residues: bool

    :arg resnames: residues name included
        Default is False
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
        angle and angleDA also work
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
    angleDHA = kwargs.pop('angleDHA', angleDHA)
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
            ligand_sel = atoms.select('all not protein and not ion')
            if ligand_sel:
                ligand_name = list(set(ligand_sel.getResnames()))[0]
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
    if atoms_KRED is None:
        LOGGER.warn('There are no side chain heavy atoms for residues K, R, E, D and H, so not salt bridges are calculated')
        return []

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

    angle_min_PS = kwargs.pop('angle_min_PS', 0)
    angle_min = kwargs.pop('angle_min', angle_min_PS)

    angle_max_PS = kwargs.pop('angle_max_PS', 360)
    angle_max = kwargs.pop('angle_max', angle_max_PS)
    
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
    with mp.Pool() as pool:
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

            for (resnum, chid) in sele2_nr:
                sele2_filter = sele2.select('resnum {0} and chid {1}'.format(resnum, chid))

                if sele1_name[0] in aromatic:
                    # avoid double counting pi stacking and don't include same residue interactions
                    sele2_filter = sele2_filter.select('all and not (resname TYR PHE TRP or resid '+str(i[0])+' and chain '+i[1]+')')
                elif sele1_name[0] not in aromatic and i in sele2_nr:
                    # don't include same residue interactions but don't worry about double counting pi stacking
                    sele2_filter = sele2_filter.select('all and not (resid '+str(i[0])+' and chain '+i[1]+')')

                if sele2_filter != None:
                    listOfAtomToCompare = cleanNumbers(findNeighbors(sele1, distA, sele2_filter))

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
    
    :arg distDB: non-zero value, maximal distance between atoms of cysteine residues.
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
    
    atoms_SG = atoms.select('protein and resname CYS and name SG')
    DisulfideBonds_list = []

    if atoms_SG is None:
        LOGGER.info('Lack of cysteines in the structure.')
    else:
        atoms_SG_res = list(set(zip(atoms_SG.getResnums(), atoms_SG.getChids())))
    
        LOGGER.info('Calculating disulfide bonds.')
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
                                                    minDistancePair[2]+'_'+str(minDistancePair[0]),
                                                    sele1_new.getChids()[0],
                                                    sele2_new.getResnames()[0]+str(sele2_new.getResnums()[0]),
                                                    minDistancePair[3]+'_'+str(minDistancePair[1]),
                                                    sele2_new.getChids()[0],
                                                    round(minDistancePair[-1],4), round(float(diheAng),4)])

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
    using default parameters or those from kwargs.

    :arg max_proc: maximum number of processes to use
        default is half of the number of CPUs
    :type max_proc: int
    """
    
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
    max_proc = kwargs.pop('max_proc', mp.cpu_count()//2)

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
        
        if isinstance(trajectory, Trajectory):
            nfi = trajectory._nfi
            trajectory.reset()
        numFrames = trajectory._n_csets
        
        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame+1]
        
        atoms_copy = atoms.copy()
        def analyseFrame(j0, frame0, interactions_all):
            LOGGER.info('Frame: {0}'.format(j0))
            atoms_copy.setCoords(frame0.getCoords())
            protein = atoms_copy.select('protein')
            interactions = interactions_dic[interaction_type](protein, **kwargs)
            interactions_all.append(interactions)

        if max_proc == 1:
            interactions_all = []
            for j0, frame0 in enumerate(traj, start=start_frame):
                analyseFrame(j0, frame0, interactions_all)
        else:
            with mp.Manager() as manager:
                interactions_all = manager.list()

                j0 = start_frame
                while j0 < traj.numConfs()+start_frame:

                    processes = []
                    for _ in range(max_proc):
                        frame0 = traj[j0-start_frame]
                        
                        p = mp.Process(target=analyseFrame, args=(j0, frame0,
                                                                 interactions_all))
                        p.start()
                        processes.append(p)

                        j0 += 1
                        if j0 >= traj.numConfs()+start_frame:
                            break

                    for p in processes:
                        p.join()

                interactions_all = interactions_all[:]

        if isinstance(trajectory, Trajectory):
            trajectory._nfi = nfi
    
    else:
        if atoms.numCoordsets() > 1:
            def analyseFrame(i, interactions_all):
                LOGGER.info('Model: {0}'.format(i+start_frame))
                atoms.setACSIndex(i+start_frame)
                protein = atoms.select('protein')
                interactions = interactions_dic[interaction_type](protein, **kwargs)
                interactions_all.append(interactions)

            if stop_frame == -1:
                stop_frame = atoms.numCoordsets()

            if max_proc == 1:
                interactions_all = []
                for i in range(len(atoms.getCoordsets()[start_frame:stop_frame+1])):
                    analyseFrame(i, interactions_all)
            else:
                with mp.Manager() as manager:
                    interactions_all = manager.list()

                    i = start_frame
                    while i < len(atoms.getCoordsets()[start_frame:stop_frame+1]):
                        processes = []
                        for _ in range(max_proc):
                            p = mp.Process(target=analyseFrame, args=(i, interactions_all))
                            p.start()
                            processes.append(p)

                            i += 1
                            if i >= len(atoms.getCoordsets()[start_frame:stop_frame]):
                                break

                        for p in processes:
                            p.join()

                    interactions_all = interactions_all[:]
        else:
            LOGGER.info('Include trajectory or use multi-model PDB file.')
    
    return interactions_all


def calcProteinInteractions(atoms, **kwargs):
    """Compute all protein interactions (shown below).
        (1) Hydrogen bonds
        (2) Salt Bridges
        (3) RepulsiveIonicBonding 
        (4) Pi stacking interactions
        (5) Pi-cation interactions
        (6) Hydrophobic interactions
        (7) Disulfide Bonds

    kwargs can be passed on to the underlying functions as described
    in their documentation. For example, distDA and angleDHA can be used
    to control hydrogen bonds, or distA and angle can be used across types.
    
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

    LOGGER.info('Calculating interactions.')
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
        distDA also works
    :type distA: int, float
    
    :arg angle: non-zero value, maximal (180 - D-H-A angle) (donor, hydrogen, acceptor).
        default is 40.
        angleDHA also works
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
        distSB also works
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
            distRB also works
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
                distPS also works
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
                distPC also works
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
        distHPh also works
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
    
    :arg distA: non-zero value, maximal distance between atoms of cysteine residues.
        default is 2.5.
        distDB also works
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
        Default is 0.1.
    :type cutoff: int, float

    :arg code: representation of the residues, 3-letter or 1-letter
        Default is 3-letter
    :type code: str

    :arg multiple_chains: display chain name for structure with many chains
        Default is False
    :type multiple_chains: bool
    
    :arg edge_cmap: color of the residue connection
        Default is plt.cm.Blues (blue color)
    :type edge_cmap: str

    :arg node_size: size of the nodes which describes residues
        Default is 300
    :type node_size: int
    
    :arg node_distance: value which will scale residue-residue interactions
        Default is 5
    :type node_distance: int

    :arg font_size: size of the font
        Default is 14
    :type font_size: int

    :arg seed: random number which affect the distribution of residues
        Default is 42
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
    

    if isinstance(statistics, int) or isinstance(statistics, str) or isinstance(statistics, Atomic):
        raise TypeError('input data must be a list, use calcStatisticsInteractions to obtain statistics for a particular interaction type')

    if isinstance(statistics, InteractionsTrajectory) or isinstance(statistics, Interactions):
        raise TypeError('use calcStatisticsInteractions to obtain statistics for a particular interaction type')
        
    else:
        if len(statistics[0]) != 5:
            raise TypeError('input data must be a list obtained from calcStatisticsInteractions')

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
    
    
def showInteractionsHist(statistics, **kwargs):
    """Return information about residue-residue interactions as a bar plot.
    
    :arg statistics: Results obtained from calcStatisticsInteractions analysis
    :type statistics: list
    
    :arg clip: maxmimal number of residue pairs on the bar plot
        Default is 20
    :type clip: int

    :arg font_size: size of the font
        Default is 18
    :type font_size: int
    """

    if isinstance(statistics, int) or isinstance(statistics, str) or isinstance(statistics, Atomic):
        raise TypeError('input data must be a list, use calcStatisticsInteractions to obtain statistics for a particular interaction type')

    if isinstance(statistics, InteractionsTrajectory) or isinstance(statistics, Interactions):
        raise TypeError('use calcStatisticsInteractions to obtain statistics for a particular interaction type')
        
    else:
        if len(statistics[0]) != 5:
            raise TypeError('input data must be a list obtained from calcStatisticsInteractions')

    import matplotlib.pyplot as plt
    import numpy as np

    clip = kwargs.pop('clip', 20)
    font_size = kwargs.pop('font_size', 18)
    
    labels = [row[0] for row in statistics]
    values = [row[2] for row in statistics]
    std_devs = [row[3] for row in statistics]
    colors = [row[1] for row in statistics]

    sorted_indices = np.argsort(colors)[-clip:]
    labels = [labels[i] for i in sorted_indices]
    values = [values[i] for i in sorted_indices]
    std_devs = [std_devs[i] for i in sorted_indices]
    colors = [colors[i] for i in sorted_indices]

    norm = plt.Normalize(min(colors), max(colors))
    cmap = plt.get_cmap("Blues")
    
    num_labels = len(labels)
    height = max(6, num_labels * 0.4)

    fig, ax = plt.subplots(figsize=(8, height))
    bars = ax.barh(labels, values, xerr=std_devs, color=cmap(norm(colors)))
    
    for bar in bars:
        bar.set_edgecolor('black')
    
    ax.set_xlabel('Average distance [\u00C5]')
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='Weight')
    ax.tick_params(axis='both', labelsize=font_size)
    plt.tight_layout()
        
    if SETTINGS['auto_show']:
        showFigure()
    
    
def calcStatisticsInteractions(data, **kwargs):
    """Return the statistics of interactions from PDB Ensemble or trajectory including:
    (1) the weight for each residue pair: corresponds to the number of counts divided by the 
    number of frames (values >1 are obtained when the residue pair creates multiple contacts); 
    (2) average distance of interactions for each pair [in Ang], 
    (3) standard deviation [Ang.],
    (4) Energy [in kcal/mol] that is not distance dependent. Energy by default is solvent-mediated
    from [OK98]_ ('IB_solv') in RT units. To use non-solvent-mediated (residue-mediated) entries ('IB_nosolv')
    from [OK98]_ in RT units or solvent-mediated values obtained from MD for InSty paper ('CS', under preparation)
    in kcal/mol, change `energy_list_type` parameter.
    If energy information is not available, please check whether the pair of residues is listed in 
    the "tabulated_energies.txt" file, which is localized in the ProDy directory.
        
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
            
            try:
                stats[element] = {
                    "stddev": np.round(np.std(values),6),
                    "mean": np.round(np.mean(values),6),
                    "weight": np.round(float(len(values))/len(data), 6),
                    "energy": get_energy([element.split('-')[0][:3], element.split('-')[1][:3]], energy_list_type)
                }
            except:
                LOGGER.warn('energy information is not available for ', element)
                stats[element] = {
                    "stddev": np.round(np.std(values),6),
                    "mean": np.round(np.mean(values),6),
                    "weight": np.round(float(len(values))/len(data), 6)
                }

    statistic = []
    unit = 'RT' if energy_list_type in ['IB_solv', 'IB_nosolv'] else 'kcal/mol'
    for key, value in stats.items():
        if float(value['weight']) > weight_cutoff:
            LOGGER.info("Statistics for {0}:".format(key))
            LOGGER.info("  Average [Ang.]: {}".format(value['mean']))
            LOGGER.info("  Standard deviation [Ang.]: {0}".format(value['stddev']))
            LOGGER.info("  Weight: {0}".format(value['weight']))
            try:
                LOGGER.info("  Energy [{0}]: {1}".format(unit, value['energy']))
                statistic.append([key, value['weight'], value['mean'], value['stddev'], value['energy']])
            except:
                statistic.append([key, value['weight'], value['mean'], value['stddev']])
        else: pass
    
    statistic.sort(key=lambda x: x[1], reverse=True)
    
    if statistic == []:
        LOGGER.info("No meaningful interactions found. Decrease weight_cutoff to obtain some results (default is 0.2).")
        
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


def saveInteractionsAsDummyAtoms(atoms, interactions, filename, **kwargs):
    '''Creates a PDB file which will contain protein structure and dummy atoms that will be placed between pairs
    of interacting residues.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg interactions: list of interactions
    :type interactions: list

    :arg filename: name of the PDB file which will contain dummy atoms and protein structure
    :type filename: str 
    
    :arg RESNAME_dummy: resname of the dummy atom, use 3-letter name
                        be default is 'DUM'
    :type RESNAME_dummy: str '''


    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
                            
    RESNAME_dummy = kwargs.pop('RESNAME_dummy', 'DUM')
    
    def calcDUMposition(coord1, coord2):
        midpoint = [
            (coord1[0] + coord2[0]) / 2,
            (coord1[1] + coord2[1]) / 2,
            (coord1[2] + coord2[2]) / 2
        ]
        return midpoint

    all_DUMs = []
    atoms_ = atoms.copy()

    for i in interactions:
        if len(i[1].split('_')) <= 3:
            res1_name = 'chain '+i[2]+' and resname '+i[0][:3]+' and resid '+i[0][3:]+' and index '+' '.join(i[1].split('_')[1:])
            res1_coords = calcCenter(atoms.select(res1_name))  
        
        if len(i[1].split('_')) > 3:
            res1_name = 'chain '+i[2]+' and resname '+i[0][:3]+' and resid '+i[0][3:]+' and index '+' '.join(i[1].split('_'))
            res1_coords = calcCenter(atoms.select(res1_name))

        if len(i[4].split('_')) <= 3:
            res2_name = 'chain '+i[5]+' and resname '+i[3][:3]+' and resid '+i[3][3:]+' and index '+' '.join(i[4].split('_')[1:])
            res2_coords = calcCenter(atoms.select(res2_name))
            
        if len(i[4].split('_')) > 3:     
            res2_name = 'chain '+i[5]+' and resname '+i[3][:3]+' and resid '+i[3][3:]+' and index '+' '.join(i[4].split('_'))
            res2_coords = calcCenter(atoms.select(res2_name))

        all_DUMs.append(calcDUMposition(res1_coords, res2_coords))
    
    if all_DUMs == []:
        LOGGER.info('Lack of interactions')
    else:
        LOGGER.info('Creating file with dummy atoms')
        dummyAtoms = AtomGroup()
        coords = array([all_DUMs], dtype=float)
        dummyAtoms.setCoords(coords)
        dummyAtoms.setNames([RESNAME_dummy]*len(dummyAtoms))
        dummyAtoms.setResnums(range(1, len(dummyAtoms)+1))
        dummyAtoms.setResnames([RESNAME_dummy]*len(dummyAtoms))

        writePDB(filename, atoms_+dummyAtoms)


def listLigandInteractions(PLIP_output, **kwargs):
    """Create a list of interactions from PLIP output created using calcLigandInteractions().
    Results can be displayed in VMD. 
    
    :arg PLIP_output: Results from PLIP for protein-ligand interactions.
    :type PLIP_output: PLIP object obtained from calcLigandInteractions()
    
    :arg output: parameter to print the interactions on the screen
                 while analyzing the structure (True | False)
                 Default is False
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


def showProteinInteractions(atoms, interactions, color='red',**kwargs):
    """Display protein interactions in an inline py3Dmol viewer.
    
    Different types of interactions can be saved separately (color can be selected) 
    or all at once for all types of interactions (hydrogen bonds - blue, salt bridges - yellow,
    pi stacking - green, cation-pi - orangem, hydrophobic - silver, and disulfide bonds - black).

    kwargs are passed on to showProtein
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg interactions: List of interactions for protein interactions.
    :type interactions: List of lists
    
    :arg color: color to draw interactions ,
                 used only for single interaction type.
                default **"red"**
    :type color: str
    """    

    import sys        
    if 'py3Dmol' not in sys.modules: 
            LOGGER.warn('py3Dmol not loaded. No visualization will be displayed.')
            return None

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
    
    view = showProtein(atoms, **kwargs);
    view.addStyle({'stick':{'radius':0.1}})

    def singleInteraction(view, interaction, color='blue'):
        """Creates cylinders for the interactions.
        """
        
        for nr_i,i in enumerate(interaction):
            try:
                at1 = atoms.select('index '+' '.join([k for k in i[1].split('_') if k.isdigit() ] ))
                at1xyz = calcCenter(at1.getCoords())
                at2 = atoms.select('index '+' '.join([kk for kk in i[4].split('_') if kk.isdigit() ] ))
                at2xyz = calcCenter(at2.getCoords())
                            
                view.addCylinder({'start':{'x':at1xyz[0],'y':at1xyz[1],'z':at1xyz[2]},
                              'end':{'x':at2xyz[0],'y':at2xyz[1],'z':at2xyz[2]},
                              'dashed': True,
                              'color':color});

            except: LOGGER.info("There was a problem.")
     
    if len(interactions) == 7 and isinstance(interactions[0][0], list):
        # For all seven types of interactions at once
        # HBs_calculations, SBs_calculations, SameChargeResidues, Pi_stacking, Pi_cation, Hydroph_calculations, Disulfide Bonds
        colors = ['blue', 'yellow', 'red', 'green', 'orange', 'silver', 'black']
        
        for nr_inter,inter in enumerate(interactions):
            singleInteraction(view, inter, color=colors[nr_inter])

    elif interactions == []:
        LOGGER.info("Lack of results")

    elif len(interactions[0]) == 0:
        LOGGER.info("Lack of results")
        
    else:
        singleInteraction(view,interactions,color)

    return view


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
     
    if len(interactions) == 7 and isinstance(interactions[0][0], list):
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

def showLigandInteraction(atoms, interactions, **kwargs):
    """Display information from PLIP for ligand-protein interactions in a py3dmol viewer.
    
    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`
    
    :arg interactions: List of interactions lists for protein-ligand interactions.
    :type interactions: list       

    To obtain protein-ligand interactions:
    >>> calculations = calcLigandInteractions(atoms)
    >>> interactions = listLigandInteractions(calculations) """

    import sys        
    if 'py3Dmol' not in sys.modules: 
            LOGGER.warn('py3Dmol not loaded. No visualization will be displayed.')
            return None

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

    view = showProtein(atoms,**kwargs)
    view.addStyle({'stick':{'radius':0.1}})

    if len(interactions[0]) >= 10: 
        dic_color = {'HBs':'blue','PiStack':'green','SBs':'yellow','PiCat':'orange',
                     'HPh':'silver','watBridge':'cyan'}
        
        for i in interactions:
            color = dic_color[i[0]]

            if i[0] == 'waterbridge':
                hoh_id = atoms.select('x `'+str(i[11][0])+'` and y `'+str(i[11][1])+'` and z `'+str(i[11][2])+'`').getResnums()[0]

                view.addCylinder({'start':{'x':i[9][0],'y':i[9][1],'z':i[9][2]},
                              'end':{'x':i[11][0],'y':i[11][1],'z':i[11][2]},
                              'color':color, 'dashed':True})
                view.addCylinder({'start':{'x':i[10][0],'y':i[10][1],'z':i[10][2]},
                              'end':{'x':i[11][0],'y':i[11][1],'z':i[11][2]},
                              'color':color, 'dashed': True})                
            else:
                view.addCylinder({'start':{'x':i[9][0],'y':i[9][1],'z':i[9][2]},
                              'end':{'x':i[10][0],'y':i[10][1],'z':i[10][2]},
                              'color':color, 'dashed': True})
    return view


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
                            Default is "protein" 
    :type protein_selection: str
    
    :arg ligand_selection: selection string for ligand,
                           e.g. "resname ADP",
                           Default is "all not protein"
    :type ligand_selection: str
    
    :arg ligand_selection: scoring function (vina or vinardo)
                           Default is "vina"
    
    :type ligand_selection: str
    
    :arg atom_terms: write per-atom interaction term values.
                     It will provide more information as dictionary for each frame/model,
                     and details for atom terms will be saved in terms_*frame_number*.txt,    
                     Default is False

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
    
        if isinstance(trajectory, Trajectory):
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
        
        if isinstance(trajectory, Trajectory):
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


def calcSminaPerAtomInteractions(atoms, list_terms, detailed=False,
                                 ligand_selection=None, match_tolerance=0.15):
    """Calculate contributions of individual ligand atoms to SMINA scoring function terms.
    SMINA package is described in [DRK13]_.

    The function reads atom-resolved scoring data generated by SMINA using the 
    *--atom_terms* option. Ligand atoms are identified by matching their coordinates 
    in the SMINA output with the coordinates of the corresponding atoms in the input 
    structure.

    :arg atoms: Atomic object containing the ligand coordinates.
        The coordinates must correspond to the structure used to generate the SMINA 
        atom-resolved scoring data.
    :type atoms: :class:`.Atomic`

    :arg list_terms: Path to a file generated using the SMINA --atom_terms option, 
        or a list of such files. When multiple files are provided, they must contain
        the same scoring terms and ligand atoms in identical order.
    :type list_terms: str, list

    :arg detailed: If True, calculate and return contributions of individual ligand 
        atoms to each scoring term, together with total contributions and atom-mapping 
        quality control. If False, return the original ProDy output format,
        containing summed atomic scoring values for each input file.
        Default is False.
    :type detailed: bool

    :arg ligand_selection: ProDy selection string identifying the ligand whose atoms 
        will be associated with the SMINA atom-resolved scoring data.
        Required when detailed=True. When detailed=False and no selection is provided,
        all atoms in the input structure are considered.
        Default is None.
    :type ligand_selection: str, None

    :arg match_tolerance: Maximum allowed difference between the coordinates of a ligand 
        atom reported in the SMINA --atom_terms output and the coordinates of the 
        corresponding atom in the input structure. The distance is expressed in Angstroms.
        Both sets of coordinates must be expressed in the same reference frame.
        Default is 0.15 Angstrom.
    :type match_tolerance: float

    :returns: Dictionary containing atom-resolved SMINA scoring information.

        When detailed=False, the dictionary contains atom
        identifiers as keys and lists of summed scoring values
        from the provided files as values.

        When detailed=True, the dictionary contains:
        - terms: names of scoring function terms.
        - atoms: atom identifiers, mapping information,
          and contributions to individual scoring terms.
        - totals: summed contributions for each scoring term.
        - global_terms: terms identified as potentially
          describing global properties of the ligand.
        - qc: atom-mapping quality control information.

        Individual atom records contain the following fields:
        - smina_row: position of the atom in the SMINA output.
        - smina_type: SMINA atom type.
        - ligand_residue: ligand residue identifier.
        - ligand_atom: ligand atom name.
        - mapping_status: coordinate-matching status.
        - mapping_distance_A: coordinate difference in Angstroms.
        - scoring-term values: contributions reported by SMINA.
        - sum_atom_terms: sum of the atom's scoring-term values.

        The qc dictionary provides the number of successfully
        matched atoms, the number of SMINA atom records,
        mapping-status counts, the absolute magnitude of
        unmatched scoring contributions, and the overall
        matching status (OK or REVIEW).
    :rtype: dict

    The sum of atom-resolved scoring contributions does not necessarily correspond 
    to the final SMINA binding affinity, because the scoring function may include torsional
    corrections and other global terms.

    Atom identification is based on the coordinates. Unmatched atoms with nonzero scoring 
    contributions require additional verification.

    The detailed analysis currently supports one static structure and one SMINA atom-terms file.  
    """

    import re
    from pathlib import Path
    import numpy as np
    from builtins import any, sum

    if isinstance(list_terms, (str, Path)):
        term_files = [Path(list_terms)]
    elif isinstance(list_terms, (list, tuple)) and list_terms:
        term_files = [Path(filename) for filename in list_terms]
    else:
        raise TypeError('list_terms must be a filename or a non-empty list of filenames')

    if not np.isfinite(match_tolerance) or match_tolerance <= 0:
        raise ValueError('match_tolerance must be a positive finite distance in Angstrom')
    if detailed and len(term_files) != 1:
        raise ValueError('detailed=True requires exactly one SMINA terms file')
    if detailed and not ligand_selection:
        raise ValueError('ligand_selection is required when detailed=True')

    def read_terms_file(filename):
        """Parse term names and ligand-atom rows from single SMINA file."""
        
        lines = filename.read_text(encoding='utf-8').splitlines()
        term_names = None

        for line in lines:
            fields = line.strip().lstrip('#').strip().split()
            if fields and fields[0].lower() == 'atomid' and 'pos' in fields:
                term_names = fields[fields.index('pos') + 1:]
                break

        if not term_names or len(term_names) != len(set(term_names)):
            raise ValueError('{}: missing or invalid SMINA atom-term header'.format(filename))

        atom_records = []
        for line_number, line in enumerate(lines, start=1):
            match = re.search(r'<([^<>]+)>', line)
            if match is None:
                continue
            try:
                position = np.asarray([float(v) for v in match.group(1).split(',')], dtype=float)
                term_values = np.asarray([float(v) for v in line[match.end():].split()], dtype=float)
            except ValueError as error:
                raise ValueError('{}:{}: invalid numeric SMINA atom record'.format(
                    filename, line_number)) from error

            if position.shape != (3,) or not np.all(np.isfinite(position)):
                raise ValueError('{}:{}: invalid 3D atom position'.format(filename, line_number))
            if term_values.shape != (len(term_names),) or not np.all(np.isfinite(term_values)):
                raise ValueError('{}:{}: expected {} finite scoring-term values'.format(
                    filename, line_number, len(term_names)))

            atom_records.append({'atom_type': line[:match.start()].strip(),
                                 'position': position,
                                 'term_values': term_values})

        if not atom_records:
            raise ValueError('{}: no atom records found'.format(filename))
        return term_names, atom_records

    term_names, reference_records = read_terms_file(term_files[0])
    records_by_file = [reference_records]
    for filename in term_files[1:]:
        other_term_names, other_records = read_terms_file(filename)
        if other_term_names != term_names:
            raise ValueError('{}: scoring-term names or order differ'.format(filename))
        if len(other_records) != len(reference_records):
            raise ValueError('{}: ligand-atom record count differs'.format(filename))
        if any(record['atom_type'] != ref['atom_type']
               for record, ref in zip(other_records, reference_records)):
            raise ValueError('{}: ligand-atom types or order differ'.format(filename))
        records_by_file.append(other_records)

    ligand = atoms.select(ligand_selection) if ligand_selection else atoms
    if ligand is None or ligand.numAtoms() == 0:
        raise ValueError('ligand_selection does not select any atoms')
    coordinates = np.asarray(ligand.getCoords(), dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or not np.all(np.isfinite(coordinates)):
        raise ValueError('Input ligand must have finite 3D coordinates')

    atom_names = ligand.getNames()
    residue_names = ligand.getResnames()
    residue_numbers = ligand.getResnums()
    chain_ids = ligand.getChids()
    used_atom_indices = set()
    output_rows = []
    legacy_output = {}

    for record_index, record in enumerate(reference_records):
        distances = np.linalg.norm(coordinates - record['position'], axis=1)
        nearest_index = int(np.argmin(distances))
        nearest_distance = float(distances[nearest_index])
        possible_indices = np.flatnonzero(distances <= match_tolerance)
        ambiguous = (len(possible_indices) > 1 and
                     abs(float(np.partition(distances, 1)[1]) - nearest_distance) < 0.01)
        duplicate = nearest_index in used_atom_indices and nearest_distance <= match_tolerance
        matched = nearest_distance <= match_tolerance and not ambiguous and not duplicate
        values = record['term_values']

        if matched:
            used_atom_indices.add(nearest_index)
            ligand_atom = str(atom_names[nearest_index])
            ligand_residue = '{}:{}{}'.format(
                chain_ids[nearest_index], residue_names[nearest_index],
                residue_numbers[nearest_index])
            mapping_status = 'MATCHED'
        else:
            ligand_atom = 'UNMATCHED'
            ligand_residue = 'UNMATCHED'
            if ambiguous:
                mapping_status = 'AMBIGUOUS'
            elif duplicate:
                mapping_status = 'DUPLICATE_MATCH'
            elif np.any(np.abs(values) > 1e-9):
                mapping_status = 'UNMATCHED_NONZERO_ENERGY'
            elif 'hydrogen' in record['atom_type'].lower():
                mapping_status = 'UNMATCHED_ZERO_ENERGY_HYDROGEN'
            else:
                mapping_status = 'UNMATCHED_ZERO_ENERGY'

        legacy_key = '{} {}'.format(record['atom_type'], ligand_atom)
        if legacy_key in legacy_output:
            legacy_key = '{} #{}'.format(legacy_key, record_index + 1)
        legacy_output[legacy_key] = [float(np.sum(file_records[record_index]['term_values']))
                                     for file_records in records_by_file]

        atom_row = {'smina_row': record_index + 1,
                    'smina_type': record['atom_type'],
                    'ligand_residue': ligand_residue,
                    'ligand_atom': ligand_atom,
                    'mapping_status': mapping_status,
                    'mapping_distance_A': nearest_distance}
                    
        atom_row.update({term: float(value) for term, value in zip(term_names, values)})
        atom_row['sum_atom_terms'] = float(np.sum(values))
        output_rows.append(atom_row)

    if not detailed:
        return legacy_output

    totals = {term: float(sum(row[term] for row in output_rows)) for term in term_names}
    mapping_status_counts = {}
    for row in output_rows:
        status = row['mapping_status']
        mapping_status_counts[status] = mapping_status_counts.get(status, 0) + 1

    absolute_unmatched_energy = float(sum(
        sum(abs(row[term]) for term in term_names)
        for row in output_rows if row['mapping_status'] != 'MATCHED'))
        
    requires_review = (absolute_unmatched_energy > 1e-9 or
                       any(row['mapping_status'] in ('AMBIGUOUS', 'DUPLICATE_MATCH')
                           for row in output_rows))
    
    quality_control = {'mapped_atoms': mapping_status_counts.get('MATCHED', 0),
                       'smina_atom_rows': len(output_rows),
                       'mapping_status_counts': mapping_status_counts,
                       'absolute_unmatched_energy': absolute_unmatched_energy,
                       'status': 'REVIEW' if requires_review else 'OK'}
        
    global_terms = [term for term in term_names
                    if term.startswith('num_') or term == 'ligand_length']
    
    return {'terms': term_names,
            'global_terms': global_terms,
            'atoms': output_rows,
            'totals': totals,
            'qc': quality_control}


def showSminaPerAtomInteractions(data, top=10, terms=None,
                                 sort_by='favorable', ax=None):
    """Show atom-resolved SMINA scoring contribution.

    The default plot compares the repulsion and hydrogen-bond scoring
    terms for ligand atoms ranked by the magnitude of their favorable
    contributions. A favorable contribution is a negative value of an
    individual scoring term. The ranking uses all available atom-resolved
    terms, excluding terms classified as global.

    :arg data: Results returned by :func:`.calcSminaPerAtomInteractions`
        with ``detailed=True``. Only uniquely matched ligand atoms are
        plotted. Unmatched zero-energy hydrogens are excluded.
    :type data: dict

    :arg top: Maximum number of ligand atoms shown. If None, show all
        eligible atoms. Default is 10.
    :type top: int, None

    :arg terms: Scoring terms to display as separate bar series. If None,
        use the repulsion and non-directional hydrogen-bond terms when
        available. Missing default terms are replaced with those having
        the largest absolute contributions. Scoring term names are read
        from ``data['terms']``, allowing non-Vina scoring functions.
        Default is None.
    :type terms: list, tuple, None

    :arg sort_by: Ranking criterion. ``'favorable'`` (default) sorts by
        the sum of negative contributions from all atom-resolved terms;
        ``'net'`` sorts by the sum of all signed contributions;
        ``'hbond'`` sorts by the most favorable hydrogen-bond contribution.
        Alternatively, provide an exact term name from ``data['terms']``
        to sort by its signed value (ascending).
    :type sort_by: str

    :arg ax: Matplotlib Axes on which to draw the plot. If None, create
        a new figure and axes. Default is None.
    :type ax: :class:`matplotlib.axes.Axes`, None

    :returns: Matplotlib Axes containing the bar plot. The associated
        figure is accessible as ``ax.figure``.
    :rtype: :class:`matplotlib.axes.Axes`

    A negative contribution to an attractive scoring term is favorable,
    whereas a positive repulsion contribution is unfavorable. The
    ranking is a diagnostic ranking of SMINA atom-term contribution. In
    particular, the favorable-only ranking does not subtract steric
    penalties.

    The function does not rerun SMINA and does not modify the input
    data or save any files. Results with unmatched nonzero-energy atom
    records cannot be plotted reliably and raise ValueError. """

    import numpy as np
    import matplotlib.pyplot as plt
    from builtins import any, sum, min, max

    if not isinstance(data, dict):
        raise TypeError('data must be a detailed SMINA per-atom result dictionary')

    term_names = data.get('terms')
    atom_rows = data.get('atoms')
    if not isinstance(term_names, (list, tuple)) or not term_names:
        raise ValueError('data must contain a non-empty list of scoring terms')
    if not isinstance(atom_rows, list):
        raise ValueError('data must contain an atoms list')
    if len(term_names) != len(set(term_names)):
        raise ValueError('data contains duplicate scoring-term names')

    if top is not None and (isinstance(top, bool) or
                            not isinstance(top, (int, np.integer)) or top <= 0):
        raise ValueError('top must be a positive integer or None')

    try:
        unmatched_energy = float(data.get('qc', {}).get('absolute_unmatched_energy', 0.0))
    except (ValueError, TypeError, AttributeError) as error:
        raise ValueError('Invalid unmatched-energy quality-control value') from error
    if not np.isfinite(unmatched_energy) or unmatched_energy > 1e-9:
        raise ValueError('Cannot plot: unmapped ligand atoms have nonzero '
                         'scoring contributions (check data["qc"])')

    global_terms = set(data.get('global_terms', []))
    atom_terms = [name for name in term_names if name not in global_terms]
    if not atom_terms:
        raise ValueError('No atom-resolved scoring terms are available')

    matched_rows = [row for row in atom_rows
                    if row.get('mapping_status') == 'MATCHED']
    if not matched_rows:
        raise ValueError('No uniquely matched ligand atoms are available')

    for row in matched_rows:
        for name in atom_terms:
            try:
                value = float(row[name])
            except (KeyError, ValueError, TypeError) as error:
                raise ValueError('Missing/invalid value for scoring term {!r}'.format(
                    name)) from error
            if not np.isfinite(value):
                raise ValueError('Non-finite scoring contribution for {!r}'.format(name))

    if terms is None:
        repulsion = next((name for name in atom_terms
                          if name.startswith('repulsion(')), None)
        hydrogen_bond = next((name for name in atom_terms
                              if 'h_bond' in name.lower() or
                              'hbond' in name.lower()), None)
        selected_terms = [name for name in (repulsion, hydrogen_bond)
                          if name is not None]

        ranked_terms = sorted(atom_terms, key=lambda name: sum(
            abs(float(row[name])) for row in matched_rows), reverse=True)
        for name in ranked_terms:
            if len(selected_terms) >= 2:
                break
            if name not in selected_terms:
                selected_terms.append(name)
    else:
        if not isinstance(terms, (list, tuple)) or not terms:
            raise TypeError('terms must be a non-empty list/tuple or None')
        selected_terms = list(terms)
        if len(selected_terms) != len(set(selected_terms)):
            raise ValueError('terms contains duplicate term names')
        invalid = [name for name in selected_terms if name not in atom_terms]
        if invalid:
            raise ValueError('Unknown atom-resolved scoring terms: {}'.format(invalid))

    rows = [row for row in matched_rows if any(
        abs(float(row[name])) > 1e-12 for name in selected_terms)]
    if not rows:
        raise ValueError('Selected scoring terms have no nonzero matched contributions')

    if sort_by == 'favorable':
        key = lambda row: sum(min(0.0, float(row[name])) for name in atom_terms)
    elif sort_by == 'net':
        key = lambda row: sum(float(row[name]) for name in atom_terms)
    elif sort_by == 'hbond':
        hbond_names = [name for name in atom_terms
                       if 'h_bond' in name.lower() or 'hbond' in name.lower()]
        if not hbond_names:
            raise ValueError("sort_by='hbond' requires a hydrogen-bond scoring term")
        key = lambda row: sum(float(row[name]) for name in hbond_names)
    elif sort_by in atom_terms:
        key = lambda row: float(row[sort_by])
    else:
        raise ValueError('sort_by must be favorable, net, hbond, or '
                         'an exact scoring-term name')

    rows = sorted(rows, key=key)
    if top is not None:
        rows = rows[:top]

    labels = [str(row['ligand_atom']) for row in rows]

    if len(labels) != len(set(labels)):
        labels = ['{}:{}'.format(row['ligand_residue'], row['ligand_atom'])
                  for row in rows]
    if len(labels) != len(set(labels)):
        labels = ['{} (#{})'.format(label, row['smina_row'])
                  for label, row in zip(labels, rows)]

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, max(4, 0.5 * len(rows) + 1.6)))
        own_figure = True
    else:
        fig = ax.figure
        own_figure = False

    positions = np.arange(len(rows))
    bar_height = min(0.38, 0.8 / len(selected_terms))
    for index, name in enumerate(selected_terms):
        offset = (index - (len(selected_terms) - 1) / 2) * bar_height
        if name.startswith('repulsion('):
            display_name = 'Repulsion'
        elif 'h_bond' in name.lower() or 'hbond' in name.lower():
            display_name = 'H-bond'
        else:
            display_name = name
        values = [float(row[name]) for row in rows]
        ax.barh(positions + offset, values, height=bar_height * 0.92,
                label=display_name)

    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.axvline(0.0, color='0.5', linewidth=0.8)
    ax.set_xlabel('SMINA scoring-term contribution [kcal/mol]')
    ax.set_title('SMINA: ligand-atom contributions (TOP {})'.format(len(rows)))
    ax.legend()
    if own_figure:
        fig.tight_layout()

    if SETTINGS['auto_show']:
        showFigure()
    return ax


def calcSminaResidueContributions(atoms, ligand_selection,
                                  protein_selection='protein', cutoff=8.0,
                                  scoring_function='vina', max_residues=None,
                                  output_tsv=None, smina_executable='smina'):
    """Measure how individual protein residues influence SMINA scoring.

    Calculate the SMINA score of a fixed protein-ligand pose and repeat
    ``--score_only`` after omitting one protein residue at a time from the
    receptor. The ligand and all remaining protein atom coordinates are
    retained. No docking or coordinate minimization is performed.

    Negative values indicate that the presence of the residue
    favors the interaction according to the SMINA scoring function,
    whereas positive values indicate an unfavorable contribution
    in this fixed-pose comparison.

    :arg atoms: Atomic object containing the protein-ligand
        complex. The structure must contain one coordinate set.
    :type atoms: :class:`.Atomic`

    :arg ligand_selection: ProDy selection string identifying
        a single ligand residue whose interactions with the
        receptor will be analyzed.
    :type ligand_selection: str

    :arg protein_selection: ProDy selection string defining
        the receptor used for SMINA scoring.
        Protein residues included in this selection are
        considered for individual removal.
        Default is 'protein'.
    :type protein_selection: str

    :arg cutoff: Maximum distance between a protein residue
        and the ligand for including the residue in the analysis.

        The distance is defined as the minimum distance between
        any heavy atom of the protein residue and any heavy atom
        of the ligand.
        Only protein residues with a minimum distance smaller
        than or equal to cutoff are considered.
        Default is 8.0 Angstrom.
    :type cutoff: float

    :arg scoring_function: Name of the SMINA scoring function
        used for the calculations.
        The available options depend on the installed SMINA version.
        Default is 'vina'.
    :type scoring_function: str

    :arg max_residues: Maximum number of protein residues
        included in the calculations.

        Residues are selected in ascending order of their
        minimum heavy-atom distance from the ligand.
        If None, all protein residues satisfying the cutoff
        criterion are included.
        Default is None.
    :type max_residues: int, None

    :arg output_tsv: Optional output filename for saving the
        calculated residue contributions as a tab-separated
        text file.

        If None, results are returned without saving a file.
        An existing output file is not overwritten.
        Default is None.
    :type output_tsv: str, None

    :arg smina_executable: Name or path of the SMINA executable
        used for scoring calculations.
        Default is 'smina'.
    :type smina_executable: str

    :returns: List of dictionaries containing the calculated
        influence of individual protein residues on SMINA scoring.

        Each dictionary contains:
            - residue: protein residue identifier.
            - segment: segment identifier.
            - chain: chain identifier.
            - resname: amino acid residue name.
            - resnum: residue number.
            - icode: insertion code.
            - min_distance_A: minimum heavy-atom distance
              between the residue and ligand in Angstroms.
            - score_full_kcal_mol: SMINA score for the
              complete protein-ligand complex.
            - score_without_kcal_mol: SMINA score obtained
              after removing the corresponding protein residue.
            - delta_score_kcal_mol: difference between the
              complete and residue-deleted receptor scores.
            - delta_terms: dictionary containing changes
              in individual weighted scoring function terms.

        The results are sorted by increasing delta_score_kcal_mol.
        Scoring term names and their number are obtained
        automatically from the SMINA output.
    :rtype: list

    The calculated residue contributions represent changes
    in the empirical SMINA scoring function for a fixed ligand
    pose. They should not be interpreted as binding free
    energies or rigorous residue-wise energy decomposition.

    Removing a protein residue can affect the preparation
    or atom typing of the remaining receptor atoms.
    Consequently, the calculated contributions may not
    be strictly additive.

    Individual scoring-term differences are reported separately.
    Global scoring terms, such as torsional corrections,
    are not assigned directly to individual protein residues. """
    
    import csv
    import subprocess
    import tempfile
    from pathlib import Path

    import numpy as np

    if not hasattr(atoms, 'select') or not hasattr(atoms, 'numCoordsets'):
        raise TypeError('atoms must be a ProDy Atomic object')
    if atoms.numCoordsets() != 1:
        raise ValueError('Exactly one coordinate set is required')
    if not isinstance(ligand_selection, str) or not ligand_selection.strip():
        raise ValueError('ligand_selection must identify one ligand residue')
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('cutoff must be a positive finite distance in Angstrom')
    if max_residues is not None and (not isinstance(max_residues, int) or max_residues < 1):
        raise ValueError('max_residues must be a positive integer or None')

    ligand = atoms.select(ligand_selection)
    receptor = atoms.select(protein_selection)
    if ligand is None or receptor is None:
        raise ValueError('ligand_selection or protein_selection returned no atoms')
    if np.intersect1d(ligand.getIndices(), receptor.getIndices()).size:
        raise ValueError('Ligand and receptor atom selections must not overlap')

    ligand_insertion_codes = ligand.getIcodes()
    if ligand_insertion_codes is None:
        ligand_insertion_codes = [''] * ligand.numAtoms()
    ligand_keys = set(zip(ligand.getChids(), ligand.getResnums(),
                          ligand.getResnames(), ligand_insertion_codes))
    if len(ligand_keys) != 1:
        raise ValueError('ligand_selection must identify exactly one ligand residue')

    protein = receptor.select('protein')
    ligand_heavy = ligand.select('noh')
    if protein is None or ligand_heavy is None:
        raise ValueError('The receptor has no protein or the ligand has no heavy atoms')
    protein_heavy = protein.select('noh')
    if protein_heavy is None:
        raise ValueError('The receptor contains no protein heavy atoms')

    def optional_strings(values, count):
        """Return residue annotations, substituting empty strings when missing."""
        if values is None:
            return [''] * count
        return [str(value).strip() for value in values]

    atom_count = protein.numAtoms()
    segments = optional_strings(protein.getSegnames(), atom_count)
    chains = optional_strings(protein.getChids(), atom_count)
    residue_numbers = [int(number) for number in protein.getResnums()]
    insertion_codes = optional_strings(protein.getIcodes(), atom_count)
    residue_names = optional_strings(protein.getResnames(), atom_count)
    atom_indices = [int(index) for index in protein.getIndices()]

    residue_groups = {}
    for position in range(atom_count):
        key = (segments[position], chains[position], residue_numbers[position],
               insertion_codes[position], residue_names[position])
        residue_groups.setdefault(key, []).append(position)

    heavy_atom_indices = set(map(int, protein_heavy.getIndices()))
    protein_coordinates = protein.getCoords()
    ligand_coordinates = ligand_heavy.getCoords()
    candidates = []
    for residue_key, residue_positions in residue_groups.items():
        heavy_positions = [position for position in residue_positions
                           if atom_indices[position] in heavy_atom_indices]
        if not heavy_positions:
            continue
        displacement = (protein_coordinates[heavy_positions, None, :] -
                        ligand_coordinates[None, :, :])
        distance = float(np.sqrt(np.min(np.einsum(
            'ijk,ijk->ij', displacement, displacement))))
        if distance <= cutoff:
            candidates.append((distance, residue_key,
                               [atom_indices[position] for position in residue_positions]))

    candidates.sort(key=lambda item: (item[0], item[1]))
    if max_residues is not None:
        candidates = candidates[:max_residues]
    if not candidates:
        raise ValueError('No protein residues were found within cutoff')

    global_term_names = {'num_tors_div', 'num_heavy_atoms_div', 'num_heavy_atoms',
                        'num_tors_add', 'num_tors_sqr', 'num_tors_sqrt',
                        'num_hydrophobic_atoms', 'ligand_length'}

    def score_structure(receptor_filename, ligand_filename):
        """Obtain the SMINA score and weighted terms for one fixed pose."""
        command = [str(smina_executable), '-r', str(receptor_filename),
                   '-l', str(ligand_filename), '--score_only',
                   '--scoring', str(scoring_function)]
        try:
            process = subprocess.run(command, capture_output=True,
                                     text=True, check=True)
        except FileNotFoundError as error:
            raise RuntimeError('SMINA executable not found: {}'.format(
                smina_executable)) from error
        except subprocess.CalledProcessError as error:
            raise RuntimeError('SMINA scoring failed (return code {}):\n{}\n{}'.format(
                error.returncode, error.stderr[-2000:], error.stdout[-2000:])) from error

        parsed = SMINA_extract_data(process.stdout)
        if 'Affinity' not in parsed or 'Term values, before weighting' not in parsed:
            raise ValueError('Incomplete SMINA score_only output')
        excluded_names = {'Affinity', 'Intramolecular energy',
                          'Term values, before weighting'}
        weights = {name: float(value) for name, value in parsed.items()
                   if name not in excluded_names}
        raw_values = parsed['Term values, before weighting']
        if len(weights) != len(raw_values):
            raise ValueError('SMINA scoring-term count differs from its weight count')
        weighted_terms = {name: weight * float(value)
                          for (name, weight), value in zip(weights.items(), raw_values)}
        return float(parsed['Affinity']), weights, weighted_terms

    rows = []
    with tempfile.TemporaryDirectory(prefix='prody_smina_residues_') as directory:
        directory = Path(directory)
        receptor_filename = directory / 'receptor.pdb'
        ligand_filename = directory / 'ligand.pdb'
        writePDB(str(ligand_filename), ligand)
        writePDB(str(receptor_filename), receptor)
        full_score, full_weights, full_terms = score_structure(
            receptor_filename, ligand_filename)

        for distance, residue_key, indices_to_remove in candidates:
            remaining = receptor.select('not (index {})'.format(
                ' '.join(map(str, indices_to_remove))))
            if remaining is None or remaining.numAtoms() == 0:
                raise ValueError('Omitting this residue would empty the receptor')

            writePDB(str(receptor_filename), remaining)
            omitted_score, omitted_weights, omitted_terms = score_structure(
                receptor_filename, ligand_filename)
                
            full_names = set(full_weights)
            omitted_names = set(omitted_weights)

            if full_names != omitted_names:
                missing = sorted(full_names - omitted_names)
                added = sorted(omitted_names - full_names)

                raise ValueError('SMINA scoring terms changed after omitting residue {}. '
                    'Missing terms: {}. Additional terms: {}.'.format(
                                                residue_key, missing, added))

            changed_weights = {
                name: (full_weights[name], omitted_weights[name])
                for name in full_weights
                if abs(full_weights[name] - omitted_weights[name]) > 1e-10}

            if changed_weights:
                raise ValueError('SMINA scoring weights changed after omitting residue {}. '
                    'Weights (full, omitted): {}'.format(residue_key, changed_weights))    
                            
            delta_terms = {name: round(full_terms[name] - omitted_terms[name], 6)
                            for name in full_terms if name not in global_term_names}
            segment, chain, residue_number, insertion_code, residue_name = residue_key
            rows.append({'residue': '{}:{}{}{}'.format(chain or '_', residue_name,
                                              residue_number, insertion_code),
                         'segment': segment,
                         'chain': chain,
                         'resname': residue_name,
                         'resnum': residue_number,
                         'icode': insertion_code,
                         'min_distance_A': round(distance, 3),
                         'score_full_kcal_mol': full_score,
                         'score_without_kcal_mol': omitted_score,
                         'delta_score_kcal_mol': round(full_score - omitted_score, 6),
                         'delta_terms': delta_terms})

    rows.sort(key=lambda row: row['delta_score_kcal_mol'])
    if output_tsv is not None:
        path = Path(output_tsv)
        term_names = list(rows[0]['delta_terms'])
        fieldnames = ([name for name in rows[0] if name != 'delta_terms'] +
                      ['delta_term: {}'.format(name) for name in term_names])
        with path.open('x', newline='', encoding='utf-8') as output_file:
            writer = csv.DictWriter(output_file, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                values = {name: value for name, value in row.items()
                          if name != 'delta_terms'}
                values.update({'delta_term: {}'.format(name): row['delta_terms'][name]
                               for name in term_names})
                writer.writerow(values)
    return rows


def listSminaResidueContributions(data, top=10, sort_by='favorable'):
    """List and rank previously calculated SMINA residue contributions.

    This function sorts and prints the results obtained with
    :func:`calcSminaResidueContributions`. All residue records remain 
    available in the returned list, even when only the first ``top`` 
    are printed.

    :arg data: Results returned by :func:`calcSminaResidueContributions`,
        or the path to a TSV file created with its ``output_tsv`` option.
        Each residue record must contain ``residue``, ``min_distance_A``,
        ``delta_score_kcal_mol`` and ``delta_terms``.
    :type data: list, str, pathlib.Path

    :arg top: Maximum number of ranked residues to print. If **None**,
        print all residues. Default is 10. This parameter does not limit
        the number of returned records.
    :type top: int, None

    :arg sort_by: Criterion used to rank the residues. Available options are:

        - ``'favorable'``: most negative total score change (default).
        - ``'unfavorable'``: most positive total score change.
        - ``'absolute'``: largest absolute total score change.
        - ``'hbond'``: most negative combined hydrogen-bond contributions.
        - ``'repulsion'``: largest positive combined repulsion contributions.
        - ``'distance'``: shortest minimum heavy-atom distance to the ligand.

        An exact term name in ``delta_terms`` can also be supplied;
        contributions are then sorted from most negative to most positive.
        Available terms are detected from the input data, without assuming
        a particular scoring function or number of terms.
    :type sort_by: str

    :returns: New list of all residue dictionaries sorted by ``sort_by``.
        The input list and its dictionaries are not modified.
    :rtype: list

    The ranking describes changes in the empirical SMINA score upon
    omission of a residue at a fixed ligand pose. It is not a binding
    free-energy decomposition. Negative score changes denote favorable
    contributions to the score. """
    
    import csv
    from pathlib import Path
    from builtins import any, sum

    if isinstance(data, (str, Path)):
        with Path(data).open('r', newline='', encoding='utf-8') as input_file:
            reader = csv.DictReader(input_file, delimiter='\t')
            fields = reader.fieldnames or []
            term_fields = [name for name in fields
                           if name.startswith('delta_term: ')]
            if not term_fields:
                raise ValueError('The TSV file has no delta_term columns')
            loaded = []
            for index, row in enumerate(reader):
                record = {name: value for name, value in row.items()
                          if name not in term_fields}
                try:
                    record['delta_terms'] = {
                        name[len('delta_term: '):]: float(row[name])
                        for name in term_fields
                    }
                    for name in ('min_distance_A', 'score_full_kcal_mol',
                                 'score_without_kcal_mol', 'delta_score_kcal_mol'):
                        if name in record:
                            record[name] = float(record[name])
                    if 'resnum' in record:
                        record['resnum'] = int(record['resnum'])
                except (ValueError, TypeError) as error:
                    raise ValueError('Invalid TSV values in row {}'.format(
                        index + 2)) from error
                loaded.append(record)
            data = loaded

    if not isinstance(data, list) or any(not isinstance(row, dict) for row in data):
        raise TypeError('data must be a list of residue dictionaries or a TSV path')
    if top is not None and (isinstance(top, bool) or
                            not isinstance(top, int) or top < 1):
        raise ValueError('top must be a positive integer or None')
    if not isinstance(sort_by, str) or not sort_by.strip():
        raise ValueError('sort_by must be a nonempty string')
    if not data:
        LOGGER.warn('SMINA residue contributions: no residues to display.')
        return []

    required = {'residue', 'min_distance_A', 'delta_score_kcal_mol', 'delta_terms'}
    term_names = list(data[0]['delta_terms']) if isinstance(
        data[0].get('delta_terms'), dict) else []
    expected_terms = set(term_names)
    if not term_names:
        raise ValueError('No scoring terms were found in delta_terms')
    for index, row in enumerate(data):
        missing = required - set(row)
        if missing:
            raise ValueError('Residue record {} is missing: {}'.format(
                index, sorted(missing)))
        if not isinstance(row['delta_terms'], dict) or set(row['delta_terms']) != expected_terms:
            raise ValueError('Inconsistent scoring terms in residue record {}'.format(index))
        try:
            float(row['min_distance_A'])
            float(row['delta_score_kcal_mol'])
            for name in term_names:
                float(row['delta_terms'][name])
        except (TypeError, ValueError, KeyError) as error:
            raise ValueError('Invalid numeric data in residue record {}'.format(index)) from error

    hbond_terms = [name for name in term_names if 'h_bond' in name.lower()
                   or 'hbond' in name.lower()]
    repulsion_terms = [name for name in term_names if 'repulsion' in name.lower()]

    def term_total(row, names):
        return sum(float(row['delta_terms'][name]) for name in names)

    def secondary(row):
        return str(row['residue'])

    if sort_by == 'favorable':
        key = lambda row: (float(row['delta_score_kcal_mol']), secondary(row))
    elif sort_by == 'unfavorable':
        key = lambda row: (-float(row['delta_score_kcal_mol']), secondary(row))
    elif sort_by == 'absolute':
        key = lambda row: (-abs(float(row['delta_score_kcal_mol'])), secondary(row))
    elif sort_by == 'hbond':
        if not hbond_terms:
            raise ValueError('No hydrogen-bond terms are available for sorting')
        key = lambda row: (term_total(row, hbond_terms), secondary(row))
    elif sort_by == 'repulsion':
        if not repulsion_terms:
            raise ValueError('No repulsion terms are available for sorting')
        key = lambda row: (-term_total(row, repulsion_terms), secondary(row))
    elif sort_by == 'distance':
        key = lambda row: (float(row['min_distance_A']), secondary(row))
    elif sort_by in term_names:
        key = lambda row: (float(row['delta_terms'][sort_by]), secondary(row))
    else:
        raise ValueError(
            'Unknown sort_by={!r}. Available options: favorable, unfavorable, '
            'absolute, hbond, repulsion, distance, or an exact term name: {}'.format(
                sort_by, term_names))

    ranked = sorted(data, key=key)
    selected = ranked if top is None else ranked[:top]

    header = '{:<13s}  {:>10s}  {:>14s}'.format('Residue', 'Dist. [A]', 'Delta score')
    if hbond_terms:
        header += '  {:>12s}'.format('H-bond')
    if repulsion_terms:
        header += '  {:>12s}'.format('Repulsion')
    if sort_by in term_names and sort_by not in hbond_terms + repulsion_terms:
        header += '  {:>12s}'.format('Selected term')

    LOGGER.info('\nSMINA residue contributions: TOP {} / {} | sort_by={}'.format(
        len(selected), len(ranked), sort_by))
    LOGGER.info(header)
    LOGGER.info('-' * len(header))
    
    for row in selected:
        line = '{:<13s}  {:>10.2f}  {:>+14.4f}'.format(
            str(row['residue']), float(row['min_distance_A']),
            float(row['delta_score_kcal_mol']))
        if hbond_terms:
            line += '  {:>+12.4f}'.format(term_total(row, hbond_terms))
        if repulsion_terms:
            line += '  {:>+12.4f}'.format(term_total(row, repulsion_terms))
        if sort_by in term_names and sort_by not in hbond_terms + repulsion_terms:
            line += '  {:>+12.4f}'.format(float(row['delta_terms'][sort_by]))
        LOGGER.info(line)

    return ranked


def showSminaResidueContributions(data, top=10, sort_by='favorable',
                                  terms=None, components=True, figsize=None):
    """Visualize SMINA residue-omission scoring for a ligand pose.

    Plot the signed difference in the final SMINA score for each selected
    protein residue, optionally alongside its weighted scoring-term differences.

    The two panels have DIFFERENT scoring conventions: the final score includes
    normalization/corrections that need not be present in the per-term scores.
    Consequently, stacked term differences should not be interpreted as a
    numerical decomposition of the final-score difference.

    :arg data: Results returned by :func:`calcSminaResidueContributions` or
        :func:`listSminaResidueContributions`, or a path to the TSV file
        created by ``calcSminaResidueContributions(output_tsv=...)``.
    :type data: list, str, pathlib.Path

    :arg top: Number of residues to plot; use **None** for all residues.
        The complete input data are preserved. Default is 10.
    :type top: int, None

    :arg sort_by: Ranking criterion: ``'favorable'`` (most negative final
        score change), ``'unfavorable'`` (most positive change),
        ``'absolute'``, ``'hbond'``, ``'repulsion'``, ``'distance'``, or
        an exact name of a scoring term present in ``delta_terms``.
        Default is ``'favorable'``.
    :type sort_by: str

    :arg terms: Exact names of weighted scoring terms to show in the right
        panel. If **None**, show all terms with a nonzero contribution among
        the selected residues. Term names are detected dynamically and are
        not limited to the Vina scoring function. Default is None.
    :type terms: list, tuple, None

    :arg components: If **True**, draw a second panel with independently
        stacked positive and negative weighted term contributions. If
        **False**, only the signed final-score differences are plotted.
        Default is True.
    :type components: bool

    :arg figsize: Optional figure size ``(width, height)`` in inches.
        If **None**, select a size based on ``top`` and ``components``.
    :type figsize: tuple, None

    :returns: The Matplotlib Axes with final-score differences. The figure
        (and, if requested, the component Axes) can be accessed through
        ``ax.figure`` (and ``ax.figure.axes``).
    :rtype: :class:`matplotlib.axes.Axes`

    A negative difference is favorable for the SMINA score, while a positive
    difference is unfavorable. These differences are defined as the score of
    the complete complex minus the score after omitting one protein residue,
    at a fixed ligand pose. """

    import builtins as _bi
    import csv
    import math
    from pathlib import Path
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    if isinstance(data, (str, Path)):
        source = Path(data)
        with source.open('r', newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            header = reader.fieldnames or []
            term_fields = [name for name in header
                           if name.startswith('delta_term: ')]
            if not term_fields:
                raise ValueError('TSV has no delta_term columns: {}'.format(source))
            missing = {'residue', 'min_distance_A', 'delta_score_kcal_mol'} - set(header)
            if missing:
                raise ValueError('TSV is missing columns: {}'.format(sorted(missing)))
            loaded = []
            for line_number, item in enumerate(reader, start=2):
                try:
                    record = {name: value for name, value in item.items()
                              if name not in term_fields}
                    record['delta_terms'] = {
                        name[len('delta_term: '):]: float(item[name])
                        for name in term_fields}
                    for key in ('min_distance_A', 'delta_score_kcal_mol',
                                'score_full_kcal_mol', 'score_without_kcal_mol'):
                        if key in record and record[key] not in ('', None):
                            record[key] = float(record[key])
                    if 'resnum' in record and record['resnum'] not in ('', None):
                        record['resnum'] = int(record['resnum'])
                except (TypeError, ValueError, KeyError) as exc:
                    raise ValueError('Bad TSV values at line {}'.format(
                        line_number)) from exc
                loaded.append(record)
            data = loaded

    if not isinstance(data, list) or not _bi.all(
            isinstance(row, dict) for row in data):
        raise TypeError('data must be a list of residue dictionaries or a TSV path')
    if not data:
        raise ValueError('No residue contributions were provided')
    if top is not None and (isinstance(top, bool) or
                            not isinstance(top, int) or top < 1):
        raise ValueError('top must be a positive integer or None')
    if not isinstance(sort_by, str) or not sort_by.strip():
        raise ValueError('sort_by must be a non-empty string')
    if not isinstance(components, bool):
        raise TypeError('components must be True or False')
    if figsize is not None:
        try:
            valid_size = (len(figsize) == 2 and
                          _bi.all(math.isfinite(float(x)) and float(x) > 0
                                  for x in figsize))
        except (TypeError, ValueError):
            valid_size = False
        if not valid_size:
            raise ValueError('figsize must be two positive finite numbers')

    required = {'residue', 'min_distance_A',
                'delta_score_kcal_mol', 'delta_terms'}
    if not isinstance(data[0].get('delta_terms'), dict) or not data[0]['delta_terms']:
        raise ValueError('No scoring-term data found in the first residue')
    term_names = list(data[0]['delta_terms'])
    expected_terms = set(term_names)
    for index, row in enumerate(data):
        if required - set(row):
            raise ValueError('Residue {} is missing: {}'.format(
                index, sorted(required - set(row))))
        if not isinstance(row['delta_terms'], dict) or set(row['delta_terms']) != expected_terms:
            raise ValueError('Inconsistent scoring terms in residue {}'.format(index))
        try:
            numbers = ([row['min_distance_A'], row['delta_score_kcal_mol']] +
                       [row['delta_terms'][name] for name in term_names])
            if not _bi.all(math.isfinite(float(x)) for x in numbers):
                raise ValueError('non-finite')
        except (TypeError, ValueError, KeyError) as exc:
            raise ValueError('Invalid numeric data in residue {}'.format(index)) from exc

    hbond_terms = [name for name in term_names
                   if 'h_bond' in name.lower() or 'hbond' in name.lower()]
    repulsion_terms = [name for name in term_names
                       if 'repulsion' in name.lower()]

    def subtotal(row, names):
        return _bi.sum(float(row['delta_terms'][name]) for name in names)

    def residue_name(row):
        return str(row['residue'])

    if sort_by == 'favorable':
        key = lambda row: (float(row['delta_score_kcal_mol']), residue_name(row))
    elif sort_by == 'unfavorable':
        key = lambda row: (-float(row['delta_score_kcal_mol']), residue_name(row))
    elif sort_by == 'absolute':
        key = lambda row: (-abs(float(row['delta_score_kcal_mol'])), residue_name(row))
    elif sort_by == 'hbond':
        if not hbond_terms:
            raise ValueError("sort_by='hbond' requires a hydrogen-bond scoring term")
        key = lambda row: (subtotal(row, hbond_terms), residue_name(row))
    elif sort_by == 'repulsion':
        if not repulsion_terms:
            raise ValueError("sort_by='repulsion' requires a repulsion scoring term")
        key = lambda row: (-subtotal(row, repulsion_terms), residue_name(row))
    elif sort_by == 'distance':
        key = lambda row: (float(row['min_distance_A']), residue_name(row))
    elif sort_by in term_names:
        key = lambda row: (float(row['delta_terms'][sort_by]), residue_name(row))
    else:
        raise ValueError('Unknown sort_by={!r}; use favorable, unfavorable, '
                         'absolute, hbond, repulsion, distance, or one of {}'.format(
                             sort_by, term_names))

    ranked = sorted(data, key=key)
    displayed = ranked if top is None else ranked[:top]
    if terms is None:
        selected_terms = [name for name in term_names
                          if _bi.any(abs(float(row['delta_terms'][name])) > 1e-12
                                     for row in displayed)]
    else:
        if not isinstance(terms, (list, tuple)) or not terms:
            raise TypeError('terms must be a nonempty list/tuple or None')
        selected_terms = list(terms)
        if len(selected_terms) != len(set(selected_terms)):
            raise ValueError('terms contains duplicate scoring terms')
        unknown = [name for name in selected_terms if name not in expected_terms]
        if unknown:
            raise ValueError('Unknown scoring terms: {}'.format(unknown))

    if components and not selected_terms:
        raise ValueError('Selected residues have no nonzero scoring terms')

    n_rows = len(displayed)
    if figsize is None:
        figsize = (13.6 if components else 8.7,
                   _bi.max(4.5, 0.48 * n_rows + 2.0))
    if components:
        fig, (ax, ax_terms) = plt.subplots(
            ncols=2, figsize=figsize, sharey=True,
            gridspec_kw={'width_ratios': [1.0, 1.45], 'wspace': 0.08})
    else:
        fig, ax = plt.subplots(figsize=figsize)
        ax_terms = None

    positions = np.arange(n_rows)
    scores = np.asarray([float(row['delta_score_kcal_mol'])
                         for row in displayed], dtype=float)
    bar_colors = ['tab:blue' if score < 0 else 'tab:orange'
                  if score > 0 else '0.65' for score in scores]
    ax.barh(positions, scores, color=bar_colors, height=0.68, alpha=0.9)
    ax.axvline(0, color='0.25', lw=1.0)
    labels = ['{}  ({:.2f} Å)'.format(row['residue'],
                                   float(row['min_distance_A']))
              for row in displayed]
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Δ final SMINA score (kcal/mol)')
    ax.set_title('Final score difference', fontsize=11)
    ax.grid(axis='x', alpha=0.18)
    ax.set_axisbelow(True)
    left = _bi.min(0, *scores)
    right = _bi.max(0, *scores)
    span = _bi.max(0.15, right - left)
    ax.set_xlim(left - 0.22 * span, right + 0.18 * span)
    for y, value in enumerate(scores):
        ax.text(value + (0.025 * span if value >= 0 else -0.025 * span),
                y, '{:+.3f}'.format(value),
                ha='left' if value >= 0 else 'right',
                va='center', fontsize=8)

    if components:
        colors = plt.get_cmap('tab20')
        negative_offsets = np.zeros(n_rows)
        positive_offsets = np.zeros(n_rows)
        legend_handles = []
        for term_index, term_name in enumerate(selected_terms):
            values = np.asarray([float(row['delta_terms'][term_name])
                                 for row in displayed], dtype=float)
            negative = np.minimum(values, 0.0)
            positive = np.maximum(values, 0.0)
            color = colors(term_index % colors.N)
            if np.any(negative != 0):
                ax_terms.barh(positions, negative, left=negative_offsets,
                              color=color, height=0.68, alpha=0.92)
            if np.any(positive != 0):
                ax_terms.barh(positions, positive, left=positive_offsets,
                              color=color, height=0.68, alpha=0.92)
            negative_offsets += negative
            positive_offsets += positive
            legend_handles.append(Patch(color=color, label=term_name))
        ax_terms.axvline(0, color='0.25', lw=1.0)
        ax_terms.tick_params(axis='y', which='both', left=False, labelleft=False)
        ax_terms.set_xlabel('Δ weighted scoring terms (before global scaling)')
        ax_terms.set_title('Weighted component differences', fontsize=11)
        ax_terms.grid(axis='x', alpha=0.18)
        ax_terms.set_axisbelow(True)
        comp_left = _bi.min(0, *negative_offsets)
        comp_right = _bi.max(0, *positive_offsets)
        comp_span = _bi.max(0.15, comp_right - comp_left)
        ax_terms.set_xlim(comp_left - 0.08 * comp_span,
                          comp_right + 0.08 * comp_span)

        fig.legend(handles=legend_handles, loc='lower center',
            bbox_to_anchor=(0.62, 0.015), frameon=False,
            ncol=_bi.min(5, len(legend_handles)),
            fontsize=7, handlelength=0.8,
            handletextpad=0.25, columnspacing=0.6)

        fig.subplots_adjust(left=0.165, right=0.985, top=0.97, bottom=0.16, wspace=0.16)

    else:
        fig.subplots_adjust(left=0.235, right=0.965, top=0.97, bottom=0.12)

    return ax


def showSminaLigandStructure(atoms, data, ligand_selection,
                             ligand_file=None,
                             protein_selection='protein',
                             residue_data=None, **kwargs):
    """Show SMINA ligand-atom contributions on a 2D chemical structure.

    Draw a chemical diagram of the selected ligand and highlight its most
    important atoms using contributions previously calculated by
    :func:`.calcSminaPerAtomInteractions` with ``detailed=True``. Blue
    highlights indicate favorable (negative) contributions; orange
    highlights indicate unfavorable (positive) contributions. Both may
    be present on the same ligand atom. By default, each ring is split
    proportionally into favorable and unfavorable sectors: colored ring
    area is proportional to magnitude, while overall radius represents
    the sum of absolute weighted scoring terms for that atom. Set
    ``atom_circle_mode='legacy'`` to reproduce equal-sized color sectors.

    Use ``show_atom_labels=False`` to hide the additional PDB atom-name labels
    (e.g., O2A, O3B, PA) while keeping RDKit chemical element symbols
    (O, N, P), the colored highlights and the terminal ranking.
    ``circle_scale`` changes the radius of the colored circles without
    changing their scoring values. Placement of the optional name labels
    accounts for the selected circle radius.

    The optional MOL, SDF or MOL2 input preserves explicitly defined
    bond orders and should be preferred to PDB coordinates when available.
    Atom identities are matched to the input ligand using 3D coordinates,
    not assumed to have the same ordering in different input files.

    The default visualization enables residue contacts and all qualifying
    residue labels; ``residue_data`` is deliberately not guessed from a file
    name. With ``residue_data`` set, ``residue_sort_by='auto'`` uses the most
    favorable Delta S values; without residue scores it ranks by distance
    and displays neutral boxes. ``protein_selection='protein'`` is kept
    structure-independent: specify a chain where appropriate. A large
    number of contacts can produce crowded figures; use a finite
    ``max_contact_residues`` in that case.
    
    RDKit is required and imported only when the function is called.
    Scoring-term signs are interpreted within the SMINA scoring model. 
    Negative and positive values for one atom can originate from different 
    scoring terms, and must not be assigned to specific protein residues.

    Missing ligand atoms with nonzero scoring contributions, duplicate
    atom-coordinate matches and mismatched chemical structures raise
    ValueError. Unmatched zero-scoring hydrogens are permitted because
    the diagram includes only heavy atoms. The input data and structure
    are not modified and no file is written automatically. 

    Delta S is displayed with TWO decimal places in labels and the contact
    terminal table; internal selection and coloring still use full precision.

    :arg atoms: ProDy Atomic object containing the protein-ligand complex
        used to generate the SMINA ``--atom_terms`` output.
    :type atoms: :class:`.Atomic`

    :arg data: Dictionary returned by :func:`.calcSminaPerAtomInteractions`
        with ``detailed=True``. Includes the per-atom scoring terms and
        coordinate-matching quality-control information.
    :type data: dict

    :arg ligand_selection: ProDy selection string identifying exactly one
        ligand residue, e.g. ``'resname UNL'``.
    :type ligand_selection: str

    :arg ligand_file: Optional ligand structure file (``.mol``, ``.sdf``
        or ``.mol2``) with the same 3D coordinates as the ligand selected
        from ``atoms``. Bond orders are read from this file. If None,
        RDKit constructs a molecular graph from the selected ligand PDB
        coordinates, where connectivity and bond orders may be uncertain.
        Default is None.
    :type ligand_file: str, None

    :arg top: Maximum number of ligand heavy atoms highlighted and printed
        in the ranking. None (default) includes all heavy atoms with nonzero
        scoring contributions.
    :type top: int, None

    :arg sort_by: Ranking criterion. ``'favorable'`` (default) ranks by
        the sum of the magnitudes of all negative scoring-term values;
        ``'net'`` by the negative of the sum of signed term values;
        ``'unfavorable'`` by the sum of positive term values;
        ``'hbond'`` by the magnitude of favorable hydrogen-bond terms.
        The ranking is used to select atoms; both favorable and
        unfavorable components remain visible for each selected atom.
    :type sort_by: str

    :arg match_tolerance: Maximum difference in Angstroms between 3D
        coordinates for the same ligand atom in the chemical structure
        and in the ProDy ligand selection. This is not a distance cutoff
        for protein-ligand interactions. Default is 0.15 Angstrom.
    :type match_tolerance: float

    :arg ax: Optional Matplotlib Axes for drawing the ligand structure.
        If None, the function creates a new figure containing just
        the ligand depiction. Default is None.
    :type ax: :class:`matplotlib.axes.Axes`, None

    :arg show_atom_labels: If True, show additional small PDB atom-name
        labels (e.g., O2A, O3B, PA) outside the highlight circles with a
        white halo. If False (default), hide these additional labels.
        Standard RDKit chemical element symbols, highlights, legend and
        terminal ranking remain unchanged.
    :type show_atom_labels: bool

    :arg circle_scale: Scale factor applied to highlight-circle radii.
        Default is 1.25 (25% larger than the original 1.0 setting).
    :type circle_scale: float

    :arg circle_linewidth: Highlight thickness control (1..10). In
        ``'fractional'`` mode, adjusts the thickness of Matplotlib
        annular sectors; in ``'legacy'`` mode, adjusts RDKit outlines.
        Default is 3.
    :type circle_linewidth: int

    :arg show_residue_contacts: If True (default), label nearby protein
        residues around the ligand diagram and connect them to their nearest
        ligand heavy atom. Connections represent 3D geometric proximity,
        NOT atom-pair SMINA scores. Set False for a ligand-only diagram.
    :type show_residue_contacts: bool

    :arg protein_selection: ProDy selection string specifying the receptor
        atoms considered for the contact overlay. Default is 'protein'.
        For a particular complex, use e.g. 'protein and chain A'.
    :type protein_selection: str

    :arg contact_cutoff: Maximum heavy-atom protein-ligand distance (Angstrom)
        for a residue to be included in the contact overlay. Default 4.0.
        This is INDEPENDENT of match_tolerance and the SMINA scoring cutoff.
        This filter is applied BEFORE residue_sort_by ranking.
    :type contact_cutoff: float

    :arg max_contact_residues: Maximum number of selected contact residues.
        None (default) displays ALL qualifying residues for the selected
        ranking after applying ``contact_cutoff``; use a positive integer to
        limit potentially crowded diagrams. A 'favorable' rank excludes
        non-negative Delta S residues regardless of this value.
    :type max_contact_residues: int, None

    :arg residue_data: Optional residue contributions previously obtained
        from calcSminaResidueContributions() or its saved TSV file. When
        supplied, append the WHOLE-RESIDUE final score difference to the
        contact label. This score is not assigned to the drawn atom pair.
        No SMINA calculations are performed here. Default None.
    :type residue_data: list, str, pathlib.Path, None

    :arg residue_fontsize: Font size of protein-residue labels. Default 16.
    :type residue_fontsize: int, float

    :arg residue_sort_by: Criterion to select residues within contact_cutoff.
        'favorable': most negative delta S; 'unfavorable': most positive
        delta S; 'absolute': largest absolute delta S; 'mixed': balanced
        favorable/unfavorable selection; 'distance': shortest heavy-atom
        distances. 'auto' (default): favorable when residue_data is given,
        otherwise distance. Energy rankings require residue_data.
    :type residue_sort_by: str

    :arg show_contact_distances: If True, include the nearest 3D heavy-atom
        contact distance in residue boxes. If False (default), omit distances
        from boxes, but keep them in the textual terminal contact table.
    :type show_contact_distances: bool

    :arg show_chain: If True, prefix residue labels with chain ID.
        If False (default), omit it unless two displayed residues would otherwise
        have identical labels, in which case chains are added to disambiguate.
    :type show_chain: bool

    :arg residue_code: 'one' (default) or 'three' for amino-acid names.
        Nonstandard/unknown residues retain their original names.
    :type residue_code: str

    :arg color_residue_boxes_by_delta: If True (default), color residue
        boxes blue (negative delta S) or orange-red (positive delta S) with
        saturation scaled by absolute magnitude. Requires residue_data to
        provide delta S; boxes without matching values stay white.
    :type color_residue_boxes_by_delta: bool

    :arg atom_circle_mode: 'fractional' (default) divides each colored
        ring by its atom-wise favorable/(favorable+unfavorable) and
        unfavorable/(favorable+unfavorable) shares, so arc length and
        colored ring AREA express the relative proportions. Overall ring
        radius represents total absolute term contribution. For clear
        proportions, fractional mode uses solid blue and orange hues.
        'legacy' restores equal-sized RDKit color sectors from v7.
    :type atom_circle_mode: str

    :arg atom_symbol_glow: If True (default), redraw original RDKit element
        glyphs (including H prefixes, subscripts and charges) ON TOP of the
        colored scoring rings, with a white halo. The chemical symbols come
        from RDKit's own SVG drawing, preserving their original geometry,
        rather than being approximated by Matplotlib text.
    :type atom_symbol_glow: bool

    :arg atom_symbol_glow_width: Width of the white outline around chemical
        element glyphs in image pixels. Default is 2.0.
    :type atom_symbol_glow_width: float

    :arg residue_text_glow: If True (default), draw a white outline around
        residue LABEL TEXT inside the blue/red boxes without changing the
        colored box background or gray contact line.
    :type residue_text_glow: bool

    :arg residue_text_glow_width: Width of the white text outline in points,
        default 1.5. Controls visual clarity independently from box color.
    :type residue_text_glow_width: float

    :arg ligand_render: 'vector' (default) rebuilds RDKit bond strokes and atom
        symbols as Matplotlib vector paths for resolution-independent SVG/PDF
        output. 'raster' uses the previous RDKit Cairo image background.
    :type ligand_render: str

    :arg residue_box_contrast: Positive multiplier of residue box color contrast.
        1.0 (default) uses the enhanced saturation range; smaller values
        create lighter boxes, larger values more saturated boxes (capped).
        Display only: does not alter the underlying Delta S or ordering.
    :type residue_box_contrast: float

    :returns: Matplotlib Axes showing the annotated chemical structure.
        The parent figure is accessible as ``ax.figure`` and may be
        saved with ``ax.figure.savefig(...)``.
    :rtype: :class:`matplotlib.axes.Axes` """
    
    import io
    import warnings
    import csv
    import re
    import xml.etree.ElementTree as ET
    from pathlib import Path
    from builtins import any, min, max, sum

    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Wedge, PathPatch
    from matplotlib.path import Path as MPLPath
    import numpy as np
    from prody import LOGGER

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
    
    except ImportError as err:
        LOGGER.warn(
            'RDKit is required for showSminaLigandStructure(). '
            'Please install RDKit using: '
            'conda install -c conda-forge rdkit')
        raise ImportError('RDKit is required for showSminaLigandStructure().') from err

    top = kwargs.pop('top', None)
    sort_by = kwargs.pop('sort_by', 'favorable')
    match_tolerance = kwargs.pop('match_tolerance', 0.15)
    ax = kwargs.pop('ax', None)
    ligand_render = kwargs.pop('ligand_render', 'vector')
    show_atom_labels = kwargs.pop('show_atom_labels', False)
    circle_scale = kwargs.pop('circle_scale', 1.25)
    circle_linewidth = kwargs.pop('circle_linewidth', 3)
    atom_circle_mode = kwargs.pop('atom_circle_mode', 'fractional')
    atom_symbol_glow = kwargs.pop('atom_symbol_glow', True)
    atom_symbol_glow_width = kwargs.pop('atom_symbol_glow_width', 2.0)
    show_residue_contacts = kwargs.pop('show_residue_contacts', True)
    contact_cutoff = kwargs.pop('contact_cutoff', 4.0)
    max_contact_residues = kwargs.pop('max_contact_residues', None)
    residue_sort_by = kwargs.pop('residue_sort_by', 'auto')
    residue_fontsize = kwargs.pop('residue_fontsize', 16)
    show_contact_distances = kwargs.pop('show_contact_distances', False)
    show_chain = kwargs.pop('show_chain', False)
    residue_code = kwargs.pop('residue_code', 'one')
    color_residue_boxes_by_delta = kwargs.pop('color_residue_boxes_by_delta', True)
    residue_box_contrast = kwargs.pop('residue_box_contrast', 1.0)
    residue_text_glow = kwargs.pop('residue_text_glow', True)
    residue_text_glow_width = kwargs.pop('residue_text_glow_width', 1.5)

    if kwargs:
        raise TypeError(
            "showSminaLigandStructure() got unexpected keyword "
            "argument(s): {}".format(', '.join(sorted(kwargs))))

    if not isinstance(data, dict) or not data.get('terms') or not isinstance(data.get('atoms'), list):
        raise ValueError('data must be a detailed calcSminaPerAtomInteractions() dictionary')
    if not isinstance(ligand_selection, str) or not ligand_selection:
        raise ValueError('ligand_selection must select exactly one ligand residue')
    if not isinstance(show_atom_labels, bool):
        raise TypeError('show_atom_labels must be True or False')
    if not isinstance(show_residue_contacts, bool):
        raise TypeError('show_residue_contacts must be True or False')
    if not isinstance(show_contact_distances, bool):
        raise TypeError('show_contact_distances must be True or False')
    if not isinstance(show_chain, bool):
        raise TypeError('show_chain must be True or False')
    if not isinstance(color_residue_boxes_by_delta, bool):
        raise TypeError('color_residue_boxes_by_delta must be True or False')
    if not isinstance(atom_symbol_glow, bool):
        raise TypeError('atom_symbol_glow must be True or False')
    if not isinstance(residue_text_glow, bool):
        raise TypeError('residue_text_glow must be True or False')
    for parameter_name, value in (('atom_symbol_glow_width', atom_symbol_glow_width),
                                  ('residue_text_glow_width', residue_text_glow_width)):
        if (isinstance(value, bool) or
                not isinstance(value, (int, float, np.integer, np.floating)) or
                not np.isfinite(value) or value < 0.0 or value > 15.0):
            raise ValueError('{} must be finite and in [0, 15]'.format(parameter_name))
    if residue_code not in ('three', 'one'):
        raise ValueError("residue_code must be 'three' or 'one'")
    if ligand_render not in ('vector', 'raster'):
        raise ValueError("ligand_render must be 'vector' or 'raster'")
    if (isinstance(residue_box_contrast, bool) or
            not isinstance(residue_box_contrast, (int, float, np.integer, np.floating)) or
            not np.isfinite(residue_box_contrast) or
            not 0.0 < residue_box_contrast <= 2.0):
        raise ValueError('residue_box_contrast must be finite and in (0, 2]')
    if atom_circle_mode not in ('fractional', 'legacy'):
        raise ValueError("atom_circle_mode must be 'fractional' or 'legacy'")
    if atom_circle_mode == 'legacy' and ligand_render == 'vector':
        warnings.warn("Legacy highlight mode uses raster ligand rendering; "
                      "switch to atom_circle_mode='fractional' for an "
                      "all-vector ligand figure.", UserWarning, stacklevel=2)
        ligand_render = 'raster'
    accepted_residue_ranks = ('auto', 'favorable', 'unfavorable', 'absolute',
                              'distance', 'mixed')
    if residue_sort_by not in accepted_residue_ranks:
        raise ValueError('residue_sort_by must be one of: {}'.format(
            ', '.join(accepted_residue_ranks)))
    if not isinstance(protein_selection, str) or not protein_selection.strip():
        raise ValueError('protein_selection must be a non-empty selection string')
    if (isinstance(contact_cutoff, bool) or
            not isinstance(contact_cutoff, (int, float, np.integer, np.floating)) or
            not np.isfinite(contact_cutoff) or contact_cutoff <= 0):
        raise ValueError('contact_cutoff must be a positive finite distance in Angstrom')
    if (max_contact_residues is not None and
            (isinstance(max_contact_residues, bool) or
             not isinstance(max_contact_residues, (int, np.integer)) or
             max_contact_residues < 1)):
        raise ValueError('max_contact_residues must be a positive integer or None')
    if (isinstance(residue_fontsize, bool) or
            not isinstance(residue_fontsize, (int, float, np.integer, np.floating)) or
            not np.isfinite(residue_fontsize) or residue_fontsize <= 0):
        raise ValueError('residue_fontsize must be positive')

    if (isinstance(circle_scale, bool) or
            not isinstance(circle_scale, (int, float, np.integer, np.floating)) or
            not np.isfinite(circle_scale) or not 0.1 <= circle_scale <= 3.0):
        raise ValueError('circle_scale must be a finite number between 0.1 and 3.0')
    if (isinstance(circle_linewidth, bool) or
            not isinstance(circle_linewidth, (int, np.integer)) or
            not 1 <= circle_linewidth <= 10):
        raise ValueError('circle_linewidth must be an integer between 1 and 10')
    if top is not None and (isinstance(top, bool) or
                            not isinstance(top, (int, np.integer)) or top <= 0):
        raise ValueError('top must be a positive integer or None')
    if not np.isfinite(match_tolerance) or match_tolerance <= 0:
        raise ValueError('match_tolerance must be a positive finite number in Angstroms')
    unmatched_energy = float(data.get('qc', {}).get('absolute_unmatched_energy', 0))
    if not np.isfinite(unmatched_energy) or unmatched_energy > 1e-9:
        raise ValueError('Unmatched SMINA atom records have nonzero contributions; review mapping')

    ligand = atoms.select(ligand_selection)
    if ligand is None or ligand.numAtoms() == 0:
        raise ValueError('ligand_selection does not select any ligand atoms')
    coords = np.asarray(ligand.getCoords(), dtype=float)
    if coords.shape != (ligand.numAtoms(), 3) or not np.isfinite(coords).all():
        raise ValueError('Invalid ligand coordinates')
    names = [str(n).strip() for n in ligand.getNames()]
    residue_keys = {(str(c), int(n), str(r)) for c, n, r in zip(
        ligand.getChids(), ligand.getResnums(), ligand.getResnames())}
    if len(residue_keys) != 1:
        raise ValueError('ligand_selection must identify exactly one residue')
    if len(names) != len(set(names)):
        raise ValueError('Ligand contains duplicate atom names; select a uniquely named residue')

    ligand_elements = ligand.getElements()
    if ligand_elements is None:
        raise ValueError('Ligand element symbols are missing from the PDB; supply an element-annotated PDB')
    ligand_elements = [str(x).strip().upper() for x in ligand_elements]
    if any(not element for element in ligand_elements):
        raise ValueError('Missing element symbol in ligand selection')

    inferred_bonds = ligand_file is None
    if ligand_file is not None:
        source_path = Path(ligand_file)
        if not source_path.is_file():
            raise FileNotFoundError('Ligand file not found: {}'.format(source_path))
        extension = source_path.suffix.lower()
        if extension == '.mol':
            mol = Chem.MolFromMolFile(str(source_path), removeHs=True, sanitize=True)
        elif extension == '.sdf':
            supplier = Chem.SDMolSupplier(str(source_path), removeHs=True, sanitize=True)
            if len(supplier) != 1:
                raise ValueError('ligand_file SDF must contain exactly one molecule')
            mol = supplier[0]
        elif extension == '.mol2':
            mol = Chem.MolFromMol2File(str(source_path), removeHs=True, sanitize=True)
        else:
            raise ValueError('ligand_file must have .mol, .sdf or .mol2 extension')
    else:
        from prody.proteins.pdbfile import writePDBStream
        pdb_buffer = io.StringIO()
        writePDBStream(pdb_buffer, ligand)
        mol = Chem.MolFromPDBBlock(pdb_buffer.getvalue(), sanitize=True,
                                   removeHs=True, proximityBonding=True)
        warnings.warn('Ligand 2D topology was inferred from PDB coordinates. '
                      'Bond orders and aromaticity may be incorrect; use '
                      'ligand_file with a validated MOL/SDF/MOL2 for '
                      'chemical interpretation.', UserWarning, stacklevel=2)

    if mol is None:
        raise ValueError('RDKit could not construct ligand chemistry; '
                         'provide a validated 3D MOL/SDF/MOL2 file')
    if mol.GetNumConformers() != 1:
        raise ValueError('The ligand chemical structure must contain one 3D conformer')
    conf = mol.GetConformer()
    if not conf.Is3D():
        raise ValueError('Ligand chemical structure must contain 3D coordinates '
                         'consistent with the ProDy input structure')

    heavy_indices = [i for i, elem in enumerate(ligand_elements) if elem not in ('H', 'D')]
    if mol.GetNumAtoms() != len(heavy_indices):
        raise ValueError('Different ligand heavy-atom counts: RDKit={} and ProDy={}'.format(
            mol.GetNumAtoms(), len(heavy_indices)))
    rdkit_to_prody = {}
    used_prody = set()
    for atom in mol.GetAtoms():
        rd_index = atom.GetIdx()
        position = conf.GetAtomPosition(rd_index)
        xyz = np.array([position.x, position.y, position.z])
        atom_symbol = atom.GetSymbol().upper()
        candidates = [(j, float(np.linalg.norm(coords[j] - xyz)))
                      for j in heavy_indices if ligand_elements[j] == atom_symbol]
        matches = [(j, d) for j, d in candidates if d <= match_tolerance]
        if len(matches) != 1:
            raise ValueError('RDKit atom {} ({}) has {} coordinate matches '
                             'within {:.2f} Angstrom; check ligand_file and '
                             'ligand_selection'.format(rd_index, atom_symbol,
                                                       len(matches), match_tolerance))
        j, _ = matches[0]
        if j in used_prody:
            raise ValueError('More than one RDKit atom maps to ligand atom {}'.format(names[j]))
        used_prody.add(j)
        rdkit_to_prody[rd_index] = j

    term_names = [term for term in data['terms']
                  if term not in set(data.get('global_terms', []))]
    if not term_names:
        raise ValueError('No atom-resolved scoring terms are available')
    rows = [row for row in data['atoms'] if row.get('mapping_status') == 'MATCHED']
    by_name = {}
    for row in rows:
        name = str(row.get('ligand_atom', '')).strip()
        if name in by_name:
            raise ValueError('Multiple mapped SMINA records for ligand atom {}'.format(name))
        by_name[name] = row
    if any(names[j] not in by_name for j in heavy_indices):
        missing = [names[j] for j in heavy_indices if names[j] not in by_name]
        raise ValueError('Missing mapped SMINA contributions for heavy ligand atoms: '
                         + ', '.join(missing))

    scores = []
    hbond_terms = [term for term in term_names
                   if 'h_bond' in term.lower() or 'hbond' in term.lower()]
    if sort_by == 'hbond' and not hbond_terms:
        raise ValueError("sort_by='hbond' requires a hydrogen-bond scoring term")
    if sort_by not in ('favorable', 'net', 'unfavorable', 'hbond'):
        raise ValueError('sort_by must be favorable, net, unfavorable or hbond')

    for rd_index, prody_index in rdkit_to_prody.items():
        atom_name = names[prody_index]
        row = by_name[atom_name]
        values = [float(row[term]) for term in term_names]
        if not np.isfinite(values).all():
            raise ValueError('Nonfinite SMINA term for ligand atom {}'.format(atom_name))
        favorable = -sum(min(0.0, val) for val in values)
        unfavorable = sum(max(0.0, val) for val in values)
        hbond = -sum(min(0.0, float(row[term])) for term in hbond_terms)
        if sort_by == 'favorable':
            score = favorable
        elif sort_by == 'net':
            score = favorable - unfavorable
        elif sort_by == 'unfavorable':
            score = unfavorable
        else:
            score = hbond
        scores.append((score, rd_index, atom_name, favorable, unfavorable))

    ranked = sorted((r for r in scores if r[3] > 1e-12 or r[4] > 1e-12),
                    key=lambda r: (-r[0], r[2]))
    if not ranked:
        raise ValueError('All ligand heavy-atom scoring contributions are zero')
    selected = ranked if top is None else ranked[:top]
    max_magnitude = max(max(row[3], row[4]) for row in ranked)
    max_total_magnitude = max(row[3] + row[4] for row in ranked)

    # Construct a new 2D depiction, preserving the original 3D structure.
    drawing_mol = Chem.Mol(mol)
    AllChem.Compute2DCoords(drawing_mol)
    blue = (0.16, 0.41, 0.73)
    orange = (0.88, 0.34, 0.19)

    def weighted_color(base, magnitude):
        """Fade smaller contributions so equal-size color sectors cannot
        misleadingly suggest equal energetic contributions."""
        opacity = 0.20 + 0.80 * np.sqrt(magnitude / max_magnitude)
        return tuple(float(1.0 - (1.0 - part) * opacity) for part in base)
    highlights = {}
    radii = {}
    atom_labels = {}
    for _, idx, name, fav, unfav in selected:
        colors = []
        if fav > 1e-12:
            colors.append(weighted_color(blue, fav))
        if unfav > 1e-12:
            colors.append(weighted_color(orange, unfav))
        highlights[idx] = colors
        if atom_circle_mode == 'fractional':
            relative_size = (fav + unfav) / max_total_magnitude
        else:
            relative_size = max(fav, unfav) / max_magnitude
        radii[idx] = float(circle_scale) * (
            0.18 + 0.22 * np.sqrt(relative_size))
        atom_labels[idx] = name

    draw_width, draw_height = 1100, 700
    painter = rdMolDraw2D.MolDraw2DCairo(draw_width, draw_height)
    options = painter.drawOptions()
    options.fillHighlights = False
    options.continuousHighlight = False
    if atom_circle_mode == 'legacy':
        painter.DrawMoleculeWithHighlights(
            drawing_mol, '', highlights, {}, radii,
            {index: int(circle_linewidth) for index in highlights})
    else:
        painter.DrawMolecule(drawing_mol)

    atom_xy = {}
    for atom in drawing_mol.GetAtoms():
        p = painter.GetDrawCoords(atom.GetIdx())
        atom_xy[atom.GetIdx()] = np.asarray((p.x, p.y), dtype=float)

    conf_2d = drawing_mol.GetConformer()
    scale_samples = []
    for bond in drawing_mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        a3 = conf_2d.GetAtomPosition(i)
        b3 = conf_2d.GetAtomPosition(j)
        molecular_distance = float(np.hypot(a3.x - b3.x, a3.y - b3.y))
        pixel_distance = float(np.linalg.norm(atom_xy[i] - atom_xy[j]))
        if molecular_distance > 1e-8:
            scale_samples.append(pixel_distance / molecular_distance)
    pixels_per_unit = float(np.median(scale_samples)) if scale_samples else 45.0
    circle_radii_pixels = {index: radius * pixels_per_unit for index, radius in radii.items()}

    painter.FinishDrawing()
    from PIL import Image
    picture = Image.open(io.BytesIO(painter.GetDrawingText())).convert('RGB')

    glyph_paths = []
    chemical_paths = []

    if atom_symbol_glow or ligand_render == 'vector':
        vector_painter = rdMolDraw2D.MolDraw2DSVG(draw_width, draw_height)
        vector_options = vector_painter.drawOptions()
        vector_options.fillHighlights = False
        vector_options.continuousHighlight = False

        if atom_circle_mode == 'legacy':
            vector_painter.DrawMoleculeWithHighlights(
                drawing_mol, '', highlights, {}, radii,
                {index: int(circle_linewidth) for index in highlights})
        else:
            vector_painter.DrawMolecule(drawing_mol)
        for atom in drawing_mol.GetAtoms():
            index = atom.GetIdx()
            xy_svg = vector_painter.GetDrawCoords(index)
            if np.linalg.norm(atom_xy[index] - (xy_svg.x, xy_svg.y)) > 0.01:
                raise RuntimeError('RDKit SVG and Cairo atom coordinates do not match')
        vector_painter.FinishDrawing()

        def parse_glyph_path(commands):
            """Convert RDKit SVG glyph outlines into a Matplotlib Path.

            RDKit's atom-symbol outlines normally consist of absolute
            M/L/Q/C/Z commands; handling these directly avoids introducing
            any extra image conversion library dependency. """
            
            tokens = re.findall(r'[MLQCZmlqcz]|[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?', commands)
            vertices, codes = [], []
            pos, command = 0, None
            current = np.zeros(2, dtype=float)
            start = None
            n_params = {'M': 2, 'L': 2, 'Q': 4, 'C': 6}
            
            while pos < len(tokens):
                if tokens[pos].isalpha():
                    command = tokens[pos]
                    pos += 1
                if command is None:
                    raise ValueError('Malformed RDKit SVG atom-symbol path')
                if command.upper() == 'Z':
                    if start is None:
                        raise ValueError('Unexpected SVG path close')
                    vertices.append(tuple(start))
                    codes.append(MPLPath.CLOSEPOLY)
                    current = start.copy()
                    command = None
                    continue
                    
                kind = command.upper()
                if kind not in n_params:
                    raise ValueError('Unsupported SVG glyph path command: {}'.format(command))
                count = n_params[kind]
                if pos + count > len(tokens):
                    raise ValueError('Truncated RDKit SVG atom-symbol path')
                values = [float(t) for t in tokens[pos:pos + count]]
                pos += count
                relative = command.islower()
                segments = [np.asarray(values[k:k + 2], dtype=float)
                            for k in range(0, count, 2)]
                if relative:
                    segments = [v + current for v in segments]
                if kind == 'M':
                    current = segments[0]
                    start = current.copy()
                    vertices.append(tuple(current))
                    codes.append(MPLPath.MOVETO)
                    command = 'l' if relative else 'L'
                elif kind == 'L':
                    current = segments[0]
                    vertices.append(tuple(current))
                    codes.append(MPLPath.LINETO)
                elif kind == 'Q':
                    vertices.extend(tuple(v) for v in segments)
                    codes.extend((MPLPath.CURVE3, MPLPath.CURVE3))
                    current = segments[-1]
                elif kind == 'C':
                    vertices.extend(tuple(v) for v in segments)
                    codes.extend((MPLPath.CURVE4,) * 3)
                    current = segments[-1]
            return MPLPath(np.asarray(vertices, dtype=float), codes)

        svg_root = ET.fromstring(vector_painter.GetDrawingText())
        for element in svg_root.iter():
            if not element.tag.endswith('path'):
                continue
            path_data = element.attrib.get('d', '')
            if not path_data:
                continue
            path = parse_glyph_path(path_data)
            classes = element.attrib.get('class', '').split()
            is_glyph = (not any(cls.startswith('bond-') for cls in classes)
                        and any(cls.startswith('atom-') for cls in classes))
            css = {}
            for statement in element.attrib.get('style', '').split(';'):
                if ':' in statement:
                    k, v = statement.split(':', 1)
                    css[k.strip()] = v.strip()
            fill = element.attrib.get('fill', css.get('fill', 'none'))
            stroke = element.attrib.get('stroke', css.get('stroke', 'none'))
            if is_glyph:
                glyph_paths.append((path, fill))
            else:
                try:
                    linewidth = float(css.get('stroke-width', '2').replace('px',''))
                except ValueError as err:
                    raise ValueError('Invalid RDKit SVG stroke-width') from err
                chemical_paths.append((path, fill, stroke, linewidth,
                                       css.get('stroke-linecap', 'butt'),
                                       css.get('stroke-linejoin', 'miter')))

    def place_atom_labels():
        """Return (name, x, y) labels with simple collision avoidance."""

        from builtins import abs

        placed = []
        occupied = []
        all_atom_xy = np.asarray(list(atom_xy.values()), dtype=float)
        ligand_center = np.mean(all_atom_xy, axis=0)
        bond_segments = [(atom_xy[b.GetBeginAtomIdx()], atom_xy[b.GetEndAtomIdx()])
                         for b in drawing_mol.GetBonds()]
        label_fontsize = 11.0
        label_height = 22.0
        angles = np.arange(0.0, 2.0 * np.pi, np.pi / 8.0)
        directions = [np.asarray((np.cos(a), np.sin(a))) for a in angles]

        for _, index, name, fav, unfav in selected:
            center = atom_xy[index]
            outward = center - ligand_center
            outward_norm = float(np.linalg.norm(outward))
            if outward_norm > 1e-9:
                outward /= outward_norm
            else:
                outward = np.asarray((0.0, -1.0))
            directions_sorted = sorted(
                directions, key=lambda direction: -float(np.dot(direction, outward)))
            label_width = max(34.0, len(name) * 9.0 + 12.0)
            best_cost, best_position = float('inf'), None

            base_offset = max(62.0,
                circle_radii_pixels.get(index, 0.0) + 0.5 * label_width + 14.0)
            offsets = [base_offset + extra for extra in
                       (0.0, 12.0, 24.0, 38.0, 52.0, 68.0, 85.0)]
            for offset in offsets:
                for direction in directions_sorted:
                    candidate = center + offset * direction
                    x, y = candidate

                    if (x < label_width / 2 + 8 or
                            x > draw_width - label_width / 2 - 8 or
                            y < label_height / 2 + 10 or
                            y > draw_height - label_height / 2 - 10):
                        continue

                    cost = 0.32 * (offset - base_offset)
                    cost += 12.0 * (1.0 - float(np.dot(direction, outward)))

                    for other_index, other in atom_xy.items():
                        if other_index == index:
                            continue
                        dist = float(np.linalg.norm(candidate - other))
                        if dist < 38.0:
                            cost += 2.0 * (38.0 - dist) ** 2
                        if (abs(x - other[0]) < label_width / 2 + 13.0 and
                                abs(y - other[1]) < label_height / 2 + 13.0):
                            cost += 1800.0
                        if other_index in atom_labels and dist < 50.0:
                            cost += 0.8 * (50.0 - dist) ** 2

                        ring_radius = circle_radii_pixels.get(other_index, 0.0)
                        if ring_radius:
                            safe_distance = ring_radius + label_width / 2 + 5.0
                            if dist < safe_distance:
                                cost += 4.0 * (safe_distance - dist) ** 2

                    for start, end in bond_segments:
                        vec = end - start
                        denom = float(np.dot(vec, vec))
                        if denom < 1e-12:
                            continue
                        fraction = np.clip(np.dot(candidate - start, vec) / denom, 0, 1)
                        near = start + fraction * vec
                        bond_distance = float(np.linalg.norm(candidate - near))
                        if bond_distance < 23.0:
                            cost += 0.8 * (23.0 - bond_distance) ** 2

                    for old_x, old_y, old_width, old_height in occupied:
                        overlap_x = (label_width + old_width) / 2 + 8.0 - abs(x - old_x)
                        overlap_y = (label_height + old_height) / 2 + 6.0 - abs(y - old_y)
                        if overlap_x > 0 and overlap_y > 0:
                            cost += 4000.0 + 20.0 * overlap_x * overlap_y
                    if cost < best_cost:
                        best_cost, best_position = cost, candidate

            if best_position is None:
                best_position = np.clip(
                    center + base_offset * outward,
                    [label_width / 2 + 8, label_height / 2 + 10],
                    [draw_width - label_width / 2 - 8,
                     draw_height - label_height / 2 - 10])
            label_x, label_y = (float(v) for v in best_position)
            placed.append((name, label_x, label_y))
            occupied.append((label_x, label_y, label_width, label_height))
        return placed, label_fontsize

    if show_atom_labels:
        placed_labels, label_fontsize = place_atom_labels()
    else:
        placed_labels, label_fontsize = [], 11.0

    from prody import LOGGER

    contact_rows = []
    if show_residue_contacts:
        protein = atoms.select(protein_selection)
        if protein is None or protein.numAtoms() == 0:
            raise ValueError('protein_selection does not select protein atoms')
        protein_heavy = protein.select('noh')
        if protein_heavy is None or protein_heavy.numAtoms() == 0:
            raise ValueError('protein_selection contains no heavy atoms')
        protein_xyz = np.asarray(protein_heavy.getCoords(), dtype=float)
        if protein_xyz.ndim != 2 or protein_xyz.shape[1] != 3 or not np.isfinite(protein_xyz).all():
            raise ValueError('Invalid protein heavy-atom coordinates')
        lig_heavy_xyz = coords[heavy_indices]

        def _field(values, n):
            return [''] * n if values is None else [str(v).strip() for v in values]

        size = protein_heavy.numAtoms()
        p_segs = _field(protein_heavy.getSegnames(), size)
        p_chains = _field(protein_heavy.getChids(), size)
        p_names = _field(protein_heavy.getResnames(), size)
        p_nums = [int(v) for v in protein_heavy.getResnums()]
        p_icodes = _field(protein_heavy.getIcodes(), size)
        best_by_residue = {}

        for start in range(0, size, 1024):
            end = min(size, start + 1024)
            block_d2 = np.sum(
                (protein_xyz[start:end, None, :] - lig_heavy_xyz[None, :, :]) ** 2,
                axis=2)
            nearest_lig = np.argmin(block_d2, axis=1)
            nearest_d2 = block_d2[np.arange(end - start), nearest_lig]
            for offset, distance2 in enumerate(nearest_d2):
                if distance2 > contact_cutoff * contact_cutoff:
                    continue
                i = start + offset
                key = (p_segs[i], p_chains[i], p_names[i], p_nums[i], p_icodes[i])
                distance = float(np.sqrt(distance2))
                if key not in best_by_residue or distance < best_by_residue[key]['distance_A']:
                    prody_lig_index = heavy_indices[int(nearest_lig[offset])]
                    best_by_residue[key] = {
                        'segment': key[0], 'chain': key[1], 'resname': key[2],
                        'resnum': key[3], 'icode': key[4],
                        'distance_A': distance,
                        'ligand_index': prody_lig_index,
                        'ligand_atom': names[prody_lig_index]}

        contact_rows = list(best_by_residue.values())
        prody_to_rdkit = {j: index for index, j in rdkit_to_prody.items()}
        for row in contact_rows:
            row['rdkit_atom'] = prody_to_rdkit[row['ligand_index']]

        if residue_data is not None:
            if isinstance(residue_data, (str, Path)):
                with Path(residue_data).open('r', newline='', encoding='utf-8') as fd:
                    stored_rows = list(csv.DictReader(fd, delimiter='\t'))
            elif isinstance(residue_data, list):
                stored_rows = residue_data
            else:
                raise TypeError('residue_data must be a list of dictionaries or a TSV path')
            score_by_key = {}
            for row in stored_rows:
                if not isinstance(row, dict):
                    raise TypeError('residue_data rows must be dictionaries')
                try:
                    key = (str(row['chain']).strip(),
                           str(row['resname']).strip(),
                           int(row['resnum']), str(row.get('icode', '')).strip())
                    delta_score = float(row['delta_score_kcal_mol'])
                except (ValueError, TypeError, KeyError) as error:
                    raise ValueError('Invalid residue_data identity/score') from error
                if not np.isfinite(delta_score):
                    raise ValueError('Nonfinite residue_data final score')
                if key in score_by_key:
                    raise ValueError('Ambiguous residue_data entries for {}'.format(key))
                score_by_key[key] = delta_score
            for row in contact_rows:
                key = (row['chain'], row['resname'], row['resnum'], row['icode'])
                if key in score_by_key:
                    row['delta_score'] = score_by_key[key]

        effective_rank = ('favorable' if residue_data is not None else 'distance') \
            if residue_sort_by == 'auto' else residue_sort_by
        if effective_rank != 'distance':
            if residue_data is None:
                raise ValueError('residue_sort_by={!r} requires residue_data '
                                 '(TSV or residue list)'.format(effective_rank))
            missing = [row for row in contact_rows if 'delta_score' not in row]
            if missing:
                raise ValueError('Missing delta S for {} of {} contact residues. '
                                 'Use matching residue_data (same protein chain) '
                                 'or residue_sort_by="distance".'.format(
                                     len(missing), len(contact_rows)))

        by_distance = lambda row: (row['distance_A'], row['chain'],
                                   row['resnum'], row['resname'])
        if effective_rank == 'distance':
            ranked_contacts = sorted(contact_rows, key=by_distance)
        elif effective_rank == 'favorable':
            ranked_contacts = sorted(
                (row for row in contact_rows if row['delta_score'] < -1e-10),
                key=lambda row: (row['delta_score'], *by_distance(row)))
        elif effective_rank == 'unfavorable':
            ranked_contacts = sorted(
                (row for row in contact_rows if row['delta_score'] > 1e-10),
                key=lambda row: (-row['delta_score'], *by_distance(row)))
        elif effective_rank == 'absolute':
            ranked_contacts = sorted(
                contact_rows,
                key=lambda row: (-abs(row['delta_score']), *by_distance(row)))
        else: 
            negatives = sorted(
                (row for row in contact_rows if row['delta_score'] < -1e-10),
                key=lambda row: (row['delta_score'], *by_distance(row)))
            positives = sorted(
                (row for row in contact_rows if row['delta_score'] > 1e-10),
                key=lambda row: (-row['delta_score'], *by_distance(row)))
            if max_contact_residues is None:
                ranked_contacts = negatives + positives
            else:
                n = int(max_contact_residues)
                target_negative = (n + 1) // 2
                target_positive = n // 2
                selected_neg = negatives[:target_negative]
                selected_pos = positives[:target_positive]
                seen = {id(r) for r in selected_neg + selected_pos}
                remainder = sorted(
                    [r for r in negatives + positives if id(r) not in seen],
                    key=lambda r: (-abs(r['delta_score']), *by_distance(r)))
                ranked_contacts = (selected_neg + selected_pos +
                                   remainder[:max(0, n - len(selected_neg) - len(selected_pos))])
        if max_contact_residues is not None:
            ranked_contacts = ranked_contacts[:max_contact_residues]
        contact_rows = ranked_contacts

        aa_to_one = {
            'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
            'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I',
            'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
            'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
            'HSD':'H', 'HSE':'H', 'HSP':'H', 'HID':'H', 'HIE':'H',
            'HIP':'H', 'CYX':'C', 'CYM':'C', 'ASH':'D', 'GLH':'E',
            'MSE':'M'}
        
        for row in contact_rows:
            resname = (aa_to_one.get(row['resname'], row['resname'])
                       if residue_code == 'one' else row['resname'])
            resid = '{}{}{}'.format(resname, row['resnum'], row['icode'])
            row['label'] = ('{}:{}'.format(row['chain'], resid)
                            if show_chain and row['chain'] else resid)

        if not show_chain:
            label_counts = {}
            for row in contact_rows:
                label_counts[row['label']] = label_counts.get(row['label'], 0) + 1
            for row in contact_rows:
                if label_counts[row['label']] > 1 and row['chain']:
                    row['label'] = '{}:{}'.format(row['chain'], row['label'])

        score_max_abs = max((abs(row['delta_score']) for row in
                             best_by_residue.values() if 'delta_score' in row),
                            default=0.0)

        LOGGER.info('\nProtein-ligand geometric contacts: {} selected / {} within '
                    '{:.2f} Å | residue_sort_by={}'.format(
                        len(contact_rows), len(best_by_residue),
                        contact_cutoff, effective_rank))
        LOGGER.info('{:<14s} {:<12s} {:>12s} {:>12s}'.format(
            'Residue', 'Ligand atom', 'Dist. [Å]', 'Delta S'))
        LOGGER.info('-' * 58)
        for row in contact_rows:
            score_text = ('{:+.2f}'.format(row['delta_score'])
                          if 'delta_score' in row else '-')
            LOGGER.info('{:<14s} {:<12s} {:>12.2f} {:>12s}'.format(
                row['label'], row['ligand_atom'], row['distance_A'], score_text))

    LOGGER.info('\nSMINA ligand atom contributions: TOP {} / {} | sort_by={}'.format(
        len(selected), len(ranked), sort_by))
    LOGGER.info('{:<8s} {:>14s} {:>14s} {:>14s}'.format(
        'Atom', 'Favorable (-)', 'Penalty (+)', 'Net'))
    LOGGER.info('-' * 55)
    for _, idx, name, fav, unfav in selected:
        LOGGER.info('{:<8s} {:>+14.3f} {:>+14.3f} {:>+14.3f}'.format(
            name, -fav, unfav, unfav - fav))

    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 8) if show_residue_contacts
                               else (11.5, 7.0))
        fig.subplots_adjust(left=0.015, right=0.985, top=0.97,
                            bottom=0.15 if show_residue_contacts else 0.13)

    if ligand_render == 'vector':
        ax.set_xlim(0, draw_width)
        ax.set_ylim(draw_height, 0)
        for path, fill, stroke, linewidth, cap, join in chemical_paths:
            ax.add_patch(PathPatch(
                path,
                facecolor=fill if fill != 'none' else 'none',
                edgecolor=stroke if stroke != 'none' else 'none',
                linewidth=linewidth,
                capstyle=cap if cap in ('butt','round','projecting') else 'butt',
                joinstyle=join if join in ('miter','round','bevel') else 'miter',
                zorder=1.3))
    else:
        ax.imshow(picture, zorder=0)
    if atom_circle_mode == 'fractional':
        for _, index, _name, fav, unfav in selected:
            total = fav + unfav
            if total <= 1e-12:
                continue
            cx, cy = (float(v) for v in atom_xy[index])
            radius_px = circle_radii_pixels[index]
            band_width = min(radius_px * 0.40,
                             3.3 * float(circle_linewidth))
            blue_fraction = float(fav / total)
            if fav > 1e-12:
                ax.add_patch(Wedge(
                    (cx, cy), radius_px,
                    90.0, 90.0 + 360.0 * blue_fraction,
                    width=band_width,
                    facecolor=blue,
                    edgecolor='none', zorder=4))
            if unfav > 1e-12:
                ax.add_patch(Wedge(
                    (cx, cy), radius_px,
                    90.0 + 360.0 * blue_fraction, 450.0,
                    width=band_width,
                    facecolor=orange,
                    edgecolor='none', zorder=4))

    for glyph_path, glyph_color in glyph_paths:
        if atom_symbol_glow and atom_symbol_glow_width > 0:
            ax.add_patch(PathPatch(
                glyph_path, facecolor='white', edgecolor='white',
                linewidth=float(atom_symbol_glow_width),
                joinstyle='round', capstyle='round',
                zorder=18))
        ax.add_patch(PathPatch(
            glyph_path, facecolor=glyph_color, edgecolor='none',
            zorder=19))

    if show_residue_contacts:
        ax.set_xlim(-350, draw_width + 350)
        ax.set_ylim(draw_height + 45, -45)
        sides = {'left': [], 'right': []}
        for row in contact_rows:
            xy = atom_xy[row['rdkit_atom']]
            side = 'left' if xy[0] < draw_width / 2.0 else 'right'
            sides[side].append((row, xy))
        for side, entries in sides.items():
            entries.sort(key=lambda item: (item[1][1], item[0]['distance_A']))
            if not entries:
                continue
            n = len(entries)
            y_min, y_max = 72.0, draw_height - 72.0
            if n == 1:
                row_y = [float(np.clip(entries[0][1][1], y_min, y_max))]
            else:
                row_y = np.linspace(y_min, y_max, n)
            dest_x = -30 if side == 'left' else draw_width + 30
            align = 'right' if side == 'left' else 'left'
            for (row, xy), target_y in zip(entries, row_y):
                label = row['label']
                if show_contact_distances:
                    label += '  {:.2f} Å'.format(row['distance_A'])
                if 'delta_score' in row:
                    label += '  ΔS={:+.2f}'.format(row['delta_score'])
                fill_color, edge_color = 'white', '0.80'
                if color_residue_boxes_by_delta and 'delta_score' in row:
                    delta = row['delta_score']
                    if abs(delta) > 1e-10:
                        base = blue if delta < 0 else orange
                        fraction = (abs(delta) / score_max_abs
                                    if score_max_abs > 0 else 0)
                        intensity = min(0.84, (0.10 +
                                         0.68 * float(fraction) ** 1.8) *
                                         float(residue_box_contrast))
                        fill_color = tuple(
                            float(1.0 - (1.0 - c) * intensity) for c in base)
                        edge_color = tuple(float(0.75 * c) for c in base)
                annotation = ax.text(
                    dest_x, float(target_y), label,
                    ha=align, va='center', fontsize=float(residue_fontsize),
                    color='0.12',
                    bbox=dict(boxstyle='round,pad=0.28', fc=fill_color,
                              ec=edge_color, lw=1.0, alpha=1.0),
                    zorder=15, clip_on=False)
                if residue_text_glow and residue_text_glow_width > 0:
                    from matplotlib.patheffects import withStroke, Normal
                    annotation.set_path_effects([
                        withStroke(linewidth=float(residue_text_glow_width),
                                   foreground='white'),
                        Normal()
                    ])

                row['_label_artist'] = annotation
                row['_contact_xy'] = np.asarray(xy, dtype=float)
                row['_contact_side'] = side

        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
        for row in contact_rows:
            annotation = row['_label_artist']
            box_patch = annotation.get_bbox_patch()
            if box_patch is None:
                continue
            bbox = box_patch.get_window_extent(renderer)
            inv = ax.transData.inverted()
            c0 = inv.transform((bbox.x0, bbox.y0))
            c1 = inv.transform((bbox.x1, bbox.y1))
            left = min(c0[0], c1[0])
            right = max(c0[0], c1[0])
            target_y = 0.5 * (c0[1] + c1[1])
            side = row['_contact_side']

            target_x = right if side == 'left' else left
            end_xy = np.asarray((target_x, target_y), dtype=float)
            center = row['_contact_xy']
            delta = end_xy - center
            length = float(np.linalg.norm(delta))
            if length <= 1e-6:
                continue
            unit = delta / length
            atom_index = row['rdkit_atom']
            start_xy = center + unit * (circle_radii_pixels.get(atom_index, 0) + 5)
            dx = end_xy[0] - start_xy[0]
            bend = abs(dx) * 0.32
            sign = -1 if side == 'left' else 1
            cfirst = (start_xy[0] + sign * bend, start_xy[1])
            csecond = (end_xy[0] - sign * bend, end_xy[1])
            path = MPLPath(
                [tuple(start_xy), cfirst, csecond, tuple(end_xy)],
                [MPLPath.MOVETO, MPLPath.CURVE4,
                 MPLPath.CURVE4, MPLPath.CURVE4])
            ax.add_patch(PathPatch(
                path, edgecolor='0.49', facecolor='none', linewidth=1.0,
                capstyle='round', zorder=0.9, clip_on=False))

            del row['_label_artist'], row['_contact_xy'], row['_contact_side']

    import matplotlib.patheffects as path_effects
    for name, label_x, label_y in placed_labels:
        annotation = ax.text(label_x, label_y, name, fontsize=label_fontsize, 
                            color='black', ha='center', va='center', zorder=23,
                            clip_on=False)
        annotation.set_path_effects([
            path_effects.withStroke(linewidth=3.0, foreground='white')])

    ax.axis('off')
    legend_handles = [Patch(color=blue, label='Favorable (-)'),
                      Patch(color=orange, label='Unfavorable (+)')]
    if show_residue_contacts:
        from matplotlib.lines import Line2D
        legend_handles.append(Line2D([0], [0], color='0.47', lw=1.4,
                                     label='3D contact ≤ {:.1f} Å'.format(contact_cutoff)))
    ax.legend(
        handles=legend_handles,
        loc='lower center', ncol=len(legend_handles), frameon=False,
        bbox_to_anchor=(0.5, -0.04), fontsize=15 if show_residue_contacts else 16,
        handlelength=1.65, handleheight=1.25, columnspacing=1.65)
    return ax


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


def createFoldseekAlignment(prot_seq, prot_foldseek, **kwargs):
    """Aligns sequences from prot_seq with homologous sequences identified in prot_foldseek, 
    generating a multiple sequence alignment.
    
    :arg prot_seq: The natural sequence extracted from the PDB (seq file)
    :type prot_seq: str
    
    :arg prot_foldseek: The results from foldseek (foldseek file)
    :type prot_foldseek: str

    :arg msa_output_name: The natural sequence extracted from the PDB (msa file)
    :type msa_output_name: str
    """

    msa_output_name = kwargs.pop('msa_output_name', 'prot_struc.msa')

    def find_match_index(tar_nogap, nat_seq):
        tar_nogap_str = ''.join(tar_nogap)
        nat_seq_str = ''.join(nat_seq)
        index = nat_seq_str.find(tar_nogap_str)
        return index

    # Read input files
    with open(prot_seq, 'r') as f:
        file1 = f.readlines()

    with open(prot_foldseek, 'r') as f:
        file2 = f.readlines()

    # Open output file
    with open(msa_output_name, 'w') as fp:
        nat_seq = list(file1[0].strip())
        
        # Write the natural sequence to the output file
        fp.write(''.join(nat_seq) + "\n")
        
        # Process each foldseek entry
        for line in file2:
            entries = line.split()
            
            if float(entries[2]) >= 0.5:
                tar_seq = list(entries[11].strip())
                mat_seq = list(entries[12].strip())
                
                tar_nogap = []
                processed_mat = []
                
                for j in range(len(tar_seq)):
                    if tar_seq[j] != '-':
                        tar_nogap.append(tar_seq[j])
                        processed_mat.append(mat_seq[j])
                
                match_index = find_match_index(tar_nogap, nat_seq)
                end_index = match_index + len(tar_nogap)
                m = 0
                
                for l in range(len(nat_seq)):
                    if l < match_index:
                        fp.write("-")
                    elif l >= match_index and l < end_index:
                        fp.write(processed_mat[m])
                        m += 1
                    else:
                        fp.write("-")
                fp.write("\n")
    
    LOGGER.info("MSA file is now created, and saved as {}.".format(msa_output_name))


def extractMultiModelPDB(multimodelPDB, **kwargs):
    """Extracts individual PDB models from multimodel PDB and places them into the pointed directory.
    If used for calculating calcSignatureInteractions align the models.

    :arg multimodelPDB: The file containing models in multi-model PDB format 
    :type multimodelPDB: str
    
    :arg folder_name: The name of the folder to which PDBs will be extracted
    :type folder_name: str
    """
    
    import os
    
    folder_name = kwargs.pop('folder_name', 'struc_homologs')
    
    with open(multimodelPDB, 'r') as f:
        file = f.readlines()
    os.makedirs(folder_name, exist_ok=True)

    fp = None
    for line in file:
        line = line.strip()
        sig1 = line[:5]
        sig2 = line[:6]
        sig3 = line[:4]

        if sig1 == 'MODEL':
            model_number = line.split()[1]
            filename = 'model{}.pdb'.format(model_number)
            fp = open(filename, 'w')
            continue

        if sig2 == 'ENDMDL':
            if fp:
                fp.close()
            os.rename(filename, './{}/{}'.format(folder_name,filename))
            continue

        if sig3 == 'ATOM' and fp:
            fp.write("{}\n".format(line))
            
    LOGGER.info("Individual models are saved in {}.".format(folder_name))


def runFoldseek(pdb_file, chain, **kwargs):
    """This script processes a PDB file to extract a specified chain's sequence, searches for 
    homologous structures using foldseek, and prepares alignment outputs for further analysis.
    
    Before using the function, install Foldseek:
    >>> conda install bioconda::foldseek
    More information can be found:
    https://github.com/steineggerlab/foldseek?tab=readme-ov-file#databasesand

    This function will not work under Windows.
    Example usage: runFoldseek('5kqm.pdb', 'A', database_folder='~/Downloads/foldseek/pdb')
    where previous a folder called 'foldseek' were created and PDB database was uploaded using:
    >>> foldseek databases PDB pdb tmp   (Linux console) 
    
    :arg pdb_file: A PDB file path
    :type pdb_file: str

    :arg chain: Chain identifier
    :type chain: str

    :arg coverage_threshold: Coverage threshold 
            Default is 0.3
    :type coverage_threshold: float

    :arg tm_threshold: TM-score threshold
            Default is 0.5
    :type tm_threshold: float
    
    :arg database_folder: Folder with the database
            Default is '~/Downloads/foldseek/pdb'
            To download the database use: 'foldseek databases PDB pdb tmp' in the console 
    :type database_folder: str

    :arg fixer: The method for fixing lack of hydrogen bonds
            Default is 'pdbfixer'
    :type fixer: 'pdbfixer' or 'openbabel'
    
    :arg subset: subsets of atoms: 'ca', 'bb', 'heavy', 'noh', 'all'  (see matchChains())
            Default is 'ca'
    :type subset: str

    :arg seqid: Minimum value of the sequence identity (see matchChains())
            Default is 5
    :type seqid: float
    
    :arg overlap: percent overlap (see matchChains())
            Default is 50
    :type overlap: float

    :arg folder_name: Folder where the results will be collected
            Default is 'struc_homologs'
    :type folder_name: str """
    
    import os
    import subprocess
    import re
    import sys

    database_folder = kwargs.pop('database_folder', '~/Downloads/foldseek/pdb')
    coverage_threshold = kwargs.pop('coverage_threshold', 0.3)
    tm_threshold = kwargs.pop('tm_threshold', 0.5)    
    folder_name = kwargs.pop('folder_name', 'struc_homologs')
    subset = kwargs.pop('subset', 'ca')
    seqid = kwargs.pop('seqid', 5)
    overlap = kwargs.pop('overlap', 50)
    
    if not isinstance(pdb_file, str):
        raise TypeError('Please provide the name of the PDB file.')
    
    full_path = os.path.expanduser(database_folder)
    if not os.path.exists(full_path.strip('pdb')):
        raise ValueError('The required database is not found in {0}. Please download it first.'.format(database_folder.strip('pdb')))
    
    # Define the amino acid conversion function
    def aa_onelet(three_letter_code):
        codes = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'MSE': 'M'
        }
        return codes.get(three_letter_code)

    # Function to extract the sequence from the PDB file
    def extract_sequence_from_pdb(pdb_file, chain, output_file):
        sequence = []
        prev_resid = -9999

        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                if line.startswith("ATOM"):
                    ch = line[21]
                    resid = int(line[22:26].strip())
                    aa = line[17:20].strip()

                    if ch == chain and resid != prev_resid:
                        one_aa = aa_onelet(aa)
                        if one_aa:
                            sequence.append(one_aa)
                            prev_resid = resid

        with open(output_file, 'w') as seq_file:
            seq_file.write(''.join(sequence))

    # Inputs
    fpath = pdb_file.strip()
    chain = chain.strip()
    cov_threshold = float(coverage_threshold)
    tm_threshold = float(tm_threshold)
    
    # Filter PDB file
    awk_command = "awk '{{if (substr($0, 22, 1) == \"{0}\") print}}'".format(chain)
    subprocess.run("cat {0} | grep '^ATOM' | {1} > inp.pdb".format(fpath, awk_command), shell=True)

    # Run foldseek and other commands
    subprocess.run([
        'foldseek', 'easy-search', 'inp.pdb', database_folder, 'prot.foldseek',
        'tmp2', '--exhaustive-search', '1', '--format-output',
        "query,target,qstart,qend,tstart,tend,qcov,tcov,qtmscore,ttmscore,rmsd,qaln,taln",
        '-c', str(cov_threshold), '--cov-mode', '0'
    ])
    
    # Extract sequence and write to prot.seq
    extract_sequence_from_pdb('inp.pdb', chain, 'prot.seq')
    createFoldseekAlignment('prot.seq', 'prot.foldseek', msa_output_name='prot_struc.msa')    
    
    # Read input files
    with open('inp.pdb', 'r') as f:
        file1 = f.readlines()

    with open('prot.foldseek', 'r') as f:
        file2 = f.readlines()

    with open('prot_struc.msa', 'r') as f:
        file3 = f.readlines()

    # Open output files
    fp1 = open("aligned_structures.pdb", 'w')
    fp2 = open("shortlisted_resind.msa", 'w')
    fp6 = open("seq_match_reconfirm.txt", 'w')
    fp7 = open("aligned_structures_extended.pdb", 'w')
    fp9 = open("shortlisted.foldseek", 'w')
    fp10 = open("shortlisted.msa", 'w')

    # Initialize variables
    mod_count = 1
    fp1.write("MODEL\t{0}\n".format(mod_count))
    fp7.write("MODEL\t{0}\n".format(mod_count))

    seq1 = list(file3[0].strip())
    ind1 = 0
    prev_id = -9999

    # Process input PDB file
    for line in file1:
        line = line.strip()
        id_ = line[0:4]
        ch = line[21]
        resid = int(line[22:26].strip())
        aa = line[17:20]

        if id_ == 'ATOM' and ch == chain and resid != prev_id:
            prev_id = resid
            one_aa = aa_onelet(aa)
            if one_aa == seq1[ind1]:
                fp1.write("{0}\n".format(line))
                fp7.write("{0}\n".format(line))
                fp2.write("{0} ".format(resid))
                ind1 += 1
            else:
                print("Mismatch in sequence and structure of Query protein at Index {0}".format(ind1))
                break
        elif id_ == 'ATOM' and ch == chain and resid == prev_id:
            fp1.write("{0}\n".format(line))
            fp7.write("{0}\n".format(line))

    fp1.write("TER\nENDMDL\n\n")
    fp7.write("TER\nENDMDL\n\n")
    fp2.write("\n")
    fp10.write("{0}\n".format(file3[0].strip()))

    # Processing foldseek results
    os.makedirs("temp", exist_ok=True)

    for i, entry in enumerate(file2):
        entries = re.split(r'\s+', entry.strip())

        if float(entries[8]) < tm_threshold:
            continue

        tstart = int(entries[4])
        tend = int(entries[5])
        pdb = entries[1][:4]
        chain = entries[1][-1]
        fname = "{0}.pdb".format(pdb)

        # Download and process the target PDB file
        subprocess.run(['wget', '-P', 'temp', "https://files.rcsb.org/download/{0}".format(fname)])

        awk_command = "awk '{{if (substr($0, 22, 1) == \"{0}\") print}}'".format(chain)
        subprocess.run("cat ./temp/{0} | grep -E '^(ATOM|HETATM)' | {1} > target.pdb".format(fname, awk_command), shell=True)

        # Check if target.pdb has fewer than 5 lines
        if sum(1 for _ in open("target.pdb")) < 5:
            LOGGER.info("target.pdb has fewer than 5 lines. Skipping further processing.")
            subprocess.run(['rm', 'target.pdb'])
            continue

        with open('target.pdb', 'r') as target_file:
            file4 = target_file.readlines()

        qseq = list(entries[11])
        tseq = list(entries[12])

        start_line = 0
        prevind = -9999
        tarind = 0

        for j, line in enumerate(file4):
            resid = int(line[22:26].strip())
            if resid != prevind:
                prevind = resid
                tarind += 1
            if tarind == tstart:
                start_line = j
                break

        prevind = -9999
        ind2 = 0
        j = start_line
        flag2 = False

        with open("temp_coord.txt", 'w') as fp4, \
             open("temp_resind.txt", 'w') as fp3, \
             open("temp_seq.txt", 'w') as fp5:

            for k in range(len(qseq)):
                if tseq[k] != '-' and qseq[k] != '-':
                    line = file4[j].strip()
                    resid = int(line[22:26].strip())
                    aa = line[17:20]
                    one_aa = aa_onelet(aa)

                    if one_aa == tseq[k]:
                        fp3.write("{0} ".format(resid))
                        fp5.write("{0}".format(one_aa))
                        prevind = resid
                    else:
                        print("Mismatch in sequence and structure of Target protein {0}{1} at line {2} Index {3}-{4} ne {5}".format(
                            pdb, chain, j, k, one_aa, tseq[k]))
                        flag2 = True
                        break

                    while resid == prevind:
                        fp4.write("{0}\n".format(line))
                        j += 1
                        if j >= len(file4):
                            break
                        line = file4[j].strip()
                        resid = int(line[22:26].strip())
                elif tseq[k] != '-' and qseq[k] == '-':
                    line = file4[j].strip()
                    resid = int(line[22:26].strip())
                    aa = line[17:20]
                    one_aa = aa_onelet(aa)

                    if one_aa == tseq[k]:
                        prevind = resid
                    else:
                        print("Mismatch in sequence and structure of Target protein {0}{1} at Line {2} Index {3}-{4} ne {5}".format(
                            pdb, chain, j, k, one_aa, tseq[k]))
                        flag2 = True
                        break

                    while resid == prevind:
                        j += 1
                        if j >= len(file4):
                            break
                        line = file4[j].strip()
                        resid = int(line[22:26].strip())

                elif tseq[k] == '-' and qseq[k] != '-':
                    fp3.write("- ")
                    fp5.write("-")

            if flag2:
                continue

        with open("temp_coord.txt", 'r') as f:
            tmpcord = f.readlines()

        with open("temp_resind.txt", 'r') as f:
            tmpresind = f.readlines()

        with open("temp_seq.txt", 'r') as f:
            tmpseq = f.readlines()

        ind3 = i + 1
        seq1 = list(file3[ind3].strip())
        seq2 = tmpresind[0].strip().split()

        fp9.write("{0}".format(file2[i]))
        fp10.write("{0}\n".format(file3[ind3].strip()))

        for m in range(len(seq1)):
            if seq1[m] == '-':
                fp2.write("{0} ".format(seq1[m]))
            else:
                break

        for n in range(len(seq2)):
            fp2.write("{0} ".format(seq2[n]))
            fp6.write("{0}".format(seq1[m]))
            m += 1

        for o in range(m, len(seq1)):
            fp2.write("{0} ".format(seq1[m]))

        fp2.write("\n")
        fp6.write("\n{0}\n\n\n".format(tmpseq[0]))

        mod_count += 1
        fp1.write("MODEL\t{0}\n".format(mod_count))
        fp7.write("MODEL\t{0}\n".format(mod_count))

        for line in file4:
            if line.strip():
                fp7.write("{0}\n".format(line.strip()))

        for line in tmpcord:
            fp1.write("{0}".format(line))

        fp1.write("TER\nENDMDL\n\n")
        fp7.write("TER\nENDMDL\n\n")

    # Cleanup
    fp1.close()
    fp2.close()
    fp6.close()
    fp7.close()
    fp9.close()
    fp10.close()

    extractMultiModelPDB('aligned_structures.pdb', folder_name=folder_name)
    subprocess.run("rm -f inp.pdb prot.seq target.pdb temp_coord.txt temp_resind.txt temp_seq.txt", shell=True)
    subprocess.run("rm -rf tmp2 temp", shell=True)

    # New part
    pwd = os.path.abspath(folder_name)
    list_pdbs = [pwd+'/'+ff for ff in os.listdir(pwd) if ff.endswith('.pdb')]
    list_pdbs.sort(key=lambda x: int(''.join(filter(str.isdigit, x)))) # all structures will be aligned on model1.pdb as in the oryginal code
    
    LOGGER.info('Adding hydrogens to the structures..')
    new_pdbids = fixStructuresMissingAtoms(list_pdbs, method='pdbfixer', model_residues=True, overwrite=True)
    
    structures = parsePDB(new_pdbids)
    target = structures[0]
    rmsds = []
    for mobile in structures[1:]:
        try:
            LOGGER.info('Aligning the structures..')
            i = mobile.getTitle()
            LOGGER.info(i)
            matches = matchChains(mobile.protein, target.protein, subset=subset, seqid=seqid, overlap=overlap)
            m = matches[0]
            m0_alg, T = superpose(m[0], m[1], weights=m[0].getFlags("mapped"))
            rmsds.append(calcRMSD(m[0], m[1], weights=m[0].getFlags("mapped")))
            source_file = pwd+'/'+'align__'+i+'.pdb'
            writePDB(source_file, mobile)
        except:
            LOGGER.warn('There is a problem with {}. Change seqid or overlap parameter to include the structure.'.format(i))
    
    
def runDali(pdb, chain, **kwargs):
    """This function calls searchDali() and downloads all the PDB files, separate by chains, add hydrogens and missing side chains,
    and finally align them and put into the newly created folder.
       
    :arg pdb: A PDB code
    :type pdb: str

    :arg chain: chain identifier
    :type chain: str

    :arg cutoff_len: Length of aligned residues < cutoff_len
            (must be an integer or a float between 0 and 1)
            See searchDali for more details 
            Default is 0.5
    :type cutoff_len: float

    :arg cutoff_rmsd: RMSD cutoff (see searchDali)
            Default is 1.0
    :type cutoff_rmsd: float
    
    :arg cutoff_Z: Z score cutoff (see searchDali)
            Default is None
    :type cutoff_Z: float

    :arg cutoff_identity: RMSD cutoff (see searchDali)
            Default is None
    :type cutoff_identity: float

    :arg stringency: stringency for Dali cutoffs (see searchDali)
            Default is False
    :type stringency: bool

    :arg subset_Dali: fullPDB, PDB25, PDB50, PDB90
            Default is 'fullPDB'
    :type subset_Dali: str    
    
    :arg fixer: The method for fixing lack of hydrogen bonds
            Default is 'pdbfixer'
    :type fixer: 'pdbfixer' or 'openbabel'
    
    :arg subset: subsets of atoms: 'ca', 'bb', 'heavy', 'noh', 'all'  (see matchChains())
            Default is 'ca'
    :type subset: str

    :arg seqid: Minimum value of the sequence identity (see matchChains())
            Default is 5
    :type seqid: float
    
    :arg overlap: percent overlap (see matchChains())
            Default is 50
    :type overlap: float

    :arg folder_name: Folder where the results will be collected
            Default is 'struc_homologs'
    :type folder_name: str """

    import os
    import shutil
    from prody.database import dali 
    
    cutoff_len = kwargs.pop('cutoff_len', 0.5)
    cutoff_rmsd = kwargs.pop('cutoff_rmsd', 1.0)
    cutoff_Z = kwargs.pop('cutoff_Z', None)
    cutoff_identity = kwargs.pop('cutoff_identity', None)
    stringency = kwargs.pop('stringency', False)

    fixer = kwargs.pop('fixer', 'pdbfixer')
    subset_Dali = kwargs.pop('subset_Dali', 'fullPDB')
    subset = kwargs.pop('subset', 'ca')
    seqid = kwargs.pop('seqid', 5)
    overlap = kwargs.pop('overlap', 50)
    folder_name = kwargs.pop('folder_name', 'struc_homologs')
    
    dali_rec = dali.searchDali(pdb, chain, subset=subset_Dali)

    while not dali_rec.isSuccess:
        dali_rec.fetch()
    
    pdb_ids = dali_rec.filter(cutoff_len=cutoff_len, cutoff_rmsd=cutoff_rmsd,
                              cutoff_Z=cutoff_Z, cutoff_identity=cutoff_identity,
                              stringency=stringency)
    pdb_hits = [ (i[:4], i[4:]) for i in pdb_ids ]
    
    list_pdbs = []
    LOGGER.info('Separating chains and saving into PDB file')
    for i in pdb_hits:
        LOGGER.info('PDB code {} and chain {}'.format(i[0], i[1]))
        p = parsePDB(i[0]).select('chain '+i[1]+' and protein')
        writePDB(i[0]+i[1]+'.pdb', p)
        list_pdbs.append(i[0]+i[1]+'.pdb')
    
    LOGGER.info('Adding hydrogens to the structures..')
    new_pdbids = fixStructuresMissingAtoms(list_pdbs, method='pdbfixer', model_residues=True, overwrite=True)

    os.makedirs(folder_name)
    structures = parsePDB(new_pdbids)
    target = structures[0]
    rmsds = []
    for mobile in structures[1:]:
        try:
            LOGGER.info('Aligning the structures..')
            i = mobile.getTitle()
            LOGGER.info(i)
            matches = matchChains(mobile.protein, target.protein, subset=subset, seqid=seqid, overlap=overlap)
            m = matches[0]
            m0_alg, T = superpose(m[0], m[1], weights=m[0].getFlags("mapped"))
            rmsds.append(calcRMSD(m[0], m[1], weights=m[0].getFlags("mapped")))
            source_file = 'align__'+i+'.pdb'
            writePDB(source_file, mobile)
            shutil.move(source_file, os.path.join(folder_name, os.path.basename(source_file)))
        except:
            LOGGER.warn('There is a problem with {}. Change seqid or overlap parameter to include the structure.'.format(i))


def runBLAST(pdb, chain, **kwargs):
    """This function calls blastPDB to find homologs and downloads all of them in PDB format to the local directory,
    separate chains that were identified by BLAST, add hydrogens and missing side chains,
    and finally align them and put into the newly created folder.
       
    :arg pdb: A PDB code
    :type pdb: str

    :arg chain: chain identifier
    :type chain: str

    :arg fixer: The method for fixing lack of hydrogen bonds
            Default is 'pdbfixer'
    :type fixer: 'pdbfixer' or 'openbabel'
    
    :arg subset: subsets of atoms: 'ca', 'bb', 'heavy', 'noh', 'all'  (see matchChains())
            Default is 'bb'
    :type subset: str

    :arg seqid: Minimum value of the sequence identity (see matchChains())
            Default is 90
    :type seqid: float
    
    :arg overlap: percent overlap (see matchChains())
            Default is 50
    :type overlap: float

    :arg folder_name: Folder where the results will be collected
            Default is 'struc_homologs'
    :type folder_name: str """
    
    import os
    import shutil
    from prody.proteins.blastpdb import blastPDB
    
    fixer = kwargs.pop('fixer', 'pdbfixer')
    seqid = kwargs.pop('seqid', 90)
    overlap = kwargs.pop('overlap', 50)
    subset = kwargs.pop('subset', 'bb')
    folder_name = kwargs.pop('folder_name', 'struc_homologs')

    ref_prot = parsePDB(pdb)
    ref_hv = ref_prot.getHierView()[chain]
    sequence = ref_hv.getSequence()

    blast_record = blastPDB(sequence)
    while not blast_record.isSuccess:
        blast_record.fetch()

    pdb_hits = []
    for key, item in blast_record.getHits(seqid).items():
        pdb_hits.append((key, item['chain_id']))

    list_pdbs = []
    LOGGER.info('Separating chains and saving into PDB file')
    for i in pdb_hits:
        LOGGER.info('PDB code {} and chain {}'.format(i[0], i[1]))
        p = parsePDB(i[0]).select('chain '+i[1]+' and protein')
        writePDB(i[0]+i[1]+'.pdb', p)
        list_pdbs.append(i[0]+i[1]+'.pdb')
    
    LOGGER.info('Adding hydrogens to the structures..')
    new_pdbids = fixStructuresMissingAtoms(list_pdbs, method='pdbfixer', model_residues=True, overwrite=True)

    os.makedirs(folder_name)
    structures = parsePDB(new_pdbids)
    target = structures[0]
    rmsds = []
    for mobile in structures[1:]:
        try:
            LOGGER.info('Aligning the structures..')
            i = mobile.getTitle()
            LOGGER.info(i)
            matches = matchChains(mobile.protein, target.protein, subset=subset, seqid=seqid, overlap=overlap)
            m = matches[0]
            m0_alg, T = superpose(m[0], m[1], weights=m[0].getFlags("mapped"))
            rmsds.append(calcRMSD(m[0], m[1], weights=m[0].getFlags("mapped")))
            source_file = 'align__'+i+'.pdb'
            writePDB(source_file, mobile)
            shutil.move(source_file, os.path.join(folder_name, os.path.basename(source_file)))
        except:
            LOGGER.warn('There is a problem with {}. Change seqid or overlap parameter to include the structure.'.format(i))


def calcSignatureInteractions(PDB_folder, **kwargs):
    """Analyzes protein structures to identify various interactions using InSty. 
    Processes data from the MSA file and folder with selected models.
    
    Example usage: 
    >>> calcSignatureInteractions('./struc_homologs')
    >>> calcSignatureInteractions('./struc_homologs', mapping_file='shortlisted_resind.msa')

    :arg PDB_folder: Directory containing PDB model files
    :type PDB_folder: str

    :arg use_prefix: Whether to use perfix to select particular file names in the PDB_folder
            Default is True
    :type use_prefix: bool

    :arg mapping_file: Aligned residue indices, MSA file type
            required in Foldseek analyisis
    :type mapping_file: str """
    
    import os
    mapping_file = kwargs.get('mapping_file')
    use_prefix = kwargs.pop('use_prefix', True)
    
    if use_prefix == True:
        align_files = [os.path.join(PDB_folder, file) for file in os.listdir(PDB_folder) if file.startswith("align_")]

    else:
        align_files = [os.path.join(PDB_folder, file) for file in os.listdir(PDB_folder)]

    functions_dict = {
        "HBs": calcHydrogenBonds,
        "SBs": calcSaltBridges,
        "RIB": calcRepulsiveIonicBonding,
        "PiStack": calcPiStacking,
        "PiCat": calcPiCation,
        "HPh": calcHydrophobic,
        "DiBs": calcDisulfideBonds
    }

    if not mapping_file:
        for bond_type, func in functions_dict.items():
            for file in align_files:
                try:
                    LOGGER.info(file)
                    atoms = parsePDB(file)
                    interactions = func(atoms.select('protein'))
                    saveInteractionsAsDummyAtoms(atoms, interactions, filename=file.rstrip(file.split('/')[-1])+'INT_'+bond_type+'_'+file.split('/')[-1])
                except: pass
        
    if mapping_file:
        # for MSA file (Foldseek)
        n_per_plot = kwargs.pop('n_per_plot', None)
        min_height = kwargs.pop('min_height', 8)
        
        # Process each bond type
        for bond_type, func in functions_dict.items():
            # Check if the consensus file already exists
            consensus_file = '{}_consensus.txt'.format(bond_type)
            if os.path.exists(consensus_file):
                log_message("Consensus file for {} already exists, skipping.".format(bond_type), "INFO")
                continue

            log_message("Processing {}".format(bond_type))
            result = process_data(mapping_file, PDB_folder, func, bond_type, files_for_analysis=align_files)

            # Check if the result is None (no valid bonds found)
            if result is None:
                log_message("No valid {} entries found, skipping further processing.".format(bond_type), "WARNING")
                continue

            result, fixed_files = result

            # Proceed with plotting
            plot_barh(result, bond_type, n_per_plot=n_per_plot, min_height=min_height)
            


class Interactions(object):

    """Class for Interaction analysis of proteins."""

    def __init__(self, title='Unknown'):
        self._title = str(title).strip()
        self._atoms = None
        self._interactions = None
        self._interactions_matrix = None
        self._interactions_matrix_en = None
        self._energy_type = None
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

        LOGGER.info('Calculating interactions.')
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
        
        :arg replace: Used with selection criteria to set the new one
                      If set to **True** the selection will be replaced by the new one.
                      Default is **False**
        :type replace: bool

        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """
        
        replace = kwargs.pop('replace', False)
                
        if len(kwargs) != 0:
            results = [filterInteractions(j, self._atoms, **kwargs) for j in self._interactions]

            if replace == True:
                LOGGER.info('New interactions are set')
                self._interactions = results           
                self._hbs = results[0]
                self._sbs = results[1]
                self._rib = results[2]
                self._piStack = results[3] 
                self._piCat = results[4]
                self._hps = results[5]
                self._dibs = results[6]

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
        
        LOGGER.info('Calculating interaction matrix')
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


    def buildInteractionMatrixEnergy(self, **kwargs):
        """Build matrix with interaction energy comming from energy of pairs of specific residues.

        :arg energy_list_type: name of the list with energies 
                            Default is 'IB_solv'
                            acceptable values are 'IB_nosolv', 'IB_solv', 'CS'
        :type energy_list_type: str

        'IB_solv' and 'IB_nosolv' are derived from empirical potentials from
        O Keskin, I Bahar and colleagues from [OK98]_ and have RT units.

        'CS' is from MD simulations of amino acid pairs from Carlos Simmerling
        and Gary Wu in the [MR25]_ and have units of kcal/mol. 
        
        .. [MR25] Mikulska-Ruminska K, Krieger JM, Cao X, Banerjee A, Wu G, 
        Bogetti AT, Zhang F, Simmerling C, Coutsias EA, Bahar I
        InSty: a new module in ProDy for evaluating the interactions 
        and stability of proteins
        *Bioinformatics* **2025** 169009
        """
        
        import numpy as np
        import matplotlib
        import matplotlib.pyplot as plt
        from prody.dynamics.plotting import pplot
        
        atoms = self._atoms   
        interactions = self._interactions
        energy_list_type = kwargs.pop('energy_list_type', 'IB_solv')

        LOGGER.info('Calculating interaction energies matrix with type {0}'.format(energy_list_type))
        InteractionsMap = np.zeros([atoms.select('name CA').numAtoms(),atoms.select('name CA').numAtoms()])
        resIDs = list(atoms.select('name CA').getResnums())
        resChIDs = list(atoms.select('name CA').getChids())
        resIDs_with_resChIDs = list(zip(resIDs, resChIDs))
            
        for i in interactions:
            if i != []:
                for ii in i: 
                    m1 = resIDs_with_resChIDs.index((int(ii[0][3:]),ii[2]))
                    m2 = resIDs_with_resChIDs.index((int(ii[3][3:]),ii[5]))
                    scoring = get_energy([ii[0][:3], ii[3][:3]], energy_list_type)
                    if InteractionsMap[m1][m2] == 0:
                        InteractionsMap[m1][m2] = InteractionsMap[m2][m1] = float(scoring)

        self._interactions_matrix_en = InteractionsMap
        self._energy_type = energy_list_type
        
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
        :type filename: str  
        
        :arg energy: sum of the energy between residues
                    default is False
        :type energy: bool 
        """
        
        energy = kwargs.pop('energy', False)
        
        if not hasattr(self, '_interactions_matrix') or self._interactions_matrix is None:
            raise ValueError('Please calculate interactions matrix first.')

        if not isinstance(energy, bool):
            raise TypeError('energy should be True or False')
        
        import numpy as np
        from collections import Counter
        
        atoms = self._atoms 
        interaction_matrix = self._interactions_matrix
        interaction_matrix_en = self._interactions_matrix_en
        
        atoms = atoms.select("protein and noh")
        lista_ext = []
        aa_counter = Counter(atoms.getResindices())
        calphas = atoms.select('name CA')
        
        for i in range(calphas.numAtoms()):
            if energy == True:
                matrix_en_sum = np.sum(interaction_matrix_en, axis=0)
                lista_ext.extend(list(aa_counter.values())[i]*[round(matrix_en_sum[i], 8)]) 
            else:
                freq_contacts_residues = np.sum(interaction_matrix, axis=0)            
                lista_ext.extend(list(aa_counter.values())[i]*[round(freq_contacts_residues[i], 8)])        

        kw = {'occupancy': lista_ext}
        if 'filename' in kwargs:
            writePDB(kwargs['filename'], atoms, **kw)
            LOGGER.info('PDB file saved.')
        else:
            writePDB('filename', atoms, **kw)
            LOGGER.info('PDB file saved.')

    
    def getInteractors(self, residue_name):
        """ Provide information about interactions for a particular residue
        
        :arg residue_name: name of a resiude
                            example: LEU234A, where A is a chain name
        :type residue_name: str
        """
        
        atoms = self._atoms   
        interactions = self._interactions
        
        InteractionsMap = np.empty([atoms.select('name CA').numAtoms(),atoms.select('name CA').numAtoms()], dtype=object)
        resIDs = list(atoms.select('name CA').getResnums())
        resChIDs = list(atoms.select('name CA').getChids())
        resIDs_with_resChIDs = list(zip(resIDs, resChIDs))
        interaction_type = ['hb','sb','rb','ps','pc','hp','dibs']
        ListOfInteractions = []
        
        for nr,i in enumerate(interactions):
            if i != []:
                for ii in i: 
                    m1 = resIDs_with_resChIDs.index((int(ii[0][3:]),ii[2]))
                    m2 = resIDs_with_resChIDs.index((int(ii[3][3:]),ii[5]))
                    ListOfInteractions.append(interaction_type[nr]+':'+ii[0]+ii[2]+'-'+ii[3]+ii[5])
        
        aa_ListOfInteractions = []
        for i in ListOfInteractions:
            inter = i.split(":")[1:][0]
            if inter.split('-')[0] == residue_name or inter.split('-')[1] == residue_name:
                LOGGER.info(i)
                aa_ListOfInteractions.append(i)
        
        return aa_ListOfInteractions
        
    
    def getFrequentInteractors(self, contacts_min=2):
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
                           Default is 2.
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

                    if InteractionsMap[m1][m2] is None:
                        InteractionsMap[m1][m2] = []
            
                    InteractionsMap[m1][m2].append(interaction_type[nr] + ':' + ii[0] + ii[2] + '-' + ii[3] + ii[5])
            
        ListOfInteractions = [list(filter(None, [row[j] for row in InteractionsMap])) for j in range(len(InteractionsMap[0]))]
        ListOfInteractions = list(filter(lambda x: x != [], ListOfInteractions))
        ListOfInteractions_flattened = [j for sublist in ListOfInteractions for j in sublist]

        swapped_ListOfInteractions_list = []
        for interaction_group in ListOfInteractions_flattened:
            swapped_group = []
            for interaction in interaction_group:
                interaction_type, pair = interaction.split(':')
                swapped_pair = '-'.join(pair.split('-')[::-1])
                swapped_group.append("{}:{}".format(interaction_type, swapped_pair))
            swapped_ListOfInteractions_list.append(swapped_group)

        doubleListOfInteractions_list = ListOfInteractions_flattened+swapped_ListOfInteractions_list
        ListOfInteractions_list = [(i[0].split('-')[-1], [j.split('-')[0] for j in i]) for i in doubleListOfInteractions_list]

        merged_dict = {}
        for aa, ii in ListOfInteractions_list:
            if aa in merged_dict:
                merged_dict[aa].extend(ii)
            else:
                merged_dict[aa] = ii

        ListOfInteractions_list = [(key, value) for key, value in merged_dict.items()] 
        ListOfInteractions_list2 = [k for k in ListOfInteractions_list if len(k[-1]) >= contacts_min]
            
        for res in ListOfInteractions_list2:
            LOGGER.info('{0}  <--->  {1}'.format(res[0], '  '.join(res[1])))

        LOGGER.info('\nLegend: hb-hydrogen bond, sb-salt bridge, rb-repulsive ionic bond, ps-Pi stacking interaction,'
                             '\npc-Cation-Pi interaction, hp-hydrophobic interaction, dibs-disulfide bonds')
        
        try:
            from toolz.curried import count
        except ImportError:
            LOGGER.warn('This function requires the module toolz')
            return
        
        return ListOfInteractions_list2
        

    def showFrequentInteractors(self, cutoff=4, **kwargs):
        """Plots regions with the most frequent interactions.
        
        :arg cutoff: minimal score per residue which will be displayed.
                     If cutoff value is too big, top 30% with the higest values will be returned.
                     Default is 4.
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
            matplotlib.rcParams['font.size'] = '12' 
            fig = plt.figure(num=None, figsize=(16,5), facecolor='w')
        
        y_pos = np.arange(len(y))
        show = plt.bar(y_pos, x, align='center', alpha=0.5, color='blue')
        plt.xticks(y_pos, y, rotation=45, fontsize=16)
        plt.ylabel('Number of interactions', fontsize=16)
        plt.tight_layout()

        if SETTINGS['auto_show']:
            showFigure()
            
        dict_counts = dict(zip(y, x))
        dict_counts_sorted = dict(sorted(dict_counts.items(), key=lambda item: item[1], reverse=True))
            
        return dict_counts_sorted


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
        :type DiBs: int, float 

        :arg selstr: selection string for focusing the plot
        :type selection: str

        :arg energy: sum of the energy between residues
                    default is False
        :type energy: bool

        :arg energy_list_type: name of the list with energies
                            default is 'IB_solv'
                            acceptable values are 'IB_nosolv', 'IB_solv', 'CS'
        :type energy_list_type: str

        :arg overwrite_energies: whether to overwrite energies
                            default is False
        :type overwrite_energies: bool

        :arg percentile: a percentile threshold to remove outliers, i.e. only showing data within *p*-th
                        to *100-p*-th percentile. Default is None, so no axis limits.
        :type percentile: float

        :arg vmin: a minimum value threshold to remove outliers, i.e. only showing data greater than vmin
                This overrides percentile. Default is None, so no axis limits and little padding at the
                bottom when energy=True.
        :type vmin: float

        :arg vmax: a maximum value threshold to remove outliers, i.e. only showing data less than vmax
                This overrides percentile. Default is None, so no axis limits and a little padding for
                interaction type labels.
        :type vmax: float

        'IB_solv' and 'IB_nosolv' are derived from empirical potentials from
        O Keskin, I Bahar and colleagues from [OK98]_ and have RT units.

        'CS' is from MD simulations of amino acid pairs from Carlos Simmerling
        and Gary Wu for [MR25]_ and have units kcal/mol.
        """

        import numpy as np
        import matplotlib
        import matplotlib.pyplot as plt
        from prody.dynamics.plotting import pplot
        
        aa_dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
               'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HSE': 'H', 'HSD': 'H'}

        atoms = self._atoms
        energy = kwargs.pop('energy', False)
        
        if not isinstance(energy, bool):
            raise TypeError('energy should be True or False')
                    
        p = kwargs.pop('percentile', None)
        vmin = vmax = None
        if p is not None:
            vmin = np.percentile(matrix, p)
            vmax = np.percentile(matrix, 100-p)

        vmin = kwargs.pop('vmin', vmin)
        vmax = kwargs.pop('vmax', vmax)

        selstr = kwargs.pop('selstr', None)
        if selstr is not None:
            atoms = atoms.select(selstr)

        ResNumb = atoms.select('protein and name CA').getResnums()
        ResName = atoms.select('protein and name CA').getResnames()
        ResChid = atoms.select('protein and name CA').getChids()
        ResList = [ i[0]+str(i[1])+i[2] for i in list(zip([ aa_dic[i] for i in ResName ], ResNumb, ResChid)) ]
        
        if energy == True:
            matrix_en = self._interactions_matrix_en
            energy_list_type = kwargs.pop('energy_list_type', 'IB_solv')
            overwrite = kwargs.pop('overwrite_energies', False)
            if matrix_en is None or overwrite:
                LOGGER.warn('The energy matrix is recalculated with type {0}'.format(energy_list_type))
                self.buildInteractionMatrixEnergy(energy_list_type=energy_list_type)
                matrix_en = self._interactions_matrix_en

            elif self._energy_type != energy_list_type:
                LOGGER.warn('The energy type is {0}, not {1}'.format(self._energy_type, energy_list_type))

            matrix_en_sum = np.sum(matrix_en, axis=0)

            width = 0.8
            fig, ax = plt.subplots(num=None, figsize=(20,6), facecolor='w')
            matplotlib.rcParams['font.size'] = '24'

            if selstr is not None:
                matrix_en_sum = sliceAtomicData(matrix_en_sum,
                                                self._atoms.ca, selstr)

            zeros_row = np.zeros(matrix_en_sum.shape)
            pplot(zeros_row, atoms=atoms.ca, **kwargs)

            ax.bar(ResList, matrix_en_sum, width, color='blue')
            
            if vmin is None:
                vmin = np.min(matrix_en_sum) * 1.2

            if vmax is None:
                vmax = np.max(matrix_en_sum)

            plt.ylim([vmin, vmax])
            plt.tight_layout()    
            plt.xlabel('Residue')

            unit = 'RT' if energy_list_type in ['IB_solv', 'IB_nosolv'] else 'kcal/mol'
            plt.ylabel('Cumulative Energy [{0}]'.format(unit))
            plt.show()
            
            return matrix_en_sum
        
        else:
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

            all_ca = self._atoms.ca
            if selstr is not None:
                matrix_hbs_sum = sliceAtomicData(matrix_hbs_sum,
                                                 all_ca, selstr)
                matrix_sbs_sum = sliceAtomicData(matrix_sbs_sum,
                                                 all_ca, selstr)
                matrix_rib_sum = sliceAtomicData(matrix_rib_sum,
                                                 all_ca, selstr)
                matrix_pistack_sum = sliceAtomicData(matrix_pistack_sum,
                                                     all_ca, selstr)
                matrix_picat_sum = sliceAtomicData(matrix_picat_sum,
                                                   all_ca, selstr)
                matrix_hph_sum = sliceAtomicData(matrix_hph_sum,
                                                 all_ca, selstr)
                matrix_dibs_sum = sliceAtomicData(matrix_dibs_sum,
                                                  all_ca, selstr)

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

            if vmin is None:
                vmin = np.min(sum_matrix)

            if vmax is None:
                vmax = np.max(sum_matrix) * 1.5

            plt.ylim([vmin, vmax])
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
    
        :arg max_proc: maximum number of processes to use
            default is half of the number of CPUs
        :type max_proc: int

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

        start_frame = kwargs.pop('start_frame', 0)
        stop_frame = kwargs.pop('stop_frame', -1)
        max_proc = kwargs.pop('max_proc', mp.cpu_count()//2)

        if trajectory is None:
            if atoms.numCoordsets() > 1:
                trajectory = atoms
            else:
                LOGGER.info('Include trajectory or use multi-model PDB file.')
                return interactions_nb_traj

        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)
    
        if isinstance(trajectory, Trajectory):
            nfi = trajectory._nfi
            trajectory.reset()
        numFrames = trajectory._n_csets

        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame+1]

        HBs_all = [[] for _ in traj]
        SBs_all = [[] for _ in traj]
        RIB_all = [[] for _ in traj]
        PiStack_all = [[] for _ in traj]
        PiCat_all = [[] for _ in traj]
        HPh_all = [[] for _ in traj]
        DiBs_all = [[] for _ in traj]

        HBs_nb = [[] for _ in traj]
        SBs_nb = [[] for _ in traj]
        RIB_nb = [[] for _ in traj]
        PiStack_nb = [[] for _ in traj]
        PiCat_nb = [[] for _ in traj]
        HPh_nb = [[] for _ in traj]
        DiBs_nb = [[] for _ in traj]

        interactions_traj = [HBs_all, SBs_all, RIB_all, PiStack_all, PiCat_all, HPh_all, DiBs_all]
        interactions_nb_traj = [HBs_nb, SBs_nb, RIB_nb, PiStack_nb, PiCat_nb, HPh_nb, DiBs_nb]

        atoms_copy = atoms.copy()
        protein = atoms_copy.protein

        def analyseFrame(j0, frame0, interactions_all, interactions_nb):
            LOGGER.info('Frame: {0}'.format(j0))
            atoms_copy.setCoords(frame0.getCoords())

            ind = j0 - start_frame
            
            hydrogen_bonds = calcHydrogenBonds(protein, **kwargs)
            salt_bridges = calcSaltBridges(protein, **kwargs)
            RepulsiveIonicBonding = calcRepulsiveIonicBonding(protein, **kwargs)
            Pi_stacking = calcPiStacking(protein, **kwargs)
            Pi_cation = calcPiCation(protein, **kwargs)
            hydrophobic = calcHydrophobic(protein, **kwargs)
            Disulfide_Bonds = calcDisulfideBonds(protein, **kwargs)

            interactions_all[0][ind].extend(hydrogen_bonds)
            interactions_all[1][ind].extend(salt_bridges)
            interactions_all[2][ind].extend(RepulsiveIonicBonding)
            interactions_all[3][ind].extend(Pi_stacking)
            interactions_all[4][ind].extend(Pi_cation)
            interactions_all[5][ind].extend(hydrophobic)
            interactions_all[6][ind].extend(Disulfide_Bonds)

            interactions_nb[0][ind].append(len(hydrogen_bonds))
            interactions_nb[1][ind].append(len(salt_bridges))
            interactions_nb[2][ind].append(len(RepulsiveIonicBonding))
            interactions_nb[3][ind].append(len(Pi_stacking))
            interactions_nb[4][ind].append(len(Pi_cation))
            interactions_nb[5][ind].append(len(hydrophobic))
            interactions_nb[6][ind].append(len(Disulfide_Bonds))

        if max_proc == 1:
            interactions_all = interactions_traj
            interactions_nb = interactions_nb_traj
            for j0, frame0 in enumerate(traj, start=start_frame):
                analyseFrame(j0, frame0, interactions_all, interactions_nb)
            interactions_nb =  [[item[0] for item in row] for row in interactions_nb]
        else:
            with mp.Manager() as manager:
                interactions_all = manager.list()
                interactions_nb = manager.list()
                for row in interactions_traj:
                    interactions_all.append([manager.list() for _ in row])
                    interactions_nb.append([manager.list() for _ in row])

                j0 = start_frame
                while j0 < traj.numConfs()+start_frame:

                    processes = []
                    for _ in range(max_proc):
                        frame0 = traj[j0-start_frame]
                        
                        p = mp.Process(target=analyseFrame, args=(j0, frame0,
                                                                  interactions_all,
                                                                  interactions_nb))
                        p.start()
                        processes.append(p)

                        j0 += 1
                        if j0 >= traj.numConfs()+start_frame:
                            break

                    for p in processes:
                        p.join()

                interactions_all = [[item[:] for item in row] for row in interactions_all]
                interactions_nb =  [[item[0] for item in row] for row in interactions_nb]
        
        self._atoms = atoms
        self._traj = trajectory
        self._interactions_traj = interactions_all
        self._interactions_nb_traj = interactions_nb
        self._hbs_traj = interactions_all[0]
        self._sbs_traj = interactions_all[1]
        self._rib_traj = interactions_all[2]
        self._piStack_traj = interactions_all[3]
        self._piCat_traj = interactions_all[4]
        self._hps_traj = interactions_all[5]
        self._dibs_traj = interactions_all[6]
        
        if filename is not None:
            import pickle
            with open(str(filename)+'.pkl', 'wb') as f:
                pickle.dump(self._interactions_traj, f)  
            LOGGER.info('File with interactions saved.')

        if isinstance(trajectory, Trajectory):
            trajectory._nfi = nfi
            
        return interactions_nb


    def getInteractions(self, **kwargs):
        """Return the list of all interactions.
        
        :arg selection: selection string
        :type selection: str
    
        :arg selection2: selection string
        :type selection2: str 
        
        :arg replace: Used with selection criteria to set the new one
                      If set to **True** the selection will be replaced by the new one.
                      Default is **False**
        :type replace: bool

        Selection:
        If we want to select interactions for the particular residue or group of residues: 
            selection='chain A and resid 1 to 50'
        If we want to study chain-chain interactions:
            selection='chain A', selection2='chain B'  """
        
        replace = kwargs.pop('replace', False)
        
        if len(kwargs) != 0:
            sele_inter = []
            for i in self._interactions_traj:
                for nr_j,j in enumerate(i):
                    sele_inter.append(filterInteractions(i[nr_j], self._atoms, **kwargs))
            
            if replace == True:
                try:
                    trajectory = self._traj
                    numFrames = trajectory._n_csets
                except:
                    # If we analyze previously saved PKL file it doesn't have dcd information
                    # We have seven type of interactions. It will give number of frames.
                    numFrames = int(len(sele_inter)/7)
                    
                self._interactions_traj = sele_inter
                self._hbs_traj = sele_inter[0:numFrames]
                self._sbs_traj = sele_inter[numFrames:2*numFrames]
                self._rib_traj = sele_inter[2*numFrames:3*numFrames]
                self._piStack_traj = sele_inter[3*numFrames:4*numFrames]
                self._piCat_traj = sele_inter[4*numFrames:5*numFrames]
                self._hps_traj = sele_inter[5*numFrames:6*numFrames]
                self._dibs_traj = sele_inter[6*numFrames:7*numFrames]
                LOGGER.info('New interactions are set')
                
                self._interactions_nb_traj = None
                self._interactions_matrix_traj = None
                
                new_interactions_nb_traj = []
                new_interactions_nb_traj.append([ len(i) for i in self._hbs_traj ])
                new_interactions_nb_traj.append([ len(i) for i in self._sbs_traj ])
                new_interactions_nb_traj.append([ len(i) for i in self._rib_traj ])
                new_interactions_nb_traj.append([ len(i) for i in self._piStack_traj ])
                new_interactions_nb_traj.append([ len(i) for i in self._piCat_traj ])
                new_interactions_nb_traj.append([ len(i) for i in self._hps_traj ])
                new_interactions_nb_traj.append([ len(i) for i in self._dibs_traj ])
                self._interactions_nb_traj = new_interactions_nb_traj
                
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
        
            if isinstance(trajectory, Trajectory):
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
                    Default is True
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
                    Default is True. 
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
                          Default is 'all not (protein or water or ion)'.
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


def calcLigandResidueDistances(atoms, trajectory=None, ligand='UNL',
                               ligand_segname=None, ligand_resnum=None,
                               ligand_chain=None, protein_selection='protein',
                               contact_cutoff=5.0, filename=None,
                               save_statistics=False, save_per_frame=False,
                               output=False):
    """Calculate distances between protein residues and a ligand.
    
    For each analyzed conformation, the minimum distance between any heavy 
    atom of a protein residue and any heavy atom of the ligand is calculated, 
    regardless of the ``contact cutoff``. A residue is considered to be in 
    contact with the ligand when this distance is less than or equal to 
    ``contact_cutoff``. Contact persistence represents the percentage of 
    analyzed conformations in which such a contact occurs. Mean distances 
    and standard deviations are calculated across all analyzed conformations, 
    including those in which the residue is outside the contact cutoff. 

    :arg atoms: Protein-ligand structure, or the filename of a PDB file. PDB
        files with multi-model PDB records are read one model at a time.
    :type atoms: :class:`.Atomic` or str

    :arg trajectory: Optional DCD filename or ProDy trajectory object. Its atom
        count and atom order must correspond to *atoms*.
    :type trajectory: :class:`.Trajectory` or str or None

    :arg ligand: Residue name of the ligand. Default is ``'UNL'``.
    :type ligand: str

    :arg ligand_segname: Optional segment name specifying the ligand.
    :type ligand_segname: str or None

    :arg ligand_resnum: Optional residue number specifying the ligand.
    :type ligand_resnum: int or None

    :arg ligand_chain: Optional chain identifier specifying the ligand.
    :type ligand_chain: str or None

    :arg protein_selection: A selection string specifying the protein
        atoms to analyze. Default is ``'protein'``.
    :type protein_selection: str

    :arg contact_cutoff: Maximum heavy-atom distance (in Angstroms) defining
        a residue-ligand contact. Default is 5.0.
    :type contact_cutoff: float

    :arg filename: Output prefix for an optional compressed ``.npz`` file.
        If None (default), no files are written.
    :type filename: str or None

    :arg save_statistics: If True, also write a readable ``_statistics.dat``
        table. Requires *filename*. Default is False.
    :type save_statistics: bool

    :arg save_per_frame: If True, also write a ``_per_frame.dat`` table of all
        residue distances. Requires *filename*. Default is False.
    :type save_per_frame: bool

    :arg output: If True, report per-residue statistics.
        Default is False.
    :type output: bool

    :returns: Dictionary containing the full frame-by-residue distance matrix,
        residue identities, frame numbers, mean distances, standard deviations,
        and contact persistence. The result can be used directly with
        :func:`.showLigandResidueDistances`.
    :rtype: dict """
    
    import os
    import gzip
    from scipy.spatial import cKDTree
    from prody.atomic import AAMAP

    if not np.isfinite(contact_cutoff) or contact_cutoff <= 0:
        raise ValueError('contact_cutoff must be a positive number.')
    if filename is None and (save_statistics or save_per_frame):
        raise ValueError('A filename is required for optional DAT output.')
    if not isinstance(protein_selection, str) or not protein_selection.strip():
        raise ValueError('protein_selection must be a nonempty selection string.')

    is_path = isinstance(atoms, (str, os.PathLike))
    pdb_path = os.fspath(atoms) if is_path else None
    has_trajectory = trajectory is not None
    
    if is_path:
        if not os.path.isfile(pdb_path):
            raise IOError('PDB file not found: ' + pdb_path)
        structure = parsePDB(pdb_path, model=1)
    else:
        structure = atoms

    if structure is None or not hasattr(structure, 'select'):
        raise TypeError('atoms must be a PDB filename or a ProDy Atomic object.')
    if structure.numAtoms() < 1:
        raise ValueError('Input structure has no atoms.')

    ligand_selection = 'resname ' + str(ligand)
    if ligand_segname is not None:
        ligand_selection += ' and segname ' + str(ligand_segname)
    if ligand_resnum is not None:
        ligand_selection += ' and resnum ' + str(int(ligand_resnum))
    if ligand_chain is not None:
        ligand_selection += ' and chain ' + str(ligand_chain)
    
    ligand_atoms = structure.select('(' + ligand_selection + ') and noh')
    protein_atoms = structure.select('(' + protein_selection + ') and noh')
    
    if ligand_atoms is None:
        raise ValueError('No ligand heavy atoms for: ' + ligand_selection)
    if protein_atoms is None:
        raise ValueError('No protein heavy atoms for: ' + protein_selection)

    ligand_ids = np.asarray(ligand_atoms.getIndices(), dtype=np.int64)
    protein_ids = np.asarray(protein_atoms.getIndices(), dtype=np.int64)
    n_atoms = structure.numAtoms()
    
    if max(int(ligand_ids.max()), int(protein_ids.max())) >= n_atoms:
        raise ValueError('Provide the full AtomGroup, not an indexed Selection.')

    def _metadata(selection, getter, default=''):
        method = getattr(selection, getter, None)
        values = method() if method is not None else None
        if values is None:
            return [default] * selection.numAtoms()
        return list(values)

    ligand_keys = set()
    for ch, seg, rn, ic in zip(_metadata(ligand_atoms, 'getChids'),
                               _metadata(ligand_atoms, 'getSegnames'),
                               _metadata(ligand_atoms, 'getResnums', 0),
                               _metadata(ligand_atoms, 'getIcodes')):
        ligand_keys.add((str(ch).strip(), str(seg).strip(),
                         int(rn), str(ic).strip()))
    
    if len(ligand_keys) != 1:
        raise ValueError('Ligand selection contains %d residues; restrict it '
                         'using ligand_chain, ligand_segname or ligand_resnum.' %
                         len(ligand_keys))

    resnames = _metadata(protein_atoms, 'getResnames')
    resnums = _metadata(protein_atoms, 'getResnums', 0)
    chids = _metadata(protein_atoms, 'getChids')
    segnames = _metadata(protein_atoms, 'getSegnames')
    icodes = _metadata(protein_atoms, 'getIcodes')

    residue_keys = []
    residue_lookup = {}
    atom_resindices = np.empty(protein_atoms.numAtoms(), dtype=np.int32)
    
    for i, (rn, num, ch, seg, ic) in enumerate(zip(resnames, resnums, chids,
                                                   segnames, icodes)):
        key = (str(ch).strip(), str(seg).strip(), str(rn).strip(),
               int(num), str(ic).strip())
        if key not in residue_lookup:
            residue_lookup[key] = len(residue_keys)
            residue_keys.append(key)
        atom_resindices[i] = residue_lookup[key]

    # Prefer PDB chain identifiers where they uniquely identify structural
    # groups. If several segments share a chain ID, distinguish them by their
    # segment names.
    chain_segments = {}
    for ch, seg, _, _, _ in residue_keys:
        chain_segments.setdefault(ch, set()).add(seg)
    group_lookup = {}
    used_labels = set()
    for ch, seg in [(key[0], key[1]) for key in residue_keys]:
        pair = (ch, seg)
        if pair in group_lookup:
            continue
        if ch and len(chain_segments[ch]) == 1:
            label, category = ch, 'Chain'
        elif seg:
            label, category = seg, 'Segment'
        elif ch:
            label, category = ch, 'Chain'
        else:
            label, category = '?', 'Group'
        if label in used_labels:
            label = (ch or '?') + '/' + (seg or '?')
        if label in used_labels:
            raise ValueError('Cannot distinguish the protein chain/segment groups.')
        used_labels.add(label)
        group_lookup[pair] = (label, category)

    groups = [group_lookup[(key[0], key[1])][0] for key in residue_keys]
    group_types = [group_lookup[(key[0], key[1])][1] for key in residue_keys]

    # The nearest ligand atom is found for each protein atom. The minimum of
    # these distances within each residue yields the same minimum heavy-atom
    # distance as a residue-by-ligand atom search.
    atom_order = np.argsort(atom_resindices, kind='stable')
    sorted_indices = atom_resindices[atom_order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_indices)) + 1]

    def _pdb_frames(path):
        """Read coordinates from PDB files, including multi-model and gzipped PDBs.
           Yield the model index and atomic coordinates for each model without
           loading the entire structure into memory. """
        opener = gzip.open if path.lower().endswith('.gz') else open
        with opener(path, 'rt') as handle:
            coords = []
            model_found = False
            model_index = 0
            for line in handle:
                record = line[:6].strip()
                if record == 'MODEL':
                    model_found = True
                    if coords:
                        raise ValueError('MODEL before ENDMDL in ' + path)
                    continue
                if record in ('ATOM', 'HETATM'):
                    if line[16:17] not in (' ', 'A'):
                        continue
                    try:
                        coords.append((float(line[30:38]), float(line[38:46]),
                                       float(line[46:54])))
                    except ValueError:
                        raise ValueError('Invalid PDB coordinates: ' + line.rstrip())
                elif record == 'ENDMDL':
                    if coords:
                        model_index += 1
                        yield model_index, np.asarray(coords, dtype=float)
                        coords = []
                elif record == 'END' and not model_found:
                    break
            if coords:
                model_index += 1
                yield model_index, np.asarray(coords, dtype=float)

    opened_traj = False
    if has_trajectory:
        if isinstance(trajectory, (str, os.PathLike)):
            traj = Trajectory(os.fspath(trajectory))
            opened_traj = True
        else:
            traj = trajectory
        if traj.numAtoms() != n_atoms:
            if opened_traj:
                traj.close()
            raise ValueError('Trajectory and PDB atom counts differ; atom order '
                             'must also be identical.')
        traj.reset()
        n_frames = int(traj.numFrames())
        frames = ((i + 1, frame.getCoords()) for i, frame in enumerate(traj))
    elif is_path:
        n_frames = None
        frames = _pdb_frames(pdb_path)
    else:
        n_frames = int(structure.numCoordsets())
        if n_frames < 1:
            raise ValueError('No conformations found in the input Atomic object.')
        acs_index = structure.getACSIndex()

        def _atomic_frames():
            try:
                for i in range(n_frames):
                    structure.setACSIndex(i)
                    yield i + 1, structure.getCoords()
            finally:
                structure.setACSIndex(acs_index)
        frames = _atomic_frames()

    distances_list = []
    frame_numbers = []
    first_coords = structure.getCoords()
    first_checked = False
    try:
        for frame_number, xyz in frames:
            xyz = np.asarray(xyz)
            if xyz.shape != (n_atoms, 3):
                raise ValueError('Frame %d: expected coordinates for %d atoms, '
                                 'received %d.' % (frame_number, n_atoms, len(xyz)))
            if is_path and not has_trajectory and not first_checked:
                if not np.allclose(xyz, first_coords, atol=0.025, rtol=0):
                    raise ValueError('PDB MODEL 1 does not match parsePDB atom order.')
                first_checked = True
            ligxyz = xyz[ligand_ids]
            proxyz = xyz[protein_ids]
            if not (np.isfinite(ligxyz).all() and np.isfinite(proxyz).all()):
                raise ValueError('Non-finite coordinates in frame %d.' % frame_number)
            per_atom = cKDTree(ligxyz).query(proxyz)[0]
            per_residue = np.minimum.reduceat(per_atom[atom_order], starts)
            distances_list.append(per_residue.astype(np.float32))
            frame_numbers.append(frame_number)
    finally:
        if opened_traj:
            traj.close()

    if not distances_list:
        raise ValueError('No conformations were analyzed.')
    if n_frames is not None and len(frame_numbers) != n_frames:
        raise ValueError('Frame count differs from the number analyzed.')

    distances = np.vstack(distances_list)
    means = np.mean(distances, axis=0)
    sds = np.std(distances, axis=0)
    persistence = 100.0 * np.mean(distances <= contact_cutoff, axis=0)

    labels = [AAMAP.get(key[2], key[2]) + str(key[3]) + key[4] for key in residue_keys]

    result = {
        'distances': distances,
        'frame_numbers': np.asarray(frame_numbers, dtype=np.int32),
        'chains': np.asarray([key[0] for key in residue_keys], dtype='U32'),
        'segments': np.asarray([key[1] for key in residue_keys], dtype='U32'),
        'groups': np.asarray(groups, dtype='U64'),
        'group_types': np.asarray(group_types, dtype='U16'),
        'resnames': np.asarray([key[2] for key in residue_keys], dtype='U8'),
        'resnums': np.asarray([key[3] for key in residue_keys], dtype=np.int32),
        'icodes': np.asarray([key[4] for key in residue_keys], dtype='U8'),
        'labels': np.asarray(labels, dtype='U32'),
        'mean': means, 'sd': sds, 'persistence': persistence,
        'contact_cutoff': np.asarray(float(contact_cutoff)),
        'ligand': np.asarray(str(ligand)),
        'format_version': np.asarray(2, dtype=np.int16)}

    if filename is not None:
        prefix = os.fspath(filename)
        parent = os.path.dirname(os.path.abspath(prefix))
        if not os.path.isdir(parent):
            os.makedirs(parent)
        np.savez_compressed(prefix + '.npz', **result)
    
        LOGGER.info('Saved residue-ligand distances to %s.npz.' % prefix)
    
        if save_statistics:
            with open(prefix + '_statistics.dat', 'w') as handle:
                handle.write('# Group Residue Mean_A SD_A Persistence_percent '
                             'Contact_frames Total_frames\n')
                for i in range(len(residue_keys)):
                    handle.write('%s %s %.5f %.5f %.3f %d %d\n' %
                                 (groups[i], labels[i], means[i], sds[i],
                                  persistence[i],
                                  np.count_nonzero(distances[:, i] <= contact_cutoff),
                                  distances.shape[0]))
            LOGGER.info('Saved %s_statistics.dat.' % prefix)
    
        if save_per_frame:
            with open(prefix + '_per_frame.dat', 'w') as handle:
                handle.write('# Frame ' + ' '.join(groups[i] + ':' + labels[i]
                                                  for i in range(len(labels))) + '\n')
                for num, row in zip(frame_numbers, distances):
                    handle.write('%d %s\n' %
                                 (num, ' '.join('%.5f' % val for val in row)))
            LOGGER.info('Saved %s_per_frame.dat.' % prefix)

    LOGGER.info('Analyzed %d conformations, %d residues and ligand %s '
                '(contact cutoff %.2f A).' %
                (distances.shape[0], len(residue_keys), ligand, contact_cutoff))
    if output:
        LOGGER.info('Group Residue Mean_A SD_A Persistence_percent')
        for i in range(len(residue_keys)):
            LOGGER.info('%s %s %.3f %.3f %.1f' %
                        (groups[i], labels[i], means[i], sds[i], persistence[i]))
    
    return result


def showLigandResidueDistances(data, min_persistence=10.0,
                               contact_cutoff=None, top_n=18, residues=None,
                               balance_groups=True, filename=None, **kwargs):
    """Show a bar plot of residue-ligand minimum distances and persistence.
    
    Each bar shows the mean of minimum heavy-atom distances across all
    analyzed conformations. Error bars show the corresponding standard
    deviation. Bar colors represent contact persistence (percentage of
    conformations with a minimum distance <= *contact_cutoff*). Only the
    plotted residues are affected by *min_persistence* and *top_n*.

    :arg data: Results of :func:`.calcLigandResidueDistances`, or a previously
        saved ``.npz`` file.
    :type data: dict or str

    :arg min_persistence: Minimum contact persistence [%] for display.
        Default is 10.0. Use 0 to include all residues.
    :type min_persistence: float

    :arg contact_cutoff: Contact threshold (Angstroms) used to calculate
        persistence. If None, use the cutoff stored in *data*. Default is None.
        Changing it does not modify the calculated mean or standard deviation.
    :type contact_cutoff: float or None

    :arg top_n: Maximum number of residues to show, ranked by contact
        persistence. Default is 18. Use 0 or None to show all eligible residues.
    :type top_n: int or None

    :arg residues: Optional comma-separated string or list of residue labels,
        e.g. ``'L420,Y185'`` or ``'A:L420,B:R75'``. Overrides *top_n*.
    :type residues: str or list or None

    :arg balance_groups: If True (default), distribute *top_n* between protein
        chains/segments, then fill any remaining positions by persistence.
    :type balance_groups: bool

    :arg filename: Optional filename prefix for PNG figures. If None (default), 
        no figure files are written.
    :type filename: str or None

    :arg save_statistics: If True, write ``_plotted_residues.dat`` alongside
        the figure files. Requires *filename*. Default is False.
    :type save_statistics: bool

    :arg label_fontsize: Font size of residue names. Default is 24.
    :type label_fontsize: int

    :arg axis_fontsize: Font size of axis and colorbar labels. Default is 25.
    :type axis_fontsize: int

    :arg tick_fontsize: Font size of axis tick labels. Default is 21.
    :type tick_fontsize: int

    :arg header_fontsize: Font size of chain/segment headings. Default is 25.
    :type header_fontsize: int

    :arg colorbar_fontsize: Font size of colorbar tick labels. Default is 21.
    :type colorbar_fontsize: int

    :arg fig_width: Figure width in inches. Default is 11.5.
    :type fig_width: float

    :arg row_height: Height in inches per bar. Default is 0.66.
    :type row_height: float

    :arg png_dpi: Resolution for PNG figures. Default is 300.
    :type png_dpi: int

    :arg xlim: Optional maximum value of the distance axis (Angstroms).
        If None, the upper limit is determined automatically.
    :type xlim: float or None

    :arg show: Whether to display the figure. If None, follow
        ``SETTINGS['auto_show']``.
    :type show: bool or None  """
    
    import os
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    from matplotlib.cm import ScalarMappable

    balance_groups = kwargs.pop('balance_groups', True)
    save_statistics = kwargs.pop('save_statistics', False)

    label_fontsize = kwargs.pop('label_fontsize', 24)
    axis_fontsize = kwargs.pop('axis_fontsize', 25)
    tick_fontsize = kwargs.pop('tick_fontsize', 21)
    header_fontsize = kwargs.pop('header_fontsize', 25)
    colorbar_fontsize = kwargs.pop('colorbar_fontsize', 21)

    fig_width = kwargs.pop('fig_width', 11.5)
    row_height = kwargs.pop('row_height', 0.66)
    png_dpi = kwargs.pop('png_dpi', 300)
    xlim = kwargs.pop('xlim', None)
    show = kwargs.pop('show', None)

    if filename is None and save_statistics:
        raise ValueError('A filename is required for optional DAT output.')
    if isinstance(data, (str, os.PathLike)):
        with np.load(os.fspath(data), allow_pickle=False) as saved:
            result = {key: saved[key] for key in saved.files}
    else:
        result = data
    dist = np.asarray(result['distances'])
   
    if dist.ndim != 2 or dist.size == 0:
        raise ValueError('Distance matrix is missing or empty.')
    if not np.isfinite(min_persistence) or not 0 <= float(min_persistence) <= 100:
        raise ValueError('min_persistence must be between 0 and 100.')
   
    cutoff = (float(np.asarray(result['contact_cutoff']).item())
              if contact_cutoff is None else float(contact_cutoff))
   
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('contact_cutoff must be a positive number.')
    if top_n is not None:
        top_n = int(top_n)
        if top_n < 0:
            raise ValueError('top_n must be nonnegative or None.')
        if top_n == 0:
            top_n = None
    if fig_width <= 0 or row_height <= 0 or png_dpi <= 0:
        raise ValueError('Figure dimensions and PNG resolution must be positive.')
    if xlim is not None and (not np.isfinite(xlim) or xlim <= 0):
        raise ValueError('xlim must be a positive number.')

    persistence = np.mean(dist <= cutoff, axis=0) * 100.0
    mean = np.mean(dist, axis=0)
    sd = np.std(dist, axis=0)
    labels = np.asarray(result['labels']).astype(str)
    resnames = np.asarray(result['resnames']).astype(str)
    resnums = np.asarray(result['resnums'], dtype=int)
    groups = np.asarray(result.get('groups', result['chains'])).astype(str)
    group_types = np.asarray(result.get('group_types',
                                       np.full(len(groups), 'Chain'))).astype(str)
    eligible = np.flatnonzero(persistence >= float(min_persistence))
    if eligible.size == 0:
        LOGGER.warn('No residues satisfy the persistence threshold %.1f%%.' %
                    min_persistence)
        return []

    def _rank(i):
        return (-persistence[i], mean[i], int(resnums[i]), labels[i])

    requested_order = None
    
    if residues is not None:
        wanted = ([item.strip() for item in residues.split(',')]
                  if isinstance(residues, str)
                  else [str(item).strip() for item in residues])
        wanted = [item for item in wanted if item]
        if not wanted:
            raise ValueError('The residue selection is empty.')
        requested_order = {}
        missing = []
        for j, item in enumerate(wanted):
            matches = [int(i) for i in eligible if item.upper() in
                       (labels[i].upper(),
                        (groups[i] + ':' + labels[i]).upper(),
                        (resnames[i] + str(resnums[i])).upper(),
                        (groups[i] + ':' + resnames[i] + str(resnums[i])).upper())]
            if not matches:
                missing.append(item)
            for i in matches:
                if i not in requested_order:
                    requested_order[i] = j
        if missing:
            raise ValueError('Residues not found or below threshold: %s.' %
                             ', '.join(missing))
        selected = list(requested_order)
        
    else:
        selected = sorted(eligible.tolist(), key=_rank)
        if top_n is not None and len(selected) > top_n:
            if balance_groups and len(set(groups[selected])) > 1:
                names = sorted(set(groups[selected]))
                chosen = []
                for j, name in enumerate(names):
                    quota = top_n // len(names) + (j < top_n % len(names))
                    candidates = [i for i in selected if groups[i] == name]
                    chosen.extend(candidates[:quota])
                chosen_set = set(chosen)
                chosen.extend([i for i in selected if i not in chosen_set]
                              [:top_n - len(chosen)])
                selected = chosen
            else:
                selected = selected[:top_n]

    names = sorted(set(groups[selected]))
    indices, positions, ticklabels, headers = [], [], [], []
    y = 0.0
    
    for name in names:
        group = [i for i in selected if groups[i] == name]
        group.sort(key=(lambda i: requested_order[i])
                   if requested_order is not None else _rank)
        head_type = group_types[group[0]]
        headers.append((head_type + ' ' + name, y))
        y += 1.08
        for i in group:
            indices.append(i)
            positions.append(y)
            ticklabels.append(labels[i])
            y += 1.0
        y += 0.85

    fig_height = max(6.4, 2.0 + row_height * y)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    blues = plt.get_cmap('Blues')
    cmap = LinearSegmentedColormap.from_list(
        'Blues_readable', blues(np.linspace(0.30, 0.95, 256)))
    norm = Normalize(vmin=0.0 if float(min_persistence) == 0 else 10.0,
                     vmax=100.0)
    colors = [cmap(norm(persistence[i])) for i in indices]
    ax.barh(positions, mean[indices], xerr=sd[indices], height=0.76,
            color=colors, edgecolor='#284764', linewidth=0.8,
            error_kw={'elinewidth': 1.9, 'ecolor': '#20252a',
                      'capsize': 4.0, 'capthick': 1.7})
    ax.set_yticks(positions)
    ax.set_yticklabels(ticklabels, fontsize=label_fontsize)
    ax.set_xlabel('Mean minimum distance [Å]', fontsize=axis_fontsize,
                  labelpad=14)
    ax.tick_params(axis='x', labelsize=tick_fontsize, width=1.3, length=6)
    ax.tick_params(axis='y', length=0, pad=12)
    
    if xlim is None:
        ax.set_xlim(left=0)
    else:
        ax.set_xlim(0, xlim)
    
    ax.set_ylim(y - 0.55, -0.7)
    
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.grid(axis='x', alpha=0.13, linewidth=0.75)
    ax.set_axisbelow(True)
    
    for title, yy in headers:
        ax.text(0.002, yy, title, transform=ax.get_yaxis_transform(),
                fontsize=header_fontsize, fontweight='bold', va='center',
                ha='left', clip_on=False)
    
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                            pad=0.035, fraction=0.052, shrink=0.92, aspect=22)
    colorbar.set_label('Contact persistence [%]', fontsize=axis_fontsize,
                       labelpad=14)
    
    if float(min_persistence) == 0:
        colorbar.set_ticks([0, 10, 25, 50, 75, 100])
    else:
        colorbar.set_ticks([10, 25, 50, 75, 100])
    colorbar.ax.tick_params(labelsize=colorbar_fontsize, length=5)
    fig.subplots_adjust(left=0.19, right=0.83, top=0.982,
                        bottom=max(0.075, min(0.16, 1.1 / fig_height)))

    plotted = [groups[i] + ':' + labels[i] for i in indices]
    
    if filename is not None:
        prefix = os.fspath(filename)
        parent = os.path.dirname(os.path.abspath(prefix))
    
        if not os.path.isdir(parent):
            os.makedirs(parent)
        
        fig.savefig(prefix + '.' + 'png', dpi=png_dpi,
                        bbox_inches='tight', pad_inches=0.22)
        LOGGER.info('Saved distance plots as %s.png.' % prefix)
    
        if save_statistics:
            with open(prefix + '_plotted_residues.dat', 'w') as handle:
                handle.write('# Group Residue Mean_A SD_A Persistence_percent\n')
                for i in indices:
                    handle.write('%s %s %.5f %.5f %.3f\n' %
                                 (groups[i], labels[i], mean[i], sd[i], persistence[i]))
            LOGGER.info('Saved %s_plotted_residues.dat.' % prefix)
    
    if show is None:
        show = SETTINGS['auto_show']
    
    if show:
        plt.show()
    
    if filename is not None and not show:
        plt.close(fig)
    
    LOGGER.info('Plotted %d residues (cutoff %.2f A, minimum persistence %.1f%%).'
                % (len(plotted), cutoff, min_persistence))
    LOGGER.info('Residues on figure: %s' % ', '.join(plotted))
    
    return plotted

