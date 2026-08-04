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
           'calcSASA', 'calcVolume','compareInteractions', 'showInteractionsGraph',
           'showInteractionsHist', 'calcLigandInteractions', 'listLigandInteractions', 
           'showProteinInteractions', 'showLigandInteraction',
           'showProteinInteractions_VMD', 'showLigandInteraction_VMD', 
           'calcHydrophobicOverlapingAreas',
           'Interactions', 'InteractionsTrajectory', 'LigandInteractionsTrajectory',
           'calcSminaBindingAffinity', 'calcSminaPerAtomInteractions', 'calcSminaTermValues',
           'showSminaTermValues', 'showPairEnergy', 'checkNonstandardResidues',
           'saveInteractionsAsDummyAtoms', 'createFoldseekAlignment', 'runFoldseek', 'runDali', 
           'runBLAST', 'extractMultiModelPDB', 'calcSignatureInteractions',
           'alignTrajectory', 'calcRMSDfromReference', 'showRMSDfromReference',
           'calcPairwiseRMSD', 'showPairwiseRMSDHeatmap', 'showRMSDHistogram',
           'calcClusterPopulations', 'getCluster', 'getClusterMedoid',
           'calcClusterStatistics', 'calcAllClusterStatistics', 'showClusterStatisticsTable',
           'showClusterRMSDComparison', 'clusterHierarchical', 'showDendrogram',
           'clusterKMedoids', 'writeClusters']


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
                    sele2_filter = sele2_filter.select(sele2.select('all and not (resid '+str(i[0])+' and chain '+i[1]+')'))

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





def alignTrajectory(atoms, trajectory, align = 'protein and backbone', select = 'all'):
    """
    Aligns each trajectory frame to the reference structure and returns a tuple of the 
    reference coordinates and the aligned coordinates of the selected atoms.
    
    The trajectory frames are aligned to the reference structure using the atoms specified by 
    ``align``. After alignment, the coordinates of the atoms specified by ``select`` are 
    extracted for every frame.
    
    
    :arg atoms: reference structure used for the alignment.
    :type atoms: :class:`prody.Atomic`
    
    :arg trajectory: trajectory containing the coordinate sets to align
    :type trajectory: :class:`prody.Trajectory`
    
    :arg align: atom selection used to calculate the alignment transformation.
                Must be a valid ProDy selection string.
                Default is ``"protein and backbone"``
    :type align: str
    
    :arg select: atom selection whose coordinates are returned.
                 Must be a valid ProDy selection string.
                 Default is ``"all"``
    :type select: str
                
    :returns: a tuple containing:
        * ref_coords (numpy.ndarray): coordinates of the selected atoms in the reference structure.
        * aligned_coords (numpy.ndarray): aligned coordinates of the selected atoms for every
        trajectory frame with shape ``(n_frames, n_atoms, 3)``.
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`) 
    
    Example usage:
    >>> pdb = prody.parsePDB("structure.pdb")
    >>> dcd = prody.Trajectory("trajectory.dcd")
    >>> ref_coords, aligned_coords = prody.alignTrajectory(pdb, dcd, select='resname IOA')
    """
    
    if trajectory.numAtoms() != atoms.numAtoms():
        raise ValueError("Trajectory atoms count does not match structure atoms count.")
    
    # Save original coordinates to restore in the end
    orig_coords = atoms.getCoords().copy()
    
    atom_align = atoms.select(align)
    if atom_align is None:
        raise ValueError(f"No atoms match '{align}' in the structure.")
    ref_align = atom_align.copy()
     
    atom_select = atoms.select(select)
    if atom_select is None:
        raise ValueError(f"No atoms match '{select}' in the structure.")
    ref_coords = atom_select.getCoords().copy()
    
    trajectory.link(atoms) # linking trajectory to update coordinates frame-by-frame
    
    n_frames = trajectory.numFrames()
    n_atoms = atom_select.numAtoms()
    
    trajectory.reset()
    aligned_coords = np.zeros((n_frames, n_atoms, 3))
    
    try:
        for i, frame in enumerate(trajectory):
            trans = calcTransformation(atom_align, ref_align)
            trans.apply(atom_select)
            aligned_coords[i] = atom_select.getCoords()
    finally:
        atoms.setCoords(orig_coords)
        trajectory.reset()
    
    return ref_coords, aligned_coords


def calcRMSDfromReference(reference_coords, aligned_coords):
    """
    Calculates the RMSD of each aligned trajectory frame from the reference coordinates.
    Uses a vectorized approach for better efficiency.
    
    
    :arg reference_coords: reference coordinates with shape ``(n_atoms, 3)``.
    :type reference_coords: :class:`numpy.ndarray`
    
    :arg aligned_coords: aligned coordinates with shape ``(n_frames, n_atoms, 3)``.
                         Recommended to generate them using :func:`alignTrajectory`.
    :type aligned_coords: :class:`numpy.ndarray`
    
    :returns: an array containing the RMSD of each frame from the reference
    :rtype: :class:`numpy.ndarray`
    
    Example usage:
    >>> pdb = prody.parsePDB("structure.pdb")
    >>> dcd = prody.Trajectory("trajectory.dcd")
    >>> ref_coords, aligned_coords = prody.alignTrajectory(pdb, dcd, select='resname IOA')
    >>> rmsd = prody.calcRMSDfromReference(ref_coords, aligned_coords)
    """
    
    aligned_coords = np.asarray(aligned_coords)
    reference_coords = np.asarray(reference_coords)
    
    if reference_coords.ndim != 2:
        raise ValueError(f"reference_coords must be a 2D array of shape (n_atoms, 3), but got {reference_coords.shape}")
    
    if aligned_coords.ndim != 3:
        raise ValueError(f"aligned_coords must have shape (n_frames, n_atoms, 3), but got {aligned_coords.shape}.")

    if reference_coords.shape != aligned_coords.shape[1:]:
        raise ValueError(f"Incompatible shapes: reference is {reference_coords.shape}, but aligned frames have {aligned_coords.shape[1:]}.")
        
    # Vectorized RMSD calculation
    sd = (aligned_coords - reference_coords)**2
    msd = np.mean(sd.sum(axis=2), axis=1)
    rmsd = np.sqrt(msd)
        
    return rmsd


def showRMSDfromReference(rmsd_array, title = "Frame-to-Reference RMSD",
                          color = "#36454F", lw = 0.75, label = None, ax = None):
    """
    Plots the RMSD of an aligned trajectory from the reference on the current axis.
    
    
    :arg rmsd_array: one-dimensional array containing the RMSD values (in Å) of each 
                     trajectory frame relative to the reference structure.
                     This array is typically generated using :func:`calcRMSDfromReference`.
    :type rmsd_array: :class:`numpy.ndarray`
    
    :arg title: the title of the generated plot.
                Default is ``"Frame-to-Reference RMSD"``
    :type title: str
    
    :arg color: the color of the plot.
                Default is ``"#36454F"``
    :type color: str
    
    :arg lw: the width of the line.
             Default is ``0.75``
    :type lw: float
    
    :arg label: the label for the plot.
                Default is ``None``
    :type label: str
    
    :arg ax: axes on which to draw the plot.
             Default is ``None`` and the current axes are used.
    :type ax: :class:`matplotlib.axes.Axes`
    
    :returns: the Matplotlib axes containing the plot.
    :rtype: :class:`matplotlib.axes.Axes`, None
    
    Example usage:
    >>> import matplotlib.pyplot as plt
    >>> rmsd_run1 = prody.calcRMSDfromReference(ref_coords, aligned_coords1)
    >>> rmsd_run2 = prody.calcRMSDfromReference(ref_coords, aligned_coords2)
    >>> plt.figure(figsize=(8, 6))
    >>> prody.showRMSDfromReference(rmsd_run1, label="Run 1")
    >>> prody.showRMSDfromReference(rmsd_run2, color="red", label="Run 2")
    >>> plt.show()
    """
    
    import matplotlib.pyplot as plt
    
    
    if ax is None:
        ax = plt.gca()
    
    ax.plot(rmsd_array, color = color, lw = lw, label = label)
    ax.set_xlabel("# Frame")
    ax.set_ylabel("RMSD from Reference [Å]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    
    if label is not None:
        ax.legend()
        
    return ax


def calcPairwiseRMSD(aligned_coords):
    """
    Calculates the frame-to-frame pairwise RMSD matrix using aligned structures
    and returns a symmetric distance matrix.
    Uses a vectorized approach for better efficiency.
    
    
    :arg aligned_coords: aligned coordinates with shape ``(n_frames, n_atoms, 3)``.
                         Recommended to generate them using :func:`alignTrajectory`.
    :type aligned_coords: :class:`numpy.ndarray`
    
    :returns: symmetric matrix containing the pairwise RMSD between all frames with shape 
              ``(n_frames, n_frames)``.
    :rtype: :class:`numpy.ndarray`
    
    Example usage:
    >>> ref_coords, aligned_coords = prody.alignTrajectory(pdb, dcd, select='resname IOA')
    >>> distance_matrix = prody.calcPairwiseRMSD(aligned_coords)
    """
    
    from scipy.spatial.distance import cdist
    
    
    aligned_coords = np.asarray(aligned_coords)
    
    if aligned_coords.ndim != 3 or aligned_coords.shape[2] != 3:
        raise ValueError(f"aligned_coords must have shape (n_frames, n_atoms, 3), but got {aligned_coords.shape}.")
        
    if aligned_coords.shape[0] == 0:
        raise ValueError("aligned_coords contains no frames.")

    if aligned_coords.shape[1] == 0:
        raise ValueError("aligned_coords contains no atoms.")
    
    
    n_frames, n_atoms, _ = aligned_coords.shape
    
    # Flatten frame coordinates (n_frames, n_atoms, 3) -> (n_frames, n_atoms * 3)
    coords_flat = aligned_coords.reshape(n_frames, 3 * n_atoms)
    
    euclidean_dists = cdist(coords_flat, coords_flat, metric = 'euclidean')
    distance_matrix = euclidean_dists / np.sqrt(n_atoms)
    
    return distance_matrix


def showPairwiseRMSDHeatmap(distance_matrix, title = "Pairwise RMSD Distance Matrix",
                            label = "RMSD [Å]", cmap = "viridis", ax = None):
    """
    Plots the heatmap of the frame-to-frame pairwise RMSDs using the distance matrix
    on the current axis.
    
    
    :arg distance_matrix: two-dimensional array constituting the pairwise distance matrix.
                          Typically generated using :func:`calcPairwiseRMSD`.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg title: the title of the generated plot.
                Default is ``"Pairwise RMSD Distance Matrix"``
    :type title: str
    
    :arg label: the label for the heatmap legend.
                Default is ``"RMSD [Å]"``
    :type label: str
    
    :arg cmap: the colormap for the heatmap.
               Default is ``"viridis"``
    :type cmap: str
    
    :arg ax: axes on which to draw the plot.
             Default is ``None`` and the current axes are used.
    :type ax: :class:`matplotlib.axes.Axes`
    
    :returns: the Matplotlib axes containing the heatmap.
    :rtype: :class:`matplotlib.axes.Axes`
    
    Example usage:
    >>> import matplotlib.pyplot as plt
    >>> distance_matrix = prody.calcPairwiseRMSD(aligned_coords)
    >>> plt.figure(figsize=(8, 6))
    >>> prody.showPairwiseRMSDHeatmap(distance_matrix)
    >>> plt.show()
    """
    
    import matplotlib.pyplot as plt
    
    
    if distance_matrix.ndim != 2:
        raise ValueError(f"distance_matrix must be a 2D array, but got shape {distance_matrix.shape}.")
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise ValueError("distance_matrix must be square.")
    
    if ax is None:
        ax = plt.gca()

    im = ax.imshow(distance_matrix, cmap = cmap)
    ax.figure.colorbar(im, ax = ax, label = label)

    ax.set_xlabel("# Frame")
    ax.set_ylabel("# Frame")
    ax.set_title(title)

    return ax
    

def showRMSDHistogram(rmsd_data, bins = 50, title = 'Distribution of RMSDs',
                              xlabel = 'RMSD [Å]', ylabel = 'Frequency', label = None,
                              kde = False, element = 'bars', stat = 'count', 
                              color = 'teal', edgecolor = 'black', lw = 0.8, alpha = 0.5, ax = None):
    """
    Plots the distribution of RMSD values on the current axis. 
    The input may be either a one-dimensional RMSD array (e.g., frame-to-reference or intra-cluster RMSDs)
    or a 2D pairwise RMSD distance matrix.
    
    
    :arg rmsd_data: 1D array of RMSDs or 2D pairwise distance matrix. 
                    If 2D, the upper triangle (excluding diagonal) is automatically extracted.
                    Recommended to generate them using :func:`calcRMSDfromReference` or
                    :func:`calcPairwiseRMSD` respectively.
    :type rmsd_data: :class:`numpy.ndarray`
    
    :arg bins: number of histogram bins. 
               Default is ``50``
    :type bins: int
    
    :arg title: the title of the generated plot. 
                Default is ``'Distribution of RMSDs'``
    :type title: str
    
    :arg xlabel: the label for the X axis.
                 Default is ``'RMSD [Å]'``
    :type xlabel: str
    
    :arg ylabel: the label for the Y axis.
                 Default is ``'Frequency'``
    :type ylabel: str
    
    :arg label: the label for the plot legend. 
                Default is ``None``
    :type label: str
    
    :arg kde: whether to plot a kernel density estimate. 
              Default is ``False``
    :type kde: bool
    
    :arg element: visual representation of the histogram bins ('bars', 'step', or 'poly'). 
                  Default is ``'bars'``
    :type element: str
    
    :arg stat: aggregate statistic to compute in each bin ('count', 'frequency', 'probability', 'percent', 'density'). 
               Default is ``'count'``
    :type stat: str
    
    :arg color: the color of the plot. 
                Default is ``'teal'``
    :type color: str
    
    :arg edgecolor: the color of the bin edges.
                    Default is ``'black'``
    :type edgecolor: str
    
    :arg lw: the line width of the bin edges. 
             Default is ``0.8``
    :type lw: float
    
    :arg alpha: the transparency of the bins. 
                Default is ``0.5``
    :type alpha: float
    
    :arg ax: axes on which to draw the plot.
             Default is ``None`` and the current axes are used.
    :type ax: :class:`matplotlib.axes.Axes`
    
    :returns: the Matplotlib axes containing the plot.
    :rtype: :class:`matplotlib.axes.Axes`
    
    Example usage:
    >>> import matplotlib.pyplot as plt
    >>> distance_matrix = prody.calcPairwiseRMSD(aligned_coords)
    >>> plt.figure(figsize = (8, 6))
    >>> prody.showRMSDHistogram(distance_matrix, bins=70, kde=True)
    >>> plt.show()
    """
    
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("The 'seaborn' package is required to display the Histogram."
                          "\nPlease install it using 'pip install seaborn'."
                          "\nAlternatively, use standard matplotlib.pyplot.hist() for basic plots.")
    
    
    data_array = np.asarray(rmsd_data)
    
    if data_array.ndim == 2:
        if data_array.shape[0] != data_array.shape[1]:
            raise ValueError(f"Pairwise RMSD matrix must be square, but got shape {data_array.shape}")
        rmsd_values = data_array[np.triu_indices_from(data_array, k=1)]
    elif data_array.ndim == 1:
        rmsd_values = data_array
    else:
        raise ValueError(f"Expected 1D or 2D array, but got shape {data_array.shape}")
    
    if ax is None:
        ax = plt.gca()
    
    sns.histplot(rmsd_values, bins = bins, element = element, stat = stat, kde = kde, 
                 alpha = alpha, color = color, edgecolor = edgecolor, linewidth = lw, 
                 label = label, ax = ax)
    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if label is not None:
        ax.legend()
        
    ax.grid(axis = 'y', alpha = 0.3)

    return ax
    

def calcClusterPopulations(cluster_ids):
    """
    Uses the cluster IDs array, assigning each frame to a cluster, to calculate the population 
    of each cluster and its corresponding percentage.
    
    
    :arg cluster_ids: a one-dimensional array matching each frame to a cluster.
    :type cluster_ids: :class:`numpy.ndarray`
    
    :returns: a dictionary mapping each cluster ID to its population statistics.
              Each value is a dictionary with the keys ``"count"`` and ``"pct"``.
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = prody.clusterHierarchical(distance_matrix)
    >>> populations = prody.calcClusterPopulations(cluster_ids)
    >>> print(f"In cluster 1: {populations[1]['pct']:.2f}% of the objects")
    """
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}")
    if cluster_ids.size == 0:
        raise ValueError("cluster_ids is empty.")
    
    total_frames = len(cluster_ids)
    clusters, frequencies = np.unique(cluster_ids, return_counts=True)
    
    population_data = {}
    for c, n in zip(clusters, frequencies):
        population_data[int(c)] = {'count': int(n), 'pct': (n / total_frames) * 100}
        
    return population_data


def getCluster(cluster_ids, cluster_number):
    """
    Returns the indices of the frames grouped in a specific cluster.
    
    
    :arg cluster_ids: a one-dimensional array assigning each frame to a cluster.
    :type cluster_ids: :class:`numpy.ndarray`
    
    :arg cluster_number: the ID of the cluster whose members we want to return.
    :type cluster_number: int
    
    :returns: an array of the trajectory indices assigned to the specified cluster.
    :rtype: :class:`numpy.ndarray`
    
    Example usage:
    >>> cluster_ids, _ = prody.clusterHierarchical(distance_matrix)
    >>> cluster1_indices = prody.getCluster(cluster_ids, 1)
    """
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}.")
    if cluster_ids.size == 0:
        raise ValueError("cluster_ids is empty.")    
    
    if isinstance(cluster_number, bool) or not isinstance(cluster_number, (int, np.integer)):
        raise TypeError(f"cluster_number must be an integer, but got {type(cluster_number).__name__}.")
    
    cluster_indices = np.where(cluster_ids == cluster_number)[0]
    
    if cluster_indices.size == 0:
        raise ValueError(f'No frames belong to Cluster {cluster_number}')
        
    return cluster_indices


def getClusterMedoid(distance_matrix, cluster_indices = None, cluster_ids = None, cluster_number = None):
    """
    Determines the medoid (representative) of the cluster from the minimum total pairwise 
    distance to all other cluster members.
    Returns a dictionary containing the medoid index within the cluster ("local") and in 
    the full trajectory ("global").
    
    The cluster can be specified either by:
        1. cluster_indices (recommended and can be generated by :func:`getCluster`) or
        2. cluster_ids and cluster_number
     
        
    :arg distance_matrix: a two-dimensional square, symmetric pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_indices: a one-dimensional array of the indices of the cluster members
                          Default is ``None``
    :type cluster_indices: :class:`numpy.ndarray`
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
                      Default is ``None``
    :type cluster_ids: :class:`numpy.ndarray`
    
    :arg cluster_number: the ID (number) of the cluster we want to investigate
                         Default is ``None``
    :type cluster_number: int
    
    :returns: a dictionary of the medoid index, with keys "local" for the index within the
              cluster and "global" for the index within the full trajectory
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = prody.clusterHierarchical(distance_matrix)
    >>> cluster1_indices = prody.getCluster(cluster_ids, 1)
    >>> cluster1_medoid = prody.getClusterMedoid(distance_matrix, cluster1_indices)
    >>> print(f"Medoid of cluster 1 is # {cluster1_medoid['local']} in the cluster.")
    """
    
    cluster_indices = _resolveClusterIndices(cluster_indices=cluster_indices, 
                                             cluster_ids=cluster_ids, 
                                             cluster_number=cluster_number)
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)

    n_frames = distance_matrix.shape[0]
    if np.any((cluster_indices < 0) | (cluster_indices >= n_frames)):
        raise ValueError("cluster_indices contain indices outside the bounds of the distance matrix.")
    
    cluster_distance_matrix = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
    
    sum_dist = cluster_distance_matrix.sum(axis=1)
    medoid_local = int(np.argmin(sum_dist))
    medoid_global = int(cluster_indices[medoid_local])
    
    return {"local": medoid_local, "global": medoid_global}


def calcClusterStatistics(distance_matrix, cluster_indices = None, cluster_medoid = None,
                          cluster_ids = None, cluster_number = None):
    """
    Calculates descriptive statistics for a single cluster.
    Returns a dictionary with the statistics quantity as a key.
    
    The cluster can be specified either by:
        1. cluster_indices (recommended and can be generated by :func:`getCluster`) or
        2. cluster_ids and cluster_number
        
    
    :arg distance_matrix: a two-dimensional square, symmetric pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_indices: a one-dimensional array of the indices of the cluster members
                          Default is ``None``
    :type cluster_indices: :class:`numpy.ndarray`
    
    :arg cluster_medoid: the dictionary of the cluster medoid.
                         Expects the 'cluster_medoid' dictionary generated by the function
                         :func:`getClusterMedoid`.
                         Default is ``None``
    :type cluster_medoid: dict
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
                      Default is ``None``
    :type cluster_ids: :class:`numpy.ndarray`
    
    :arg cluster_number: the ID (number) of the cluster we want to investigate
                         Default is ``None``
    :type cluster_number: int
    
    :returns: a dictionary of statistics quantities with their names as keys
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = prody.clusterHierarchical(distance_matrix)
    >>> cluster1_indices = prody.getCluster(cluster_ids, 1)
    >>> cluster_stats1 = prody.calcClusterStatistics(distance_matrix, cluster1_indices)
    >>> print(f"Mean distance from cluster Medoid: {cluster_stats1['mean']:.2f} ± {cluster_stats1['std']:.2f} [Å]")
    """   

    cluster_indices = _resolveClusterIndices(cluster_indices, cluster_ids, cluster_number)
    population = len(cluster_indices)
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    if cluster_medoid is None:
        cluster_medoid = getClusterMedoid(distance_matrix, cluster_indices=cluster_indices)
    else:
        if not isinstance(cluster_medoid, dict):
            raise TypeError("cluster_medoid must be a dictionary with keys 'global' and 'local'.")
    
        required_keys = {"global", "local"}
        if not required_keys.issubset(cluster_medoid.keys()):
            missing = required_keys - cluster_medoid.keys()
            raise ValueError(f"Expected keys are 'global' and 'local'.\ncluster_medoid is missing required key(s): {missing}.")   
        
        if cluster_medoid["global"] not in cluster_indices:
            raise ValueError("Index not found.\ncluster_medoid does not belong to the specified cluster.")
            
        if not (0 <= cluster_medoid["local"] < population):
            raise ValueError("Local index is out of bounds.")
            
        if cluster_indices[cluster_medoid["local"]] != cluster_medoid["global"]:
            raise ValueError("cluster_medoid['local'] and cluster_medoid['global'] are inconsistent.")
        
    distance_to_medoid = distance_matrix[cluster_indices, cluster_medoid["global"]]
    total_frames = distance_matrix.shape[0]
    
    stats = {
        "population"    : population,
        "pct"           : (population / total_frames) * 100,
        "medoid_global" : cluster_medoid["global"],
        "medoid_local"  : cluster_medoid["local"],
        "distances"     : distance_to_medoid,
        "mean"          : np.mean(distance_to_medoid),
        "std"           : np.std(distance_to_medoid),
        "median"        : np.median(distance_to_medoid),
        "iqr"           : np.percentile(distance_to_medoid, 75) - np.percentile(distance_to_medoid, 25),
        "p95"           : np.percentile(distance_to_medoid, 95),
        "max"           : np.max(distance_to_medoid)
    }
    
    if cluster_number is not None:
        stats["cluster"] = cluster_number
        
    return stats


def calcAllClusterStatistics(distance_matrix, cluster_ids):
    """
    Calculates descriptive statistics for all clusters.
    Returns a list of dictionaries, each corresponding to a cluster with the statistic quantities
    as keys.
    
    
    :arg distance_matrix: a two-dimensional square, symmetric pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
    :type cluster_ids: :class:`numpy.ndarray`
    
    :returns: a list of dictionaries, one dictionary for each cluster with the statistic 
              quantities as keys.
    :rtype: list of dict
    
    Example usage:
    >>> cluster_ids, _= prody.clusterHierarchical(distance_matrix, 4)
    >>> all_stats = prody.calcAllClusterStatistics(distance_matrix, cluster_ids)
    """
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}.")
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    if len(cluster_ids) != distance_matrix.shape[0]:
        raise ValueError(f"Dimension mismatch: cluster_ids length ({len(cluster_ids)}) "
                         f"does not match distance_matrix frames ({distance_matrix.shape[0]}).")
    
    clusters = np.unique(cluster_ids)
    all_stats = []
    
    for c in clusters:
        cluster_stats = calcClusterStatistics(distance_matrix, cluster_ids = cluster_ids, cluster_number = int(c))
        all_stats.append(cluster_stats) 
        
    return all_stats


def _resolveClusterIndices(cluster_indices, cluster_ids, cluster_number):
    """
    Checks whether the input is cluster_indices or cluster_ids and cluster_number, and for the latter returns
    the cluster indices.
    """
   
    has_indices = cluster_indices is not None
    has_ids_info = (cluster_ids is not None) and (cluster_number is not None)
    
    if has_indices == has_ids_info:
        raise ValueError("Provide exactly one of 'cluster_indices' or both 'cluster_ids' and 'cluster_number'.")
    
    if has_indices:
        cluster_indices = np.asarray(cluster_indices)
        
        if cluster_indices.ndim != 1:
            raise ValueError(f"cluster_indices must be a 1D array, but shape {cluster_indices.shape} was given.")
        if cluster_indices.size == 0:
            raise ValueError("cluster_indices must not be empty.")
    else:
        cluster_ids = np.asarray(cluster_ids)
        cluster_indices = getCluster(cluster_ids, cluster_number)  
        
    return cluster_indices


def _validateDistanceMatrix(distance_matrix):
    """
    Validates that a distance matrix is a non-empty square 2D NumPy array.
    """
    
    distance_matrix = np.asarray(distance_matrix)
    if distance_matrix.ndim != 2:
        raise ValueError(f"The distance matrix must be a 2D array, but got shape {distance_matrix.shape}.")
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise ValueError(f"distance_matrix must be square, but got shape {distance_matrix.shape}.")
    if distance_matrix.size == 0:
        raise ValueError("distance_matrix cannot be empty.")
        
    return distance_matrix


def showClusterStatisticsTable(all_stats, dissimilarity = "RMSD", units = "Å", show = True):
    """
    Prints the cluster statistics in a table, where each column corresponds to a cluster.
    Expects either the 'stats' dictionary generated by :func:`calcClusterStatistics` or
    the 'all_stats' list of dictionaries generated by :func:`calcAllClusterStatistics`.
    
    
    :arg all_stats: a dictionary, or list of dictionaries, with the cluster descriptive statistics
                    Recommended to generate the dictionaries from the functions :func:`calcClusterStatistics`
                    of a single cluster or :func:`calcAllClusterStatistics`.
    :type all_stats: dict, or list of dict
                    
    :arg dissimilarity: the dissimilarity measure
                        Default is `"RMSD"`
    :type dissimilarity: str
    
    :arg units: the units of the dissimilarity measure
                Default is `"Å"`
    :type units: str
    
    :arg show: whether to print the table
               Default is `True`
    :type show: bool
    
    :returns: a dictionary containing the table, the row labels and the headers, with their names as keys
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = prody.clusterHierarchical(distance_matrix)
    >>> all_stats = prody.calcAllClusterStatistics(distance_matrix, cluster_ids)
    >>> prody.showClusterStatisticsTable(all_stats);
    """
    
    try:
        from tabulate import tabulate
    except ImportError:
        raise ImportError("The 'tabulate' package is required to display the table. "
                          "Please install it using 'pip install tabulate'.")
    
    
    if isinstance(all_stats, dict):
        all_stats = [all_stats]
    elif not isinstance(all_stats, list):
        raise TypeError(f"all_stats must be a dict or list of dicts, but got {type(all_stats).__name__}.")
        
    if not all_stats:
        raise ValueError("all_stats cannot be empty.")
    
    required = {"population", "pct", "medoid_global", "medoid_local", "mean", "std",
                "median", "iqr", "p95", "max"}
    
    headers = []
    for i, stat in enumerate(all_stats):
        if not isinstance(stat, dict):
            raise TypeError(f"Expected a dictionary in all_stats, but got {type(stat).__name__} at index {i}.")
        
        cluster_id = stat.get('cluster', i + 1)
        headers.append(f"Cluster {int(cluster_id)}")
        
        missing = required - stat.keys()
        if missing:
            raise ValueError(f"Statistics dictionary is missing required keys: {sorted(missing)}.")
    
    metrics = [
        ("population",    "Total Frames"),
        ("pct",           "Population Percentage (%)"),
        ("medoid_global", "Medoid Frame (Global)"),
        ("medoid_local",  "Medoid Frame (Within Cluster)"),
        ("mean",          f"Mean {dissimilarity} [{units}]"),
        ("std",           f"Std {dissimilarity} [{units}]"),
        ("median",        f"Median {dissimilarity} [{units}]"),
        ("iqr",           f"IQR [{units}]"),
        ("p95",           f"95th Percentile [{units}]"),
        ("max",           f"Max {dissimilarity} [{units}]")
    ]
    
    table = []
    row_labels = []

    for key, label in metrics:
        row_labels.append(label)
        table.append([stat[key] for stat in all_stats])
        
    if not isinstance(show, bool):
        raise TypeError("show must be a bool.")
    
    if show:
        print(tabulate(table, headers=headers, showindex=row_labels, 
                       tablefmt="fancy_grid", floatfmt=".4f", stralign='center'))
    
    return {"table": table, "row_labels": row_labels, "headers": headers}


def showClusterRMSDComparison(all_stats, bins = 50, title = 'RMSD Distributions',
                              xlabel = 'RMSD to Medoid [Å]', ylabel = 'Frequency', 
                              kde = True, element = 'step', stat = 'count', 
                              alpha = 0.5, lw = 1.5, ax = None):
    """
    Overlays the internal RMSD distributions of all clusters onto a single plot 
    on the current axis.
    
    Expects the list of cluster statistics dictionaries generated by :func:`calcAllClusterStatistics`.
    
    
    :arg all_stats: list of cluster statistics dictionaries. Each dictionary
                    must contain the keys ``"cluster"`` and ``"distances"``.
    :type all_stats: list of dict
    
    :arg bins: number of histogram bins.
               Default is ``50``
    :type bins: int
    
    :arg title: the title of the generated plot.
                Default is ``'Intra-Cluster RMSD Distributions'``
    :type title: str
    
    :arg xlabel: the label for the X axis.
                 Default is ``'RMSD to Medoid [Å]'``
    :type xlabel: str
    
    :arg ylabel: the label for the Y axis.
                 Default is ``'Frequency'``
    :type ylabel: str
    
    :arg kde: whether to plot kernel density estimates.
              Default is ``True``
    :type kde: bool
    
    :arg element: visual representation of the histogram bins ('step', 'bars', or 'poly').
                  Default is ``'step'``
    :type element: str
    
    :arg stat: aggregate statistic to compute in each bin.
               Default is ``'count'``
    :type stat: str
    
    :arg alpha: the transparency of the bins.
                Default is ``0.5``
    :type alpha: float
    
    :arg lw: line width of the distribution lines.
             Default is ``1.5``
    :type lw: float
    
    :arg ax: axes on which to draw the plot.
             Default is ``None`` and the current axes are used.
    :type ax: :class:`matplotlib.axes.Axes`
    
    :returns: the Matplotlib axes containing the plot.
    :rtype: :class:`matplotlib.axes.Axes`
    
    Example usage:
    >>> import matplotlib.pyplot as plt
    >>> all_stats = prody.calcAllClusterStatistics(distance_matrix, clusterIDs_array)
    >>> plt.figure(figsize = (8, 6))
    >>> prody.showClusterRMSDComparison(all_stats)
    >>> plt.show()
    """
    
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("The 'seaborn' package is required to display the Histogram."
                          "\nPlease install it using 'pip install seaborn'."
                          "\nAlternatively, use standard matplotlib.pyplot.hist() for basic plots.")
    
    
    if len(all_stats) == 0:
        raise ValueError("all_stats cannot be empty.")
    
    if ax is None:
        ax = plt.gca()
    
    for cluster_stats in all_stats:
        sns.histplot(cluster_stats["distances"], bins = bins, element = element, stat = stat, kde = kde,
                     alpha = alpha, linewidth = lw, label = f"Cluster {cluster_stats['cluster']}", 
                     ax = ax)
        
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.legend(title = "Clusters")
    ax.grid(axis = 'y', alpha = 0.3)

    return ax


def clusterHierarchical(distance_matrix, method='average', cutoff=None):
    """
    Performs bottom-up hierarchical clustering from a pairwise distance matrix.
    
    Note that the resulting cluster IDs are 1-based (starting at 1, not 0).
    
    If cutoff is ``None``, it is automatically chosen as the midpoint of the largest
    gap between consecutive linkage distances.
    
    Recommendation: Plot the dendrogram using showDendrogram() and choose the cutoff 
                    manually whenever possible.
    
        
    :arg distance_matrix: either a one-dimensional condensed distance matrix or 
                          a two-dimensional pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg method: linkage criterion used when constructing the linkage matrix.
                 Default is ``'average'``
    :type method: str
    
    :arg cutoff: the cutoff distance that determines the number of clusters.
                 Default is ``None`` and the midpoint of the largest linkage gap
                 is used automatically.
    :type cutoff: float
    
    :returns:
        * a one-dimensional array containing the cluster ID for each object.
          IDs are 1-indexed (1 to number of clusters for the cutoff)
        * the hierarchical clustering linkage matrix
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
    
    Example usage:
    >>> distance_matrix = prody.calcPairwiseRMSD(aligned_coords)
    >>> cluster_ids, linkage_matrix = prody.clusterHierarchical(distance_matrix)
    """

    from scipy.cluster.hierarchy import linkage, fcluster


    condensed_distance_matrix = _condenseDistanceMatrix(distance_matrix)    
    linkage_matrix = linkage(condensed_distance_matrix, method=method)
    
    if cutoff is None:
        cutoff = _calcAutoCutoff(linkage_matrix)
    elif isinstance(cutoff, str):
        raise ValueError(f"Invalid input for cutoff: '{cutoff}'. Use numeric value or None.")
        
    cluster_ids = fcluster(linkage_matrix, t=cutoff, criterion='distance')
    
    return cluster_ids, linkage_matrix


def showDendrogram(distance_matrix = None, linkage_matrix = None, method = 'average', cutoff = None,
                   truncate_mode = None, p = 30, title = "Hierarchical Clustering Dendrogram", 
                   ylabel = "RMSD [Å]", ax = None):
    """
    Plots a hierarchical clustering dendrogram using either a pairwise RMSD distance matrix 
    or a pre-computed linkage matrix.

    By default, the complete dendrogram is shown. To display a truncated dendrogram, 
    use ``truncate_mode`` together with ``p``.

    Exactly one of ``distance_matrix`` or ``linkage_matrix`` must be provided.
    
    
    :arg distance_matrix: one-dimensional condensed or two-dimensional square pairwise 
                          distance matrix. Typically generated using :func:`calcPairwiseRMSD`.
                          Default is ``None``
    :type distance_matrix: :class:`numpy.ndarray`

    :arg linkage_matrix: pre-computed hierarchical clustering linkage matrix.
                         Default is ``None``
    :type linkage_matrix: :class:`numpy.ndarray`

    :arg method: linkage criterion used when constructing the linkage matrix from 
                 ``distance_matrix``.
                 Default is ``"average"``
    :type method: str

    :arg cutoff: dendrogram color threshold. If ``"auto"``, the cutoff is placed at the midpoint
                 of the largest linkage gap.
                 Default is ``None``
    :type cutoff: float, str

    :arg truncate_mode: dendrogram truncation mode passed to :func:`scipy.cluster.hierarchy.dendrogram`.
                        Default is ``None``
    :type truncate_mode: str, None

    :arg p: truncation parameter used together with ``truncate_mode``.
            Default is ``30``
    :type p: int

    :arg title: title of the generated plot.
                Default is ``"Hierarchical Clustering Dendrogram"``
    :type title: str
    
    :arg ylabel: Distance matric label for the y-axis
                 Default is ``"RMSD [Å]"``
    :type ylabel: str

    :arg ax: axes on which to draw the plot.
             Default is ``None`` and the current axes are used.
    :type ax: :class:`matplotlib.axes.Axes`

    :returns: the matplotlib axes and the linkage matrix.
    :rtype: tuple(:class:`matplotlib.axes.Axes`, :class:`numpy.ndarray`)

    Example usage:
    >>> import matplotlib.pyplot as plt
    >>> distance_matrix = prody.calcPairwiseRMSD(aligned_coords)
    >>> plt.figure(figsize=(12, 8))
    >>> ax, linkage = prody.showDendrogram(distance_matrix, cutoff='auto', truncate_mode='lastp', p=30)
    >>> plt.show()
    """
    
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram, linkage


    if (distance_matrix is None) == (linkage_matrix is None):
        raise ValueError("Provide exactly one of 'distance_matrix' or 'linkage_matrix'.")

    if linkage_matrix is None:
        condensed_distance_matrix = _condenseDistanceMatrix(distance_matrix)
        linkage_matrix = linkage(condensed_distance_matrix, method=method)
    else:
        linkage_matrix = np.asarray(linkage_matrix)

        if linkage_matrix.ndim != 2 or linkage_matrix.shape[1] != 4 or linkage_matrix.shape[0] == 0:
            raise ValueError("linkage_matrix must be a non-empty array with shape (n_samples-1, 4).")

    if cutoff == "auto":
        cutoff = _calcAutoCutoff(linkage_matrix)
    elif isinstance(cutoff, str):
        raise ValueError(f"Invalid string for cutoff: '{cutoff}'. Use 'auto' or a numeric value.")

    if ax is None:
        ax = plt.gca()

    with plt.rc_context({"lines.linewidth": 0.6}):
        dendrogram(linkage_matrix, truncate_mode = truncate_mode, p = p, color_threshold = cutoff,
                   no_labels = True, ax = ax)

    ax.set_xlabel("# Frame")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if cutoff is not None:
        ax.axhline(cutoff, color = "black", linestyle = "--", linewidth = 1.2, label = f"Cutoff ({cutoff:.2f} Å)")
        ax.legend()

    return ax, linkage_matrix


def _calcAutoCutoff(linkage_matrix):
    """
    Calculates the optimal distance cutoff based on the midpoint of the largest gap
    in the linkage matrix.
    """
    
    linkage_matrix = np.asarray(linkage_matrix)
    merge_dist = linkage_matrix[:, 2] # locations of the merges
    
    if len(merge_dist) == 0:
        raise ValueError("At least two frames are required for clustering")
    elif len(merge_dist) == 1:
        return merge_dist[0] + 1.0

    gaps = np.diff(merge_dist)
    largest_gap_idx = np.argmax(gaps)
    
    return (merge_dist[largest_gap_idx] + merge_dist[largest_gap_idx + 1]) / 2.0


def _condenseDistanceMatrix(distance_matrix):
    """
    Converts a pairwise distance matrix into condensed form.

    Accepts either:
        * a condensed one-dimensional distance matrix, or
        * a square two-dimensional distance matrix.
    
    Returns the condensed distance matrix suitable for scipy.cluster.hierarchy.linkage.
    """
    from scipy.spatial.distance import squareform
    
    
    distance_matrix = np.asarray(distance_matrix)
    
    if distance_matrix.size == 0:
        raise ValueError("distance_matrix cannot be empty.")
    
    if distance_matrix.ndim == 2:
        if distance_matrix.shape[0] != distance_matrix.shape[1]:
            raise ValueError(f"distance_matrix must be square, but got shape {distance_matrix.shape}.")
        condensed_distance_matrix = squareform(distance_matrix)
    elif distance_matrix.ndim == 1:
        condensed_distance_matrix = distance_matrix
    else:
        raise ValueError("distance_matrix must be either a condensed 1D array or a square 2D array.")
        
    return condensed_distance_matrix


def clusterKMedoids(distance_matrix, k, method = 'alternate', method_sklearn = 'pam', 
                    initial_medoids = None, seed = None, max_iter = 100, n_init = 10):
    """
    Performs K-Medoids clustering using various algorithms.

    This function acts as a facade, routing the clustering task to the specified backend 
    ('alternate', 'pam', or 'sklearn'). 'alternate' is generally faster, while 'pam' typically
    finds solutions with lower cost. 'sklearn' needs the sklearn_extra.cluster module to be installed.
    
    Note that the resulting cluster IDs are 1-based (starting at 1, not 0).


    :arg distance_matrix: square, symmetric matrix of pairwise distances between all objects
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg k: prespecified number of clusters to form
    :type k: int
    
    :arg method: the clustering algorithm to use.
                 Options are 'alternate' (custom Alternating Medoids), 'pam' (custom PAM), 'sklearn' (sklearn_extra)
                 Default is 'alternate'
    :type method: str
    
    :arg method_sklearn: the specific method to pass to sklearn_extra if method = 'sklearn' is chosen.
                         Options are 'pam', 'alternate'
                         Default is 'pam'
    :type method_sklearn: str
    
    :arg initial_medoids: one-dimensional array of indices to use as starting medoids.
                          Default is ``None`` and the starting medoids are picked randomly.
                          It is not supported for the 'sklearn' method.
    :type initial_medoids: :class:`numpy.ndarray`
    
    :arg seed: random seed for reproducibility
               Default is ``None``
    :type seed: int
    
    :arg max_iter: maximum number of iterations per a single run
                   Default is ``100``
    :type max_iter: int
    
    :arg n_init: number of times the algorithm will be run with different initial medoids
                 and the best result with the lowest cost is returned
    
    :returns: a tuple of:
        * an one-dimensional array containing the cluster ID for each object.
          IDs are 1-indexed (1 to k)
        * a one-dimensional array of shape (k,) containing the indices of the final cluster
          medoids
        * the final sum of distances from each point to its nearest medoid
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, float)
    
    Example usage:
    >>> distance_matrix = prody.calcPairwiseRMSD(aligned_coords)
    >>> cluster_ids, medoids, _ = prody.clusterKMedoids(distance_matrix, 4, method = 'alternate', 
                                                        seed = 42, n_init = 30)
    """
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    n = distance_matrix.shape[0]
    
    if not isinstance(k, (int, np.integer)):
        raise TypeError(f"k must be an integer, got {type(k).__name__}")
    if k <= 0 or k > n:
        raise ValueError(f"k must be between 1 and {n} (the number of points).")
        
    if not isinstance(max_iter, (int, np.integer)) or max_iter <= 0:
        raise ValueError("max_iter must be a positive integer.")

    if not isinstance(n_init, (int, np.integer)) or n_init <= 0:
        raise ValueError("n_init must be a positive integer.")
        
    if seed is not None and not isinstance(seed, (int, np.integer)):
        raise TypeError("seed must be an integer or None.")
        
    if initial_medoids is not None:
        if method == 'sklearn':
            raise ValueError("The 'sklearn' method does not support custom 'initial_medoids'. Use 'alternate' or 'pam' instead.")
                             
        initial_medoids = np.asarray(initial_medoids, dtype=int)
        if initial_medoids.ndim != 1:
            raise ValueError("initial_medoids must be a 1D array.")
        if len(initial_medoids) != k:
            raise ValueError(f"Number of initial medoids ({len(initial_medoids)}) must equal k ({k}).")
        if len(np.unique(initial_medoids)) != k:
            raise ValueError("initial_medoids must be unique.")
        if np.any((initial_medoids < 0) | (initial_medoids >= n)):
            raise ValueError(f"initial_medoids contain indices out of bounds (must be 0 to {n-1}).")
            
        n_init = 1

    if method == 'alternate':
        return _clusterKMedoidsAlternating(distance_matrix, k, initial_medoids, seed, max_iter, n_init)
    elif method == 'pam':
        return _clusterKMedoidsPAM(distance_matrix, k, initial_medoids, seed, max_iter, n_init)
    elif method == 'sklearn':
        if method_sklearn not in ['alternate', 'pam']:
            raise ValueError(f"method_sklearn '{method_sklearn}' is not valid. Options: 'alternate', 'pam'")
        return _clusterKMedoidsSklearn(distance_matrix, k, method_sklearn, seed, max_iter, n_init)
    else:
        raise ValueError(f"Method '{method}' is not valid. Options: 'alternate', 'pam', 'sklearn'")


def _cost(distance_matrix, medoids):
    """Calculates the sum of distances from all points to their nearest medoid."""
    distances = distance_matrix[:, medoids] 
    closest_distance = np.min(distances, axis=1) 
    return np.sum(closest_distance)


def _clusterKMedoidsAlternating(distance_matrix, k, initial_medoids, seed, max_iter, n_init):
    """K-Medoids clustering with an Alternating Medoids algorithm"""
    
    n = distance_matrix.shape[0]
    rng = np.random.default_rng(seed)

    best_cost = np.inf
    best_medoids = None
    
    for run in range(n_init): 
        if initial_medoids is None:
            medoids = np.asarray(rng.choice(n, size = k, replace = False))
        else:
            medoids = initial_medoids.copy()

        for iteration in range(max_iter):
            distances_to_medoids = distance_matrix[:, medoids] 
            cluster_ids = np.argmin(distances_to_medoids, axis = 1) 

            new_medoids = np.zeros_like(medoids)
            
            for i in range(len(medoids)):
                cluster_members = np.where(cluster_ids == i)[0]
                
                # Handles empty clusters by retaining the old medoid
                if len(cluster_members) == 0:
                    new_medoids[i] = medoids[i]
                    continue
                
                intra_cluster_distances = distance_matrix[np.ix_(cluster_members, cluster_members)]
                sum_distances = intra_cluster_distances.sum(axis=1)
                best_medoid_idx_in_cluster = np.argmin(sum_distances)
                new_medoids[i] = cluster_members[best_medoid_idx_in_cluster]
            
            current_cost = _cost(distance_matrix, new_medoids)

            if np.array_equal(medoids, new_medoids):
                break
                
            medoids = new_medoids
    
        if current_cost < best_cost:
            best_cost = current_cost
            best_medoids = medoids.copy()

    distances = distance_matrix[:, best_medoids]
    cluster_ids = np.argmin(distances, axis=1) + 1 

    return cluster_ids, best_medoids, best_cost


def _clusterKMedoidsPAM(distance_matrix, k, initial_medoids, seed, max_iter, n_init):
    """K-Medoids clustering with PAM (Partitioning Around Medoids)"""
    
    n = distance_matrix.shape[0]
    rng = np.random.default_rng(seed)

    best_cost = np.inf
    best_medoids = None

    for run in range(n_init): 
        if initial_medoids is None:
            medoids = np.asarray(rng.choice(n, size=k, replace=False))
        else:
            medoids = initial_medoids.copy()
            
        # Using a set for removal/addition operations
        non_medoids = set(i for i in range(n) if i not in medoids)
        current_cost = _cost(distance_matrix, medoids)

        for iteration in range(max_iter):
            best_cost_swap = current_cost
            best_swap = None
            
            for medoid_idx, old_medoid in enumerate(medoids):
                for new_medoid in non_medoids:
                    candidate = medoids.copy()
                    candidate[medoid_idx] = new_medoid
                    candidate_cost = _cost(distance_matrix, candidate)

                    if candidate_cost < best_cost_swap:
                        best_cost_swap = candidate_cost
                        best_swap = (medoid_idx, old_medoid, new_medoid)

            if best_swap is None:
                break

            idx, old_medoid, new_medoid = best_swap
            medoids[idx] = new_medoid
            non_medoids.remove(new_medoid)
            non_medoids.add(old_medoid)
            current_cost = best_cost_swap

        if current_cost < best_cost:
            best_cost = current_cost
            best_medoids = medoids.copy()

    distances = distance_matrix[:, best_medoids]
    cluster_ids = np.argmin(distances, axis=1) + 1 

    return cluster_ids, best_medoids, best_cost


def _clusterKMedoidsSklearn(distance_matrix, k, method_sklearn, seed, max_iter, n_init):
    """K-Medoids clustering with sklearn_extra"""
    try:
        from sklearn_extra.cluster import KMedoids
    except ImportError:
        raise ImportError("The 'sklearn_extra' package is required for this K-Medoids approach. "
                          "Please install it using 'pip install scikit-learn-extra'.")
    
    
    best_cost = np.inf
    best_labels = None
    best_medoids = None

    rng = np.random.default_rng(seed)

    for run in range(n_init): 
        kmedoids = KMedoids(n_clusters=k, metric='precomputed', method=method_sklearn, init='random', 
                            max_iter=max_iter, random_state=int(rng.integers(0, 1000000)))

        kmedoids.fit(distance_matrix)
        medoids = kmedoids.medoid_indices_

        cost = np.sum(np.min(distance_matrix[:, medoids], axis=1))

        if cost < best_cost:
            best_cost = cost
            best_labels = kmedoids.labels_ + 1 
            best_medoids = medoids.copy()

    return best_labels, best_medoids, best_cost


def writeClusters(atoms, trajectory, distance_matrix, cluster_ids, write_dcd = True,
                  align ="protein and backbone", system = "system", tag ="cluster"):
    """
    Aligns the trajectory once, then exports representative medoid structures as PDB files and 
    cluster-specific DCD trajectories.
    Returns the name of all exported .pdb and .dcd files.
    
    Note: This function loads all aligned coordinates into memory. It is highly optimized 
    for speed, provided the trajectory fits within available system RAM.
    
    :arg atoms: reference structure used for the alignment.
    :type atoms: :class:`prody.Atomic`
    
    :arg trajectory: trajectory containing the coordinate sets to align
    :type trajectory: :class:`prody.Trajectory`
    
    :arg distance_matrix: square, symmetric matrix of pairwise distances between all objects
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
    :type cluster_ids: :class:`numpy.ndarray`   
    
    :arg write_dcd: determines whether to save a .dcd file of the cluster frames
                    Default is ``True``
    :type write_dcd: bool
    
    :arg align: atom selection used to calculate the alignment transformation.
                Must be a valid ProDy selection string.
                Default is ``"protein and backbone"``
    :type align: str
    
    :arg system: the name of the system under investigation
                 Default is ``"system"``
    :type system: str
    
    :arg tag: the name of the clustering method
              Default is ``"cluster"``
    :type tag: str
    
    :returns: list of exported filenames
    :rtype: list[str]
    
    Example usage:
    >>> pdb = prody.parsePDB("structure.pdb")
    >>> dcd = prody.Trajectory("trajectory.dcd")
    >>> _, aligned_coords = prody.alignTrajectory(pdb, dcd, select='resname IOA')
    >>> distance_matrix = prody.calcPairwiseRMSD(aligned_coords)
    >>> cluster_ids, _ = prody.clusterHierarchical(distance_matrix)
    >>> prody. writeClusters(pdb, dcd, distance_matrix, cluster_ids,
                      system="type1_RUN23", tag="hier")
    """
    
    from prody import writeDCD
    
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}.")
    if cluster_ids.size == 0:
        raise ValueError("cluster_ids is empty.")
    if len(cluster_ids) != trajectory.numFrames():
        raise ValueError("cluster_ids must have one entry per trajectory frame.")
        
    if not isinstance(write_dcd, bool):
        raise TypeError(f"write_dcd must be a bool, but got {type(write_dcd).__name__}")

    clusters = np.unique(cluster_ids)
    clusters = clusters[clusters > 0]
    num_clusters = clusters.size
    if num_clusters == 0:
        raise ValueError("No clusters were found. All frames are labeled as noise")
    
    exported_files = []
    
    # NOTE: Loads the entire aligned trajectory into memory
    _, aligned_coords = alignTrajectory(atoms, trajectory, align=align, select="all")
    
    for cluster in clusters:
        cluster_indices = getCluster(cluster_ids, cluster)
        medoid = getClusterMedoid(distance_matrix, cluster_indices = cluster_indices)

        cluster_coords = aligned_coords[cluster_indices]

        cluster_atoms = atoms.copy()
        cluster_atoms.setCoords(cluster_coords[0])

        if len(cluster_coords) > 1:
            cluster_atoms.addCoordset(cluster_coords[1:])

        if write_dcd:
            dcd_filename = f"{system}_{tag}_n{num_clusters}_cluster{cluster}.dcd"
            writeDCD(dcd_filename, cluster_atoms)
            exported_files.append(dcd_filename)

        medoid_atoms = atoms.copy()
        medoid_atoms.setCoords(aligned_coords[medoid["global"]])

        pdb_filename = f"{system}_{tag}_n{num_clusters}_cluster{cluster}_medoid.pdb"
        writePDB(pdb_filename, medoid_atoms)
        exported_files.append(pdb_filename)

    return exported_files





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

