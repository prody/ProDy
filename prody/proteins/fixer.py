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

from prody import LOGGER
from numbers import Integral, Number

__all__ = ['addMissingAtoms', 'fixStructuresMissingAtoms']

def addMissingAtoms(infile, method='openbabel', pH=7.0, outfile=None, **kwargs):
    """This function will add hydrogens to the protein and ligand structure using Openbabel [NO11]_
    or PDBFixer with OpenMM. 
    
    There are also options whether to *model_residues* (default False), *remove_heterogens* 
    (default False), *keep_waters* (default True), *overwrite* (default False), *keep_ids* (default True).
    
    :arg infile: PDB file name
    :type infile: str

    :arg method: Name of program which will be use to fix protein structure.
            Two alternative options are available: 'openbabel' and 'pdbfixer'.
            For either option additional software need to be installed:
            'openbabel': OpenBabel
            'pdbfixer': PDBFixer and OpenMM
            default is 'openbabel'
    :type method: str
    
    :arg pH: pH value applied only for PDBfixer.
    :type pH: int, float
    
    :arg model_residues: add all missing atoms from residues, applied only for PDBfixer.
                    default is False
    :type model_residues: bool
    
    :arg keep_ids: keep the original residue number, applied only for PDBfixer.
                    default is True
    :type keep_ids: bool

    Instalation of Openbabel:
    conda install -c conda-forge openbabel    

    Find more information here: https://anaconda.org/conda-forge/openbabel
                                https://github.com/openmm/pdbfixer
    Program will create new file in the same directory with 'addH_' prefix.
    
    .. [NO11] O'Boyle, N. M., Banck M., James C. A., Morley C., Vandermeersch T., Hutchison G. R. 
    Open Babel: An open chemical toolbox *Journal of cheminformatics* **2011** 3:1-14. """
    
    model_residues = kwargs.get("model_residues", False)
    remove_heterogens = kwargs.get("remove_heterogens", False)
    keep_water = kwargs.get("keep_water", True)
    overwrite = kwargs.get("overwrite", False)
    keep_ids = kwargs.get("keep_ids", True)

    import os

    if not isinstance(model_residues, bool):
        raise TypeError('model_residues should be True or False')

    if not isinstance(remove_heterogens, bool):
        raise TypeError('remove_heterogens should be True or False')    

    if not isinstance(keep_water, bool):
        raise TypeError('keep_water should be True or False')

    if not isinstance(keep_ids, bool):
        raise TypeError('keep_ids should be True or False')

    if not isinstance(overwrite, bool):
        raise TypeError('overwrite should be True or False')
    
    if not isinstance(infile, str):
        raise TypeError('infile should be a string pointing to a file')
    
    if not os.path.exists(infile):
        raise ValueError('infile {0} does not exist'.format(infile))
    
    if not isinstance(pH, Number):
        raise TypeError('pH should be a number')

    if outfile == None:
        outfile = os.path.join(os.path.split(infile)[0],
                               "addH_" + os.path.split(infile)[1])

    if os.path.exists(outfile) and not overwrite:
        LOGGER.warn('outfile {0} already exists, so returning it. \
Set overwrite=True to overwrite it'.format(outfile))
        return outfile
        
    if outfile == infile:
        raise ValueError('outfile cannot be the same as infile')

    if method == 'openbabel':
        if model_residues:
            LOGGER.warn("Openbabel cannot add missing residues, skipping this step")

        if infile.endswith('cif'):
            raise ValueError('Openbabel cannot handle cif files')

        try:
            from openbabel import openbabel
        except ImportError:
            raise ImportError("Install Openbabel to add hydrogens to the structure or use PDBFixer/OpenMM.")
    
        obconversion = openbabel.OBConversion()
        obconversion.SetInFormat("pdb")
        mol = openbabel.OBMol()
        obconversion.ReadFile(mol, infile)
        mol.AddHydrogens()
        obconversion.WriteFile(mol, outfile)
        LOGGER.info("Hydrogens were added to the structure. Structure {0} is saved in the local directry.".format(outfile))

            
    elif method == 'pdbfixer':
        try:
            from pdbfixer import PDBFixer
            try:
                from openmm.app import PDBFile
            except ImportError:
                from simtk.openmm.app import PDBFile
            
            fixer = PDBFixer(filename=infile)

            if model_residues:
                fixer.findMissingResidues()
            else:
                fixer.missingResidues = {}

            if remove_heterogens:
                fixer.removeHeterogens(keepWater=keep_water)

            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(pH)
            PDBFile.writeFile(fixer.topology, fixer.positions, open(outfile, 'w'), keepIds=keep_ids)
            LOGGER.info("Hydrogens were added to the structure. New structure is saved as {0}.".format(outfile))

        except ImportError:
            raise ImportError('Install PDBFixer and OpenMM in order to fix the protein structure.')

    else:
        raise TypeError('Method should be openbabel or pdbfixer')
    
    return outfile


def fixStructuresMissingAtoms(infiles, method='openbabel', pH=7.0, outfiles=None, **kwargs):
    """This function will add hydrogens to the protein and ligand structure from a set of files
    using Openbabel [NO11]_ or PDBFixer with OpenMM.
    
    There are also options whether to *model_residues* (default False), *remove_heterogens* 
    (default False) and *keep_waters* (default True).
    
    :arg infiles: a list of PDB file names
    :type infile: list

    :arg method: Name of program which will be use to fix protein structure.
            Two alternative options are available: 'openbabel' and 'pdbfixer'.
            For either option additional software need to be installed:
            'openbabel': OpenBabel
            'pdbfixer': PDBFixer and OpenMM
            default is 'openbabel'
    :type method: str

    :arg model_residues: add all missing atoms from residues, applied only for PDBfixer.
                    default is False
    :type model_residues: bool
    
    :arg keep_ids: keep the original residue number, applied only for PDBfixer.
                    default is True
    :type keep_ids: bool
    
    :arg pH: pH value applied only for PDBfixer.
    :type pH: int, float
    
    Instalation of Openbabel:
    conda install -c conda-forge openbabel    

    Find more information here: https://anaconda.org/conda-forge/openbabel
                                https://github.com/openmm/pdbfixer
    Program will create new file in the same directory with 'addH_' prefix.
    
    .. [NO11] O'Boyle, N. M., Banck M., James C. A., Morley C., Vandermeersch T., Hutchison G. R. 
    Open Babel: An open chemical toolbox *Journal of cheminformatics* **2011** 3:1-14. """

    if not isinstance(infiles, list):
        raise TypeError('infiles should be a list')
    
    if outfiles is None:
        outfiles = [None for infile in infiles]

    if not isinstance(outfiles, list):
        raise TypeError('outfiles should be None or a list')
    if len(outfiles) != len(infiles):
        raise ValueError('outfiles should have the same length as infiles')
    
    results = []
    for i, infile in enumerate(infiles):
        results.append(addMissingAtoms(infile, method, pH, 
                                       outfiles[i], **kwargs))
    return results
