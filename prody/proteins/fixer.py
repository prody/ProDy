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

__all__ = ['addMissingAtoms']

def addMissingAtoms(infile, method='openbabel', pH=7.0, outfile=None, **kwargs):
    """Function will add hydrogens to the protein and ligand structure using Openbabel [NO11]_
    or PDBFixer with OpenMM.
    
    :arg infile: PDB file name
    :type infile: str

    :arg method: Name of program which will be use to fix protein structure.
            Two alternative options are available: 'openbabel' and 'pdbfixer'.
            For either option additional software need to be installed:
            'openbabel': OpenBabel
            'pdbfixer': PDBFixer and OpenMM
            default is 'openbabel'
    :type method: str
    
    :arg pH: pH value applyed only for PDBfixer.
    :type pH: int, float
    
    Instalation of Openbabel:
    conda install -c conda-forge openbabel    

    Find more information here: https://anaconda.org/conda-forge/openbabel
                                https://github.com/openmm/pdbfixer
    Program will create new file in the same directory with 'addH_' prefix.
    
    .. [NO11] O'Boyle, N. M., Banck M., James C. A., Morley C., Vandermeersch T., Hutchison G. R. 
    Open Babel: An open chemical toolbox *Journal of cheminformatics* **2011** 3:1-14. """
    
    model_residues = kwargs.get("model_residues", False)

    import os

    if outfile == None:
        outfile = os.path.join(os.path.split(infile)[0], "addH_" + os.path.split(infile)[1])
        
    if outfile == infile:
        raise ValueError('outfile cannot be the same as infile')

    if method == 'openbabel':
        if model_residues:
            LOGGER.warn("Openbabel cannot add missing residues, skipping this step")

        if infile.endswith('cif'):
            raise ValueError('Openbabel cannot handle cif files')

        try:
            #import openbabel
            from openbabel import openbabel
            obconversion = openbabel.OBConversion()
            obconversion.SetInFormat("pdb")
            mol = openbabel.OBMol()
            obconversion.ReadFile(mol, infile)
            mol.AddHydrogens()
            obconversion.WriteFile(mol, outfile)
            LOGGER.info("Hydrogens were added to the structure. Structure {0} is saved in the local directry.".format(outfile))
        except ImportError:
            raise ImportError("Install Openbabel to add hydrogens to the structure or use PDBFixer/OpenMM.")
            
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

            fixer.removeHeterogens(True)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(pH)
            PDBFile.writeFile(fixer.topology, fixer.positions, open(outfile, 'w'))
            LOGGER.info("Hydrogens were added to the structure. New structure is saved as {0}.".format(outfile))

        except ImportError:
            raise ImportError('Install PDBFixer and OpenMM in order to fix the protein structure.')

    else:
        raise TypeError('Method should be openbabel or pdbfixer')
    
    return outfile

