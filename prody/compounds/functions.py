"""This module defines functions for using compounds from the PDB and elsewhere."""

from .pdbligands import PDBLigandRecord
import numpy as np

__all__ = ['calc2DSimilarity', 'calc2DSimilarityMatrix']


def calc2DSimilarity(smiles1, smiles2):
    """Calculate 2D similarity using Morgan Fingerprints
    
    :arg smiles1: first SMILES string or PDBLigandRecord containing one
    :type smiles1: str, :class:`.PDBLigandRecord`

    :arg smiles2: second SMILES string or PDBLigandRecord containing one
    :type smiles2: str, :class:`.PDBLigandRecord`
    """
    try:
        from rdkit import Chem
        from rdkit import DataStructs
        from rdkit.Chem.Fingerprints import FingerprintMols
        from rdkit.Chem import AllChem
    except ImportError:
        raise ImportError('rdkit is a required package for calc2DSimilarity')

    if not isinstance(smiles1, str):
        if hasattr(smiles1, 'getCanonicalSMILES'):
            smiles1 = smiles1.getCanonicalSMILES()
        else:
            raise TypeError('smiles1 should be a string or an object with method getCanonicalSMILES')

    if not isinstance(smiles2, str):
        if hasattr(smiles2, 'getCanonicalSMILES'):
            smiles2 = smiles2.getCanonicalSMILES()
        else:
            raise TypeError('smiles2 should be a string or an object with method getCanonicalSMILES')

    m1 = Chem.MolFromSmiles(smiles1)
    m2 = Chem.MolFromSmiles(smiles2)
    if m1 is not None and m2 is not None:
        fp1 = AllChem.GetMorganFingerprint(m1, 2, useFeatures=True)
        fp2 = AllChem.GetMorganFingerprint(m2, 2, useFeatures=True)
        simi_score = DataStructs.TanimotoSimilarity(fp1, fp2)

    return simi_score


def calc2DSimilarityMatrix(smiles_set):
    num_smiles = len(smiles_set)
    sim_mat = np.ones((num_smiles, num_smiles))

    for i in range(num_smiles):
        for j in range(i+1, num_smiles):
            sim_mat[i, j] = sim_mat[j, i] = calc2DSimilarity(smiles_set[i], smiles_set[j])

    return sim_mat
    