# -*- coding: utf-8 -*-
"""This module defines functions for performing perturbation response scanning
from PCA and normal modes."""

import time

import numpy as np

from prody import LOGGER
from prody.proteins import parsePDB
from prody.atomic import AtomGroup, Selection
from prody.ensemble import Ensemble, Conformation
from prody.trajectory import TrajBase
from prody.utilities import importLA
from numpy import sqrt, arange, log, polyfit, array

from .nma import NMA
from .modeset import ModeSet
from .mode import VectorBase, Mode, Vector
from .gnm import GNMBase
from .analysis import calcCovariance

__all__ = ['calcPerturbResponse', 'parsePerturbResponseMatrix',
           'calcPerturbResponseProfiles', 'writePerturbResponsePDB']

class PRSMatrixParseError(Exception):
    pass


def calcPerturbResponse(model, **kwargs):

    """Returns a matrix of profiles from scanning the response of the
    structure to random perturbations at specific atom (or node) positions.
    The function implements the perturbation response scanning (PRS) method
    described in [CA09]_.  Rows of the matrix are the average magnitude of the
    responses obtained by perturbing the atom/node position at that row index,
    i.e. ``prs_profile[i,j]`` will give the response of residue/node *j* to
    perturbations in residue/node *i*.  PRS is performed using the covariance
    matrix from *model*, e.g. :class:`.ANM` instance.

    When an *atoms* instance is given, the PRS matrix will be added as data, 
    which can be retrieved with ``atoms.getData('prs_matrix')``.  

    *model* and *atoms* must have the same number of atoms. *atoms* must be an
    :class:`.AtomGroup` instance. 

    .. [CA09] Atilgan C, Atilgan AR, Perturbation-Response Scanning
       Reveals Ligand Entry-Exit Mechanisms of Ferric Binding Protein.
       *PLoS Comput Biol* **2009** 5(10):e1000544.

    The PRS matrix can be calculated and saved as follows::

      prs_matrix = calcPerturbResponse(p38_anm, saveMatrix=True)
      
    The PRS matrix can also be save later as follows::
    
      writeArray('prs_matrix.txt', prs_matrix, format='%8.6f', delimiter='\t')

    :arg saveMatrix: whether to save the last matrix generated to a text file.
        Default is False
    :type saveMatrix: bool

    :arg saveName: The file name for saved matrices
        Default is 'response_matrix.txt'.
    :type saveName: str
    """

    if not isinstance(model, (NMA, ModeSet, Mode)):
        raise TypeError('model must be an NMA, ModeSet, or Mode instance')

    if isinstance(model, NMA) and len(model) == 0:
        raise ValueError('model must have normal modes calculated')

    atoms = kwargs.get('atoms',None)
    if atoms is not None:
        if isinstance(atoms, Selection):
            atoms = atoms.copy()
        if not isinstance(atoms, AtomGroup):
            raise TypeError('atoms must be an AtomGroup instance')
        elif atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')

    n_atoms = model.numAtoms()
    LOGGER.timeit('_prody_prs_all')
    LOGGER.info('Calculating covariance matrix')
    LOGGER.timeit('_prody_cov')

    cov = calcCovariance(model)
    if cov is None:
        raise ValueError('model did not return a covariance matrix')

    LOGGER.clear()
    LOGGER.report('Covariance matrix calculated in %.1fs.',
                  '_prody_cov')

    LOGGER.progress('Calculating perturbation response', n_atoms, '_prody_prs_mat')

    if not model.is3d():
        prs_matrix = cov**2

    else:
        cov_squared = cov**2
        n_by_3n_cov_squared = np.zeros((n_atoms, 3 * n_atoms))
        prs_matrix = np.zeros((n_atoms, n_atoms))
        i3 = -3
        i3p3 = 0
        for i in range(n_atoms):
            i3 += 3
            i3p3 += 3
            n_by_3n_cov_squared[i,:] = (cov_squared[i3:i3p3,:]).sum(0)

        j3 = -3
        j3p3 = 0
        for j in range(n_atoms):
            j3 += 3
            j3p3 += 3                
            prs_matrix[:,j] = (n_by_3n_cov_squared[:,j3:j3p3]).sum(1)

    LOGGER.clear()
    LOGGER.report('Perturbation response matrix calculated in %.1fs.',
                      '_prody_prs_mat')

    saveMatrix = kwargs.get('saveMatrix',False)
    suppressDiag = kwargs.get('suppressDiag',False)
    saveName = kwargs.get('saveName','response_matrix.txt')

    norm_prs_matrix = np.zeros((n_atoms, n_atoms))
    self_dp = np.diag(prs_matrix)  
    self_dp = self_dp.reshape(n_atoms, 1)
    norm_prs_matrix = prs_matrix / np.repeat(self_dp, n_atoms, axis=1)

    if suppressDiag == True:
       # suppress the diagonal (self displacement) to facilitate
       # visualizing the response profile
       norm_prs_matrix = norm_prs_matrix - np.diag(np.diag(norm_prs_matrix))

    if saveMatrix == True:
        np.savetxt(saveName, norm_prs_matrix, delimiter='\t', fmt='%8.6f')

    LOGGER.report('Perturbation response scanning completed in %.1fs.',
                  '_prody_prs_all')

    if atoms is not None:
        atoms.setData('prs_matrix',norm_prs_matrix)
        return atoms, norm_prs_matrix
    else:
        return norm_prs_matrix

def parsePerturbResponseMatrix(prs_matrix_file='prs_matrix.txt',normMatrix=False):
    """Parses a perturbation response matrix from a file into a numpy ndarray.

    :arg prs_matrix_file: name of the file containing a PRS matrix, default is
        'prs_matrix.txt' as is used in the example under calcPerturbResponse.
    :type prs_matrix_file: str

    :arg normMatrix: whether to normalize the PRS matrix after parsing it.
        Default is False. If you used an old version of the script 
        and didn't normalize before saving, set this to True.
    :type norm: bool

    """
    fmat = open(prs_matrix_file,'r')
    matlines = fmat.readlines()
    fmat.close()

    prs_matrix = []
    for line in matlines:
       prs_matrix.append(line.split())

    for i in range(len(prs_matrix)):
     for j in range(len(prs_matrix)):
        prs_matrix[i][j] = float(prs_matrix[i][j])

    prs_matrix = np.array(prs_matrix)

    if normMatrix == True:
       # normalize the PRS matrix
       self_dp = np.diag(prs_matrix)  # using self displacement (diagonal of
                              # the original matrix) as a
                              # normalization factor
       self_dp = self_dp.reshape(len(prs_matrix), 1)
       norm_PRS_mat = prs_matrix / np.repeat(self_dp, len(prs_matrix), axis=1)
       return norm_PRS_mat

    else:
       return prs_matrix

def calcPerturbResponseProfiles(prs_matrix,atoms=None):
    """ Calculate the effectiveness and sensitivity
    profiles, which are the averages over the rows
    and columns of the PRS matrix.

    :arg prs_matrix: a perturbation response matrix
    :type prs_matrix: ndarray 

    When an *atoms* instance is given, the profiles will be added as data, 
    which can be retrieved with ``atoms.getData('effectiveness')`` and 
    ``atoms.getData('sensitivity')``. 
    """

    effectiveness = []
    sensitivity = []
    for i in range(len(prs_matrix)):
        effectiveness.append(np.mean(prs_matrix[i]))
        sensitivity.append(np.mean(prs_matrix.T[i]))

    effectiveness = np.array(effectiveness)
    sensitivity = np.array(sensitivity)

    if atoms is not None:
        atoms.setData('effectiveness',effectiveness)
        atoms.setData('sensitivity',sensitivity)
        return atoms, effectiveness, sensitivity
    else:
        return effectiveness, sensitivity

def writePerturbResponsePDB(prs_matrix,pdbIn=None,**kwargs):
    """ Write the average response to perturbation of
    a particular residue (a row of a perturbation response matrix)
    or the average effect of perturbation of a particular residue
    (a column of a normalized perturbation response matrix)
    into the b-factor field of a PDB file for visualisation in a
    molecular graphics program.
    If no chain is given this will be done for that residue in all chains.

    If no residue number is given then the effectiveness and sensitivity
    profiles will be written out instead. These two profiles are also returned
    as arrays for further analysis if they aren't already provided.

    :arg prs_matrix: a perturbation response matrix 
        or a :class:`.AtomGroup` object with a PRS matrix associated as data
    :type prs_matrix: array or :class:`.AtomGroup`

    :arg pdbIn: file name for the input PDB file where you would like the PRS
        data mapped
    :type pdbIn: str

    :arg pdbOut: a list of file names (enclosed in square
        brackets) for the output PDB file, default is to append
        the chain and residue info (name and number) onto the pdbIn stem.
        The input for pdbOut can also be used as a stem if you enter a 
        single string enclosed in quotes.
        If no residue number is supplied, chain is ignored and the default 
        is to append '_effectiveness' and '_sensitivity' onto the stem.
    :type pdbOut: list

    :arg chain: chain identifier for the residue of interest, default is all chains
        If you want to analyse residues in a subset of chains, concatentate them
        together e.g. 'AC'
    :type chain: str

    :arg resnum: residue number for the residue of interest
    :type resnum: int

    :arg direction: the direction you want to use to read data out
        of the PRS matrix for plotting: the options are 'effect' or 'response'.
        Default is 'effect'.
        A row gives the effect on each residue of peturbing the specified 
        residue.
        A column gives the response of the specified residue to perturbing 
        each residue.
        If no residue number is provided then this option will be ignored
    :type direction: str

    :arg returnData: whether to return effectiveness and sensitivity for analysis
        default is False
    :type returnProfiles: bool

    :arg effectiveness: effectiveness profile
    :type array

    :arg sensitivity: sensitivity profile
    :type array
    """

    if not isinstance(prs_matrix,np.ndarray):
        try:
            prs_matrix = prs_matrix.getData('prs_matrix')
        except:
            raise TypeError('Please provide a valid PRS matrix in numpy ndarray format.')

    try:
        fi = open(pdbIn,'r')
        lines = fi.readlines()
        fi.close()
    except:
        raise PRSMatrixParseError('Please provide a valid file name for the input PDB.')
 
    chain = kwargs.get('chain', None)

    structure = parsePDB(pdbIn,subset='ca')
    structure.setData('prs_matrix',prs_matrix)

    hv = structure.getHierView()
    chains = []
    for i in range(len(list(hv))):
        chainAg = list(hv)[i]
        chains.append(chainAg.getChids()[0])

    chains = np.array(chains)
    if chain is None:
        chain = ''.join(chains)

    resnum = kwargs.get('resnum', None)
    pdbOut = kwargs.get('pdbOut', None)
    if pdbOut is None:
        out_stem = pdbIn.split('.')[0]
    elif type(pdbOut) is str:
        out_stem = pdbOut.split('.')[0]
        pdbOut = None

    if resnum is None:
        effectiveness = kwargs.get('effectiveness',None)
        sensitivity = kwargs.get('sensitivity',None)
        if effectiveness is None or sensitivity is None:
            effectiveness, sensitivity = calcPerturbResponseProfiles(prs_matrix)

        structure.setData('effectiveness',effectiveness)
        structure.setData('sensitivity',sensitivity)

        file_effs_name = '{0}_effectiveness.pdb'.format(out_stem)
        file_sens_name = '{0}_sensitivity.pdb'.format(out_stem)
        fileEffs = open(file_effs_name,'w')
        fileSens = open(file_sens_name,'w')

        for line in lines:            
            if line.find('ATOM') != 0 and line.find('HETATM') != 0 and line.find('ANISOU') != 0:
                fileEffs.write(line)                    
                fileSens.write(line)
            elif line.find('ATOM') == 0:
                fileEffs.write(line[:60] + '{:6.2f}'.format(float(structure.select( \
                               'chain {0} and resnum {1}'.format(line[21],line[22:26])) \
                               .getData('effectiveness')) * 100/np.max( \
                               structure.getData('effectiveness'))) + line[66:])
                fileSens.write(line[:60] + '{:6.2f}'.format(float(structure.select( \
                               'chain {0} and resnum {1}'.format(line[21],line[22:26])) \
                               .getData('sensitivity')) * 100/np.max( \
                               structure.getData('sensitivity'))) + line[66:])
            elif line.find('HETATM') == 0:
                fileEffs.write(line[:60] + '  0.00' + line[66:])
                fileSens.write(line[:60] + '  0.00' + line[66:])
                      
        fileEffs.close()
        fileSens.close()
        LOGGER.info('The effectiveness and sensitivity profiles were written' \
                    ' to {0} and {1}.'.format(file_effs_name,file_sens_name))

        returnData = kwargs.get('returnData',False)
        if returnData:
            return structure, effectiveness, sensitivity
        else:
            return
 
    direction = kwargs.get('direction','effect')
    for n in range(len(chain)):
        if not chain[n] in chains:
            raise PRSMatrixParseError('Chain {0} was not found in {1}'.format(chain[n], pdbIn))

    if pdbOut is None:
        pdbOut = []
        for c in chain:
            pdbOut.append('{0}_{1}_{2}{3}_{4}.pdb' \
                          .format(out_stem, c, \
                                  str(structure.select('chain {0} and resnum {1}' \
                                      .format(c, resnum)).getResnames()), \
                                  resnum, direction))

    for c in chain:
        fo = open(pdbOut[n],'w')
        for line in lines:
            if line.find('ATOM') != 0 and line.find('HETATM') != 0 and line.find('ANISOU') != 0:
                fo.write(line)
            elif line.find('ATOM') == 0:
                if direction is 'effect':
                    fo.write(line[:60] + '{:6.2f}'.format(float(structure.getData('prs_matrix') \
                                         [structure.select('chain {0} and resnum {1}' \
                                          .format(c, resnum)).getResindices(), \
                                          structure.select('chain {0} and resnum {1}' \
                                          .format(line[21], line[22:26])).getResindices()])*100) \
                             + line[66:])
                else:
                    fo.write(line[:60] + '{:6.2f}'.format(float(structure.getData('prs_matrix') \
                                         [structure.select('chain {0} and resnum {1}' \
                                          .format(line[21], line[22:26])).getResindices(), \
                                          structure.select('chain {0} and resnum {1}' \
                                          .format(c, resnum)).getResindices()])*100) \
                             + line[66:])
            elif line.find('HETATM') == 0:
                fo.write(line[:60] + '  0.00' + line[66:])

        LOGGER.info('Perturbation responses for specific residues were written' \
                    ' to {0}.'.format(', '.join(pdbOut)))

