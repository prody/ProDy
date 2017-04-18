# -*- coding: utf-8 -*-
"""This module defines functions for calculating physical properties from normal
modes."""

import time

import numpy as np

from prody import LOGGER
from prody.proteins import parsePDB
from prody.atomic import AtomGroup
from prody.ensemble import Ensemble, Conformation
from prody.trajectory import TrajBase
from prody.utilities import importLA
from numpy import sqrt, arange, log, polyfit

from .nma import NMA
from .modeset import ModeSet
from .mode import VectorBase, Mode, Vector
from .gnm import GNMBase

__all__ = ['calcCollectivity', 'calcCovariance', 'calcCrossCorr',
           'calcFractVariance', 'calcSqFlucts', 'calcTempFactors',
           'calcProjection', 'calcCrossProjection', 'calcPerturbResponse', 
           'parsePerturbResponseMatrix', 'writePerturbResponsePDB',
           'calcSpecDimension', 'calcPairDeformationDist',]

class PRSMatrixParseError(Exception):
    pass

def calcCollectivity(mode, masses=None):
    """Returns collectivity of the mode.  This function implements collectivity
    as defined in equation 5 of [BR95]_.  If *masses* are provided, they will
    be incorporated in the calculation.  Otherwise, atoms are assumed to have
    uniform masses.

    .. [BR95] Bruschweiler R. Collective protein dynamics and nuclear
       spin relaxation. *J Chem Phys* **1995** 102:3396-3403.

    :arg mode: mode or vector
    :type mode: :class:`.Mode` or :class:`.Vector`

    :arg masses: atomic masses
    :type masses: :class:`numpy.ndarray`"""

    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance')

    is3d = mode.is3d()
    if masses is not None:
        if len(masses) != mode.numAtoms():
            raise ValueError('length of masses must be equal to number of atoms')
        if is3d:
            u2in = (mode.getArrayNx3() ** 2).sum(1) / masses
    else:
        if is3d:
            u2in = (mode.getArrayNx3() ** 2).sum(1)
        else:
            u2in = (mode.getArrayNx3() ** 2)
    u2in = u2in * (1 / u2in.sum() ** 0.5)
    coll = np.exp(-(u2in * np.log(u2in)).sum()) / mode.numAtoms()
    return coll

def calcSpecDimension(mode):

    """
    :arg mode: mode or vector
    :type mode: :class:`.Mode` or :class:`.Vector`

    """
    # if not isinstance(mode, Mode):
    #     raise TypeError('mode must be a Mode instance')
    
    length = mode.shape[0]
    numbers = arange(2,length+1)
    ds,p=polyfit(log(sqrt(mode[0:int(length*0.25)])),log(numbers[0:int(length*0.25)]),1)
    
    return ds

def calcFracDimension(mode):
    """
    :arg mode: mode or vector
    :type mode: mode or vector """




def calcFractVariance(mode):
    """Returns fraction of variance explained by the *mode*.  Fraction of
    variance is the ratio of the variance along a mode to the trace of the
    covariance matrix of the model."""

    if isinstance(mode, Mode):
        var = mode.getVariance()
        trace = mode.getModel()._getTrace()
    elif isinstance(mode, (ModeSet, NMA)):
        var = mode.getVariances()
        if isinstance(mode, ModeSet):
            trace = mode.getModel()._getTrace()
        else:
            trace = mode._getTrace()
    else:
        raise TypeError('mode must be a Mode instance')
    if trace is None:
        raise ValueError('modes are not calculated')

    return var / trace


def calcProjection(ensemble, modes, rmsd=True):
    """Returns projection of conformational deviations onto given modes.
    *ensemble* coordinates are used to calculate the deviations that are
    projected onto *modes*.  For K conformations and M modes, a (K,M)
    matrix is returned.

    :arg ensemble: an ensemble, trajectory or a conformation for which
        deviation(s) will be projected, or a deformation vector
    :type ensemble: :class:`.Ensemble`, :class:`.Conformation`,
        :class:`.Vector`, :class:`.Trajectory`
    :arg modes: up to three normal modes
    :type modes: :class:`.Mode`, :class:`.ModeSet`, :class:`.NMA`

    By default root-mean-square deviation (RMSD) along the normal mode is
    calculated. To calculate the projection pass ``rmsd=True``.
    :class:`.Vector` instances are accepted as *ensemble* argument to allow
    for projecting a deformation vector onto normal modes."""

    if not isinstance(ensemble, (Ensemble, Conformation, Vector, TrajBase)):
        raise TypeError('ensemble must be Ensemble, Conformation, Vector, '
                        'or a TrajBase, not {0}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, VectorBase)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0}'
                        .format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be 3-dimensional')
    if isinstance(ensemble, Vector):
        n_atoms = ensemble.numAtoms()
    else:
        n_atoms = ensemble.numSelected()
    if n_atoms != modes.numAtoms():
        raise ValueError('number of atoms are not the same')
    if isinstance(ensemble, Vector):
        if not ensemble.is3d():
            raise ValueError('ensemble must be a 3d vector instance')
        deviations = ensemble._getArray()
    elif isinstance(ensemble, (Ensemble, Conformation)):
        deviations = ensemble.getDeviations()
    else:
        nfi = ensemble.nextIndex()
        ensemble.goto(0)
        deviations = np.array([frame.getDeviations() for frame in ensemble])
        ensemble.goto(nfi)
    if deviations.ndim == 3:
        deviations = deviations.reshape((deviations.shape[0],
                                         deviations.shape[1] * 3))
    elif deviations.ndim == 2:
        deviations = deviations.reshape((1, deviations.shape[0] * 3))
    else:
        deviations = deviations.reshape((1, deviations.shape[0]))
    projection = np.dot(deviations, modes._getArray())
    if rmsd:
        projection = (1 / (n_atoms ** 0.5)) * projection
    return projection


def calcCrossProjection(ensemble, mode1, mode2, scale=None, **kwargs):
    """Returns projection of conformational deviations onto modes from
    different models.

    :arg ensemble: ensemble for which deviations will be projected
    :type ensemble: :class:`.Ensemble`
    :arg mode1: normal mode to project conformations onto
    :type mode1: :class:`.Mode`, :class:`.Vector`
    :arg mode2: normal mode to project conformations onto
    :type mode2: :class:`.Mode`, :class:`.Vector`
    :arg scale: scale width of the projection onto mode1 (``x``) or mode2(``y``),
        an optimized scaling factor (scalar) will be calculated by default 
        or a value of scalar can be passed."""

    if not isinstance(ensemble, (Ensemble, Conformation, Vector, TrajBase)):
        raise TypeError('ensemble must be Ensemble, Conformation, Vector, '
                        'or a Trajectory, not {0}'.format(type(ensemble)))
    if not isinstance(mode1, VectorBase):
        raise TypeError('mode1 must be a Mode instance, not {0}'
                        .format(type(mode1)))
    if not mode1.is3d():
        raise ValueError('mode1 must be 3-dimensional')
    if not isinstance(mode2, VectorBase):
        raise TypeError('mode2 must be a Mode instance, not {0}'
                        .format(type(mode2)))
    if not mode2.is3d():
        raise ValueError('mode2 must be 3-dimensional')

    if scale is not None:
        assert isinstance(scale, str), 'scale must be a string'
        scale = scale.lower()
        assert scale in ('x', 'y'), 'scale must be x or y'

    xcoords = calcProjection(ensemble, mode1, kwargs.get('rmsd', True))
    ycoords = calcProjection(ensemble, mode2, kwargs.pop('rmsd', True))
    if scale:
        scalar = kwargs.get('scalar', None)
        if scalar:
            assert isinstance(scalar, (float, int)), 'scalar must be a number'
        else:
            scalar = ((ycoords.max() - ycoords.min()) /
                      (xcoords.max() - xcoords.min())
                      ) * np.sign(np.dot(xcoords, ycoords))
            if scale == 'x':
                LOGGER.info('Projection onto {0} is scaled by {1:.2f}'
                            .format(mode1, scalar))
            else:
                scalar = 1 / scalar
                LOGGER.info('Projection onto {0} is scaled by {1:.2f}'
                            .format(mode2, scalar))

        if scale == 'x':
            xcoords = xcoords * scalar
        else:
            ycoords = ycoords * scalar

    return xcoords, ycoords


def calcSqFlucts(modes):
    """Returns sum of square-fluctuations for given set of normal *modes*.
    Square fluctuations for a single mode is obtained by multiplying the
    square of the mode array with the variance (:meth:`.Mode.getVariance`)
    along the mode.  For :class:`.PCA` and :class:`.EDA` models built using
    coordinate data in Å, unit of square-fluctuations is |A2|, for
    :class:`.ANM` and :class:`.GNM`, on the other hand, it is arbitrary or
    relative units."""

    if not isinstance(modes, (VectorBase, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0}'.format(type(modes)))
    is3d = modes.is3d()
    if isinstance(modes, Vector):
        if is3d:
            return (modes._getArrayNx3()**2).sum(axis=1)
        else:
            return (modes._getArray() ** 2)
    else:
        sq_flucts = np.zeros(modes.numAtoms())
        if isinstance(modes, VectorBase):
            modes = [modes]
        for mode in modes:
            if is3d:
                sq_flucts += ((mode._getArrayNx3()**2).sum(axis=1) *
                              mode.getVariance())
            else:
                sq_flucts += (mode._getArray() ** 2) * mode.getVariance()
        return sq_flucts


def calcCrossCorr(modes, n_cpu=1):
    """Returns cross-correlations matrix.  For a 3-d model, cross-correlations
    matrix is an NxN matrix, where N is the number of atoms.  Each element of
    this matrix is the trace of the submatrix corresponding to a pair of atoms.
    Covariance matrix may be calculated using all modes or a subset of modes
    of an NMA instance.  For large systems, calculation of cross-correlations
    matrix may be time consuming.  Optionally, multiple processors may be
    employed to perform calculations by passing ``n_cpu=2`` or more."""

    if not isinstance(n_cpu, int):
        raise TypeError('n_cpu must be an integer')
    elif n_cpu < 1:
        raise ValueError('n_cpu must be equal to or greater than 1')

    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0}'.format(type(modes)))

    if modes.is3d():
        model = modes
        if isinstance(modes, (Mode, ModeSet)):
            model = modes._model
            if isinstance(modes, (Mode)):
                indices = [modes.getIndex()]
                n_modes = 1
            else:
                indices = modes.getIndices()
                n_modes = len(modes)
        else:
            n_modes = len(modes)
            indices = np.arange(n_modes)
        array = model._array
        n_atoms = model._n_atoms
        variances = model._vars
        if n_cpu == 1:
            s = (n_modes, n_atoms, 3)
            arvar = (array[:, indices]*variances[indices]).T.reshape(s)
            array = array[:, indices].T.reshape(s)
            covariance = np.tensordot(array.transpose(2, 0, 1),
                                      arvar.transpose(0, 2, 1),
                                      axes=([0, 1], [1, 0]))
        else:
            import multiprocessing
            n_cpu = min(multiprocessing.cpu_count(), n_cpu)
            queue = multiprocessing.Queue()
            size = n_modes / n_cpu
            for i in range(n_cpu):
                if n_cpu - i == 1:
                    indices = modes.indices[i*size:]
                else:
                    indices = modes.indices[i*size:(i+1)*size]
                process = multiprocessing.Process(
                    target=_crossCorrelations,
                    args=(queue, n_atoms, array, variances, indices))
                process.start()
            while queue.qsize() < n_cpu:
                time.sleep(0.05)
            covariance = queue.get()
            while queue.qsize() > 0:
                covariance += queue.get()
    else:
        covariance = calcCovariance(modes)
    diag = np.power(covariance.diagonal(), 0.5)
    return covariance / np.outer(diag, diag)


def _crossCorrelations(queue, n_atoms, array, variances, indices):
    """Calculate covariance-matrix for a subset of modes."""

    n_modes = len(indices)
    arvar = (array[:, indices] * variances[indices]).T.reshape((n_modes,
                                                                n_atoms, 3))
    array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
    covariance = np.tensordot(array.transpose(2, 0, 1),
                              arvar.transpose(0, 2, 1),
                              axes=([0, 1], [1, 0]))
    queue.put(covariance)


def calcTempFactors(modes, atoms):
    """Returns temperature (β) factors calculated using *modes* from a
    :class:`.ANM` or :class:`.GNM` instance scaled according to the
    experimental β-factors from *atoms*."""

    model = modes.getModel()
    if not isinstance(model, GNMBase):
        raise TypeError('modes must come from GNM or ANM')
    if model.numAtoms() != atoms.numAtoms():
        raise ValueError('modes and atoms must have same number of nodes')
    sqf = calcSqFlucts(modes)
    return sqf / ((sqf**2).sum()**0.5) * (atoms.getBetas()**2).sum()**0.5


def calcCovariance(modes):
    """Returns covariance matrix calculated for given *modes*."""

    if isinstance(modes, Mode):
        array = modes._getArray()
        return np.outer(array, array) * modes.getVariance()
    elif isinstance(modes, ModeSet):
        array = modes._getArray()
        return np.dot(array, np.dot(np.diag(modes.getVariances()), array.T))
    elif isinstance(modes, NMA):
        return modes.getCovariance()
    else:
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance')


def calcPerturbResponse(model, atoms=None, repeats=100, saveOrig=False, \
                        normMatrix=False, suppressDiag=False, saveNorm=False, \
                        baseSaveName='response_matrix', operation='mean'):
    """Returns a matrix of profiles from scanning of the response of the
    structure to random perturbations at specific atom (or node) positions.
    The function implements the perturbation response scanning (PRS) method
    described in [CA09]_.  Rows of the matrix are the average magnitude of the
    responses obtained by perturbing the atom/node position at that row index,
    i.e. ``prs_profile[i,j]`` will give the response of residue/node *j* to
    perturbations in residue/node *i*.  PRS is performed using the covariance
    matrix from *model*, e.t. :class:`.ANM` instance.  Each residue/node is
    perturbed *repeats* times with a random unit force vector.  When *atoms*
    instance is given, PRS profile for residues will be added as an attribute
    which then can be retrieved as ``atoms.getData('prs_profile')``.  *model*
    and *atoms* must have the same number of atoms. *atoms* must be an
    :class:`.AtomGroup` instance.

    :arg operation: which operation to perform to get a single response matrix::
        the mean, variance or max across the set of repeats. Default is mean. 
        To obtain all response matrices, set operation=None without quotes.
        You can also ask for 'all 3' operations or provide a list containing
        any set of them.
    :type operation: str or list

    .. [CA09] Atilgan C, Atilgan AR, Perturbation-Response Scanning
       Reveals Ligand Entry-Exit Mechanisms of Ferric Binding Protein.
       *PLoS Comput Biol* **2009** 5(10):e1000544.

    The PRS matrix can be calculated and saved as follows::

      prs_matrix = calcPerturbationResponse(p38_anm, saveMatrix=True)
      
    The PRS matrix can also be save later as follows::
    
      writeArray('prs_matrix.txt', prs_matrix, format='%8.6f', delimiter='\t')
      
    """

    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance')
    elif not model.is3d():
        raise TypeError('model must be a 3-dimensional NMA instance')
    elif len(model) == 0:
        raise ValueError('model must have normal modes calculated')
    if atoms is not None:
        if not isinstance(atoms, AtomGroup):
            raise TypeError('atoms must be an AtomGroup instance')
        elif atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')

    assert isinstance(repeats, int), 'repeats must be an integer'
    cov = calcCovariance(model)
    if cov is None:
        raise ValueError('model did not return a covariance matrix')

    n_atoms = model.numAtoms()
    response_matrix = np.zeros((repeats, n_atoms, n_atoms))
    LOGGER.progress('Calculating perturbation response', n_atoms, '_prody_prs')
    i3 = -3
    i3p3 = 0
    for i in range(n_atoms):
        i3 += 3
        i3p3 += 3
        forces = np.random.rand(repeats * 3).reshape((repeats, 3))
        forces /= ((forces**2).sum(1)**0.5).reshape((repeats, 1))
        for n in range(repeats):
            force = forces[n]
            response_matrix[n,i,:] = (
                np.dot(cov[:, i3:i3p3], force)
                ** 2).reshape((n_atoms, 3)).sum(1)
        LOGGER.update(i, '_prody_prs')

    if operation is not None:
        if type(operation) is str:
            if operation == 'all 3' or operation == 'all operations':
                operationList = ['var','max','mea']
            else:
                operationList = []
                operationList.append(operation.lower()[:3])
        elif type(operation) is list:
            operationList = operation
            for i in range(len(operationList)):
                operationList[i] = operationList[i].lower()[:3]

        operationList = np.array(operationList) 
        matrix_set = np.zeros((len(operationList),n_atoms,n_atoms))
        found_valid_operation = False

        if 'var' in operationList:
            found_valid_operation = True
            var_response_matrix = np.zeros((n_atoms, n_atoms))
            for i in range(n_atoms):
                for j in range(n_atoms):
                    var_response_matrix[i,j] = np.var(response_matrix[:,i,j])
            matrix_set[np.where(operationList == 'var')[0][0]] = var_response_matrix

        if 'max' in operationList:
            found_valid_operation = True
            max_response_matrix = np.zeros((n_atoms, n_atoms))
            for i in range(n_atoms):
                for j in range(n_atoms):
                    max_response_matrix[i,j] = np.max(response_matrix[:,i,j])
            matrix_set[np.where(operationList == 'max')[0][0]] = max_response_matrix

        if 'mea' in operationList:
            found_valid_operation = True
            mean_response_matrix = np.zeros((n_atoms, n_atoms))
            for i in range(n_atoms):
                for j in range(n_atoms):
                    mean_response_matrix[i,j] = np.mean(response_matrix[:,i,j]) 
            matrix_set[np.where(operationList == 'mea')[0][0]] = mean_response_matrix

        if not found_valid_operation:
            raise ValueError('Operation should be mean, variance, max in quotes ' \
                             'or a list of containing a set of these or None.')

    LOGGER.clear()
    LOGGER.report('Perturbation response scanning completed in %.1fs.',
                  '_prody_prs')


    if operation is None:
        LOGGER.info('Operation is None so all {0} repeats are output.' \
                    ' This is not compatible with saving, normalizing' \
                    ' or mapping to atoms at present.'.format(repeats))
        return response_matrix

    if atoms is not None:
        atoms.setData('prs_profile', response_matrix)

    if saveOrig == True:
       # save the original PRS matrix for each operation
       for i in range(len(operationList)):
           np.savetxt('orig_{0}_{1}'.format(baseSaveName,operationList[i]), \
                      matrix_set[i], delimiter='\t', fmt='%8.6f')
           
    if normMatrix == True:
        norm_PRS_mat = np.zeros((len(operationList),n_atoms,n_atoms))
        # calculate the normalized PRS matrix for each operation
        for i in range(len(operationList)):
            self_dp = np.diag(matrix_set[i])  # using self displacement (diagonal of
                                              # the original matrix) as a
                                              # normalization factor
            self_dp = self_dp.reshape(n_atoms, 1)
            norm_PRS_mat[i] = matrix_set[i] / np.repeat(self_dp, n_atoms, axis=1)

            if suppressDiag == True:
                # suppress the diagonal (self displacement) to facilitate
                # visualizing the response profile
                norm_PRS_mat[i] = norm_PRS_mat[i] - np.diag(np.diag(norm_PRS_mat[i]))

            if saveNorm == True:
                np.savetxt('norm_{0}_{1}'.format(baseSaveName,operationList[i]), \
                           norm_PRS_mat[i], delimiter='\t', fmt='%8.6f')
           
    if normMatrix == True:
        return norm_PRS_mat
    else:
        return matrix_set

def parsePerturbResponseMatrix(prs_matrix_file='prs_matrix.txt',normMatrix=True):
    """Parses a perturbation response matrix from a file into a numpy ndarray.

    :arg prs_matrix_file: name of the file containing a PRS matrix, default is
        'prs_matrix.txt' as is used in the example under calcPerturbResponse.
    :type prs_matrix_file: str

    :arg normMatrix: whether to normalise the PRS matrix after parsing it.
        Default is True. If you have already normalised this before saving,
        or you do not want it normalised then set this to False.
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

def writePerturbResponsePDB(prs_matrix,pdbIn,**kwargs):
    """ This function writes the average response to perturbation of
    a particular residue (a row of a perturbation response matrix)
    into the b-factor field of a PDB file for visualisation in PyMOL.
    If no chain is given this will be done for that residue in all chains.
    
    If no residue number is given then the effectiveness and sensitivity
    profiles will be written out instead. These two profiles are also output 
    as arrays for further analysis.

    :arg prs_matrix: name of the variable containing a perturbation response 
        matrix
    :type prs_matrix: ndarray

    :arg pdbInFile: file name for the input PDB file where you would like the PRS
        data mapped
    :type pdbIn: str

    :arg pdbOutFiles: a list of file names (enclosed in square
        brackets) for the output PDB file, default is to append
        the chain and residue info (name and number) onto the pdbIn stem.
        When multiple chains are to be used, a single file name can be 
        entered as string (encloded in quotes) and the chain IDs will be
        appended onto the pdbOut stem.
        If no residue number is supplied, pdbOut and chain are ignored and
        it appends 'effectiveness' and 'sensitivity' onto the pdbIn stem.
    :type pdbOut: list

    :arg chain: chain identifier for the residue of interest, default is all chains
        If you want to analyse residues in a subset of chains, concatentate them
        together e.g. 'AC'
    :type chain: str

    :arg resnum: residue number for the residue of interest
    :type resnum: int
    """

    if not type(prs_matrix) is np.ndarray:
        raise TypeError('Please provide a valid PRS matrix in numpy ndarray format.')
           
    try:
        fi = open(pdbIn,'r')
        lines = fi.readlines()
        fi.close()
    except:
        raise PRSMatrixParseError('Please provide a valid file name for the input PDB.')

    chain = kwargs.get('chain', None)
    structure = parsePDB(pdbIn).calpha
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
        pdbOut = []
        for c in chain:
            pdbOut.append('{0}_{1}_{2}{3}.pdb'.format(pdbIn.split('.')[0], c, \
                              structure.getResnames()[i], resnum))
    elif type(pdbOut) is str:
        pdbOut2 = []
        for c in chain:
            pdbOut2.append(pdbOut.split('.')[0] + '_' + c + pdbOut.split('.')[1])
        pdbOut = pdbOut2

    if resnum is None:
        effectiveness = []
        sensitivity = []
        for i in range(len(prs_matrix)):
            effectiveness.append(np.mean(prs_matrix[i]))
            sensitivity.append(np.mean(prs_matrix.T[i]))

        file_effs_name = '{0}_effectiveness.pdb'.format(pdbIn.split('.')[0])
        file_sens_name = '{0}_sensitivity.pdb'.format(pdbIn.split('.')[0])
        fileEffs = open(file_effs_name,'w')
        fileSens = open(file_sens_name,'w')

        for line in lines:            
            if line.find('ATOM') != 0:  
                fileEffs.write(line)                    
                fileSens.write(line)
            else:
                resnum_matrix_offset = np.where(structure.getChids() == line[21]) \
                                       [0][0] - structure.getResnums() \
                                       [np.where(structure.getChids() == line[21])[0][0]]
                i = int(line.split()[5]) + resnum_matrix_offset
                fileEffs.write(line[:60] + ' '*(6-len('{:3.2f}'.format((effectiveness[i])*100))) \
                         + '{:3.2f}'.format((effectiveness[i])*100) + line[66:])
                fileSens.write(line[:60] + ' '*(6-len('{:3.2f}'.format((sensitivity[i])*10))) \
                         + '{:3.2f}'.format((sensitivity[i])*10) + line[66:])
                      
        fileEffs.close()
        fileSens.close()
        LOGGER.info('The effectiveness and sensitivity profiles were written', \
                       ' to {0} and {1}.'.format(file_effs_name,file_sens_name))
        return effectiveness, sensitivity

    outFiles = []
    for n in range(len(chain)):
        if not chain[n] in chains:
            raise PRSMatrixParseError('Chain {0} was not found in {1}'.format(chain[n], pdbIn))

        chainNum = int(np.where(chains == chain[n])[0])
        chainAg = list(hv)[chainNum]
        if not resnum in chainAg.getResnums():
            raise PRSMatrixParseError('A residue with number {0} was not found', 
                                      ' in chain {1}'.format(resnum, chain[n]))

        resnum_matrix_offset = np.where(structure.getResnums() == \
                                        chainAg.getResnums()[0])[0][chainNum] \
                               - chainAg.getResnums()[0]
        i = resnum + resnum_matrix_offset

        fo = open(pdbOut[n],'w')
        for line in lines:
            if line.find('ATOM') != 0:
                fo.write(line)
            else:
                resnum_matrix_offset = np.where(structure.getChids() == line[21]) \
                                       [0][0] - structure.getResnums() \
                                       [np.where(structure.getChids() == line[21])[0][0]]
                j = int(line.split()[5]) + resnum_matrix_offset
            fo.write(line[:60] + ' '*(6-len('{:3.2f}'.format((prs_matrix[i][j])*10))) \
                     + '{:3.2f}'.format((prs_matrix[i][j])*10) + line[66:])
        fo.close()
        outFiles.append(fo)
        LOGGER.report('Perturbation responses for specific residues were written', 
                       ' to {0} and {1}.'.format(' '.join(outFiles)))
    return outFiles


def calcPairDeformationDist(model, coords, ind1, ind2, kbt=1.):
                                                
    """Returns distribution of the deformations in the distance contributed by each mode 
    for selected pair of residues *ind1* *ind2* using *model* from a :class:`.ANM`.
    Method described in [EB08]_ equation (10) and figure (2).     
    
    .. [EB08] Eyal E., Bahar I. Toward a Molecular Understanding of 
        the Anisotropic Response of Proteins to External Forces:
        Insights from Elastic Network Models. *Biophys J* **2008** 94:3424-34355. 
    
    :arg model: this is an 3-dimensional NMA instance from a :class:`.ANM
    calculations.
    :type model: :class:`.ANM`  
    :arg coords: a coordinate set or an object with ``getCoords`` method.
      Recommended: coords = parsePDB('pdbfile').select('protein and name CA').
    :type coords: :class:`numpy.ndarray`.
    :arg ind1: first residue number.
    :type ind1: int 
    :arg ind2: secound residue number.
    :type ind2: int 
    """

    try:
        resnum_list = coords.getResnums()
        resnam_list = coords.getResnames()
        coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                coords.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be a Numpy array or an object '
                            'with `getCoords` method')
    
    if not isinstance(model, NMA):
        raise TypeError('model must be a NMA instance')
    elif not model.is3d():
        raise TypeError('model must be a 3-dimensional NMA instance')
    elif len(model) == 0:
        raise ValueError('model must have normal modes calculated')
    elif model.getStiffness() is None:
        raise ValueError('model must have stiffness matrix calculated')
    
    linalg = importLA()
    n_atoms = model.numAtoms()
    n_modes = model.numModes()
    LOGGER.timeit('_pairdef')

    r_ij = np.zeros((n_atoms,n_atoms,3))
    r_ij_norm = np.zeros((n_atoms,n_atoms,3))

    for i in range(n_atoms):
        for j in range(i+1,n_atoms):
            r_ij[i][j] = coords[j,:] - coords[i,:]
            r_ij[j][i] = r_ij[i][j]
            r_ij_norm[i][j] = r_ij[i][j]/linalg.norm(r_ij[i][j])
            r_ij_norm[j][i] = r_ij_norm[i][j]

    eigvecs = model.getEigvecs()
    eigvals = model.getEigvals()
    
    D_pair_k = []
    mode_nr = []
    ind1 = ind1 - resnum_list[0]
    ind2 = ind2 - resnum_list[0]

    for m in xrange(6,n_modes):
        U_ij_k = [(eigvecs[ind1*3][m] - eigvecs[ind2*3][m]), (eigvecs[ind1*3+1][m] \
            - eigvecs[ind2*3+1][m]), (eigvecs[ind1*3+2][m] - eigvecs[ind2*3+2][m])] 
        D_ij_k = abs(np.sqrt(kbt/eigvals[m])*(np.vdot(r_ij_norm[ind1][ind2], U_ij_k)))  
        D_pair_k.append(D_ij_k)
        mode_nr.append(m)

    LOGGER.report('Deformation was calculated in %.2lfs.', label='_pairdef')
    
    return mode_nr, D_pair_k

