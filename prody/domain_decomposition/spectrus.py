import numpy as np
from scipy import sparse
from prody import LOGGER, SETTINGS
from prody import buildDistMatrix, calcDistFlucts

__all__ = ['calcSpectrusSims', 'MBSPointMutation', 'calcMBS', 'calcMBSfromSim']


def calcSpectrusSims(distFlucts, pdb, cutoff=10., sigma='MRSDF', **kwargs):
    
    coords = pdb.getCoords()
    n = coords.shape[0]
    
    if distFlucts.shape != (n, n):
        raise ValueError('distFlucts and atoms must have same linear '
                         'size (now %d and %d)' %(distFlucts.shape[0], n))

    # identify atom pairs within cutoff and store relative dist. flucts
    nearestNeighs = np.full((n, n), True, dtype=bool)
    np.fill_diagonal(nearestNeighs, False)
    if isinstance(cutoff, (int, float)):
        # compute inter-atomic distances
        dist = buildDistMatrix(coords)
        nearestNeighs &= (dist <= cutoff)
    elif cutoff is not None:
        raise ValueError('cutoff must be either a number or None. '
                         'Got: {0}'.format(type(cutoff)))
    nnDistFlucts = distFlucts[nearestNeighs]

    # set the sigma parameter for the Gaussian weights
    if sigma == 'MRSDF':
        # sigma is computed as the average of the root distance fluctuations
        # between residues within the distance cutoff, as defined in the 
        # SPECTRUS algorithm
        sigma = np.mean(np.sqrt(nnDistFlucts))
    elif sigma == 'RMSDF':
        # sigma is computed as the root mean squared dist. fluctuations
        # (faster to compute than MRSDF)
        sigma = np.sqrt(np.mean(nnDistFlucts))

    # check if sigma is a number
    try:
        ss = 2.*sigma**2
    except:
        raise ValueError('sigma must be \'MRSDF\', \'RMSDF\' or a number.')

    # compute the Gaussian weights only for residue pairs
    # within the distance cutoff
    reducedSims = np.where(nearestNeighs, np.exp(-distFlucts/ss), 0)
    np.fill_diagonal(reducedSims, 1.)
    sparseSims  = sparse.csr_matrix(reducedSims)
    sparse.csr_matrix.eliminate_zeros(sparseSims)
    
    return sparseSims, sigma


def MBSPointMutation(simMatrix, index, **kwargs):

    if sparse.issparse(simMatrix):
        # slightly faster
        newSim = simMatrix.tolil()
    else:
        newSim = simMatrix.copy()
    n = simMatrix.shape[0]

    # cut non-adjacent links around atom 'index'
    nonNearestNeighs = list(range(0,index-1)) + list(range(index+2,n))
    for j in nonNearestNeighs:
        newSim[index, j] = 0
        newSim[j, index] = 0
    return newSim


def _removeOutliers(data, Delta=100., **kwargs):
    assert Delta > 0.
    med = np.median(data[~np.isnan(data)])
    if not 0 <= med <= 1e+4:
        raise RuntimeError('MBS profile is not well defined: Check for ' +\
        'possible disconnected components in the PDB structure and remove them.')
    # compute median of |distances from median|
    dists = np.abs(data - med)
    mdev = np.median(dists[~np.isnan(dists)])
    # replace entries with 'nan' if outside the "safe"
    # interval (median - Delta*mdev, median + Delta*mdev)
    Delta_mdev = Delta * mdev
    for i, dist in enumerate(dists):
        if np.isnan(dist) or dist>Delta_mdev:
            data[i] = np.nan
    return data 


def calcMBSfromSim(simMatrix, nEvals=20, remove_outliers=True,
                   remove_offset=True, **kwargs):

    LOGGER.timeit('_MBS')
    n = simMatrix.shape[0]
    mbs = np.zeros(n) 
    for i in range(n):
        try:
            # cut "non-covalent" bonds around atom 'i'
            modSim = MBSPointMutation(simMatrix, i)
            # compute laplacian's spectrum of eigvals
            laplacian = sparse.csgraph.laplacian(modSim, normed=True)
            evals = sparse.linalg.eigsh(laplacian, k=min(nEvals, n-1), 
                                        which='SM', return_eigenvectors=False)
            # sort eigvals in ascending order
            evals = np.sort(evals)
            # compute MBS at site i
            mbs[i] = np.sum(1./evals[1:])
        except Exception as err:
            LOGGER.warn('Unable to compute MBS at position '
                        '{0}. {1}'.format(i, err))
            mbs[i] = np.nan
    if any(~np.isnan(mbs)):
        # remove outliers
        if remove_outliers is True:
            mbs = _removeOutliers(mbs, **kwargs)
        # remove offset
        if remove_offset is True:
            offset = min(mbs[~np.isnan(mbs)])
            mbs = mbs - offset 
    LOGGER.report('MBS computed in %.1fs.', '_MBS')

    return mbs


def calcMBS(anm, atomGroup, remove_outliers=True, remove_offset=True, 
            **kwargs):

    distFlucts = calcDistFlucts(anm, norm=False)
    sparseSim, sigma = calcSpectrusSims(distFlucts, atomGroup, **kwargs)
    mbs = calcMBSfromSim(sparseSim, **kwargs)

    return mbs



