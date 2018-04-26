import numpy as np
from scipy import sparse
from prody import *


def calcSpectrusSims(distFlucts, pdb, cutoff=15., sigma='MRSDF'):
    
    coords = pdb.getCoords()
    n = coords.shape[0]
    
    if distFlucts.shape != (n, n):
        raise ValueError('distFlucts and atoms must have same linear '
                         'size (now %d and %d)' %(distFlucts.shape[0], n))
    # identify atom pairs within cutoff and store relative dist. flucts
    
    nearestNeighs = np.full((n, n), True)
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
    
    return sparseSims.toarray(), sigma
