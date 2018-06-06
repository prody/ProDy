import numpy as np

from prody.chromatin.functions import div0
from prody.utilities import importLA
from prody import LOGGER

__all__ = ['VCnorm', 'SQRTVCnorm', 'Filenorm', 'SCN']

def VCnorm(M, **kwargs):
    """ Performs vanilla coverage normalization on matrix *M*."""

    total_count = kwargs.get('total_count', 'original')

    C = np.diag(div0(1., np.sum(M, axis=0)))
    R = np.diag(div0(1., np.sum(M, axis=1)))

    # N = R * M * C
    N = np.dot(np.dot(R,M),C)

    if total_count is 'original':
        total_count = np.sum(M)

    if total_count is not None:
        sum_N = np.sum(N)
        k = total_count / sum_N
        N = N * k
    return N

def SQRTVCnorm(M, **kwargs):
    """ Performs square-root vanilla coverage normalization on matrix *M*."""

    total_count = kwargs.get('total_count', 'original')

    C = np.diag(np.sqrt(div0(1., np.sum(M, axis=0))))
    R = np.diag(np.sqrt(div0(1., np.sum(M, axis=1))))

    # N = R * M * C
    N = np.dot(np.dot(R,M),C)

    if total_count is 'original':
        total_count = np.sum(M)

    if total_count is not None:
        sum_N = np.sum(N)
        k = total_count / sum_N
        N = N * k
    return N

def SCN(M, **kwargs):
    la = importLA()
    total_count = kwargs.pop('total_count', None)
    max_loops = kwargs.pop('max_loops', 100)
    tol = kwargs.pop('tol', 1e-5)

    N = M.copy()
    n = 0
    d0 = None
    p = 1
    last_p = None

    while True:
        C = np.diag(div0(1., np.sum(N, axis=0)))
        N = np.dot(N, C)

        R = np.diag(div0(1., np.sum(N, axis=1)))
        N = np.dot(R, N)

        n += 1

        # check convergence of symmetry
        d = np.mean(np.abs(N - N.T))
        
        if d0 is not None:
            p = div0(d, d0)
            dp = np.abs(p - last_p)
            if dp < tol:
                break
        else:
            d0 = d
        LOGGER.debug('Iteration {0}: d = {1}, p = {2}'.format(str(n), str(d), str(p)))
        last_p = p
        
        if max_loops is not None:
            if n >= max_loops:
                LOGGER.warn('The SCN algorithm did not converge after {0} '
                            'iterations.'.format(max_loops))
                break
    # guarantee symmetry
    N = (N + N.T) / 2.
    if total_count is 'original':
        total_count = np.sum(M)

    if total_count is not None:
        sum_N = np.sum(N)
        k = total_count / sum_N
        N = N * k
    return N

def Filenorm(M, **kwargs):
    """ Performs normalization on matrix *M* given a file. *filename* specifies 
    the path to the file. The file should be one-column, and ideally has the 
    same number of entries with the size of *M* (extra entries will be ignored). 
    Say *F* is vector of the normalization factors, *N* is the normalized matrix, 
    if *expected* is **True**, ``N[i,j] = M[i,j]/F[i]/F[j]``. If *expected* is
    **True**, ``N[i,j] = M[i,j]/F[|i-j|]``."""

    filename = kwargs.get('filename')
    expected = kwargs.get('expected', False)

    if filename is None:
        raise IOError("'filename' is not specified.")
    factors = np.loadtxt(filename)
    L = M.shape[0]
    if not expected:
        factors.resize(L)
        norm_mat = np.outer(factors, factors)
        
        N = div0(M, norm_mat)
        return N
    else:
        N = np.zeros(M.shape)
        for i in range(L):
            for j in range(L):
                N[i,j] = div0(M[i,j], factors[abs(i-j)])
        return N
