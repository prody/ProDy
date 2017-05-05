import numpy as np

__all__ = ['VCnorm', 'SQRTVCnorm', 'Filenorm']

## normalization methods ##
def div0(a, b):
    """ Performs ``true_divide`` but ignores the error when division by zero 
    (result is set to zero instead). """

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        if np.isscalar(c):
            if not np.isfinite(c):
                c = 0
        else:
            c[~np.isfinite(c)] = 0.  # -inf inf NaN
    return c

def VCnorm(M, **kwargs):
    """ Performs vanilla coverage normalization on matrix *M*."""

    C = np.diag(div0(1., np.sum(M, axis=0)))
    R = np.diag(div0(1., np.sum(M, axis=1)))

    # N = R * M * C
    N = np.dot(np.dot(R,M),C)

    sum_M = np.sum(M)
    sum_N = np.sum(N)
    k = sum_M / sum_N
    N = N * k
    return N

def SQRTVCnorm(M, **kwargs):
    """ Performs square-root vanilla coverage normalization on matrix *M*."""

    C = np.diag(np.sqrt(div0(1., np.sum(M, axis=0))))
    R = np.diag(np.sqrt(div0(1., np.sum(M, axis=1))))

    # N = R * M * C
    N = np.dot(np.dot(R,M),C)

    sum_M = np.sum(M)
    sum_N = np.sum(N)
    k = sum_M / sum_N
    N = N * k
    return N

def Filenorm(M, **kwargs):
    """ Performs normalization on matrix *M* given a file. *filename* specifies 
    the path to the file. The file should be one-column, and ideally has the 
    same number of entries with the size of *M* (extra entries will be ignored). 
    Say *F* is vector of the normalization factors, *N* is the normalized matrix, 
    if *expected* is ``False``, ``N[i,j] = M[i,j]/F[i]/F[j]``. If *expected* is
    ``True``, ``N[i,j] = M[i,j]/F[|i-j|]``."""

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
