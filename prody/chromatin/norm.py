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

    if total_count == 'original':
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

    if total_count == 'original':
        total_count = np.sum(M)

    if total_count is not None:
        sum_N = np.sum(N)
        k = total_count / sum_N
        N = N * k
    return N

def SCN(M, **kwargs):
    """ Performs Sequential Component Normalization on matrix *M*.
    
    .. [AC12] Cournac A, Marie-Nelly H, Marbouty M, Koszul R, Mozziconacci J. 
       Normalization of a chromosomal contact map. *BMC Genomics* **2012**. 
    """

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
    if total_count == 'original':
        total_count = np.sum(M)

    if total_count is not None:
        sum_N = np.sum(N)
        k = total_count / sum_N
        N = N * k
    return N

def bnewt(A, mask=[], tol = 1e-6, delta_lower = 0.1, delta_upper = 3, fl = 0, check = 1, largemem = 0, chunk_size = 10000):
    """
    BNEWT A balancing algorithm for symmetric matrices
    X = BNEWT(A) attempts to find a vector X such that
    diag(X)*A*diag(X) is close to doubly stochastic. A must
    be symmetric and nonnegative.
    
    X0: initial guess. TOL: error tolerance.
    delta/Delta: how close/far balancing vectors can get
    to/from the edge of the positive cone.
    We use a relative measure on the size of elements.
    FL: intermediate convergence statistics on/off.
    RES: residual error, measured by norm(diag(x)*A*x - e).
    """
    # see details in Knight and Ruiz (2012)
    (n,m) = A.shape
    #np.seterr(divide='ignore')
    print('Verifying Matrix')
    if (n != m):
        print('Matrix must be symmetric to converge')
        return 'NaN'
    if (check):
        for i in range(0,n):
            for j in range(i,n):
                if (A[i][j] != A[j][i])or(A[i][j] < 0):
                    print('Matrix must be symmetric and nonnegative to converge')
                    return 'NaN'
        print('Check OK\n')
    else:
        print('Check escaped\n')
  
    e        = np.ones((n,1))
    e[mask]  = 0
    #res      = np.empty((n,1))
    
    g        = 0.9
    etamax   = 0.1
    eta      = etamax
    stop_tol = tol*0.5
    x        = e #initial guess
    rt       = tol*tol
    if largemem:
        v      = x * chunking_dot(A,x,chunk_size=chunk_size)
    else:
        v      = x*np.dot(A,x)
    rk       = 1 - v
    rk[mask] = 0
    rho_km1  = np.dot(np.transpose(rk),rk)
    rout     = rho_km1
    rold     = rout
    
    MVP = 0 #matrix vector products
    i = 0
  
    while rout > rt:
        i = i+1
        k=0
        y=e
        innertol = max(eta*eta*rout,rt)
    
        while rho_km1 > innertol: #inner iteration by CG
            k = k+1
            if k==1:
                with np.errstate(invalid='ignore'):
                    Z       = rk/v
                Z[mask] = 0
                p       = Z
                rho_km1 = np.dot(np.transpose(rk),Z)
            else:
                beta = rho_km1/rho_km2
                p    =  Z + beta*p
      
            #update search direction 
            if largemem:
                w   = x*chunking_dot(A,x*p,chunk_size=chunk_size) + v*p
            else:
                w   = x*np.dot(A,x*p) + v*p
      
            alpha = rho_km1/np.dot(np.transpose(p),w)
            ap = alpha*p
      
            #test distance to boundary of cone
            ynew = y + ap
            if min(np.delete(ynew,mask)) <= delta_lower:
                if delta_lower == 0:
                    break
                else:
                    ind = np.nonzero(ap < 0)
                    gamma = min((delta_lower - y[ind])/ap[ind])
                    y = y + gamma*ap
                    break
                if max(ynew) >= delta_upper:
                    ind = np.nonzero(ynew > delta_upper)
                    gamma = min((delta_upper-y[ind])/ap[ind])
                    y = y + gamma*ap
                    break
      
            y       = ynew
            rk      = rk - alpha*w
            rho_km2 = rho_km1
            with np.errstate(invalid='ignore'):
                Z       = rk/v
            Z[mask] = 0
            rho_km1 = np.dot(np.transpose(rk),Z)
        #end inner iteration
    
        x        = x*y
        if largemem:
            v      = x * chunking_dot(A,x,chunk_size=chunk_size)
        else:
            v      = x*np.dot(A,x)
        rk       = 1-v
        rk[mask] = 0
        rho_km1  = np.dot(np.transpose(rk),rk)
        rout     = rho_km1
        MVP      = MVP + k + 1
        #print MVP,res
        #update inner iteration stopping criterion
        rat      = rout/rold
        rold     = rout
        res_norm = math.sqrt(rout)
        eta_o    = eta
        eta      = g*rat
    
        if g*eta_o*eta_o > 0.1:
            eta = max(eta,g*eta_o*eta_o)
        eta = max(min(eta,etamax),stop_tol/res_norm)
    
        if fl == 1:
            print('%3d %6d %.3f \n' % (i,k,r_norm))
      
        if MVP > 50000:
            break
    #end outer
  
    print('Matrix vector products = %6d' % MVP)
    #x = np.array(x)
    #x[mask] = 0
    return x

def KRnorm(M, mask=None, **kwargs):
    """
    Uses Knight-Ruiz normalization to balance the matrix.
    """

    x = bnewt(M, mask=mask, check=0, **kwargs)*100
    M *= X 
    M *= X.T

    return M

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
