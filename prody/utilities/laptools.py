"""This module defines functions for solving linear assignment problems."""

from numpy import arange, insert, delete, zeros, unique, vstack, ceil, argmin
from scipy.optimize import linear_sum_assignment as lap

__all__ = ['multilap']

def expand_nodes(node):
    includes0, excludes0, solution = node

    pairs = []
    R, C = solution
    for r, c in zip(R, C):
        pairs.append((r, c))

    nodes = []; previous_pairs = []
    for i, pair in enumerate(pairs):
        if i == len(pairs) - 1:
            continue
        if pair in includes0:
            continue
        if pair in excludes0: # unlikely
            raise ValueError('%s should be excluded'%str(pair))
        
        includes = list(includes0)
        includes.extend(previous_pairs)

        excludes = list(excludes0)
        excludes.append(pair)

        nodes.append([includes, excludes, None])
        previous_pairs.append(pair)

    return nodes
    
def multilap(cost_matrix, nodes=[], BIG_NUMBER=1e6):
    """ Finds the (next) optimal solution to the linear assignment problem. 
    The function can handle the cases where each row can be assigned to multiple 
    columns. 

    .. [KM68] Murty KG. Letter to the editor—An algorithm for ranking 
       all the assignments in order of increasing cost.
       *Operations research* **1968** 16(3):682-687.
    """

    if len(nodes):
        M, N = cost_matrix.shape

        costs = []
        for node in nodes:
            includes, excludes, assignments = node
            if assignments is None:
                D = cost_matrix.copy()
                R = arange(M)
                C = arange(N)

                # mask excludes
                for r, c in excludes:
                    D[r, c] = BIG_NUMBER
                
                # remove includes
                delR = []; delC = []
                for r, c in includes:
                    delR.append(r)
                    delC.append(c)

                if delC:
                    D = delete(D, delC, axis=1)
                    C = delete(C, delC)
                m, n = D.shape
                if m > n:
                    if delR:
                        D = delete(D, delR, axis=0)
                        R = delete(R, delR)
                
                # solve lap
                I, J = multilap_solver(D)
                cost_ = D[I, J].sum()
                # check for existance of excludes. -0.01 is to avoid precision problem
                if cost_ - BIG_NUMBER >= -0.01:  
                    return None
                row_ind = R[I]; col_ind = C[J]
                # add back includes
                row_ind = insert(row_ind, 0, delR)
                col_ind = insert(col_ind, 0, delC)

                assignments = node[-1] = (row_ind, col_ind)
            
            R, C = assignments
            cost = cost_matrix[R, C].sum()
            costs.append(cost)

        i = argmin(costs)
        cost = costs[i]
        node = nodes.pop(i)
        assignments = node[-1]
    else:
        assignments = multilap_solver(cost_matrix)
        node = ([], [], assignments)

        R, C = assignments
        cost = cost_matrix[R, C].sum()

    new_nodes = expand_nodes(node)
    nodes.extend(new_nodes)

    return gen_mappings(assignments)

def gen_mappings(assignments):
    from itertools import product as iproduct

    I, J = assignments
    m = len(unique(I))

    if m == len(I):
        mappings = [(I, J)]
    else:
        pool = [[] for _ in range(m)]
        for i, j in zip(I, J):
            pool[i].append((i, j))
        
        mappings = []
        for mapping in iproduct(*pool):
            r = zeros(m, dtype=int); c = zeros(m, dtype=int)
            for i, pair in enumerate(mapping):
                r[i] = pair[0]
                c[i] = pair[1]
            mappings.append((r, c))
        
    return mappings

def multilap_solver(cost_matrix): 
    from scipy.optimize import linear_sum_assignment as lap

    m, n = cost_matrix.shape
    
    # make duplicates for alternative assignments
    n_copies = int(ceil(float(n) / m))
    if n_copies > 1:
        D = vstack((cost_matrix,)*n_copies)
    else:
        D = cost_matrix.copy()
    
    # solve lap
    I, J = lap(D)

    # format outputs
    if n_copies > 1:
        # map combinations of assignments
        for idx, i_ in enumerate(I):
            i = i_ % m
            I[idx] = i

    return I, J
