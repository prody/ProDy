"""This module defines functions for solving linear assignment problems."""

from numpy import arange, insert, delete, zeros, unique, vstack, ceil, argmin, tile
from scipy.optimize import linear_sum_assignment as lap

__all__ = ['multilap', 'SolutionDepletionException']

class SolutionDepletionException(Exception):
    pass

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

    .. [KM68] Murty KG. Letter to the editor-An algorithm for ranking 
       all the assignments in order of increasing cost.
       *Operations research* **1968** 16(3):682-687.
    """

    m, n = cost_matrix.shape
    row_labels = arange(m)
    
    # make duplicates for alternative assignments
    n_copies = int(ceil(float(n) / m))
    if n_copies > 1:
        cost_matrix = vstack((cost_matrix,)*n_copies)
        row_labels = tile(row_labels, n_copies)
    else:
        cost_matrix = cost_matrix.copy()

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
                    r_ = row_labels[r]
                    D[row_labels==r_, c] = BIG_NUMBER
                
                # remove includes
                delR = []; delC = []
                for r, c in includes:
                    delR.append(r)
                    delC.append(c)

                if delC:
                    D = delete(D, delC, axis=1)
                    C = delete(C, delC)
                
                if delR:
                    D = delete(D, delR, axis=0)
                    R = delete(R, delR)
                
                # solve lap
                I, J = lap(D)
                cost = D[I, J].sum()
                
                if includes:
                    cost += cost_matrix[delR, delC].sum()
                row_ind = R[I]; col_ind = C[J]
                # add back includes
                if includes:
                    row_ind = insert(row_ind, 0, delR)
                    col_ind = insert(col_ind, 0, delC)

                assignments = node[-1] = (row_ind, col_ind)
            else:
                R, C = assignments
                cost = cost_matrix[R, C].sum()
            
            costs.append(cost)
        
        # remove assignments with violations. -0.01 is to avoid precision problem
        for i in reversed(range(len(costs))):
            if costs[i] - BIG_NUMBER >= -0.01:
                costs.pop(i)
                nodes.pop(i)

        if len(costs) == 0:
            raise SolutionDepletionException('solution depleted')

        i = argmin(costs)
        node = nodes.pop(i)
        R, C = node[-1]
    else:
        R, C = lap(cost_matrix)
        node = ([], [], (R, C))
        cost = cost_matrix[R, C].sum()

    new_nodes = expand_nodes(node)
    nodes.extend(new_nodes)

    R_ = row_labels[R]

    return (R_, C), gen_mappings((R_, C))

def gen_mappings(assignments):
    from itertools import product as iproduct

    I, J = assignments
    m = len(unique(I))
    M = I.max() + 1

    if m == len(I):
        mappings = [(I, J)]
    else:
        pool = [[] for _ in range(M)]
        for i, j in zip(I, J):
            pool[i].append((i, j))

        mappings = []
        for mapping in iproduct(*pool):
            r = zeros(len(mapping), dtype=int)
            c = zeros(len(mapping), dtype=int)
            for i, pair in enumerate(mapping):
                r[i] = pair[0]
                c[i] = pair[1]
            mappings.append((r, c))
        
    return mappings
