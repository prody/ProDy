# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
#
# This file is originally part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
# 
# The file is adapted by She Zhang (shezhang620@gmail.com) on 04/13/2020 to fix a 
# bug in the construction of the UPGMA tree and made suitable for being used by 
# ProDy.

"""Classes and methods for tree construction."""

import itertools
import copy
import numbers
from Bio.Phylo import BaseTree
from Bio.Phylo.TreeConstruction import _DistanceMatrix as DistanceMatrix

__all__ = ['_Matrix', 'DistanceMatrix', 'TreeConstructor', 'DistanceTreeConstructor']

class _Matrix:
    """Base class for distance matrix or scoring matrix.

    Accepts a list of names and a lower triangular matrix.::

        matrix = [[0],
                  [1, 0],
                  [2, 3, 0],
                  [4, 5, 6, 0]]
        represents the symmetric matrix of
        [0,1,2,4]
        [1,0,3,5]
        [2,3,0,6]
        [4,5,6,0]

    :Parameters:
        names : list
            names of elements, used for indexing
        matrix : list
            nested list of numerical lists in lower triangular format

    Examples
    --------
    >>> from Bio.Phylo.TreeConstruction import _Matrix
    >>> names = ['Alpha', 'Beta', 'Gamma', 'Delta']
    >>> matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]
    >>> m = _Matrix(names, matrix)
    >>> m
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]])

    You can use two indices to get or assign an element in the matrix.

    >>> m[1,2]
    3
    >>> m['Beta','Gamma']
    3
    >>> m['Beta','Gamma'] = 4
    >>> m['Beta','Gamma']
    4

    Further more, you can use one index to get or assign a list of elements related to that index.

    >>> m[0]
    [0, 1, 2, 4]
    >>> m['Alpha']
    [0, 1, 2, 4]
    >>> m['Alpha'] = [0, 7, 8, 9]
    >>> m[0]
    [0, 7, 8, 9]
    >>> m[0,1]
    7

    Also you can delete or insert a column&row of elemets by index.

    >>> m
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])
    >>> del m['Alpha']
    >>> m
    _Matrix(names=['Beta', 'Gamma', 'Delta'], matrix=[[0], [4, 0], [5, 6, 0]])
    >>> m.insert('Alpha', [0, 7, 8, 9] , 0)
    >>> m
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])

    """

    def __init__(self, names, matrix=None):
        """Initialize matrix.

        Arguments are a list of names, and optionally a list of lower
        triangular matrix data (zero matrix used by default).
        """
        # check names
        if isinstance(names, list) and all(isinstance(s, str) for s in names):
            if len(set(names)) == len(names):
                self.names = names
            else:
                raise ValueError("Duplicate names found")
        else:
            raise TypeError("'names' should be a list of strings")

        # check matrix
        if matrix is None:
            # create a new one with 0 if matrix is not assigned
            matrix = [[0] * i for i in range(1, len(self) + 1)]
            self.matrix = matrix
        else:
            # check if all elements are numbers
            if (
                isinstance(matrix, list)
                and all(isinstance(l, list) for l in matrix)
                and all(
                    isinstance(n, numbers.Number)
                    for n in [item for sublist in matrix for item in sublist]
                )
            ):
                # check if the same length with names
                if len(matrix) == len(names):
                    # check if is lower triangle format
                    if [len(m) for m in matrix] == list(range(1, len(self) + 1)):
                        self.matrix = matrix
                    else:
                        raise ValueError("'matrix' should be in lower triangle format")
                else:
                    raise ValueError("'names' and 'matrix' should be the same size")
            else:
                raise TypeError("'matrix' should be a list of numerical lists")

    def __getitem__(self, item):
        """Access value(s) by the index(s) or name(s).

        For a _Matrix object 'dm'::

            dm[i]                   get a value list from the given 'i' to others;
            dm[i, j]                get the value between 'i' and 'j';
            dm['name']              map name to index first
            dm['name1', 'name2']    map name to index first

        """
        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            return [self.matrix[index][i] for i in range(0, index)] + [
                self.matrix[i][index] for i in range(index, len(self))
            ]
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            if row_index > col_index:
                return self.matrix[row_index][col_index]
            else:
                return self.matrix[col_index][row_index]
        else:
            raise TypeError("Invalid index type.")

    def __setitem__(self, item, value):
        """Set value by the index(s) or name(s).

        Similar to __getitem__::

            dm[1] = [1, 0, 3, 4]    set values from '1' to others;
            dm[i, j] = 2            set the value from 'i' to 'j'

        """
        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, list) and all(
                isinstance(n, numbers.Number) for n in value
            ):
                if len(value) == len(self):
                    for i in range(0, index):
                        self.matrix[index][i] = value[i]
                    for i in range(index, len(self)):
                        self.matrix[i][index] = value[i]
                else:
                    raise ValueError("Value not the same size.")
            else:
                raise TypeError("Invalid value type.")
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, numbers.Number):
                if row_index > col_index:
                    self.matrix[row_index][col_index] = value
                else:
                    self.matrix[col_index][row_index] = value
            else:
                raise TypeError("Invalid value type.")
        else:
            raise TypeError("Invalid index type.")

    def __delitem__(self, item):
        """Delete related distances by the index or name."""
        index = None
        if isinstance(item, int):
            index = item
        elif isinstance(item, str):
            index = self.names.index(item)
        else:
            raise TypeError("Invalid index type.")
        # remove distances related to index
        for i in range(index + 1, len(self)):
            del self.matrix[i][index]
        del self.matrix[index]
        # remove name
        del self.names[index]

    def insert(self, name, value, index=None):
        """Insert distances given the name and value.

        :Parameters:
            name : str
                name of a row/col to be inserted
            value : list
                a row/col of values to be inserted

        """
        if isinstance(name, str):
            # insert at the given index or at the end
            if index is None:
                index = len(self)
            if not isinstance(index, int):
                raise TypeError("Invalid index type.")
            # insert name
            self.names.insert(index, name)
            # insert elements of 0, to be assigned
            self.matrix.insert(index, [0] * index)
            for i in range(index, len(self)):
                self.matrix[i].insert(index, 0)
            # assign value
            self[index] = value
        else:
            raise TypeError("Invalid name type.")

    def __len__(self):
        """Matrix length."""
        return len(self.names)

    def __repr__(self):
        """Return Matrix as a string."""
        return self.__class__.__name__ + "(names=%s, matrix=%s)" % tuple(
            map(repr, (self.names, self.matrix))
        )

    def __str__(self):
        """Get a lower triangular matrix string."""
        matrix_string = "\n".join(
            [
                self.names[i] + "\t" + "\t".join([str(n) for n in self.matrix[i]])
                for i in range(0, len(self))
            ]
        )
        matrix_string = matrix_string + "\n\t" + "\t".join(self.names)
        return matrix_string

class TreeConstructor:
    """Base class for all tree constructor."""

    def build_tree(self, msa):
        """Caller to built the tree from a MultipleSeqAlignment object.

        This should be implemented in subclass.
        """
        raise NotImplementedError("Method not implemented!")


class DistanceTreeConstructor(TreeConstructor):
    """Distance based tree constructor.

    :Parameters:
        method : str
            Distance tree construction method, 'nj'(default) or 'upgma'.

    Examples
    --------
    Loading a small PHYLIP alignment from which to compute distances, and then
    build a upgma Tree::

        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
        from Bio.Phylo.TreeConstruction import DistanceCalculator
        from Bio import AlignIO
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        constructor = DistanceTreeConstructor()
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        upgmatree = constructor.upgma(dm)
        print(upgmatree)

    Output::

        Tree(rooted=True)
            Clade(branch_length=0, name='Inner4')
                Clade(branch_length=0.18749999999999994, name='Inner1')
                    Clade(branch_length=0.07692307692307693, name='Epsilon')
                    Clade(branch_length=0.07692307692307693, name='Delta')
                Clade(branch_length=0.11057692307692304, name='Inner3')
                    Clade(branch_length=0.038461538461538464, name='Inner2')
                        Clade(branch_length=0.11538461538461536, name='Gamma')
                        Clade(branch_length=0.11538461538461536, name='Beta')
                    Clade(branch_length=0.15384615384615383, name='Alpha')

    Build a NJ Tree::

        njtree = constructor.nj(dm)
        print(njtree)

    Output::

        Tree(rooted=False)
            Clade(branch_length=0, name='Inner3')
                Clade(branch_length=0.18269230769230765, name='Alpha')
                Clade(branch_length=0.04807692307692307, name='Beta')
                Clade(branch_length=0.04807692307692307, name='Inner2')
                    Clade(branch_length=0.27884615384615385, name='Inner1')
                        Clade(branch_length=0.051282051282051266, name='Epsilon')
                        Clade(branch_length=0.10256410256410259, name='Delta')
                    Clade(branch_length=0.14423076923076922, name='Gamma')

    """

    methods = ["nj", "upgma"]

    def __init__(self, method="nj"):
        """Initialize the class."""
        if isinstance(method, str) and method in self.methods:
            self.method = method
        else:
            raise TypeError(
                "Bad method: "
                + method
                + ". Available methods: "
                + ", ".join(self.methods)
            )

    def upgma(self, distance_matrix):
        """Construct and return an UPGMA tree.

        Constructs and returns an Unweighted Pair Group Method
        with Arithmetic mean (UPGMA) tree.

        :Parameters:
            distance_matrix : DistanceMatrix
                The distance matrix for tree construction.

        """
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Must provide a DistanceMatrix object.")

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        while len(dm) > 1:
            min_dist = dm[1, 0]
            # find minimum index
            for i in range(1, len(dm)):
                for j in range(0, i):
                    if min_dist >= dm[i, j]:
                        min_dist = dm[i, j]
                        min_i = i
                        min_j = j

            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            clade1.branch_length = min_dist * 1.0 / 2 - self._height_of(clade1)
            clade2.branch_length = min_dist * 1.0 / 2 - self._height_of(clade2)

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]) * 1.0 / 2

            dm.names[min_j] = "Inner" + str(inner_count)

            del dm[min_i]
        inner_clade.branch_length = 0
        return BaseTree.Tree(inner_clade)

    def nj(self, distance_matrix):
        """Construct and return a Neighbor Joining tree.

        :Parameters:
            distance_matrix : DistanceMatrix
                The distance matrix for tree construction.

        """
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Must provide a DistanceMatrix object.")

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init node distance
        node_dist = [0] * len(dm)
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        # special cases for Minimum Alignment Matrices
        if len(dm) == 1:
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)
        elif len(dm) == 2:
            # minimum distance will always be [1,0]
            min_i = 1
            min_j = 0
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            clade1.branch_length = dm[min_i, min_j] / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
            inner_clade = BaseTree.Clade(None, "Inner")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            clades[0] = inner_clade
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)
        while len(dm) > 2:
            # calculate nodeDist
            for i in range(0, len(dm)):
                node_dist[i] = 0
                for j in range(0, len(dm)):
                    node_dist[i] += dm[i, j]
                node_dist[i] = node_dist[i] / (len(dm) - 2)

            # find minimum distance pair
            min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
            min_i = 0
            min_j = 1
            for i in range(1, len(dm)):
                for j in range(0, i):
                    temp = dm[i, j] - node_dist[i] - node_dist[j]
                    if min_dist > temp:
                        min_dist = temp
                        min_i = i
                        min_j = j
            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            clade1.branch_length = (
                dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]
            ) / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (
                        dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]
                    ) / 2.0

            dm.names[min_j] = "Inner" + str(inner_count)
            del dm[min_i]

        # set the last clade as one of the child of the inner_clade
        root = None
        if clades[0] == inner_clade:
            clades[0].branch_length = 0
            clades[1].branch_length = dm[1, 0]
            clades[0].clades.append(clades[1])
            root = clades[0]
        else:
            clades[0].branch_length = dm[1, 0]
            clades[1].branch_length = 0
            clades[1].clades.append(clades[0])
            root = clades[1]

        return BaseTree.Tree(root, rooted=False)

    def _height_of(self, clade):
        """Calculate clade height -- the longest path to any terminal (PRIVATE)."""
        if clade.is_terminal():
            height = 0 
        else:
            height = max(self._height_of(c) + c.branch_length for c in clade.clades)

        return height

