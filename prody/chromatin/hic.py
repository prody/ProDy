from numpy import ma
import numpy as np
from scipy.sparse import coo_matrix
from prody.chromatin.norm import VCnorm, SQRTVCnorm,Filenorm

from prody.dynamics import GNM

__all__ = ['HiC', 'parseHiC', 'parseHiCStream', 'plotmat']

class HiC(object):

    """This class is used to store and preprocess Hi-C contact map. A :class:`.GNM`
    instance for analyzing the contact map can be also created by using this class.
    """

    __slots__ = ['_map', 'mask', 'useTrimed', '_title', 'bin']

    def __init__(self, title='Unknown', map=None, bin=None):
        self._title = title
        self._map = None
        self.mask = False
        self.useTrimed = True
        self.bin = bin
        self.Map = map
    
    @property
    def Map(self):
        if self.useTrimed:
            return self.getTrimedMap()
        else:
            return self._map

    @Map.setter
    def Map(self, map):
        if map is None: 
            self._map = None
        else:
            self._map = np.array(map)
            self._makeSymmetric()
            self._maskUnmappedRegions()

    def __repr__(self):

        return '<HiC: {0} ({1} mapped loci; {2} in total)>'.format(self._title, len(self.getTrimedMap()), len(self._map))

    def __str__(self):

        return 'HiC ' + self._title

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.Map.flatten()[index]
        else:
            i, j = index
            return self.Map[i,j]

    def __len__(self):

        return len(self.Map)
    
    def getTitle(self):
        """Returns title of the instance."""

        return self._title

    def setTitle(self, title):
        """Sets title of the instance."""

        self._title = str(title)

    def getCompleteMap(self):
        """Obtains the complete contact map with unmapped regions."""
        
        return self._map
        
    def getTrimedMap(self):
        """Obtains the contact map without unmapped regions."""

        if self._map is None: 
            return None
        if self.mask is False:
            return self._map

        M = ma.array(self._map)
        M.mask = np.diag(self.mask)
        return ma.compress_rowcols(M)
    
    def getKirchhoff(self):
        """Builds a Kirchhoff matrix based on the contact map."""

        if self.Map is None:
            return None
        else:
            M = self.Map
            
            I = np.eye(M.shape[0], dtype=bool)
            A = M.copy()
            A[I] = 0.
            D = np.diag(sum(A, axis=0))
            K = D - A
            return K

    def _maskUnmappedRegions(self):
        """Finds and masks unmapped regions in the contact map."""

        M = self._map
        if M is None: return
        # Obtain the diagonal values, need to make sure d is an array 
        # instead of a matrix, otherwise diag() later will not work as 
        # intended.
        d = np.array(np.diag(M))
        # mask if a diagonal value is zero
        mask_zero = np.array(d==0)
        # mask if a diagonal value is NAN
        mask_nan = np.isnan(d)
        # combine two masks
        mask = np.logical_or(mask_nan, mask_zero)
        self.mask = mask
        return mask

    def _makeSymmetric(self):
        """Ensures the symmetricity of the contact map."""

        M = self._map
        if M is None: return
        
        # determine which part of the matrix has values
        U = np.triu(M, k=1)
        L = np.tril(M, k=-1)

        if np.sum(U) == 0:
            M += L.T
        elif np.sum(L) == 0:
            M += U.T
        else:
            M = (M + M.T) / 2.
        return self._map
    
    def calcGNM(self, n_modes=None):
        """Calculates GNM on the current Hi-C map."""

        gnm = GNM(self._title)
        gnm.setKirchhoff(self.getKirchhoff())
        gnm.calcModes(n_modes=n_modes)

        return gnm
    
    def normalize(self, method=VCnorm, **kwargs):
        """Applies chosen normalization on the current Hi-C map."""

        M = self._map
        N = method(M, **kwargs)
        self.Map = N
        return N

def parseHiC(filename, **kwargs):
    """Returns an :class:`.HiC` from a Hi-C data file.

    This function extends :func:`.parseHiCStream`.

    :arg filename: the filename to the Hi-C data file.
    :type filename: str
    """

    import os
    title = kwargs.get('title')
    if title is None:
        title = os.path.basename(filename)
    with open(filename, 'rb') as filestream:
        hic = parseHiCStream(filestream, title=title, **kwargs)
    return hic

def parseHiCStream(stream, **kwargs):
    """Returns an :class:`.HiC` from a stream of Hi-C data lines.

    :arg stream: Anything that implements the method ``read``, ``seek``
        (e.g. :class:`file`, buffer, stdin)
    """

    title = kwargs.get('title', 'Unknown')

    import csv
    dialect = csv.Sniffer().sniff(stream.read(1024))
    stream.seek(0)
    reader = csv.reader(stream, dialect)
    D = list()
    for row in reader:
        d = list()
        for element in row:
            d.append(np.double(element))
        D.append(d)
    D = np.array(D)

    bin = None
    size = D.shape
    if len(D.shape) <= 1:
        raise ValueError("Cannot parse the file: input file only contains one column.")
    if size[0] == size[1]:
        bin = None
        M = D
    else:
        i, j, value = D.T
        # determine the bin size by the first and second locus
        bins = np.unique(np.sort(i))
        bin = bins[1] - bins[0]
        # convert coordinate from basepair to locus index
        i = i//bin
        j = j//bin
        # make sure that the matrix is square
        if max(i) != max(j):
            b = max(append(i, j))
            i = append(i, b)
            j = append(j, b)
            value = append(value, 0.)
        # Convert to sparse matrix format, then full matrix format
        # and finally array type. Matrix format is avoided because
        # diag() won't work as intended for Matrix instances.
        M = np.array(coo_matrix((value, (i,j))).todense())
    return HiC(title=title, map=M, bin=bin)

def plotmat(mat, spec='', p = 5, **kwargs):
    """A convenient function to visualize Hi-C contact map. *kwargs* will be passed to 
    :func:`matplotlib.pyplot.imshow`.

    :arg mat: a matrix.
    :type mat: :class:`numpy.ndarray`

    :arg spec: a string specifies how to preprocess the matrix. Blank for no preprocessing,
    'p' for showing only data from *p*-th to *100-p*-th percentile.
    :type spec: str

    :arg p: specifies the percentile threshold.
    :type p: double
    """

    from matplotlib.pyplot import figure, imshow

    if not '_' in spec:
        figure()
    if 'p' in spec:
        vmin = np.percentile(mat, p)
        vmax = np.percentile(mat, 100-p)
    else:
        vmin = vmax = None
    
    return imshow(mat, vmin=vmin, vmax=vmax, **kwargs)
