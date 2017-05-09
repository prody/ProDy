from numpy import ma
import numpy as np
from scipy.sparse import coo_matrix
from prody.chromatin.norm import VCnorm, SQRTVCnorm,Filenorm
from prody.chromatin.cluster import KMeans

from prody.dynamics import GNM
from prody.dynamics.functions import writeArray
from prody.dynamics.mode import Mode
from prody.dynamics.modeset import ModeSet

__all__ = ['HiC', 'parseHiC', 'parseHiCStream', 'plotmap', 'writeHiC', 'writeMap']

class HiC(object):

    """This class is used to store and preprocess Hi-C contact map. A :class:`.GNM`
    instance for analyzing the contact map can be also created by using this class.
    """

    __slots__ = ['_map', 'mask', 'useTrimed', '_title', 'bin', '_labels']

    def __init__(self, title='Unknown', map=None, bin=None):
        self._title = title
        self._map = None
        self.mask = False
        self._labels = 0
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
            self._labels = np.zeros(len(self._map))

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
        return len(self._map)
    
    def numAtoms(self):
        return len(self._map)

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
            D = np.diag(np.sum(A, axis=0))
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
    
    def segment(self, modes, method=KMeans, **kwargs):

        if isinstance(modes, ModeSet):
            V = modes.getEigvecs()
        elif isinstance(modes, Mode):
            V = modes.getEigvec()
        elif isinstance(modes, np.ndarray):
            V = modes
        else:
            try:
                mode0 = modes[0]
                if isinstance(mode0, Mode):
                    V = np.empty((len(mode0),0))
                    for mode in modes:
                        assert isinstance(mode, Mode), 'Modes should be a list of modes.'
                        v = mode.getEigvec()
                        v = np.expand_dims(v, axis=1)
                        V = np.hstack((V, v))
                else:
                    V = np.array(modes)
            except TypeError:
                TypeError('Modes should be a list of modes.')
        if V.ndim == 1:
            V = np.expand_dims(V, axis=1)

        if len(self.Map) != V.shape[0]:
            raise ValueError('Modes (%d) and the Hi-C map (%d) should have the same number'
                             ' of atoms. Turn off "useTrimed" if you intended to apply the'
                             ' modes to the full map.'
                             %(V.shape[0], len(self.Map)))

        labels = method(V, **kwargs)

        if self.useTrimed:
            full_length = len(self._map)
            if full_length != len(labels):
                _labels = np.empty(full_length)
                _labels.fill(np.nan)
                _labels[~self.mask] = labels

                currlbl = labels[np.argmax(~np.isnan(labels))]

                for i in range(len(_labels)):
                    l = _labels[i]
                    if np.isnan(l):
                        _labels[i] = currlbl
                    elif currlbl != l:
                        currlbl = l
                labels = _labels
            
        return labels

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
        if np.max(i) != np.max(j):
            b = np.max(np.append(i, j))
            i = np.append(i, b)
            j = np.append(j, b)
            value = np.append(value, 0.)
        # Convert to sparse matrix format, then full matrix format
        # and finally array type. Matrix format is avoided because
        # diag() won't work as intended for Matrix instances.
        M = np.array(coo_matrix((value, (i,j))).todense())
    return HiC(title=title, map=M, bin=bin)

def writeMap(filename, map, bin=None, format='%f'):
    """Writes *map* to the file designated by *filename*.

    :arg filename: the file to be written.
    :type filename: str

    :arg map: a Hi-C contact map.
    :type map: :class:`numpy.ndarray`

    :arg bin: bin size of the *map*. If bin is `None`, *map* will be 
    written in full matrix format.
    :type bin: int

    :arg format: output format for map elements.
    :type format: str
    """

    assert isinstance(map, np.ndarray), 'map must be a numpy.ndarray.'

    if bin is None:
        return writeArray(filename, map, format=format)
    else:
        L = int(map.size - np.diag(map).size)/2 + np.diag(map).size
        spmat = np.zeros((L,3))
        m,n = map.shape
        l = 0
        for i in range(m):
            for j in range(i,n):
                spmat[l,0] = i * bin
                spmat[l,1] = j * bin
                spmat[l,2] = map[i,j]
                l += 1
        fmt = ['%d', '%d', format]
        return writeArray(filename, spmat, format=fmt)

def writeHiC(filename, hic, format='%f'):
    """Writes *hic.Map* instance to the file designated by *filename*.

    :arg filename: the file to be written.
    :type filename: str

    :arg hic: a HiC instance.
    :type hic: :class:`.HiC`

    :arg format: output format for map elements.
    :type format: str
    """
    assert isinstance(hic, HiC), 'hic must be a HiC instance.'
    return writeMap(filename, hic.Map, hic.bin, format)

def plotmap(map, spec='', p=5, **kwargs):
    """A convenient function to visualize Hi-C contact map. *kwargs* will be passed to 
    :func:`matplotlib.pyplot.imshow`.

    :arg map: a Hi-C contact map.
    :type map: :class:`numpy.ndarray`

    :arg spec: a string specifies how to preprocess the matrix. Blank for no preprocessing,
    'p' for showing only data from *p*-th to *100-p*-th percentile.
    :type spec: str

    :arg p: specifies the percentile threshold.
    :type p: double
    """

    assert isinstance(map, np.ndarray), 'map must be a numpy.ndarray.'
    
    from matplotlib.pyplot import figure, imshow

    if not '_' in spec:
        figure()
    if 'p' in spec:
        vmin = np.percentile(map, p)
        vmax = np.percentile(map, 100-p)
    else:
        vmin = vmax = None
    
    return imshow(map, vmin=vmin, vmax=vmax, **kwargs)
