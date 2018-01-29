from numpy import ma
import numpy as np
from scipy.sparse import coo_matrix
from scipy.stats import mode
from prody.chromatin.norm import VCnorm, SQRTVCnorm, Filenorm
from prody.chromatin.cluster import KMeans, Hierarchy
from prody.chromatin.functions import div0, showMap, showDomains, _getEigvecs

from prody.dynamics import GNM, TrimedGNM
from prody.dynamics.functions import writeArray
from prody.dynamics.mode import Mode
from prody.dynamics.modeset import ModeSet

from prody.utilities import openFile, importLA

__all__ = ['HiC', 'TrimedGNM', 'parseHiC', 'parseHiCStream', 'saveHiC', 'loadHiC', 'writeMap']

class HiC(object):

    """This class is used to store and preprocess Hi-C contact map. A :class:`.GNM`
    instance for analyzing the contact map can be also created by using this class.
    """

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
        if np.isscalar(self.mask):
            return self._map

        M = ma.array(self._map)
        M.mask = np.diag(~self.mask)
        return ma.compress_rowcols(M)
    
    def align(self, array, axis=None):
        if not isinstance(map, np.ndarray):
            array = np.array(array)

        ret = array = array.copy()

        if np.isscalar(self.mask):
            return ret

        mask = self.mask.copy()

        l_full = self.getCompleteMap().shape[0]
        l_trim = self.getTrimedMap().shape[0]
        
        if len(array.shape) == 0:
            raise ValueError('Aligned array cannot be empty.')
        elif len(array.shape) == 1:
            l = array.shape[0]
            if l == l_trim:
                N = len(mask)
                ret = np.zeros(N)
                ret[mask] = array
            elif l == l_full:
                ret = array[mask]
            else:
                raise ValueError('The length of the array (%d) does not '
                                'match that of either the full (%d) '
                                'or trimed (%d).'
                                %(l, l_full, l_trim))
        elif len(array.shape) == 2:
            s = array.shape

            if axis is None:
                if s[0] != s[1]:
                    raise ValueError('The array must be a square matrix '
                                     'if axis is set to None.')
                if s[0] == l_trim:
                    N = len(mask)
                    whole_mat = np.zeros((N,N))
                    mask = np.outer(mask, mask)
                    whole_mat[mask] = array.flatten()
                    ret = whole_mat
                elif s[0] == l_full:
                    M = ma.array(array)
                    M.mask = np.diag(mask)
                    ret = ma.compress_rowcols(M)
                else:
                    raise ValueError('The size of the array (%d) does not '
                                    'match that of either the full (%d) '
                                    'or trimed (%d).'
                                    %(s[0], l_full, l_trim))
            else:
                new_shape = list(s)
                otheraxis = 0 if axis!=0 else 1
                if s[axis] == l_trim:
                    N = len(mask)
                    new_shape[axis] = N
                    whole_mat = np.zeros(new_shape)
                    mask = np.expand_dims(mask, axis=otheraxis)
                    mask = mask.repeat(s[otheraxis], axis=otheraxis)
                    whole_mat[mask] = array.flatten()
                    ret = whole_mat
                elif s[axis] == l_full:
                    mask = np.expand_dims(mask, axis=otheraxis)
                    mask = mask.repeat(s[otheraxis])
                    ret = map[mask]
                else:
                    raise ValueError('The size of the array (%d) does not '
                                    'match that of either the full (%d) '
                                    'or trimed (%d).'
                                    %(sh[0], l_full, l_trim))
        
        return ret

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
        self.mask = ~mask
        return self.mask

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

        if self.useTrimed:
            gnm = TrimedGNM(self._title, self.mask)
        else:
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
    
    def segment(self, modes, method=Hierarchy, **kwargs):
        """Uses spectral clustering to identify structural domains on the chromosome.
        
        :arg modes: GNM modes used for segmentation
        :type modes: :class:`ModeSet`

        :arg method: Label assignment algorithm used after Laplacian embedding.
        :type method: func
        """

        V = _getEigvecs(modes, True)

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
                _labels[self.mask] = labels

                currlbl = labels[np.argmax(~np.isnan(labels))]

                for i in range(len(_labels)):
                    l = _labels[i]
                    if np.isnan(l):
                        _labels[i] = currlbl
                    elif currlbl != l:
                        currlbl = l
                labels = _labels
        
        self._labels = labels
        return self.getDomains()
    
    def getDomains(self):
        """Returns an 1D :class:`numpy.ndarray` whose length is the number of loci. Each 
        element is an index denotes to which domain the locus belongs."""

        lbl = self._labels
        mask = self.mask
        if self.useTrimed:
            lbl = lbl[mask]
        return lbl

    def getDomainList(self):
        """Returns a list of domain separations. The list has two columns: the first is for 
        the domain starts and the second is for the domain ends."""

        indicators = np.diff(self.getDomains())
        indicators = np.append(1., indicators)
        indicators[-1] = 1
        sites = np.where(indicators != 0)[0]
        starts = sites[:-1]
        ends = sites[1:]
        domains = np.array([starts, ends]).T

        return domains

    def view(self, spec='p', **kwargs):
        """Visualization of the Hi-C map and domains (if present). The function makes use 
        of :func:`.showMap`."""

        dm_kwargs = {}
        keys = kwargs.keys()
        for k in keys:
            if k.startswith('dm_'):
                dm_kwargs[k[3:]] = kwargs.pop(k)
            elif k.startswith('domain_'):
                dm_kwargs[k[7:]] = kwargs.pop(k)

        im = showMap(self.Map, spec, **kwargs)

        domains = self.getDomainList()
        if len(domains) > 1:
            showDomains(domains, **dm_kwargs)

        return im

    def copy(self):
        new = type(self)()
        new.__dict__.update(self.__dict__)
        return new
    
    __copy__ = copy


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
    else:
        title = kwargs.pop('title')
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

    bin = kwargs.get('bin', None)
    size = D.shape
    if len(D.shape) <= 1:
        raise ValueError("Cannot parse the file: input file only contains one column.")
    if size[0] == size[1]:
        M = D
    else:
        i, j, value = D.T
        # determine the bin size by the most frequent interval
        if bin is None:
            loci = np.unique(np.sort(i))
            bins = np.diff(loci)
            bin = mode(bins)[0][0]
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

def saveHiC(hic, filename=None, map=True, **kwargs):
    """Saves *HiC* model data as :file:`filename.hic.npz`. If *map* is ``False``, 
    Hi-C contact map will not be saved and it can be loaded from raw data file 
    later. If *filename* is ``None``, name of the Hi-C instance will be used as 
    the filename, after ``" "`` (white spaces) in the name are replaced with 
    ``"_"`` (underscores). Upon successful completion of saving, filename is 
    returned. This function makes use of :func:`numpy.savez` function."""

    assert isinstance(hic, HiC), 'hic must be a HiC instance.'
    
    if filename is None:
        filename = hic.getTitle().replace(' ', '_')
    
    if filename.endswith('.hic'):
        filename += '.npz'
    elif not filename.endswith('.hic.npz'):
        filename += '.hic.npz'

    attr_dict = hic.__dict__.copy()
    if not map:
        attr_dict.pop('_map')

    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()

    return filename

def loadHiC(filename):
    """Returns HiC instance after loading it from file (*filename*).
    This function makes use of :func:`numpy.load` function. See also 
    :func:`saveHiC`."""
    
    attr_dict = np.load(filename)
    hic = HiC()

    keys = attr_dict.keys()

    for k in keys:
        val = attr_dict[k]
        if len(val.shape) == 0:
            val = np.asscalar(val)
        setattr(hic, k, val)
    return hic

