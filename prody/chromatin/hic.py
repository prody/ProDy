from numbers import Integral

from numpy import ma
import numpy as np
from scipy.sparse import coo_matrix
from scipy.stats import mode
from prody.chromatin.norm import VCnorm, SQRTVCnorm, Filenorm
from prody.chromatin.functions import div0, showDomains, _getEigvecs

from prody import PY2K
from prody.dynamics import GNM, MaskedGNM
from prody.dynamics.functions import writeArray
from prody.dynamics.mode import Mode
from prody.dynamics.modeset import ModeSet

from prody.utilities import openFile, importLA, showMatrix, isURL, fixArraySize, makeSymmetric

__all__ = ['HiC', 'parseHiC', 'parseHiCStream', 'parseHiCBinary', 'saveHiC', 'loadHiC', 'writeMap']

class HiC(object):

    """This class is used to store and preprocess Hi-C contact map. A :class:`.GNM`
    instance for analyzing the contact map can be also created by using this class.
    """

    def __init__(self, title='Unknown', map=None, bin=None):
        self._title = title
        self._map = None
        self.mask = False
        self._labels = 0
        self.masked = True
        self.bin = bin
        self.map = map
    
    @property
    def map(self):
        if self.masked:
            return self.getTrimedMap()
        else:
            return self._map

    @map.setter
    def map(self, value):
        if value is None: 
            self._map = None
        else:
            self._map = np.asarray(value)
            self._map = makeSymmetric(self._map)
            self._maskUnmappedRegions()
            self._labels = np.zeros(len(self._map), dtype=int)

    def __repr__(self):
        mask = self.mask
        
        if np.isscalar(mask):
            return '<HiC: {0} ({1} loci)>'.format(self._title, len(self._map))
        else:
            return '<HiC: {0} ({1} mapped loci; {2} in total)>'.format(self._title, np.count_nonzero(mask), len(self._map))

    def __str__(self):

        return 'HiC ' + self._title

    def __getitem__(self, index):
        if isinstance(index, Integral):
            return self.map.flatten()[index]
        else:
            i, j = index
            return self.map[i,j]

    def __len__(self):
        mask = self.mask 
        
        if np.isscalar(mask):
            return len(self._map)
        else:
            return np.count_nonzero(mask)
    
    def numAtoms(self):
        return len(self.map)

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
        if not isinstance(array, np.ndarray):
            array = np.array(array)

        ret = array = array.copy()

        if np.isscalar(self.mask):
            return ret

        mask = self.mask.copy()

        l_full = self.getCompleteMap().shape[0]
        l_trim = self.getTrimedMap().shape[0]
        
        if len(array.shape) == 0:
            raise ValueError('array cannot be empty')
        elif len(array.shape) == 1:
            l = array.shape[0]
            if l == l_trim:
                N = len(mask)
                ret = np.zeros(N, dtype=array.dtype)
                ret[mask] = array
            elif l == l_full:
                ret = array[mask]
            else:
                raise ValueError('The length of array (%d) does not '
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
                    whole_mat = np.zeros((N,N), dtype=array.dtype)
                    mask = np.outer(mask, mask)
                    whole_mat[mask] = array.flatten()
                    ret = whole_mat
                elif s[0] == l_full:
                    M = ma.array(array)
                    M.mask = np.diag(mask)
                    ret = ma.compress_rowcols(M)
                else:
                    raise ValueError('The size of array (%d) does not '
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
                    ret = self._map[mask]
                else:
                    raise ValueError('The size of array (%d) does not '
                                    'match that of either the full (%d) '
                                    'or trimed (%d).'
                                    %(s[0], l_full, l_trim))
        
        return ret

    def getKirchhoff(self):
        """Builds a Kirchhoff matrix based on the contact map."""

        if self._map is None:
            return None
        else:
            M = self.map
            
            I = np.eye(M.shape[0], dtype=bool)
            A = M.copy()
            A[I] = 0.
            D = np.diag(np.sum(A, axis=0))
            K = D - A
            return K

    def _maskUnmappedRegions(self, diag=False):
        """Finds and masks unmapped regions in the contact map."""

        M = self._map
        if M is None: return

        if diag:
            # Obtain the diagonal values, need to make sure d is an array 
            # instead of a matrix, otherwise diag() later will not work as 
            # intended.
            d = np.array(np.diag(M))
        else:
            d = np.array(M.sum(0))

        # mask if a diagonal value is zero
        mask_zero = np.array(d==0)
        # mask if a diagonal value is NAN
        mask_nan = np.isnan(d)
        # combine two masks
        mask = np.logical_or(mask_nan, mask_zero)
        self.mask = ~mask

        return self.mask
    
    def calcGNM(self, n_modes=None, **kwargs):
        """Calculates GNM on the current Hi-C map. By default, ``n_modes`` is 
        set to **None** and ``zeros`` to **True**."""

        if 'zeros' not in kwargs:
            kwargs['zeros'] = True
            
        if self.masked:
            gnm = MaskedGNM(self._title, self.mask)
        else:
            gnm = GNM(self._title)
        gnm.setKirchhoff(self.getKirchhoff())
        gnm.calcModes(n_modes=n_modes, **kwargs)
        return gnm
    
    def normalize(self, method=VCnorm, **kwargs):
        """Applies chosen normalization on the current Hi-C map."""

        M = self._map
        N = method(M, **kwargs)
        self.map = N
        return N
    
    def setDomains(self, labels, **kwargs):
        """Uses spectral clustering to identify structural domains on the chromosome.
        
        :arg labels: domain labels
        :type labels: :class:`~numpy.ndarray`, list

        :arg method: Label assignment algorithm used after Laplacian embedding.
        :type method: func
        """
        wastrimmed = self.masked

        self.masked = True
        if len(labels) == self.numAtoms():
            full_length = self.numAtoms()
            if full_length != len(labels):
                _labels = np.empty(full_length)
                _labels.fill(np.nan)
                _labels[self.mask] = labels

                currlbl = labels[0]

                for i in range(len(_labels)):
                    l = _labels[i]
                    if np.isnan(l):
                        _labels[i] = currlbl
                    elif currlbl != l:
                        currlbl = l
                labels = _labels
        else:
            self.masked = False
            if len(labels) != self.numAtoms():
                raise ValueError('The length of the labels should match either the length '
                                 'of masked or complete Hi-C map. Turn off "masked" if '
                                 'you intended to set the labels to the full map.')
        
        self.masked = wastrimmed
        self._labels = labels
        return self.getDomains()
    
    def getDomains(self):
        """Returns an 1D :class:`numpy.ndarray` whose length is the number of loci. Each 
        element is an index denotes to which domain the locus belongs."""

        lbl = self._labels
        mask = self.mask
        if self.masked:
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
        of :func:`.showMatrix`.
        
        :arg spec: a string specifies how to preprocess the matrix. Blank for no preprocessing,
                'p' for showing only data from *p*-th to *100-p*-th percentile. '_' is to suppress 
                creating a new figure and paint to the current one instead. The letter specifications 
                can be applied sequentially, e.g. 'p_'.
        :type spec: str

        :arg p: specifies the percentile threshold.
        :type p: double
        """

        dm_kwargs = {}
        keys = list(kwargs.keys())
        for k in keys:
            if k.startswith('dm_'):
                dm_kwargs[k[3:]] = kwargs.pop(k)
            elif k.startswith('domain_'):
                dm_kwargs[k[7:]] = kwargs.pop(k)

        M = self.map
        if 'p' in spec:
            p = kwargs.pop('p', 5)
            lp = kwargs.pop('lp', p)
            hp = kwargs.pop('hp', 100-p)
            vmin = np.percentile(M, lp)
            vmax = np.percentile(M, hp)
        else:
            vmin = vmax = None

        if not 'vmin' in kwargs:
            kwargs['vmin'] = vmin
        if not 'vmax' in kwargs:
            kwargs['vmax'] = vmax

        im = showMatrix(M, **kwargs)

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

    import os, struct
    title = kwargs.get('title')
    if title is None:
        title = os.path.basename(filename)
    else:
        title = kwargs.pop('title')

    if isURL(filename):
        M, res = parseHiCBinary(filename, title=title, **kwargs)
    else:
        with open(filename,'rb') as req:
            magic_number = struct.unpack('<3s',req.read(3))[0]
        if magic_number == b"HIC":
            M, res = parseHiCBinary(filename, title=title, **kwargs)
        else:
            with open(filename, 'r') as filestream:
                M, res = parseHiCStream(filestream, title=title, **kwargs)
    
    hic = HiC(title=title, map=M, bin=res)

    return hic

def _sparse2dense(I, J, values, bin=None):
    I = np.asarray(I, dtype=int)
    J = np.asarray(J, dtype=int)
    values = np.asarray(values, dtype=float)
    # determine the bin size by the most frequent interval
    if bin is None:
        loci = np.unique(np.sort(I))
        bins = np.diff(loci)
        bin = mode(bins)[0][0]
    # convert coordinate from basepair to locus index
    bin = int(bin)
    I = I // bin
    J = J // bin
    # make sure that the matrix is square
    # if np.max(I) != np.max(J):
    #     b = np.max(np.append(I, J))
    #     I = np.append(I, b)
    #     J = np.append(J, b)
    #     values = np.append(values, 0.)
    # Convert to sparse matrix format, then full matrix format
    # and finally array type. Matrix format is avoided because
    # diag() won't work as intended for Matrix instances.
    M = np.array(coo_matrix((values, (I, J))).todense())
    return M, bin

def parseHiCStream(stream, **kwargs):
    """Returns an :class:`.HiC` from a stream of Hi-C data lines.

    :arg stream: Anything that implements the method ``read``, ``seek``
        (e.g. :class:`file`, buffer, stdin)
    """

    issparse = kwargs.get('sparse', None)

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

    res = kwargs.get('bin', None)
    if res is not None:
        res = int(res)
    size = D.shape
    if len(D.shape) <= 1:
        raise ValueError("cannot parse the file: input file only contains one column.")
    
    if issparse is None:
        issparse = size[1] == 3

    if not issparse:
        M = D
    else:
        try:
            I, J, values = D.T[:3]
        except ValueError:
            raise ValueError('the sparse matrix format should have three columns')
        
        M, res = _sparse2dense(I, J, values, bin=res)
    return M, res

def parseHiCBinary(filename, **kwargs):

    chrloc = kwargs.get('chrom', None)
    if chrloc is None:
        raise ValueError('chrom needs to be specified when parsing .hic format')
    chrloc1 = kwargs.get('chrom1', chrloc)
    chrloc2 = kwargs.get('chrom2', chrloc)
    norm = kwargs.get('norm', 'NONE')
    unit = kwargs.get('unit', 'BP')
    res = kwargs.get('binsize', None)
    res = kwargs.get('bin', res)
    if res is None:
        raise ValueError('bin needs to be specified when parsing .hic format')
    res = int(res)

    from .straw import straw
    result = straw(norm, filename, chrloc1, chrloc2, unit, res)

    M, res = _sparse2dense(*result, bin=res)
    return M, res

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
        L = int(map.size - np.diag(map).size)//2 + np.diag(map).size
        spmat = np.zeros((L, 3))
        m,n = map.shape
        l = 0
        for i in range(m):
            for j in range(i,n):
                spmat[l, 0] = i * bin
                spmat[l, 1] = j * bin
                spmat[l, 2] = map[i, j]
                l += 1
        fmt = ['%d', '%d', format]
        return writeArray(filename, spmat, format=fmt)

def saveHiC(hic, filename=None, map=True, **kwargs):
    """Saves *HiC* model data as :file:`filename.hic.npz`. If *map* is **True**, 
    Hi-C contact map will not be saved and it can be loaded from raw data file 
    later. If *filename* is **None**, name of the Hi-C instance will be used as 
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
    np.savez_compressed(ostream, **attr_dict)
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

def saveHiC_h5(hic, filename=None, **kwargs):
    """Saves *HiC* model data as :file:`filename.hic.npz`. If *filename* is 
    **None**, name of the Hi-C instance will be used as 
    the filename, after ``" "`` (white spaces) in the name are replaced with 
    ``"_"`` (underscores). Upon successful completion of saving, filename is 
    returned. This function makes use of :func:`numpy.savez` function."""

    try:
        import h5py
    except:
        raise ImportError('h5py needs to be installed for using this function')

    assert isinstance(hic, HiC), 'hic must be a HiC instance.'
    
    if filename is None:
        filename = hic.getTitle().replace(' ', '_')
    
    if filename.endswith('.hic'):
        filename += '.hic'
    elif not filename.endswith('.hic.h5'):
        filename += '.hic.h5'

    attr_dict = hic.__dict__.copy()

    with h5py.File(filename, 'w') as f:
        for key in attr_dict:
            value = attr_dict[key]
            compression = None if np.isscalar(value) else 'gzip'
            f.create_dataset(key, data=value, compression=compression)

    return filename

def loadHiC_h5(filename):
    """Returns HiC instance after loading it from file (*filename*).
    This function makes use of :func:`numpy.load` function. See also 
    :func:`saveHiC`."""
    
    try:
        import h5py
    except:
        raise ImportError('h5py needs to be installed for using this function')

    hic = HiC()
    with h5py.File(filename, 'r') as f:
        for key in f.keys():
            try:
                value = f[key][:]
            except:
                value = f[key][()]
            setattr(hic, key, value)

    return hic
