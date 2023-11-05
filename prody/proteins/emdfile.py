# -*- coding: utf-8 -*-
"""This module defines functions for parsing and writing `EMD map files`_.

.. _EMD map files: http://emdatabank.org/mapformat.html"""

from numbers import Number
import os.path

from prody.atomic import AtomGroup
from prody.atomic import flags
from prody.atomic import ATOMIC_FIELDS

from prody.utilities import openFile, isListLike, copy
from prody import LOGGER

from .localpdb import fetchPDB

import struct as st
import numpy as np

__all__ = ['parseEMDStream', 'parseEMD', 'writeEMD', 'TRNET', 'EMDMAP']


class EMDParseError(Exception):
    pass


def parseEMD(emd, **kwargs):
    """Parses an EM density map in EMD/MRC2014 format and 
    optionally returns an :class:`.AtomGroup` containing  
    beads built in the density using the TRN algorithm [_TM94]. 

    This function extends :func:`.parseEMDStream`.

    See :ref:`cryoem_analysis` for a usage example. 

    :arg emd: an EMD identifier or a file name. A 4-digit 
              EMDataBank identifier can be provided to download 
              it via FTP.
    :type emd: str

    :arg min_cutoff: minimum density cutoff to read EMD map. The regions with 
                 lower density than this cutoff are discarded.
                 This corresponds to the previous cutoff and take values from it.
    :type min_cutoff: float

    :arg max_cutoff: maximum density cutoff to read EMD map. The regions with 
                 higher density than this cutoff are discarded.
    :type max_cutoff: float

    :arg n_nodes: A bead based network will be constructed into the provided density map. 
                  This parameter will set the number of beads to fit to density map. 
                  Default is 0. Please change it to some number to run the TRN algorithm.
                  Other parameters are passed through as kwargs to :meth:`.TRNET.run`
                  as described in its docs.
    :type n_nodes: int

    :arg map: Return the density map itself. Default is **False** in line with previous behaviour.
        This value is reset to **True** if n_nodes is 0 or less.
    :type map: bool
    """

    title = kwargs.get('title', None)
    if not os.path.isfile(emd):
        if emd.startswith('EMD-') and len(emd[4:]) in [4, 5]:
            emd = emd[4:]

        if len(emd) in [4, 5] and emd.isdigit():
            if title is None:
                title = emd
                kwargs['title'] = title

            if os.path.isfile(emd + '.map'):
                filename = emd + '.map'
            elif os.path.isfile(emd + '.map.gz'):
                filename = emd + '.map.gz'
            else:
                filename = fetchPDB(emd, report=True,
                                    format='emd', compressed=False)
                if filename is None:
                    raise IOError('EMD map file for {0} could not be downloaded.'
                                  .format(emd))
            emd = filename
        else:
            raise IOError('EMD file {0} is not available in the directory {1}'
                          .format(emd, os.getcwd()))
    if title is None:
        kwargs['title'], ext = os.path.splitext(os.path.split(emd)[1])

    emdStream = openFile(emd, 'rb')
    result = parseEMDStream(emdStream, **kwargs)
    emdStream.close()

    return result


def parseEMDStream(stream, **kwargs):
    """Parse lines of data stream from an EMD/MRC2014 file and 
    optionally return an :class:`.AtomGroup` containing TRN 
    nodes based on it.

    :arg stream: Any object with the method ``readlines``
                (e.g. :class:`file`, buffer, stdin)
    """
    cutoff = kwargs.get('cutoff', None)
    min_cutoff = kwargs.get('min_cutoff', cutoff)
    if min_cutoff is not None:
        if isinstance(min_cutoff, Number):
            min_cutoff = float(min_cutoff)
        else:
            raise TypeError('min_cutoff should be a number or None')

    max_cutoff = kwargs.get('max_cutoff', None)
    if max_cutoff is not None:
        if isinstance(max_cutoff, Number):
            max_cutoff = float(max_cutoff)
        else:
            raise TypeError('max_cutoff should be a number or None')

    n_nodes = kwargs.get('n_nodes', 0)
    map = kwargs.get('map', False)

    if not isinstance(n_nodes, int):
        raise TypeError('n_nodes should be an integer')

    if n_nodes > 0:
        make_nodes = True
    else:
        make_nodes = False
        map = True
        LOGGER.info('As n_nodes is less than or equal to 0, no nodes will be'
                    ' made and the raw map will be returned')

    emd = EMDMAP(stream, min_cutoff, max_cutoff)

    if make_nodes:
        title_suffix = kwargs.get('title_suffix', '')
        atomgroup = AtomGroup(str(kwargs.get('title', 'Unknown')) + title_suffix)
        atomgroup._n_atoms = n_nodes

        coordinates = np.zeros((n_nodes, 3), dtype=float)
        atomnames = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['name'].dtype)
        resnames = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['resname'].dtype)
        resnums = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['resnum'].dtype)
        chainids = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['chain'].dtype)

        trn = TRNET(n_nodes=n_nodes)
        trn.inputMap(emd, sample='density')

        trn.run(**kwargs)
        for i in range(n_nodes):
            coordinates[i, :] = trn.W[i, :]
            atomnames[i] = 'B'
            resnames[i] = 'CGB'
            resnums[i] = i+1
            chainids[i] = 'X'

        atomgroup.setCoords(coordinates)
        atomgroup.setNames(atomnames)
        atomgroup.setResnames(resnames)
        atomgroup.setResnums(resnums)
        atomgroup.setChids(chainids)

    if make_nodes:
        if map:
            return atomgroup, emd
        else:
            return atomgroup
    else:
        return emd


def writeEMD(filename, emd):
    '''
    Writes a map file in MRC2014 format (counting words 25 to 49 as 'extra').

    :arg filename: name for output file
    :type filename: str

    :arg emd: an EMD object containing data to be written to file
    :type emd: :class:`.EMD`
    '''

    f = open(filename, "wb")
    f.write(st.pack('<L', emd.NC))
    f.write(st.pack('<L', emd.NR))
    f.write(st.pack('<L', emd.NS))
    f.write(st.pack('<L', emd.mode))
    f.write(st.pack('<l', emd.ncstart))
    f.write(st.pack('<l', emd.nrstart))
    f.write(st.pack('<l', emd.nsstart))
    f.write(st.pack('<L', emd.Nx))
    f.write(st.pack('<L', emd.Ny))
    f.write(st.pack('<L', emd.Nz))
    f.write(st.pack('<f', emd.Lx))
    f.write(st.pack('<f', emd.Ly))
    f.write(st.pack('<f', emd.Lz))
    f.write(st.pack('<f', emd.a))
    f.write(st.pack('<f', emd.b))
    f.write(st.pack('<f', emd.c))
    f.write(st.pack('<L', emd.mapc))
    f.write(st.pack('<L', emd.mapr))
    f.write(st.pack('<L', emd.maps))
    f.write(st.pack('<f', emd.dmin))
    f.write(st.pack('<f', emd.dmax))
    f.write(st.pack('<f', emd.dmean))
    f.write(st.pack('<L', emd.ispg))
    f.write(st.pack('<L', emd.nsymbt))
    f.write(st.pack('<100s', emd.extra))
    f.write(st.pack('<f', emd.x0))
    f.write(st.pack('<f', emd.y0))
    f.write(st.pack('<f', emd.z0))
    f.write(st.pack('<4s', emd.wordMAP))
    f.write(st.pack('<4s', emd.machst))
    f.write(st.pack('<f', emd.rms))
    f.write(st.pack('<L', emd.nlabels))
    f.write(st.pack('<800s', emd.labels))

    for s in range(0, emd.NS):
        for r in range(0, emd.NR):
            for c in range(0, emd.NC):
                f.write(st.pack('<f', emd.density[s, r, c]))

    f.close()


class EMDMAP(object):
    """Class for handling EM density maps in EMD/MRC2014 format.
    
    :arg stream: a file stream containing data from an EMD/MRC file.
    
    :arg min_cutoff: minimum cutoff for thresholding
    :type min_cutoff: None, float

    :arg max_cutoff: maximum cutoff for thresholding
    :type max_cutoff: None, float
    """
    def __init__(self, stream, min_cutoff, max_cutoff):
        if min_cutoff is not None and not isinstance(min_cutoff, Number):
            raise TypeError('min_cutoff should be a number or None')

        if max_cutoff is not None and not isinstance(max_cutoff, Number):
            raise TypeError('max_cutoff should be a number or None')
        
        self._filename = stream.name

        # Number of columns, rows, and sections (3 words, 12 bytes, 1-12)
        self.NC = st.unpack('<L', stream.read(4))[0]
        self.NR = st.unpack('<L', stream.read(4))[0]
        self.NS = st.unpack('<L', stream.read(4))[0]
        self.Ntot = self.NC * self.NR * self.NS

        # Mode (1 word, 4 bytes, 13-16)
        self.mode = st.unpack('<L', stream.read(4))[0]

        # Number of first column, row, section (3 words, 12 bytes, 17-28)
        self.ncstart = st.unpack('<l', stream.read(4))[0]
        self.nrstart = st.unpack('<l', stream.read(4))[0]
        self.nsstart = st.unpack('<l', stream.read(4))[0]

        # Number of intervals along x, y, z (3 words, 12 bytes, 29-40)
        self.Nx = st.unpack('<L', stream.read(4))[0]
        self.Ny = st.unpack('<L', stream.read(4))[0]
        self.Nz = st.unpack('<L', stream.read(4))[0]

        # Cell dimensions (Angstroms) (3 words, 12 bytes, 41-52)
        self.Lx = st.unpack('<f', stream.read(4))[0]
        self.Ly = st.unpack('<f', stream.read(4))[0]
        self.Lz = st.unpack('<f', stream.read(4))[0]

        # Cell angles (Degrees) (3 words, 12 bytes, 53-64)
        self.a = st.unpack('<f', stream.read(4))[0]
        self.b = st.unpack('<f', stream.read(4))[0]
        self.c = st.unpack('<f', stream.read(4))[0]

        # Which axis corresponds to column, row, and sections (1, 2, 3 for x, y ,z)
        # (3 words, 12 bytes, 65-76)
        self.mapc = st.unpack('<L', stream.read(4))[0]
        self.mapr = st.unpack('<L', stream.read(4))[0]
        self.maps = st.unpack('<L', stream.read(4))[0]

        # Density values (min, max, mean) (3 words, 12 bytes, 77-88)
        self.dmin = st.unpack('<f', stream.read(4))[0]
        self.dmax = st.unpack('<f', stream.read(4))[0]
        self.dmean = st.unpack('<f', stream.read(4))[0]

        # Space group number (1 word, 4 bytes, 89-92)
        # For EM/ET, this encodes the type of data:
        # 0 for 2D images and image stacks, 1 for 3D volumes, 401 for volume stacks
        self.ispg = st.unpack('<L', stream.read(4))[0]

        # size of extended header (1 word, 4 bytes, 93-96)
        # contained symmetry records in original format definition
        self.nsymbt = st.unpack('<L', stream.read(4))[0]

        # we treat this all as extra stuff like MRC2014 format (25 word, 100 bytes, 97-196)
        self.extra = st.unpack('<100s', stream.read(100))[0]

        # origins for x, y, z (3 words, 12 bytes, 197-208)
        self.x0 = st.unpack('<f', stream.read(4))[0]
        self.y0 = st.unpack('<f', stream.read(4))[0]
        self.z0 = st.unpack('<f', stream.read(4))[0]

        # the character string 'MAP' to identify file type (1 word, 4 bytes, 209-212)
        self.wordMAP = st.unpack('<4s', stream.read(4))[0]

        # machine stamp encoding byte ordering of data (1 word, 4 bytes, 213-216)
        self.machst = st.unpack('<4s', stream.read(4))[0]

        # rms deviation of map from mean density (1 word, 4 bytes, 217-220)
        self.rms = st.unpack('<f', stream.read(4))[0]

        # number of labels being used (1 word, 4 bytes, 221-224)
        self.nlabels = st.unpack('<L', stream.read(4))[0]

        # 10 80-character text labels, which we leave concatenated (200 words, 800 bytes, 225-1024)
        self.labels = st.unpack('<800s', stream.read(800))[0]

        # Data blocks (1024-end)
        self.density = np.empty([self.NS, self.NR, self.NC])
        for s in range(0, self.NS):
            for r in range(0, self.NR):
                for c in range(0, self.NC):
                    d = st.unpack('<f', stream.read(4))[0]
                    if min_cutoff is not None and d < min_cutoff:
                        d = 0
                    if max_cutoff is not None and d > min_cutoff:
                        d = 0
                    self.density[s, r, c] = d

        self.sampled = False

    def thresholdMap(self, min_cutoff=None, max_cutoff=None):
        """Thresholds a map and returns a new map like the equivalent function in TEMPy"""
        newDensity = self.density.copy()

        if max_cutoff is not None:
            if isinstance(max_cutoff, Number):
                min_cutoff = float(max_cutoff)
            else:
                raise TypeError('max_cutoff should be a number or None')

            newDensity = newDensity * (newDensity < max_cutoff)
            
        if min_cutoff is not None:
            if isinstance(min_cutoff, Number):
                min_cutoff = float(min_cutoff)
            else:
                raise TypeError('min_cutoff should be a number or None')

            newDensity = newDensity * (newDensity > min_cutoff)

        newMap = self.copyMap()
        newMap.density = newDensity
        return newMap

    def numidx2matidx(self, numidx):
        """ Given index of the position, it will return the numbers of section, row and column. """
        # calculate section idx
        s = int(numidx / (self.NC * self.NR))
        numidx = numidx - s * self.NC * self.NR
        # calculate row idx
        r = int(numidx / self.NC)
        # calculate column idx
        c = int(numidx - r * self.NC)
        return s, r, c

    def drawsample(self):
        if not self.sampled:
            self.cumsumdens = np.cumsum(self.density)
            self.sampled = True
        summ = self.cumsumdens[-1]
        r = np.random.rand() * summ
        j = np.searchsorted(self.cumsumdens, r)
        return self.numidx2matidx(j)

    def drawsample_uniform(self):
        r = int(np.random.rand() * self.Ntot)
        return self.numidx2matidx(r)

    def center(self):
        return int(self.NS / 2), int(self.NR / 2), int(self.NC / 2)

    def getOrigin(self):
        return self.x0, self.y0, self.z0

    def setOrigin(self, x0, y0, z0):
        self.x0, self.y0, self.z0 = x0, y0, z0

    origin = property(getOrigin, setOrigin)

    def getTitle(self):
        return self._filename

    def setTitle(self, title):
        self._filename = title

    filename = property(getTitle, setTitle)

    def getApix(self):
        return np.array((self.Lx / self.NS,
                         self.Ly / self.NR,
                         self.Lz / self.NC))

    def setApix(self, apix):
        if not isListLike(apix):
            try:
                apix = [apix, apix, apix]
            except:
                raise TypeError('apix must be a single value or list-like')

        if len(apix) != 3:
            raise ValueError('apix must be a single value or 3 values')
        
        self._apix = apix
        self.Lx = apix[0] * self.NS
        self.Ly = apix[1] * self.NR
        self.Lz = apix[2] * self.NC

    apix = property(getApix, setApix)

    def coordinate(self, sec, row, col):
        """Given a position as *sec*, *row* and *col*, 
        it will return its coordinate in Angstroms. """
        # calculate resolution
        res = np.empty(3)
        res[self.mapc - 1] = self.NC
        res[self.mapr - 1] = self.NR
        res[self.maps - 1] = self.NS
        res = np.divide(np.array([self.Lx, self.Ly, self.Lz]), res)

        # find coordinates in voxels relative to start
        ret = np.empty(3)
        ret[self.mapc - 1] = col + self.ncstart
        ret[self.mapr - 1] = row + self.nrstart
        ret[self.maps - 1] = sec + self.nsstart

        # convert to Angstroms
        ret = np.multiply(ret, res)
        return ret

    def toTEMPyMap(self):
        """Convert to a TEMPy Map."""
        try:
            from TEMPy.maps.em_map import Map
            from TEMPy.maps.map_parser import MapParser
        except ImportError:
            raise ImportError('TEMPy needs to be installed for this functionality')
        
        header = MapParser.readMRCHeader(self.filename)
        newOrigin = np.array((self.ncstart, self.nrstart, self.nsstart)) * self.apix
        return Map(self.density, newOrigin, self.apix, self.filename, header)

    def copyMap(self):
        """
        Copy to a new object.
        """
        return copy(self)



class TRNET(object):
    """Class for building topology representing networks using 
    EM density maps. It uses the algorithm described in [TM94]_.

    .. [TM94] Martinetz T, Schulten K, Topology Representing Networks.
       *Neural Networks* **1994** 7(3):507-552."""

    def __init__(self, n_nodes):
        self.N = n_nodes
        self.W = np.empty([n_nodes, 3])
        self.C = np.eye(n_nodes, n_nodes)
        # test
        self.V = np.array([])

    def inputMap(self, emdmap, sample='density'):
        self.map = emdmap
        # initialize the positions of nodes
        for i in range(self.N):
            if sample == 'density':
                p = self.map.drawsample()
            elif sample == 'uniform':
                p = self.map.drawsample_uniform()
            elif sample == 'center':
                p = self.map.center()
            else:
                p = (0, 0, 0)
            self.W[i, :] = self.map.coordinate(p[0], p[1], p[2])

    def runOnce(self, t, l, ep, T, c=0):
        # draw a point from the map
        p = self.map.drawsample()
        v = self.map.coordinate(p[0], p[1], p[2])
        if len(self.V) == 0:
            self.V = v
        else:
            self.V = np.vstack((self.V, v))

        # calc the squared distances \\ws - v\\^2
        D = v - self.W
        sD = np.empty(self.N)
        for i in range(self.N):
            d = D[i, :]
            sD[i] = np.dot(d, d)

        # calc the closeness rank k's
        I = np.argsort(sD)
        K = np.empty(I.shape)
        K[I] = range(len(I))

        # move the nodes
        if c == 0:
            K = K[:, np.newaxis]
            self.W += ep * np.exp(-K/l) * D
        else:
            kc = - l * np.log(c/ep)
            idx = K < kc
            K = K[:, np.newaxis]
            self.W[idx, :] += ep * np.exp(-K[idx]/l) * D[idx, :]

        if T >= 0:
            # search for i0 and i1
            i0 = I[0]
            i1 = I[1]

            # refresh connections
            for i in range(self.N):
                if i == i1:
                    self.C[i0, i] = 1
                    self.C[i, i0] = 1
                elif i != i0 and self.C[i0, i] > 0:
                    self.C[i0, i] = self.C[i0, i] + 1
                    if self.C[i0, i] > T:
                        self.C[i0, i] = 0
                    self.C[i, i0] = self.C[i0, i]

    def run(self, **kwargs):
        """
        :arg tmax: multiplicative factor such that the maximum total number of 
            iterations is tmax times the number of beads
            default 200
        :type tmax: int

        :arg li: initial Gaussian bandwidth for determining how much each node is moved
            As the iterations progress, the bandwidth increases from li to lf.
            default 0.2
        :type li: float

        :arg lf: final Gaussian bandwidth for determining how much each node is moved
            As the iterations progress, the bandwidth increases from li to lf.
            default 0.01
        :type lf: float

        :arg ei: initial value of the adaptive step size
            As the iterations progress, the step size increases from ei to ef.
            default 0.3
        :type ei: float

        :arg ef: final value of the adaptive step size
            As the iterations progress, the step size increases from ei to ef.
            default 0.05
        :type ef: float

        :arg c: cutoff for moving the nodes. When c=0, all nodes are moved in each iteration. 
            When c>0, only the nearest c/#nodes nodes are moved. This parameter is used for optimization.
            default 0
        :type c: float

        :arg calcC: whether to calculate the connectivity matrix from TRN. This is **False** by default 
            because the connectivity is usually built by ANM or GNM.
            default **False**
        :type calcC: bool

        :arg Ti: initial value of the adaptive threshold for building the connectivity. Not used if calcC is False.
            default 0.1
        :type Ti: float

        :arg Tf: final value of the adaptive threshold for building the connectivity. Not used if calcC is False.
            default 2
        :type Tf: float
        """

        tmax = kwargs.get('tmax', 200)
        li = kwargs.get('li', 0.2)
        lf = kwargs.get('lf', 0.01)
        ei = kwargs.get('ei', 0.3)
        ef = kwargs.get('ef', 0.05)
        Ti = kwargs.get('Ti', 0.1)
        Tf = kwargs.get('Tf', 2)
        c = kwargs.get('c', 0)
        calcC = kwargs.get('calcC', False)

        LOGGER.info('Building coordinates from electron density map. This may take a while.')
        LOGGER.timeit('_prody_make_nodes')
        tmax = int(tmax * self.N)
        li = li * self.N
        if calcC:
            Ti = Ti * self.N
            Tf = Tf * self.N
        for t in range(1, tmax + 1):
            # calc the parameters
            tt = float(t) / tmax
            l = li * np.power(lf / li, tt)
            ep = ei * np.power(ef / ei, tt)
            if calcC:
                T = Ti * np.power(Tf / Ti, tt)
            else:
                T = -1
            self.runOnce(t, l, ep, T, c)
        LOGGER.report('{0} pseudoatoms were fitted in %.2fs.'.format(
            self.N), '_prody_make_nodes')
        return

    def run_n_pause(self, k0, k, tmax=200, li=0.2, lf=0.01, ei=0.3,
                    ef=0.05, Ti=0.1, Tf=2):
        tmax = int(tmax * self.N)
        li = li * self.N
        Ti = Ti * self.N
        Tf = Tf * self.N
        for t in range(k0, k + 1):
            # calc the parameters
            tt = float(t) / tmax
            l = li * np.power(lf / li, tt)
            ep = ei * np.power(ef / ei, tt)
            T = Ti * np.power(Tf / Ti, tt)
            # run once
            self.runOnce(t, l, ep, T)
        return

    def outputEdges(self):
        return self.C > 0
