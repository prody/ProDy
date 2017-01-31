# -*- coding: utf-8 -*-
"""This module defines functions for parsing and writing `PDB files`_.

.. _EMD maps: http://emdatabank.org/mapformat.html"""

from collections import defaultdict
import os.path

from prody.atomic import AtomGroup
from prody.atomic import flags
from prody.atomic import ATOMIC_FIELDS

from prody.utilities import openFile
from prody import LOGGER, SETTINGS

import struct as st
import numpy as np

__all__ = ['parseEMDStream', 'parseEMD']

class EMDParseError(Exception):
    pass

""" For  documentation"""

def parseEMD(emd, **kwargs):
    """Returns an :class:`.AtomGroup` containing the information parsed from EMD file. 

    This function extends :func:`.parseEMDStream`.

    See :ref:`parseEMD` for a detailed usage example. 

    :arg emd: an EMD identifier or a file name, EMD files should be locally available. 

    :arg cutoff: density cutoff to read EMD map. The regions with lower density than given cutoff 
    are discarded.
    :type cutoff: float

    :arg n_nodes: A bead based network will be constructed over provided density map. n_nodes parameter
    will show the number of beads that will fit to density map. 
    :type n_nodes: integer  

    :arg num_iter: After reading density map, coordinates are predicted with topological domain reconstruction
    method. This parameter is the total number of iterations of this algorithm: 
    :type num_iter: integer
    """

    title = kwargs.get('title', None)
    cutoff = float(kwargs.get('cutoff', 1.20))
    n_nodes = int(kwargs.get('n_nodes', 1000))
    num_iter = int(kwargs.get('num_iter', 20))
    if not os.path.isfile(emd):
        raise IOError('EMD file {0} is not available in the directory {1}'
                        .format(emd),os.getcwd())
    if title is None:
        title, ext = os.path.splitext(os.path.split(emd)[1])
        kwargs['title'] = title
    kwargs['cutoff'] = cutoff
    kwargs['n_nodes'] = n_nodes
    emd = openFile(emd, 'rt')
    result = parseEMDStream(emd, **kwargs)
    emd.close()
    return result

def _parseEMDLines(atomgroup, stream, cutoff=1.2, n_nodes=1000, num_iter=20, format='EMD'):
    """ Returns an AtomGroup. see also :func:`.parseEMDStream()`.

    :arg stream: stream from parser.
    """

    format = format.upper()
    if format == 'EMD' or format == 'MAP':
        isEMD = True
    else:
        isEMD = False

    n_atoms = n_nodes
    if n_atoms > 0:
        asize = n_atoms

    alength = asize
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    hetero = np.zeros(asize, dtype=bool)
    termini = np.zeros(asize, dtype=bool)
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
    elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)
    bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
    occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)
    anisou = None
    siguij = None

    emd = EMDMAP(stream, cutoff)
    trn = TRNET(n_nodes = asize)
    trn.inputMap(emd, sample='density')

    trn.run(tmax = num_iter)
    for i in range(asize):
        coordinates[i,:] = trn.W[i,:]
        atomnames[i] = 'B'
        resnames[i] = 'CGB'
        resnums[i] = i+1
        chainids[i] = 'X'

    atomgroup.setCoords(coordinates)
    atomgroup.setNames(atomnames)
    atomgroup.setResnames(resnames)
    atomgroup.setResnums(resnums)
    atomgroup.setChids(chainids)

    return atomgroup


def parseEMDStream(stream, **kwargs):
    """ Returns an :class:`.AtomGroup` containing EMD data parsed from a stream of EMD file.

    :arg stream: Anything that implements the method ``readlines``
        (e.g. :class:`file`, buffer, stdin)"""

    cutoff = float(kwargs.get('cutoff', 1.20))
    n_nodes = int(kwargs.get('n_nodes', 1000))
    num_iter = int(kwargs.get('num_iter', 20))

    ag = None
    title_suffix = ''
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    else: 
        ag = AtomGroup(str(kwargs.get('title', 'Unknown')) + title_suffix)
        n_csets = 0

    biomol = kwargs.get('biomol', False)
    hd = None
    LOGGER.warn('Building coordinates from electron density map. This may take a while.')
    LOGGER.timeit()
    _parseEMDLines(ag, stream, cutoff=cutoff, n_nodes=n_nodes, num_iter=num_iter, format='EMD')
    LOGGER.report('{0} atoms and {1} coordinate sets were '
                      'parsed in %.2fs.'.format(ag.numAtoms(),
                         ag.numCoordsets() - n_csets))

class EMDMAP:
    def __init__(self, stream, cutoff):
        # Number of columns, rows, and sections (3 words, 12 bytes, 1-12)
        self.NC = st.unpack('<L', stream.read(4))[0]
        self.NR = st.unpack('<L', stream.read(4))[0]
        self.NS = st.unpack('<L', stream.read(4))[0]
        self.Ntot = self.NC * self.NR * self.NS

        # Mode (1 word, 4 bytes, 13-16)
        self.mode = st.unpack('<L', stream.read(4))[0]

        # Number of first column, row, section (3 words, 12 bytes, 17-28)
        self.ncstart = st.unpack('<L', stream.read(4))[0]
        self.nrstart = st.unpack('<L', stream.read(4))[0]
        self.nsstart = st.unpack('<L', stream.read(4))[0]

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

        # Not interested (1 word, 4 bytes, 89-92)
        stream.read(1*4)

        # Not interested (1 word, 4 bytes, 93-96)
        self.nsym = st.unpack('<f', stream.read(4))[0]

        # Not interested (25 word, 4 bytes, 97-196)
        stream.read(25*4)

        # origins for x, y, z (3 words, 12 bytes, 197-208)
        self.x0 = st.unpack('<f', stream.read(4))[0]
        self.y0 = st.unpack('<f', stream.read(4))[0]
        self.z0 = st.unpack('<f', stream.read(4))[0]

        # Not interested (204 words, 816 bytes, 209-1024)
        stream.read(204*4)

        # Data blocks (1024-end)
        self.density = np.empty([self.NS, self.NR, self.NC])
        for s in xrange(0, self.NS):
            for r in xrange(0, self.NR):
                for c in xrange(0, self.NC):
                    d = st.unpack('<f', stream.read(4))[0]
                    if not np.isnan(cutoff) and d < cutoff:
                        d = 0
                    self.density[s, r, c] = d

        self.sampled = False

    def numidx2matidx(self, numidx):
        """ Given index of the position, it will return the numbers of section, row and column. """
        # calculate section idx
        s = numidx / (self.NC * self.NR)
        numidx = numidx - s * self.NC * self.NR
        # calculate row idx
        r = numidx / self.NC
        # calculate column idx
        c = numidx - r * self.NC
        return s, r, c

    def drawsample(self):
        if not self.sampled:
            self.cumsumdens = np.cumsum(self.density)
            self.sampled = True
        summ = self.cumsumdens[-1]
        r = np.random.rand() * summ
        j = np.searchsorted(self.cumsumdens,r)
        return self.numidx2matidx(j)

    def drawsample_uniform(self):
        r = int(np.random.rand() * self.Ntot)
        return self.numidx2matidx(r)

    def center(self):
        return self.NS / 2, self.NR / 2, self.NC / 2

    def coordinate(self, sec, row, col ):
        # calculate resolution
        res = np.empty(3)
        res[self.mapc - 1] = self.NC
        res[self.mapr - 1] = self.NR
        res[self.maps - 1] = self.NS
        res = np.divide(np.array([self.Lx, self.Ly, self.Lz]), res)
        
        ret = np.empty(3)
        ret[self.mapc - 1] = col + self.ncstart
        ret[self.mapr - 1] = row + self.nrstart
        ret[self.maps - 1] = sec + self.nsstart
        
        ret = np.multiply(ret, res)
        return ret

class TRNET:
    def __init__(self, n_nodes):
        self.N = n_nodes
        self.W = np.empty([n_nodes, 3])
        self.C = np.eye(n_nodes, n_nodes)
        # test
        self.V = np.array([])
        
    def inputMap(self, emdmap, sample = 'density'):
        self.map = emdmap
        # initialize the positions of nodes
        for i in xrange(self.N):
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
        for i in xrange(self.N):
            d = D[i, :]
            sD[i] = np.dot(d,d)
        
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
            
        
        if T>=0:
            # search for i0 and i1
            i0 = I[0]
            i1 = I[1]
            
            # refresh connections
            for i in xrange(self.N):
                if i == i1:
                    self.C[i0, i] = 1
                    self.C[i, i0] = 1
                elif i!=i0 and self.C[i0, i] > 0:
                    self.C[i0, i] = self.C[i0, i] + 1
                    if self.C[i0, i] > T:
                        self.C[i0, i] = 0
                    self.C[i, i0] = self.C[i0, i]
                
    def run(self, tmax = 200, li = 0.2, lf = 0.01, ei = 0.3,
            ef = 0.05, Ti = 0.1, Tf = 2, c = 0, calcC = False):
        tmax = int(tmax * self.N)
        li = li * self.N
        if calcC:
            Ti = Ti * self.N
            Tf = Tf * self.N        
        for t in xrange(1, tmax + 1):
            # calc the parameters
            tt = float(t) / tmax
            l = li * np.power(lf / li, tt)
            ep = ei * np.power(ef / ei, tt)
            if calcC:
                T = Ti * np.power(Tf / Ti, tt)
            else:
                T = -1
            # run once
            self.runOnce(t, l, ep, T, c)
            
            #if t % 1000 == 0:
            #    print str(t) + " steps have been run"
    
    def run_n_pause(self, k0, k, tmax = 200, li = 0.2, lf = 0.01, ei = 0.3,
            ef = 0.05, Ti = 0.1, Tf = 2):
        tmax = int(tmax * self.N)
        li = li * self.N
        Ti = Ti * self.N
        Tf = Tf * self.N        
        for t in xrange(k0, k + 1):
            # calc the parameters
            tt = float(t) / tmax
            l = li * np.power(lf / li, tt)
            ep = ei * np.power(ef / ei, tt)
            T = Ti * np.power(Tf / Ti, tt)            
            # run once
            self.runOnce(t, l, ep, T)
            
            #if t % 1000 == 0:
            #    print str(t) + " steps have been run"
        
    def outputEdges(self):
        return self.C > 0
