# -*- coding: utf-8 -*-
"""This module defines functions for parsing and writing `EMD map files`_.

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

__all__ = ['parseEMDStream', 'parseEMD', 'writeEMD', 'TRNET']

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
    :type cutoff: float or None

    :arg n_nodes: A bead based network will be constructed over provided density map. n_nodes parameter
    will show the number of beads that will fit to density map. 
    :type n_nodes: integer  

    :arg num_iter: After reading density map, coordinates are predicted with topological domain reconstruction
    method. This parameter is the total number of iterations of this algorithm: 
    :type num_iter: integer

    :arg return_map: Return the density map itself. Default is False in line with previous behaviour.
    :type return_map: bool

    :arg make_nodes: Use the topology representing network algorithm to fit pseudoatom nodes to the map.
        Default is True. Changing this to False is equivalent to asking for return_map to be True.
    :type make_nodes: bool
    """

    title = kwargs.get('title', None)
    if not os.path.isfile(emd):
        raise IOError('EMD file {0} is not available in the directory {1}'
                        .format(emd),os.getcwd())
    if title is None:
        kwargs['title'], ext = os.path.splitext(os.path.split(emd)[1])

    emdStream = openFile(emd, 'rb')
    result = parseEMDStream(emdStream, **kwargs)
    emdStream.close()
    return result

def _parseEMDLines(atomgroup, stream, cutoff=None, n_nodes=1000, num_iter=20,return_map=False,make_nodes=True):
    """ Returns an AtomGroup. see also :func:`.parseEMDStream()`.

    :arg stream: stream from parser.
    """

    if not n_nodes > 0:
        raise ValueError('n_nodes should be larger than 0')

    emd = EMDMAP(stream, cutoff)

    if make_nodes:
        coordinates = np.zeros((n_nodes, 3), dtype=float)
        atomnames = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['name'].dtype)
        resnames = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['resname'].dtype)
        resnums = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['resnum'].dtype)
        chainids = np.zeros(n_nodes, dtype=ATOMIC_FIELDS['chain'].dtype)

        trn = TRNET(n_nodes = n_nodes)
        trn.inputMap(emd, sample='density')

        trn.run(tmax = num_iter)
        for i in range(n_nodes):
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
 
    if make_nodes:
        if return_map:
            return emd, atomgroup
        else:
            return atomgroup
    else:
        return emd


def parseEMDStream(stream, **kwargs):
    """ Returns an :class:`.AtomGroup` containing EMD data parsed from a stream of EMD file.

    :arg stream: Anything that implements the method ``readlines``
        (e.g. :class:`file`, buffer, stdin)"""

    cutoff = kwargs.get('cutoff', None)
    if cutoff is not None:
        cutoff = float(cutoff)

    n_nodes = int(kwargs.get('n_nodes', 1000))
    num_iter = int(kwargs.get('num_iter', 20))
    return_map = kwargs.get('return_map',False)
    make_nodes = kwargs.get('make_nodes',True)

    if return_map is False and make_nodes is False:
        LOGGER.warn('At least one of return_map and make_nodes should be True. '
                    'Setting make_nodes to False was an intentional change from the default '
                    'so return_map has been set to True.')
        kwargs['return_map'] = True

    title_suffix = kwargs.get('title_suffix','')
    atomgroup = AtomGroup(str(kwargs.get('title', 'Unknown')) + title_suffix)

    if make_nodes:
        LOGGER.info('Building coordinates from electron density map. This may take a while.')
        LOGGER.timeit()

        if return_map:
            emd, atomgroup = _parseEMDLines(atomgroup, stream, cutoff=cutoff, n_nodes=n_nodes, \
                                            num_iter=num_iter, return_map=return_map, \
                                            make_nodes=make_nodes)
        else:
            atomgroup = _parseEMDLines(atomgroup, stream, cutoff=cutoff, n_nodes=n_nodes, \
                                       num_iter=num_iter, return_map=return_map, \
                                       make_nodes=make_nodes)
        LOGGER.report('{0} atoms and {1} coordinate sets were '
                      'parsed in %.2fs.'.format(atomgroup.numAtoms(), atomgroup.numCoordsets()))
    else: 
        emd = _parseEMDLines(atomgroup, stream, cutoff=cutoff, n_nodes=n_nodes, \
                             num_iter=num_iter, return_map=return_map, \
                             make_nodes=make_nodes)

    if make_nodes:
        if return_map:
            return emd, atomgroup
        else:
            return atomgroup
    else:
        return emd

def writeEMD(filename,emd):
    '''
    Write a map file in MRC2014 format (counting words 25 to 49 as 'extra').

    :arg emd: an EMD object containing data to be written to file
    :type emd: :class:`.EMD`
    '''

    f = open(filename, "wb")
    f.write(st.pack('<L',emd.NC))
    f.write(st.pack('<L',emd.NR))
    f.write(st.pack('<L',emd.NS))
    f.write(st.pack('<L',emd.mode))
    f.write(st.pack('<L',emd.ncstart))
    f.write(st.pack('<L',emd.nrstart))
    f.write(st.pack('<L',emd.nsstart))
    f.write(st.pack('<L',emd.Nx))
    f.write(st.pack('<L',emd.Ny))
    f.write(st.pack('<L',emd.Nz))
    f.write(st.pack('<f',emd.Lx))
    f.write(st.pack('<f',emd.Ly))
    f.write(st.pack('<f',emd.Lz))
    f.write(st.pack('<L',emd.a))
    f.write(st.pack('<L',emd.b))
    f.write(st.pack('<L',emd.c))
    f.write(st.pack('<L',emd.mapc))
    f.write(st.pack('<L',emd.mapr))
    f.write(st.pack('<L',emd.maps))
    f.write(st.pack('<f',emd.dmin))
    f.write(st.pack('<f',emd.dmax))
    f.write(st.pack('<f',emd.dmean))
    f.write(st.pack('<L',emd.ispg))
    f.write(st.pack('<L',emd.nsymbt))
    f.write(st.pack('<100s',emd.extra))
    f.write(st.pack('<f',emd.x0))
    f.write(st.pack('<f',emd.y0))
    f.write(st.pack('<f',emd.z0))
    f.write(st.pack('<4s',emd.wordMAP))
    f.write(st.pack('<4s',emd.machst))
    f.write(st.pack('<f',emd.rms))
    f.write(st.pack('<L',emd.nlabels))
    f.write(st.pack('<800s',emd.labels))

    for s in xrange(0, emd.NS):
        for r in xrange(0, emd.NR):
            for c in xrange(0, emd.NC):
                f.write(st.pack('<f',emd.density[s,r,c]))

    f.close()

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
        self.ncstart = st.unpack('<l', stream.read(4))[0]
        self.nrstart = st.unpack('<l', stream.read(4))[0]
        self.nsstart = st.unpack('<l', stream.read(4))[0
]
        # Number of intervals along x, y, z (3 words, 12 bytes, 29-40)
        self.Nx = st.unpack('<L', stream.read(4))[0]
        self.Ny = st.unpack('<L', stream.read(4))[0]
        self.Nz = st.unpack('<L', stream.read(4))[0]

        # Cell dimensions (Angstroms) (3 words, 12 bytes, 41-52)
        self.Lx = st.unpack('<f', stream.read(4))[0]
        self.Ly = st.unpack('<f', stream.read(4))[0]
        self.Lz = st.unpack('<f', stream.read(4))[0]
    
        # Cell angles (Degrees) (3 words, 12 bytes, 53-64)
        self.a = st.unpack('<L', stream.read(4))[0]
        self.b = st.unpack('<L', stream.read(4))[0]
        self.c = st.unpack('<L', stream.read(4))[0]

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
        for s in xrange(0, self.NS):
            for r in xrange(0, self.NR):
                for c in xrange(0, self.NC):
                    d = st.unpack('<f', stream.read(4))[0]
                    if cutoff is not None and d < cutoff:
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
