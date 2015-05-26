# -*- coding: utf-8 -*-
"""This module defines a class and a function for rotating translating blocks
(RTB) calculations."""

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB
from prody.utilities import importLA, checkCoords
from numpy import sqrt, zeros, linalg, min, max
from subprocess import call

from .anm import ANMBase, calcANM
from .editing import reduceModel

__all__ = ['RTB', 'imANM', 'exANM']

class Increment(object):

    def __init__(self, s=0):

        self._i = s

    def __call__(self, i=1):

        self._i += i
        return self._i


class RTB(ANMBase):

    """Class for Rotations and Translations of Blocks (RTB) method ([FT00]_).
    Optional arguments permit imposing constrains along Z-direction as in
    *imANM* method described in [TL12]_.

    .. [FT00] Tama F, Gadea FJ, Marques O, Sanejouand YH. Building-block
       approach for determining low-frequency normal modes of macromolecules.
       *Proteins* **2000** 41:1-7.

    .. [TL12] Lezon TR, Bahar I, Constraints Imposed by the Membrane
       Selectively Guide the Alternating Access Dynamics of the Glutamate
       Transporter GltPh

    """

    def __init__(self, name='Unknown'):

        super(RTB, self).__init__(name)
        self._project = None


    def buildHessian(self, coords, blocks, cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg blocks: a list or array of block identifiers
        :type blocks: list, :class:`numpy.ndarray`

        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float

        :arg scale: scaling factor for force constant along Z-direction,
            default is 1.0
        :type scale: float

        """


        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')

        LOGGER.timeit('_rtb')
        self._n_atoms = natoms = int(coords.shape[0])
        if natoms != len(blocks):
            raise ValueError('len(blocks) must match number of atoms')
        from collections import defaultdict
        i = Increment()
        d = defaultdict(i)
        blocks = np.array([d[b] for b in blocks], np.int64)

        try:
            from collections import Counter
        except ImportError:
            counter = defaultdict(int)
            for b in blocks:
                counter[b] += 1
        else:
            counter = Counter(blocks)

        nblocks = len(counter)
        maxsize = 1
        nones = 0
        while counter:
            _, size = counter.popitem()
            if size == 1:
                nones += 1
            if size > maxsize:
                maxsize = size
        LOGGER.info('System has {0} blocks largest with {1} of {2} units.'
                    .format(nblocks, maxsize, natoms))
        nb6 = nblocks * 6 - nones * 3

        coords = coords.T.copy()

        self._hessian = hessian = np.zeros((nb6, nb6), float)
        self._project = project = np.zeros((natoms * 3, nb6), float)

        from .rtbtools import buildhessian
        buildhessian(coords, blocks, hessian, project,
                     natoms, nblocks, maxsize,
                     float(cutoff), float(gamma),
                     scale=float(kwargs.get('scale', 1.0)),
                     memlo=float(kwargs.get('membrane_low', 1.0)),
                     memhi=float(kwargs.get('membrane_high', -1.0)),)

        LOGGER.report('Hessian was built in %.2fs.', label='_rtb')

    def getProjection(self):
        """Return a copy of the projection matrix."""

        if self._project is not None:
            return self._project.copy()

    def _getProjection(self):

        return self._project

    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh`
        function to diagonalize the Hessian matrix. When Scipy is not found,
        :func:`numpy.linalg.eigh` is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
            If ``None`` is given, all modes will be calculated.
        :type n_modes: int or None, default is 20

        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``

        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        """

        super(RTB, self).calcModes(n_modes, zeros, turbo)

        self._array = np.dot(self._project, self._array)

def imANM(pdb='2nwl-mem.pdb', blk='2nwl.blk', scale=1.):
    
    from prody import parsePDB
    from numpy import zeros, dot

    pdb = parsePDB(pdb, subset='ca')
    pdb.setData('block', zeros(len(pdb), int))
    with open(blk) as inp:
        for line in inp:
            if line.startswith('BLOCK'):
                _, b, n1, c1, r1, n2, c2, r2 = line.split()
                sel = pdb.select('chain {} and resnum {} to {}'.format(c1, r1, r2))
                if sel:
                    sel.setData('block', int(b))
    pdb.setBetas(pdb.getData('block'))
    rtb = RTB(pdb)
    rtb.buildHessian(pdb, pdb.getData('block'), scale)
    h_prime = rtb.getHessian()
    p = rtb.getProjection()
    values, vectors = linalg.eigh(h_prime)
    vv = dot(p, vectors)
    return vv

def assign_lpvs(lat):
    lpv = zeros((3,3))
    if lat=='FCC':
        lpv[0,1]=1./sqrt(2)
        lpv[0,2]=1./sqrt(2)
        lpv[1,0]=1./sqrt(2)
        lpv[1,2]=1./sqrt(2)
        lpv[2,0]=1./sqrt(2)
        lpv[2,1]=1./sqrt(2)
    elif lat=='SC':
        lpv[0,0]=1
        lpv[1,1]=1
        lpv[2,2]=1
    elif lat=='SH':
        lpv[0,0]=1./2
        lpv[0,1]=-sqrt(3)/2
        lpv[1,0]=1./2
        lpv[1,1]=sqrt(3)/2
        lpv[2,2]=1.
    return lpv

def checkClash(coordinates, pdb, radius):
    for i in range(pdb.ca.getCoords().shape[0]):
        if linalg.norm(coordinates-pdb.ca.getCoords()[i])<radius:
            return False
    return True
def cgmembrane_wpdb(pdb_id=None, membrane_hi=None, membrane_lo=None, R=None, r=None, lat=None, outputName=None):
    
    if pdb_id==None:
        pdb_id = '2nwl-mem.pdb'
    if lat == None:
        lat = 'FCC'
    if outputName == None:
        outputName = 'membrane.pdb'
    if membrane_hi ==None:
        membrane_hi = 13.
    if membrane_lo ==None:
        membrane_lo = -13.
    if R == None:
        R=80
    if r == None:
        r=2.5
        
    pdb = parsePDB(pdb_id)
    coords = pdb.ca.getCoords()
    
    pxlo = min(coords[:,0])
    pxhi = max(coords[:,0])
    pylo = min(coords[:,1])
    pyhi = max(coords[:,1])
    pzlo = min(coords[:,2])
    pzhi = max(coords[:,2])
    pxlo = min([pxlo, 10000])
    pylo = min([pylo, 10000])
    pzlo = min([pzlo, 10000])
    pxhi = max([pxhi, -10000])
    pyhi = max([pyhi, -10000])
    pzhi = max([pzhi, -10000])
    if lat == None:
        lat = 'FCC'
    lpv = assign_lpvs(lat)
    imax = (R + lpv[0,2] * (membrane_hi - membrane_lo)/2.)/r
    jmax = (R + lpv[1,2] * (membrane_hi - membrane_lo)/2.)/r
    kmax = (R + lpv[2,2] * (membrane_hi - membrane_lo)/2.)/r
    f = open(outputName, 'w')
    atm = 0
    for i in range(-int(imax),int(imax+1)):
        for j in range(-int(jmax),int(jmax+1)):
            for k in range(-int(kmax),int(kmax+1)):
                X = zeros(3)
                for p in range(3):
                    X[p]=2.*r*(i*lpv[0,p]+j*lpv[1,p]+k*lpv[2,p])
                dd=0
                for p in range(3):
                    dd += X[p] ** 2
                if dd<R**2 and X[2]>membrane_lo and X[2]<membrane_hi:
                    if X[0]>pxlo and X[0]<pxhi and X[1]>pylo and X[1]<pyhi and X[2]>pzlo and X[2]<pzhi:
                        if checkClash(X, pdb, radius=5):
                            atm = atm + 1
                            f.write('ATOM%7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n' % (atm, atm, X[0], X[1], X[2]))
    f.close()
                


                
def cgmembrane_wopdb(membrane_hi=None, membrane_lo=None, R=None, r=None, lat=None, outputName=None):
    
    if lat == None:
        lat = 'FCC'
    if outputName == None:
        outputName = 'membrane.pdb'
    if membrane_hi ==None:
        membrane_hi = 13.
    if membrane_lo ==None:
        membrane_lo = -13.
    if R == None:
        R=80
    if r == None:
        r=2.5
        
    lpv = assign_lpvs(lat)
    imax = (R + lpv[0,2] * (membrane_hi - membrane_lo)/2.)/r
    jmax = (R + lpv[1,2] * (membrane_hi - membrane_lo)/2.)/r
    kmax = (R + lpv[2,2] * (membrane_hi - membrane_lo)/2.)/r
    f = open(outputName, 'w')
    atm = 0
    for i in range(-int(imax),int(imax+1)):
        for j in range(-int(jmax),int(jmax+1)):
            for k in range(-int(kmax),int(kmax+1)):
                X = zeros(3)
                for p in range(3):
                    X[p]=2.*r*(i*lpv[0,p]+j*lpv[1,p]+k*lpv[2,p])
                dd=0
                for p in range(3):
                    dd += X[p] ** 2
                if dd<R**2 and X[2]>membrane_lo and X[2]<membrane_hi:
                    if dd<R**2+4*r**2-4*R*r:
                        atm = atm + 1
                        f.write('ATOM%7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n' % (atm, atm, X[0], X[1], X[2]))
                        
    f.close()              
    

def exANM(pdb_id, membrane_hi=None, membrane_lo=None, R=None, r=None, lat=None, outputName=None):
    if outputName == None:
        outputName = 'membrane.pdb'
    cgmembrane_wpdb(pdb_id)
    f = open('final.pdb', 'w')
    call(["cat",pdb_id,outputName], stdout=f)
    f.close()
    rt = parsePDB('final.pdb')
    anm, sel = calcANM(rt, selstr="name CA or name Q1")
    anm_red_membrane, sel_membrane = reduceModel(anm, rt, 'not chain Q')
    modes = anm_red_membrane.calcModes()
    call(["rm", "final.pdb"])
    call(["rm", outputName])
    return modes
    

def test(pdb='2nwl-mem.pdb', blk='2nwl.blk'):

    from prody import parsePDB
    from numpy import zeros, linalg

    pdb = parsePDB(pdb, subset='ca')
    pdb.setData('block', zeros(len(pdb), int))
    with open(blk) as inp:
        for line in inp:
            if line.startswith('BLOCK'):
                _, b, n1, c1, r1, n2, c2, r2 = line.split()
                sel = pdb.select('chain {} and resnum {} to {}'
                                 .format(c1, r1, r2))
                if sel:
                    sel.setData('block', int(b))
    pdb.setBetas(pdb.getData('block'))
    from prody import writePDB
    writePDB('pdb2gb1_truncated.pdb', pdb)
    rtb = RTB('2nwl')
    return rtb

