# -*- coding: utf-8 -*-
"""This module defines a class and a function for rotating translating blocks
(RTB) calculations."""

import numpy as np
import scipy as sp

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB
from prody.utilities import importLA, checkCoords, sqrtm
from numpy import sqrt, zeros, linalg, min, max, unique, mean, eye, outer, dot
from scipy import sparse
from subprocess import call

from .anm import ANMBase, calcANM, ANM
from .gnm import checkENMParameters
from .editing import reduceModel

__all__ = ['RTB']

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

	:arg membrane_low: minimum z-coordinate at which membrane scaling
            is applied
            default is 1.0
	:type membrane_low: float

	:arg membrane_high: maximum z-coordinate at which membrane scaling
            is applied.  If membrane_high < membrane_low, scaling will be 
	    applied to the entire structure
            default is -1.0
         :type membrane_high: float
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

        #coords = coords.T.copy()

        #self._hessian = hessian = np.zeros((nb6, nb6), float)
        #self._project = project = np.zeros((natoms * 3, nb6), float)

        hessian, project = buildBlockHessian(coords, blocks,
                     natoms, nblocks, nones, cutoff, gamma)

        self._hessian = hessian
        self._project = project
        self._dof = self._hessian.shape[0]
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
        n_modes = self._dof
        super(RTB, self).calcModes(n_modes, zeros, turbo)
        print self._project.shape
        print self._array.shape
        #self._project = self._project.T.copy()
        print self._project.shape
        self._array = np.dot(self._project, self._array)

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
    coords = pdb.getCoords() 
    blocks = pdb.getBetas()
    from prody import writePDB
    writePDB('pdb2gb1_truncated.pdb', pdb)
    rtb = RTB('2nwl')
    rtb.buildHessian(coords, blocks, scale=64)
    #rtb.calcModes()
    return rtb

def buildBlockHessian(coords, blocks, natoms, nblocks, nones, cutoff, gamma):
    project = zeros((3*natoms,6*nblocks-3*nones))
        
    block_uniq = unique(blocks)
    block_res = []
    block_sizes = []
    for i in range(nblocks):
        block_res.append(blocks==block_uniq[i])
    cm = []
    for i in range(nblocks):         
        cm.append(mean(coords[block_res[i]],0))
    block_coords = []
    for i in range(nblocks):
        block_coords.append(coords[block_res[i]])
        block_sizes.append(block_coords[i].shape[0])
    block_coords_tr = block_coords[:]
    for i in range(nblocks):
        block_coords_tr[i]=block_coords[i]-cm[i]
    I=[]
    Isqt=[]
    for i in range(nblocks):
        if block_sizes[i]==1:
            I.append(eye(3))
        else:
            dum = zeros((3,3))
            for k in range(block_coords[i].shape[0]):
                dd = linalg.norm(block_coords_tr[i][k])
                for j in range(3):
                    dum[j,j]=(dd-block_coords_tr[i][k,j]**2)
                    for l in range(j+1,3):
                        dum[j,l]-=block_coords_tr[i][k,j]*block_coords_tr[i][k,l]
                        dum[l,j]=dum[j,l]
            I.append(dum)
        Isqt.append(linalg.inv(sqrtm(I[i])))

    none_passed=0
    for k in range(natoms):
        for l in range(nblocks):
            if blocks[k]==block_uniq[l]:
                if block_sizes[l]==1:
                    project[3*k:3*k+3,l*6-none_passed*3:l*6-none_passed*3+3]=eye(3)
                    none_passed+=1
                else:
                    project[3*k:3*k+3,l*6-none_passed*3:l*6-none_passed*3+3]=eye(3)/block_sizes[l]
                    project[3*k,l*6-none_passed*3+1]=Isqt[l][0,1]*(coords[k,2]-cm[l][2])-Isqt[l][0,2]*(coords[k,1]-cm[l][1])
                    project[3*k,l*6-none_passed*3+2]=Isqt[l][0,2]*(coords[k,1]-cm[l][1])-Isqt[l][0,1]*(coords[k,2]-cm[l][2])
                    project[3*k+1,l*6-none_passed*3]=-Isqt[l][1,0]*(coords[k,2]-cm[l][2])+Isqt[l][1,2]*(coords[k,0]-cm[l][0])
                    project[3*k+1,l*6-none_passed*3+2]=Isqt[l][1,0]*(coords[k,2]-cm[l][2])-Isqt[l][1,2]*(coords[k,0]-cm[l][0])
                    project[3*k+2,l*6-none_passed*3]=Isqt[l][2,0]*(coords[k,1]-cm[l][1])-Isqt[l][2,1]*(coords[k,0]-cm[l][0])
                    project[3*k+2,l*6-none_passed*3+1]=-Isqt[l][2,0]*(coords[k,1]-cm[l][1])+Isqt[l][2,1]*(coords[k,0]-cm[l][0])

    cutoff, g, gamma = checkENMParameters(cutoff, gamma)
    total_hessian = zeros((natoms * 3, natoms * 3))               
    cutoff2 = cutoff * cutoff
    for i in range(natoms):
        res_i3 = i*3
        res_i33 = res_i3+3
        i_p1 = i+1
        i2j_all = coords[i_p1:, :] - coords[i]
        for j, dist2 in enumerate((i2j_all ** 2).sum(1)):
            if dist2 > cutoff2:
                continue
            i2j = i2j_all[j]
            j += i_p1
            g = gamma(dist2, i, j)
            res_j3 = j*3
            res_j33 = res_j3+3
            super_element = np.outer(i2j, i2j) * (- g / dist2)
            total_hessian[res_i3:res_i33, res_j3:res_j33] = super_element
            total_hessian[res_j3:res_j33, res_i3:res_i33] = super_element
            total_hessian[res_i3:res_i33, res_i3:res_i33] = \
                total_hessian[res_i3:res_i33, res_i3:res_i33] - super_element
            total_hessian[res_j3:res_j33, res_j3:res_j33] = \
                total_hessian[res_j3:res_j33, res_j3:res_j33] - super_element

    project = sp.sparse.coo_matrix(project)
    total_hessian = sp.sparse.coo_matrix(total_hessian)
    hessian = dot(dot(project.transpose(),total_hessian), project)
    hessian = hessian.todense()
    project = project.todense()
    return (hessian, project)
