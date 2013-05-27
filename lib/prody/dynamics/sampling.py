# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
#
# Copyright (C) 2010-2012 Ahmet Bakan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines functions for generating alternate conformations along
normal modes


.. ipython:: python

   from prody import *
   from pylab import *
   ion()
   p38_pca = loadModel('p38_xray.pca.npz')
   p38_anm = loadModel('1p38.anm.npz')
   p38_ensemble = loadEnsemble('p38_X-ray.ens.npz')
   p38_structure = parsePDB('p38_ref_chain.pdb')"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.ensemble import Ensemble

from .nma import NMA
from .mode import Mode, VectorBase
from .modeset import ModeSet

__all__ = ['deformAtoms', 'sampleModes', 'traverseMode']


def sampleModes(modes, atoms=None, n_confs=1000, rmsd=1.0):
    """Return an ensemble of randomly sampled conformations along given
    *modes*.  If *atoms* are provided, sampling will be around its active
    coordinate set.  Otherwise, sampling is around the 0 coordinate set.

    :arg modes: Modes along which sampling will be performed.
    :type modes: :class:`.Mode`, :class:`.ModeSet`, :class:`.PCA`,
                 :class:`.ANM` or :class:`.NMA`

    :arg atoms: Atoms whose active coordinate set will be used as the initial
        conformation.
    :type atoms: :class:`.Atomic`

    :arg n_confs: Number of conformations to generate. Default is 1000.
    :type n_steps: int

    :arg rmsd: The average RMSD that the conformations will have with
        respect to the initial conformation. Default is 1.0 A.
    :type rmsd: float

    :returns: :class:`.Ensemble`

    For given normal modes :math:`[u_1 u_2 ... u_m]` and their eigenvalues
    :math:`[\lambda_1 \lambda_2 ... \lambda_m]`, a new conformation
    is sampled using the relation:

    .. math::

       R_k = R_0 + s \sum_{i=1}^{m} r_i^k \lambda^{-0.5}_i u_i

    :math:`R_0` is the active coordinate set of *atoms*.
    :math:`[r_1^k r_2^k ... r_m^k]` are normally distributed random numbers
    generated for conformation :math:`k` using :func:`numpy.random.randn`.

    RMSD of the new conformation from :math:`R_0` can be calculated as

    .. math::

      RMSD^k = \sqrt{ {\\left( s \sum_{i=1}^{m} r_i^k \lambda^{-0.5}_i u_i  \\right)}^{2} / N } = \\frac{s}{ \sqrt{N}} \sqrt{ \sum_{i=1}^{m} (r_i^k)^2 \lambda^{-1}_i  }


    Average :math:`RMSD` of the generated conformations from the initial conformation is:

    .. math::

      \\left< RMSD^k \\right> = \\frac{s}{ \sqrt{N}} \\left< \sqrt{ \sum_{i=1}^{m} (r_i^k)^2 \lambda^{-1}_i } \\right>


    From this relation :math:`s` scaling factor obtained using the relation

    .. math::

       s =  \\left< RMSD^k \\right> \sqrt{N} {\\left< \sqrt{ \sum_{i=1}^{m} (r_i)^2 \lambda^{-1}_i} \\right>}^{-1}


    Note that random numbers are generated before conformations are
    sampled, hence exact value of :math:`s` is known from this relation to
    ensure that the generated ensemble will have user given average *rmsd*
    value.

    Note that if modes are from a :class:`.PCA`, variances are used instead of
    inverse eigenvalues, i.e. :math:`\sigma_i \sim \lambda^{-1}_i`.

    See also :func:`.showEllipsoid`.

    .. ipython:: python

       # Generate 300 conformations using ANM modes 1-3
       ensemble = sampleModes(p38_anm[:3], n_confs=500)
       # Project these conformations onto the space spanned by these modes
       @savefig reference_dynamics_sampling_ensemble.png width=4in
       showProjection(ensemble, p38_anm[:3], rmsd=True);"""

    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a NMA or ModeSet instance, '
                        'not {0}'.format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be from a 3-dimensional model')
    n_confs = int(n_confs)
    n_atoms = modes.numAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, (Atomic)):
            raise TypeError('{0} is not correct type for atoms'
                            .format(type(atoms)))
        if atoms.numAtoms() != n_atoms:
            raise ValueError('number of atoms do not match')
        initial = atoms.getCoords()

    rmsd = float(rmsd)
    LOGGER.info('Parameter: rmsd = {0:.2f} A'.format(rmsd))
    n_confs = int(n_confs)
    LOGGER.info('Parameter: n_confs = {0}'.format(n_confs))

    if isinstance(modes, Mode):
        n_modes = 1
        variances = np.array([modes.getVariance()])
    else:
        n_modes = len(modes)
        variances = modes.getVariances()
    if np.any(variances == 0):
        raise ValueError('one or more modes has zero variance')
    randn = np.random.standard_normal((n_confs, n_modes))
    coef = ((randn ** 2 * variances).sum(1) ** 0.5).mean()
    scale = n_atoms**0.5 * rmsd / coef

    LOGGER.info('Modes are scaled by {0}.'.format(scale))

    confs = []
    append = confs.append
    scale = scale * variances ** 0.5
    array = modes._getArray()
    if array.ndim > 1:
        for i in range(n_confs):
            append((array * scale * randn[i]).sum(1).reshape((n_atoms, 3)))
    else:
        for i in range(n_confs):
            append((array * scale * randn[i]).reshape((n_atoms, 3)))

    ensemble = Ensemble('Conformations along {0}'.format(modes))
    if initial is None:
        ensemble.setCoords(np.zeros((n_atoms, 3)))
        ensemble.addCoordset(np.array(confs))
    else:
        ensemble.setCoords(initial)
        ensemble.addCoordset(np.array(confs) + initial)
    return ensemble


def traverseMode(mode, atoms, n_steps=10, rmsd=1.5):
    """Generates a trajectory along a given *mode*, which can be used to
    animate fluctuations in an external program.

    :arg mode: Mode along which a trajectory will be generated.
    :type mode: :class:`.Mode`

    :arg atoms: Atoms whose active coordinate set will be used as the initial
        conformation.
    :type atoms: :class:`.Atomic`

    :arg n_steps: Number of steps to take along each direction.
        For example, for ``n_steps=10``, 20 conformations will be
        generated along the first mode. Default is 10.
    :type n_steps: int

    :arg rmsd: The maximum RMSD that the conformations will have with
        respect to the initial conformation. Default is 1.5 A.
    :type rmsd: float

    :returns: :class:`.Ensemble`

    For given normal mode :math:`u_i`, its eigenvalue
    :math:`\lambda_i`, number of steps :math:`n`, and maximum :math:`RMSD`
    conformations :math:`[R_{-n} R_{-n+1} ... R_{-1} R_0 R_1 ... R_n]` are
    generated.

    :math:`R_0` is the active coordinate set of *atoms*.
    :math:`R_k = R_0 + sk\lambda_iu_i`, where :math:`s` is found using
    :math:`s = ((N (\\frac{RMSD}{n})^2) / \lambda_i^{-1}) ^{0.5}`, where
    :math:`N` is the number of atoms.


    .. ipython:: python

       trajectory = traverseMode(p38_anm[0], p38_structure.select('calpha'),
                                 n_steps=8, rmsd=1.4)
       rmsd = calcRMSD(trajectory)
       plot(rmsd, '-o');
       xlabel('Frame index');
       @savefig reference_dynamics_sampling_traverse.png width=4in
       ylabel('RMSD (A)');"""

    if not isinstance(mode, VectorBase):
        raise TypeError('mode must be a Mode or Vector instance, '
                        'not {0}'.format(type(mode)))
    if not mode.is3d():
        raise ValueError('mode must be from a 3-dimensional model.')
    n_atoms = mode.numAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, Atomic):
            raise TypeError('{0} is not correct type for atoms'
                            .format(type(atoms)))
        if atoms.numAtoms() != n_atoms:
            raise ValueError('number of atoms do not match')
        initial = atoms.getCoords()

    name = str(mode)

    rmsd = float(rmsd) + 0.000004
    LOGGER.info('Parameter: rmsd = {0:.2f} A'.format(rmsd))
    n_steps = int(n_steps)
    LOGGER.info('Parameter: n_steps = {0}'.format(n_steps))
    step = rmsd / n_steps
    LOGGER.info('Step size is {0:.2f} A RMSD'.format(step))
    arr = mode.getArrayNx3()
    var = mode.getVariance()
    scale = ((n_atoms * step**2) / var) ** 0.5
    LOGGER.info('Mode is scaled by {0}.'.format(scale))

    array = arr * var**0.5 * scale
    confs_add = [initial + array]
    for s in range(1, n_steps):
        confs_add.append(confs_add[-1] + array)
    confs_sub = [initial - array]
    for s in range(1, n_steps):
        confs_sub.append(confs_sub[-1] - array)
    confs_sub.reverse()
    ensemble = Ensemble('Conformations along {0}'.format(name))
    ensemble.setCoords(initial)
    ensemble.addCoordset(np.array(confs_sub + [initial] + confs_add))
    return ensemble


def deformAtoms(atoms, mode, rmsd=None):
    """Generate a new coordinate set for *atoms* along the *mode*.  *atoms*
    must be a :class:`.AtomGroup` instance.  New coordinate set will be
    appended to *atoms*. If *rmsd* is provided, *mode* will be scaled to
    generate a coordinate set with given RMSD distance to the active coordinate
    set.  Below example shows how to deform a structure along a normal mode
    or linear combinations of normal modes:

    .. ipython:: python

       deformAtoms(p38_structure, p38_pca[0] * p38_pca[0].getVariance()**0.5)
       deformAtoms(p38_structure, -p38_pca[1] * p38_pca[1].getVariance()**0.5)
       deformAtoms(p38_structure, p38_pca[0] * p38_pca[0].getVariance()**0.5 +
                   p38_pca[1] * p38_pca[1].getVariance()**0.5)
       deformAtoms(p38_structure, p38_pca[0], rmsd=1.0)
       calcRMSD(p38_structure)"""

    if not isinstance(atoms, AtomGroup):
        raise TypeError('atoms must be an AtomGroup, not {0}'
                        .format(type(atoms)))
    if not isinstance(mode, VectorBase):
        raise TypeError('mode must be a Mode or Vector instance, '
                        'not {0}'.format(type(mode)))
    if not mode.is3d():
        raise ValueError('mode must be from a 3-dimensional model.')
    if atoms.numAtoms() != mode.numAtoms():
        raise ValueError('number of atoms do not match')

    array = mode.getArrayNx3()

    if rmsd is not None:
        rmsd = float(rmsd)
        # rmsd = ( ((scalar * array)**2).sum() / n_atoms )**0.5
        scalar = (atoms.numAtoms() * rmsd**2 / (array**2).sum())**0.5
        LOGGER.info('Mode is scaled by {0}.'.format(scalar))
        atoms.addCoordset(atoms.getCoords() + array * scalar)
    else:
        atoms.addCoordset(atoms.getCoords() + array)
