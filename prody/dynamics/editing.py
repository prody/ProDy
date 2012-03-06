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

"""This module defines functions for editing normal mode data."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody.atomic import Atomic, AtomGroup, AtomMap, AtomSubset
from prody.atomic import Selection, SELECT
from prody.tools import importLA

from nma import NMA
from modeset import ModeSet
from mode import VectorBase, Mode, Vector
from gnm import GNM
from anm import ANM
from pca import PCA

__all__ = ['extendModel', 
           'sliceMode', 'sliceModel', 'sliceVector',
           'reduceModel',]

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

def extendModel(model, nodes, atoms):
    """Extend a coarse grained *model* built for *nodes* to *atoms*.  *model*
    may be :class:`~.ANM`, :class:`~.GNM`, :class:`~.PCA`, or :class:`~.NMA` 
    instance.  This will take part of the normal modes for each node (i.e. Cα 
    atoms) and extend it to all other atoms in the same residue.  For each atom
    in *nodes* argument *atoms* argument must contain a corresponding residue.  
    Note that modes in the extended model will not be normalized.  For a usage
    example see :ref:`extendmodel`."""
    
    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance')
    if not isinstance(nodes, Atomic):
        raise TypeError('nodes must be an Atomic instance')
    if model.numAtoms() != nodes.numAtoms():
        raise ValueError('model and nodes must have same number of atoms')
    
    if isinstance(atoms, Atomic):
        is3d = model.is3d()            
        atom_indices = []
        indices = []
        hierview = atoms.getHierView()
        for i, node in enumerate(nodes):
            res = hierview[node.getChid(), node.getResnum(), node.getIcode()]
            if res is None:
                raise ValueError('hierview must contain a residue for all atoms')
            atom_indices.append(res.getIndices())
            if is3d:
                indices.append(range(i*3, (i+1)*3) * len(res))
            else:
                indices.append([i] * len(res))
        atom_indices = np.concatenate(atom_indices)
        indices = np.concatenate(indices)
        
        array = model.getArray()[indices,:]
        extended = NMA('Extended ' + str(model))
        extended.setEigens(array, model.getEigenvalues())
        if isinstance(atoms, AtomGroup):
            ag = atoms
        else: 
            ag = atoms.getAtomGroup()
        atommap = AtomMap(ag, atom_indices, np.arange(len(atom_indices)), 
                          np.array([]), str(atoms), atoms.getACSIndex())
        return extended, atommap
    else:
        raise TypeError('atoms must be an Atomic instance')
    

def sliceVector(vector, atoms, selstr):
    """Return part of the *vector* for *atoms* matching *selstr*.  Note that 
    returned :class:`~.Vector` instance is not normalized.
    
    :arg vector: vector instance to be sliced
    :type vector: :class:`~.VectorBase`
    
    :arg atoms: atoms for which *vector* describes a deformation, motion, etc.
    :type atoms: :class:`~.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: (:class:`~.Vector`, :class:`~.Selection`)"""
    
    if not isinstance(vector, VectorBase):
        raise TypeError('vector must be a VectorBase instance, not {0:s}'
                        .format(type(vector)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.numAtoms() != vector.numAtoms(): 
        raise ValueError('number of atoms in *vector* and *atoms* must be '
                         'equal')
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getACSIndex())
    vec = Vector(vector.getArrayNx3()[
                 which, :].flatten(),
                 '{0:s} slice {1:s}'.format(str(vector), repr(selstr)), 
                 vector.is3d())
    return (vec, sel)

def sliceMode(mode, atoms, selstr):
    """Return part of the *mode* for *atoms* matching *selstr*.  This works 
    slightly different from :func:`~.sliceVector`. Mode array (eigenvector) is 
    multiplied by square-root of the variance along the mode.  If mode is from
    an elastic network model, variance is defined as the inverse of the 
    eigenvalue.  Note that returned :class:`~.Vector` instance is not 
    normalized.
    
    :arg mode: mode instance to be sliced
    :type mode: :class:`~.Mode`
    
    :arg atoms: atoms for which *mode* describes a deformation, motion, etc.
    :type atoms: :class:`~.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: (:class:`~.Vector`, :class:`~.Selection`)"""
    
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, not {0:s}'
                        .format(type(mode)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.numAtoms() != mode.numAtoms(): 
        raise ValueError('number of atoms in *mode* and *atoms* must be equal')
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getACSIndex())
    vec = Vector(mode.getArrayNx3()[
                 which,:].flatten() * mode.getVariance()**0.5,
                 '{0:s} slice {1:s}'.format(str(mode), repr(selstr)), 
                 mode.is3d()) 
    return (vec, sel)

def sliceModel(model, atoms, selstr):
    """Return a part of the *model* for *atoms* matching *selstr*.  Note that 
    normal modes (eigenvectors) are not normalized.
    
    :arg mode: NMA model instance to be sliced
    :type mode: :class:`~.NMA`
    
    :arg atoms: atoms for which the *model* was built
    :type atoms: :class:`~.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: (:class:`~.NMA`, :class:`~.Selection`)"""
    
    if not isinstance(model, NMA):
        raise TypeError('mode must be a NMA instance, not {0:s}'
                        .format(type(model)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.numAtoms() != model.numAtoms(): 
        raise ValueError('number of atoms in *model* and *atoms* must be '
                         'equal')
    
    array = model._getArray()
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getACSIndex())

    nma = type(model)('{0:s} slice {1:s}'
                      .format(model.getTitle(), repr(selstr)))
    if model.is3d():
        which = [which.reshape((len(which),1))*3]
        which.append(which[0]+1)
        which.append(which[0]+2)
        which = np.concatenate(which, 1).flatten()
    nma.setEigens( array[which, :], model.getEigenvalues() )
    return (nma, sel)
    
def reduceModel(model, atoms, selstr):
    """Return reduced NMA model.  Reduces a :class:`~.NMA` model to a subset of 
    *atoms* matching a selection *selstr*.  This function behaves differently 
    depending on the type of the *model* argument.  For :class:`~.ANM` and 
    :class:`~.GNM` or other :class:`~.NMA` models, this functions derives the 
    force constant matrix for system of interest (specified by the *selstr*) 
    from the force constant matrix for the *model* by assuming that for any 
    given displacement of the system of interest, the other atoms move along in
    such a way as to minimize the potential energy.  This is based on the 
    formulation in in [KH00]_.  For :class:`~.PCA` models, this function simply
    takes the sub-covariance matrix for the selected atoms.

    :arg model: dynamics model
    :type model: :class:`~.ANM`, :class:`~.GNM`, or :class:`~.PCA`
    :arg atoms: atoms that were used to build the model
    :arg selstr: a selection string specifying subset of atoms"""
    
    linalg = importLA()

    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance, not {0:s}'.format(type(model)))
    if not isinstance(atoms, (AtomGroup, AtomSubset, AtomMap)):
        raise TypeError('atoms type is not valid')
    if len(atoms) <= 1:
        raise TypeError('atoms must contain more than 1 atoms')

    if isinstance(model, GNM):
        matrix = model._kirchhoff
    elif isinstance(model, ANM):
        matrix = model._hessian
    elif isinstance(model, PCA):
        matrix = model._cov
    else:
        raise TypeError('model does not have a valid type derived from NMA')
    if matrix is None:
        raise ValueError('model matrix (Hessian/Kirchhoff/Covariance) is not '
                         'built')

    system = SELECT.getBoolArray(atoms, selstr)
    other = np.invert(system)
    n_sel = sum(system) 
    if n_sel == 0:
        LOGGER.warning('selection has 0 atoms')
        return None
    if len(atoms) == n_sel:
        LOGGER.warning('selection results in same number of atoms, '
                       'model is not reduced')
        return None

    if model.is3d():
        system = np.tile(system, (3,1)).transpose().flatten()
        other = np.tile(other, (3,1)).transpose().flatten()
    ss = matrix[system,:][:,system]
    if isinstance(model, PCA):
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(ss)
        return eda, system
    so = matrix[system,:][:,other]
    os = matrix[other,:][:,system]
    oo = matrix[other,:][:,other]
    matrix = ss - np.dot(so, np.dot(linalg.inv(oo), os))
    
    if isinstance(model, GNM):
        gnm = GNM(model.getTitle() + ' reduced')
        gnm.setKirchhoff(matrix)
        return gnm, system
    elif isinstance(model, ANM):
        anm = ANM(model.getTitle() + ' reduced')
        anm.setHessian(matrix)
        return anm, system
    elif isinstance(model, PCA):
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(matrix)
        return eda, system
