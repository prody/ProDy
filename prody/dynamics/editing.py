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

from prody.atomic import Atomic, AtomGroup, AtomMap, AtomSubset, HierView
from prody.atomic import Selection, SELECT
from prody.utilities import importLA
from prody import LOGGER, _PY3K

from nma import NMA
from modeset import ModeSet
from mode import VectorBase, Mode, Vector
from gnm import GNM
from anm import ANM
from pca import PCA

if not _PY3K:
    range = xrange

__all__ = ['extendModel', 'extendMode', 'extendVector',
           'sliceMode', 'sliceModel', 'sliceVector',
           'reduceModel',]

def extend(model, nodes, atoms):
    """Return mapping indices and an :class:`.AtomMap`."""

    try:        
        n_atoms = model.numAtoms()
        is3d = model.is3d()
    except AttributeError:
        raise ValueError('model must be an NMA instance')
    
    try:        
        n_nodes = nodes.numAtoms()
        i_nodes = nodes.iterAtoms()
    except AttributeError:
        raise ValueError('nodes must be an Atomic instance')
    
    if n_atoms != n_nodes:
        raise ValueError('atom numbers must be the same')

    if not nodes in atoms:
        raise ValueError('nodes must be a subset of atoms')
    
                
    atom_indices = []
    indices = []
    hierview = HierView(atoms)._dict
    
    for i, node in enumerate(i_nodes):
        res = hierview.pop((node.getSegname() or None, node.getChid() or None, 
                           node.getResnum(), node.getIcode() or None), None)
        if res is None:
            raise ValueError('atoms must contain a residue for all atoms')
        atom_indices.append(res._getIndices())
        if is3d:
            indices.append(list(range(i*3, (i+1)*3)) * len(res))
        else:
            indices.append([i] * len(res))
    atom_indices = np.concatenate(atom_indices)
    indices = np.concatenate(indices)
    
    try:
        ag = atoms.getAtomGroup()
    except AttributeError: 
        ag = atoms
    atommap = AtomMap(ag, atom_indices, atoms.getACSIndex(), 
                      title=str(atoms), intarrays=True)
    return indices, atommap


def extendModel(model, nodes, atoms):
    """Extend a coarse grained *model* built for *nodes* to *atoms*.  *model*
    may be :class:`~.ANM`, :class:`~.GNM`, :class:`~.PCA`, or :class:`~.NMA` 
    instance.  This function will take part of the normal modes for each node 
    (i.e. Cα atoms) and extend it to all other atoms in the same residue.  For 
    each atom in *nodes* argument *atoms* argument must contain a corresponding
    residue.  Note that modes in the extended model will not be normalized.  
    For a usage example see :ref:`extendmodel`."""
    
    try:        
        evecs = model._getArray()
        evals = model.getEigvals()
    except AttributeError:
        raise ValueError('model must be an NMA instance')
    
    indices, atommap = extend(model, nodes, atoms)
    
    evecs = evecs[indices, :]
    if model.is3d():
        extended = NMA('Extended ' + str(model))
    else:
        extended = GNM('Extended ' + str(model))
    extended.setEigens(evecs, evals)
    return extended, atommap

    
def extendMode(mode, nodes, atoms):
    """Extend a coarse grained normal *mode* built for *nodes* to *atoms*.  
    This function will take part of the normal modes for each node (i.e. Cα 
    atoms) and extend it to all other atoms in the same residue.  For each atom
    in *nodes* argument *atoms* argument must contain a corresponding residue.
    """
    
    try:        
        vec = mode._getArray()
        std = mode.getVariance() ** 0.5
    except AttributeError:
        raise ValueError('mode must be a normal Mode instance')
    
    indices, atommap = extend(mode, nodes, atoms)
    extended = Vector(vec[indices] * std, 'Extended ' + str(mode), mode.is3d())
    return extended, atommap


def extendVector(vector, nodes, atoms):
    """Extend a coarse grained *vector* for *nodes* to *atoms*.  This function
    will take part of the normal modes for each node (i.e. Cα atoms) and extend
    it to all other atoms in the same residue.  For each atom in *nodes*,
    *atoms* argument must contain a corresponding residue."""
    
    try:        
        vec = vector._getArray()
    except AttributeError:
        raise ValueError('vector must be a Vector instance')
    
    indices, atommap = extend(vector, nodes, atoms)
    extended = Vector(vec[indices], 'Extended ' + str(vector), vector.is3d())
    return extended, atommap


def sliceVector(vector, atoms, select):
    """Return part of the *vector* for *atoms* matching *select*.  Note that 
    returned :class:`.Vector` instance is not normalized.
    
    :arg vector: vector instance to be sliced
    :type vector: :class:`.VectorBase`
    
    :arg atoms: atoms for which *vector* describes a deformation, motion, etc.
    :type atoms: :class:`.Atomic`
    
    :arg select: an atom selection or a selection string 
    :type select: :class:`.Selection`, str 
    
    :returns: (:class:`.Vector`, :class:`.Selection`)"""
    
    if not isinstance(vector, VectorBase):
        raise TypeError('vector must be a VectorBase instance, not {0:s}'
                        .format(type(vector)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.numAtoms() != vector.numAtoms(): 
        raise ValueError('number of atoms in model and atoms must be equal')
    
    if isinstance(select, str):
        selstr = select
        if isinstance(atoms, AtomGroup):
            sel = atoms.select(selstr)
            which = sel._getIndices()
        else:
            which = SELECT.getIndices(atoms, selstr)
            sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                            selstr, atoms.getACSIndex())
            
    elif isinstance(select, AtomSubset):
        sel = select
        if isinstance(atoms, AtomGroup):
            if sel.getAtomGroup() != atoms:
                raise ValueError('select and atoms do not match')
            which = sel._getIndices()
        else:
            if atoms.getAtomGroup() != sel.getAtomGroup():
                raise ValueError('select and atoms do not match')
            elif not sel in atoms:
                raise ValueError('select is not a subset of atoms')
            idxset = set(atoms._getIndices())
            which = np.array([idx in idxset for idx in sel._getIndices()])
            which = which.nonzero()[0]
        selstr = sel.getSelstr()
    
    else:
        raise TypeError('select must be a string or a Selection instance')
        
    vec = Vector(vector.getArrayNx3()[
                 which, :].flatten(),
                 '{0:s} slice {1:s}'.format(str(vector), repr(selstr)), 
                 vector.is3d())
    return (vec, sel)


def sliceMode(mode, atoms, select):
    """Return part of the *mode* for *atoms* matching *select*.  This works 
    slightly different from :func:`.sliceVector`. Mode array (eigenvector) is 
    multiplied by square-root of the variance along the mode.  If mode is from
    an elastic network model, variance is defined as the inverse of the 
    eigenvalue.  Note that returned :class:`~.Vector` instance is not 
    normalized.
    
    :arg mode: mode instance to be sliced
    :type mode: :class:`.Mode`
    
    :arg atoms: atoms for which *mode* describes a deformation, motion, etc.
    :type atoms: :class:`.Atomic`
    
    :arg select: an atom selection or a selection string 
    :type select: :class:`.Selection`, str 
    
    :returns: (:class:`~.Vector`, :class:`~.Selection`)"""
    
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, not {0:s}'
                        .format(type(mode)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.numAtoms() != mode.numAtoms(): 
        raise ValueError('number of atoms in model and atoms must be equal')
    
    if isinstance(select, str):
        selstr = select
        if isinstance(atoms, AtomGroup):
            sel = atoms.select(selstr)
            which = sel._getIndices()
        else:
            which = SELECT.getIndices(atoms, selstr)
            sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                            selstr, atoms.getACSIndex())
        
    elif isinstance(select, AtomSubset):
        sel = select
        if isinstance(atoms, AtomGroup):
            if sel.getAtomGroup() != atoms:
                raise ValueError('select and atoms do not match')
            which = sel._getIndices()
        else:
            if atoms.getAtomGroup() != sel.getAtomGroup():
                raise ValueError('select and atoms do not match')
            elif not sel in atoms:
                raise ValueError('select is not a subset of atoms')
            idxset = set(atoms._getIndices())
            which = np.array([idx in idxset for idx in sel._getIndices()])
            which = which.nonzero()[0]
        selstr = sel.getSelstr()
    
    else:
        raise TypeError('select must be a string or a Selection instance')
    
    vec = Vector(mode.getArrayNx3()[
                 which,:].flatten() * mode.getVariance()**0.5,
                 '{0:s} slice {1:s}'.format(str(mode), repr(selstr)), 
                 mode.is3d()) 
    return (vec, sel)


def sliceModel(model, atoms, select):
    """Return a part of the *model* for *atoms* matching *select*.  Note that 
    normal modes (eigenvectors) are not normalized.
    
    :arg mode: NMA model instance to be sliced
    :type mode: :class:`.NMA`
    
    :arg atoms: atoms for which the *model* was built
    :type atoms: :class:`.Atomic`
    
    :arg select: an atom selection or a selection string 
    :type select: :class:`.Selection`, str 
    
    :returns: (:class:`.NMA`, :class:`.Selection`)"""
    
    if not isinstance(model, NMA):
        raise TypeError('mode must be a NMA instance, not {0:s}'
                        .format(type(model)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.numAtoms() != model.numAtoms(): 
        raise ValueError('number of atoms in model and atoms must be equal')
    
    array = model._getArray()
    
    if isinstance(select, str):
        selstr = select
        if isinstance(atoms, AtomGroup):
            sel = atoms.select(selstr)
            which = sel.getIndices()
        else:
            which = SELECT.getIndices(atoms, selstr)
            sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                            selstr, atoms.getACSIndex())
        
    elif isinstance(select, AtomSubset):
        sel = select
        if isinstance(atoms, AtomGroup):
            if sel.getAtomGroup() != atoms:
                raise ValueError('select and atoms do not match')
            which = sel._getIndices()
        else:
            if atoms.getAtomGroup() != sel.getAtomGroup():
                raise ValueError('select and atoms do not match')
            elif not sel in atoms:
                raise ValueError('select is not a subset of atoms')
            idxset = set(atoms._getIndices())
            which = np.array([idx in idxset for idx in sel._getIndices()])
            which = which.nonzero()[0]
        selstr = sel.getSelstr()
    
    else:
        raise TypeError('select must be a string or a Selection instance')
        
    nma = type(model)('{0:s} slice {1:s}'
                      .format(model.getTitle(), repr(selstr)))
    if model.is3d():
        which = [which.reshape((len(which),1))*3]
        which.append(which[0]+1)
        which.append(which[0]+2)
        which = np.concatenate(which, 1).flatten()
    nma.setEigens( array[which, :], model.getEigvals() )
    return (nma, sel)
    
    
def reduceModel(model, atoms, select):
    """Return reduced NMA model.  Reduces a :class:`~.NMA` model to a subset of 
    *atoms* matching *select*.  This function behaves differently depending on 
    the type of the *model* argument.  For :class:`.ANM` and :class:`.GNM` or 
    other :class:`.NMA` models, force constant matrix for system of interest 
    (specified by the *select*) is derived from the force constant matrix for 
    the *model* by assuming that for any given displacement of the system of 
    interest, other atoms move along in such a way as to minimize the potential
    energy.  This is based on the formulation in [KH00]_.  For :class:`.PCA` 
    models, this function simply takes the sub-covariance matrix for selection.

    :arg model: dynamics model
    :type model: :class:`.ANM`, :class:`.GNM`, or :class:`.PCA`
    
    :arg atoms: atoms that were used to build the model
    :type atoms: :class:`.Atomic`
    
    :arg select: an atom selection or a selection string 
    :type select: :class:`.Selection`, str 
    
    :returns: (:class:`.NMA`, :class:`.Selection`)"""
    
    linalg = importLA()

    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance, not {0:s}'
                        .format(type(model)))
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

    if isinstance(select, str):
        system = SELECT.getBoolArray(atoms, select)
        n_sel = sum(system)
        if n_sel == 0:
            raise ValueError('select matches 0 atoms')
        if len(atoms) == n_sel:
            raise ValueError('select matches all atoms')

        if isinstance(atoms, AtomGroup):
            ag = atoms
            which = np.arange(len(atoms))[system]
        else:
            ag = atoms.getAtomGroup()
            which = atoms._getIndices()[system]
        sel = Selection(ag, which, select, atoms.getACSIndex())
        
    elif isinstance(select, AtomSubset):
        sel = select
        if isinstance(atoms, AtomGroup):
            if sel.getAtomGroup() != atoms:
                raise ValueError('select and atoms do not match')
            system = np.zeros(len(atoms), bool)
            system[sel._getIndices()] = True 
        else:
            if atoms.getAtomGroup() != sel.getAtomGroup():
                raise ValueError('select and atoms do not match')
            elif not sel in atoms:
                raise ValueError('select is not a subset of atoms')
            idxset = set(atoms._getIndices())
            system = np.array([idx in idxset for idx in sel._getIndices()])
    
    else:
        raise TypeError('select must be a string or a Selection instance')
    
    other = np.invert(system)

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
        return gnm, sel
    elif isinstance(model, ANM):
        anm = ANM(model.getTitle() + ' reduced')
        anm.setHessian(matrix)
        return anm, sel
    elif isinstance(model, PCA):
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(matrix)
        return eda, sel
