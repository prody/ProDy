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

"""This module defines a functions for handling conformational ensembles."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

from prody.proteins import fetchPDB, parsePDB, writePDB
from prody.utilities import openFile
from prody import LOGGER, SETTINGS

from ensemble import *
from pdbensemble import *
from conformation import *

__all__ = ['saveEnsemble', 'loadEnsemble', 'trimPDBEnsemble', 
           'calcOccupancies', 'showOccupancies', 'alignPDBEnsemble']


def saveEnsemble(ensemble, filename=None, **kwargs):
    """Save *ensemble* model data as :file:`filename.ens.npz`.  If *filename* 
    is ``None``, title of the *ensemble* will be used as the filename, after 
    white spaces in the title are replaced with underscores.  Extension is 
    :file:`.ens.npz`. Upon successful completion of saving, filename is 
    returned. This function makes use of :func:`numpy.savez` function."""
    
    if not isinstance(ensemble, Ensemble):
        raise TypeError('invalid type for ensemble, {0:s}'
                        .format(type(ensemble)))
    if len(ensemble) == 0:
        raise ValueError('ensemble instance does not contain data')
    
    dict_ = ensemble.__dict__
    attr_list = ['_title', '_confs', '_weights', '_coords']
    if isinstance(ensemble, PDBEnsemble):
        attr_list.append('_labels')
        attr_list.append('_trans')
    if filename is None:
        filename = ensemble.getTitle().replace(' ', '_')
    attr_dict = {}
    for attr in attr_list:
        value = dict_[attr]
        if value is not None:
            attr_dict[attr] = value
    filename += '.ens.npz'
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename


def loadEnsemble(filename):
    """Return ensemble instance loaded from *filename*.  This function makes 
    use of :func:`numpy.load` function.  See also :func:`saveEnsemble`"""
    
    attr_dict = np.load(filename)
    weights = attr_dict['_weights']
    isPDBEnsemble = False
    try:
        title = str(attr_dict['_title'])
    except KeyError:
        title = str(attr_dict['_name'])
    if weights.ndim == 3:
        isPDBEnsemble = True
        ensemble = PDBEnsemble(title)
    else:
        ensemble = Ensemble(title)
    ensemble.setCoords(attr_dict['_coords'])
    if isPDBEnsemble:
        ensemble.addCoordset(attr_dict['_confs'], weights)
        if '_identifiers' in attr_dict.files:
            ensemble._labels = list(attr_dict['_identifiers'])
        if '_labels' in attr_dict.files:
            ensemble._labels = list(attr_dict['_labels'])
        if '_trans' in attr_dict.files:
            ensemble._trans = attr_dict['_trans']
    else:
        ensemble.addCoordset(attr_dict['_confs'])
        if weights != np.array(None): 
            ensemble.addCoordset(weights)
    return ensemble


def trimPDBEnsemble(pdb_ensemble, **kwargs):
    """Return a new PDB ensemble obtained by trimming given *pdb_ensemble*.
    This function helps selecting atoms in a pdb ensemble based on one of the 
    following criteria, and returns them in a new :class:`~.PDBEnsemble` 
    instance.
        
    **Occupancy**
    
    Resulting PDB ensemble will contain atoms whose occupancies are greater 
    or equal to *occupancy* keyword argument.  Occupancies for atoms will be
    calculated using ``calcOccupancies(pdb_ensemble, normed=True)``.
    
    :arg occupancy: occupancy for selecting atoms, must satisfy
        ``0 < occupancy <= 1``
    :type occupancy: float
    
    """
    
    if not isinstance(pdb_ensemble, PDBEnsemble):
        raise TypeError('pdb_ensemble argument must be a PDBEnsemble')
    if pdb_ensemble.numConfs() == 0 or pdb_ensemble.numAtoms() == 0:
        raise ValueError('pdb_ensemble must have conformations')
    
    if 'occupancy' in kwargs:
        occupancy = float(kwargs['occupancy'])
        assert 0 < occupancy <=1, ('occupancy is not > 0 and <= 1: '
                                   '{0:s}'.format(repr(occupancy)))
        n_confs = pdb_ensemble.numConfs()
        assert n_confs > 0, 'pdb_ensemble does not contain any conformations'
        occupancies = calcOccupancies(pdb_ensemble, normed=True)
        #assert weights is not None, 'weights must be set for pdb_ensemble'
        #weights = weights.flatten()
        #mean_weights = weights / n_confs
        torf = occupancies >= occupancy
    else:
        return None
    
    trimmed = PDBEnsemble(pdb_ensemble.getTitle())
    coords = pdb_ensemble.getCoords()
    if coords is not None:
        trimmed.setCoords( coords[torf] )
    confs = pdb_ensemble.getCoordsets()
    if confs is not None:
        weights = pdb_ensemble.getWeights()
        trimmed.addCoordset( confs[:, torf], weights[:, torf] )
    return trimmed


def calcOccupancies(pdb_ensemble, normed=False):
    """Return occupancy calculated from weights of a :class:`~.PDBEnsemble`.
    Any non-zero weight will be considered equal to one.  Occupancies are 
    calculated by binary weights for each atom over the conformations in 
    the ensemble. When *normed* is ``True``, total weights will be divided 
    by the number of atoms.
    
    When analyzing an ensemble of X-ray structures, this function can be used 
    to see how many times a residue is resolved.
    
    >>> from prody import *
    >>> ens = loadEnsemble('p38_X-ray.ens.npz')
    >>> print( calcOccupancies(ens) ) # doctest: +ELLIPSIS
    [ 74.  75.  75.  75.  75.  75.  75.  75.  75.  73.  73.  74.  75.  75.  75.
      ...
      75.  75.  75.  75.  75.  75.]
            
    Each number in the above example corresponds to a residue (or atoms) and 
    shows the number of structures in which the corresponding residue is 
    resolved."""
    
    if not isinstance(pdb_ensemble, PDBEnsemble):
        raise TypeError('pdb_ensemble must be a PDBEnsemble instance')
    if len(pdb_ensemble) == 0:
        raise ValueError('pdb_ensemble does not contain any conformations')
    assert isinstance(normed, bool), 'normed must be a boolean'
    weights = pdb_ensemble.getWeights()
    if weights is None:
        raise ValueError('pdb_ensemble weights are not set')

    occupancies = weights.astype(bool).sum(0).astype(float).flatten()
    if normed:
        return occupancies / len(pdb_ensemble)
    else:
        return occupancies
    
    
def showOccupancies(pdbensemble, *args, **kwargs):
    """Show occupancies for the PDB ensemble using :func:`~matplotlib.pyplot.
    plot`.  Occupancies are calculated using :meth:`calcOccupancies`."""
    
    import matplotlib.pyplot as plt

    if not isinstance(pdbensemble, PDBEnsemble):
        raise TypeError('pdbensemble must be a PDBEnsemble instance')
    weights = calcOccupancies(pdbensemble)
    if weights is None:
        return None
    show = plt.plot(weights, *args, **kwargs)
    axis = list(plt.axis())
    axis[2] = 0
    axis[3] += 1
    plt.axis(axis)
    plt.xlabel('Atom index')
    plt.ylabel('Sum of weights')
    if SETTINGS['auto_show']:
        plt.show(block=False)
    return show


def checkWeights(weights, n_atoms, n_csets=None):
    """Return weights if checks pass, otherwise raise an exception."""
    
    assert isinstance(n_atoms, int) and n_atoms > 0, \
        'n_atoms must be a positive integer'
    assert n_csets is None or isinstance(n_csets, int) and n_csets > 0, \
        'n_csets must be a positive integer'
    
    if not isinstance(weights, np.ndarray):
        raise TypeError('weights must be a Numpy array')
    elif n_csets is None and weights.ndim not in (1, 2):
        raise ValueError('weights.dim must be 1 or 2')
    elif n_csets is not None:
        if weights.ndim not in (1, 2, 3):
            raise ValueError('weights.dim must be 1, 2, or 3')
        elif weights.ndim == 3 and weights.shape[0] != n_csets:
            raise ValueError('weights.shape must be (n_csets, n_atoms, 1)')
    elif weights.ndim in (2, 3) and weights.shape[-1] != 1:
        raise ValueError('shape of weights must be ([n_csets,] n_atoms, 1)')
    elif weights.dtype != float:
        try:
            weights = weights.astype(float)
        except ValueError:
            raise ValueError('weights.astype(float) failed, float type could '
                             'not be assigned')
    if np.any(weights < 0):
        raise ValueError('all weights must greater or equal to 0')
    if weights.ndim == 1:
        weights = weights.reshape((n_atoms, 1))
    if n_csets is not None and weights.ndim == 2:
        weights = np.tile(weights.reshape((1, n_atoms, 1)), (n_csets, 1, 1))
    return weights


def alignPDBEnsemble(ensemble, suffix='_aligned', outdir='.', gzip=False):
    """Align PDB files using transformations from *ensemble*, which may be 
    a :class:`.PDBEnsemble` or a :class:`.PDBConformation` instance. Label of 
    the conformation (see :meth:`~.PDBConformation.getLabel`) will be used to 
    determine the PDB structure and model number.  First four characters of 
    the label is expected to be the PDB identifier and ending numbers to be the 
    model number.  For example, the :class:`.Transformation` from conformation 
    with label *2k39_ca_selection_'resnum_<_71'_m116* will be applied to 116th 
    model of structure **2k39**.  After applicable transformations are made, 
    structure will be written into *outputdir* as :file:`2k39_aligned.pdb`.  
    If *gzip* is **True**, output files will be compressed.  Return value is
    the output filename or list of filenames, in the order files are processed.
    Note that if multiple models from a file are aligned, that filename will
    appear in the list multiple times."""
    
    if not isinstance(ensemble, (PDBEnsemble, PDBConformation)):
        raise TypeError('ensemble must be a PDBEnsemble or PDBConformation')
    if isinstance(ensemble, PDBConformation):
        ensemble = [ensemble]
    if gzip:
        gzip = '.gz'
    else:
        gzip = ''
    output = []
    pdbdict = {}    
    for conf in ensemble:
        trans = conf.getTransformation()
        if trans is None:
            raise ValueError('transformations are not calculated, call '
                             '`superpose` or `iterpose`')
        label = conf.getLabel()

        pdb = label[:4]
        filename = pdbdict.get(pdb, fetchPDB(pdb))
        if filename is None:
            LOGGER.warning('PDB file for conformation {0:s} is not found.'
                           .format(label))
            output.append(None)
            continue
        LOGGER.info('Parsing PDB file {0:s} for conformation {1:s}.'
                    .format(pdb, label))

        acsi = None
        model = label.rfind('m')
        if model > 3:
            model = label[model+1:]
            if model.isdigit():
                acsi = int(model) - 1
            LOGGER.info('Applying transformation to model {0:s}.'
                        .format(model))

        if isinstance(filename, str):
            ag = parsePDB(filename)
        else:
            ag = filename

        if acsi is not None:
            if acsi >= ag.numCoordsets():
                LOGGER.warn('Model number {0:s} for {1:s} is out of range.'
                            .format(model, pdb))
                output.append(None)
                continue
            ag.setACSIndex(acsi)
        trans.apply(ag)
        outfn = os.path.join(outdir, pdb + suffix + '.pdb' + gzip)
        if ag.numCoordsets() > 1:
            pdbdict[pdb] = ag
        else:            
            writePDB(outfn, ag)
        output.append(outfn)
        
    for pdb, ag in pdbdict.iteritems():
        writePDB(os.path.join(outdir, pdb + suffix + '.pdb' + gzip), ag)
    if len(output) == 1:
        return output[0]
    else:
        return output
