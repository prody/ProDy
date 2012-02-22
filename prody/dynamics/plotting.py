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

"""This module defines plotting functions for protein dynamics analysis.

Plotting functions are called by the name of the plotted data/property 
and are prefixed with ``show``.  Function documentations refers to the 
:mod:`matplotlib.pyplot` function utilized for actual plotting. 
Arguments and keyword arguments are passed to the Matplotlib functions.


.. plot::
   :nofigs: 
   :context: 
    
   from prody import *
   import matplotlib.pyplot as plt
   import numpy as np

   p38_pca = loadModel('p38_xray.pca.npz')
   p38_anm = loadModel('1p38.anm.npz') 
   p38_ensemble = loadEnsemble('p38_X-ray.ens.npz')
   p38_structure = parsePDB('p38_ref_chain.pdb')   
   plt.close('all')""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody.ensemble import Ensemble, Conformation

from nma import NMA
from gnm import GNMBase
from mode import Mode, VectorBase, Vector
from modeset import ModeSet
from analysis import calcSqFlucts, calcProjection, calcCrossCorr
from compare import calcOverlap

__all__ = ['showContactMap', 'showCrossCorr',  
           'showCumOverlap', 'showFractOfVar',  
           'showCumFractOfVar', 'showMode', 
           'showOverlap', 'showOverlapTable', 'showProjection', 
           'showCrossProjection', 'showEllipsoid', 'showSqFlucts', 
           'showScaledSqFlucts', 'showNormedSqFlucts', 'resetTicks', ]

pkg = __import__(__package__)
LOGGER = pkg.LOGGER
           
def showEllipsoid(modes, onto=None, n_std=2, scale=1., *args, **kwargs):
    """Show an ellipsoid using  :meth:`~mpl_toolkits.mplot3d.Axes3D.plot_wireframe`.
    
    Ellipsoid volume gives an analytical view of the conformational space that
    given modes describe.
    
    :arg modes: 3 modes for which ellipsoid will be drawn.
    :type modes: :class:`~.ModeSet`, :class:`~.PCA`, :class:`~.ANM`, 
         :class:`~.NMA`
    
    :arg onto: 3 modes onto which ellipsoid will be projected.
    :type modes: :class:`~.ModeSet`, :class:`~.PCA`, :class:`~.ANM`,
        :class:`~.NMA`
       
    :arg n_std: Number of standard deviations to scale the ellipsoid.
    :type n_std: float
    
    :arg scale: Used for scaling the volume of ellipsoid. This can be
        obtained from :func:`~.sampleModes`.
    :type scale: float


    .. plot::
       :context:
       :include-source:
        
       # Show projection of subspace spanned by ANM 1-3 onto subspace of PC 1-3 
       plt.figure(figsize=(5,4))
       showEllipsoid(p38_anm[:3], p38_pca[:3])
       
       # Let's compare this with that of ANM modes 18-20
       showEllipsoid(p38_anm[17:], p38_pca[:3], color='red')
       # This ANM subspace appears as a tiny volume at the center
       # since faster ANM modes does not correspond to top ranking PCA modes
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be a NMA or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be from a 3-dimensional model')
    if len(modes) != 3:
        raise ValueError('length of modes is not equal to 3')
    if onto is not None:
        if not isinstance(onto, (NMA, ModeSet)):
            raise TypeError('onto must be a NMA or ModeSet instance, '
                            'not {0:s}'.format(type(onto)))
        if not onto.is3d():
            raise ValueError('onto must be from a 3-dimensional model')
        if len(onto) != 3:
            raise ValueError('length of onto is not equal to 3')
        if onto.numAtoms() != modes.numAtoms():
            raise ValueError('modes and onto must have same number of atoms')
        
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    
    
    var = modes.getVariances()
    #randn = np.random.standard_normal((1000, 3))
    #coef = ((randn ** 2 * var).sum(1) ** 0.5).mean()
    #scale=float(n_std)*modes.numAtoms()**.5 * float(rmsd) / coef * var **.5
    scale = float(n_std) * scale * var ** 0.5
    #scale=float(n_std)*modes.numAtoms()**.5*float(rmsd)/var.sum()**.5*var**.5   

    x = scale[0] * np.outer(np.cos(u), np.sin(v))
    y = scale[1] * np.outer(np.sin(u), np.sin(v))
    z = scale[2] * np.outer(np.ones(np.size(u)), np.cos(v))
    if onto is not None:
        change_of_basis = np.dot(modes._getArray().T, onto._getArray())

        xyz = np.array([x.flatten(), y.flatten(), z.flatten()])
        xyz = np.dot(xyz.T, change_of_basis)
        x = xyz[:,0].reshape((100,100))
        y = xyz[:,1].reshape((100,100))
        z = xyz[:,2].reshape((100,100))

    cf = plt.gcf()
    show = None
    for child in cf.get_children():
        if isinstance(child, Axes3D):
            show = child
            break 
    if show is None:
        show = Axes3D(cf)
    show.plot_wireframe(x, y, z, rstride=6, cstride=6, *args, **kwargs)
    if onto is not None:
        onto = onto.getModes()
        show.set_xlabel('Mode {0:d} coordinate'.format(onto[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(onto[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(onto[2].getIndex()+1))
    else:
        modes = modes.getModes()
        show.set_xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(modes[2].getIndex()+1))
    return show

def showFractOfVar(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`~matplotlib.pyplot.
    bar`.  Note that mode indices are incremented by 1.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showFractOfVar(p38_pca) 
       showCumFractOfVar(p38_pca)
      
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    if not isinstance(modes, (ModeSet, NMA)):
        raise TypeError('modes must be NMA, or ModeSet, not {0:s}'
                        .format(type(modes)))
    
    fracts = [(mode.getIndex(), mode.getFractOfVariance()) for mode in modes]
    fracts = np.array(fracts)
    show = plt.bar(fracts[:,0]+0.5, fracts[:,1], *args, **kwargs)
    axis = list(plt.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    plt.axis(axis)
    plt.xlabel('Mode index')
    plt.ylabel('Fraction of variance')
    return show

def showCumFractOfVar(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`~matplotlib.pyplot.
    plot`.
    
    Note that mode indices are incremented by 1.
    See :func:`~.showFractOfVar` for an example."""
    
    import matplotlib.pyplot as plt
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if isinstance(modes, Mode):
        indices = modes.getIndices() + 0.5
        modes = [modes]
    elif isinstance(modes, ModeSet):
        indices = modes.getIndices() + 0.5
    else:
        indices = np.arange(len(modes)) + 0.5
    fracts = np.array([mode.getFractOfVariance() for mode in modes]).cumsum()
    show = plt.plot(indices, fracts, *args, **kwargs)
    axis = list(plt.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    plt.axis(axis)
    plt.xlabel('Mode index')
    plt.ylabel('Fraction of variance')
    return show

def showProjection(ensemble, modes, *args, **kwargs):
    """Show a projection of conformational deviations onto up to three normal 
    modes from the same model.
    
    :arg ensemble: a :class:`~.Ensemble` instance
    :arg modes: a :class:`~.Mode`, :class:`~.ModeSet`, :class:`~.NMA`
    
    The projected values are by default converted to RMSD.  Pass 
    ``rmsd=False`` to use projection itself. :class:`~.Vector` instances 
    are accepted as *ensemble* argument to allow for projecting a 
    deformation vector onto normal modes.  
    
    Matplotlib function used for plotting depends on the number of modes:
        
      * 1 mode: :func:`~matplotlib.pyplot.hist`
      * 2 modes: :func:`~matplotlib.pyplot.plot`
      * 3 modes: :meth:`~mpl_toolkits.mplot3d.Axes3D.plot`
          
    By default ``marker='o', ls='None'`` is passed to the plotting function 
    to disable lines in projections onto 2 or 3-d spaces.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showProjection(p38_ensemble, p38_pca[0]) 
       plt.title('Projection onto PC1')

       plt.figure(figsize=(5,4))
       showProjection(p38_ensemble, p38_pca[:2])
       plt.title('Projection onto PC1-2')
       
       plt.figure(figsize=(5,4))
       showProjection(p38_ensemble, p38_pca[:3]) # onto top 3 PCs
       plt.title('Projection onto PC1-3')

       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    
    if not isinstance(ensemble, (Ensemble, Conformation, Vector)):
        raise TypeError('ensemble must be Ensemble, Conformation, or Vector, '
                        'not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    if not modes.is3d(): 
        raise Exception('modes must be 3-dimensional')
    
    if isinstance(modes, Mode) or (isinstance(modes, (ModeSet, NMA)) and 
                                   len(modes)==1):
        if not isinstance(modes, Mode):
            modes = modes[0]
        projection = calcProjection(ensemble, modes, kwargs.pop('rmsd', True))
        show = plt.hist(projection.flatten(), *args, **kwargs)
        plt.xlabel('Mode {0:d} coordinate'.format(modes.getIndex()+1))
        plt.ylabel('Number of conformations')
    elif len(modes) == 2:
        if 'ls' not in kwargs:
            kwargs['ls'] = 'None'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        projection = calcProjection(ensemble, modes, kwargs.pop('rmsd', True))
        show = plt.plot(projection[:, 0], projection[:, 1], *args, **kwargs)
        modes = [m for m in modes]
        plt.xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        plt.ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
    elif len(modes) == 3:
        if 'ls' not in kwargs:
            kwargs['ls'] = 'None'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        from mpl_toolkits.mplot3d import Axes3D
        projection = calcProjection(ensemble, modes, kwargs.pop('rmsd', True)) 
        modes = [m for m in modes]
        cf = plt.gcf()
        show = None
        for child in cf.get_children():
            if isinstance(child, Axes3D):
                show = child
                break 
        if show is None:
            show = Axes3D(cf)
        show.plot(projection[:, 0], projection[:, 1], projection[:, 2], 
                  *args, **kwargs)
        show.set_xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(modes[2].getIndex()+1))
    else:
        raise ValueError('Projection onto upto 3 modes can be shown. '
                         'You have given {0:d} mode.'.format(len(modes)))
    return show

def showCrossProjection(ensemble, mode_x, mode_y, scale=None, scalar=None, 
                        *args, **kwargs):
    """Show a projection of conformational deviations onto modes from
    different models using :func:`~matplotlib.pyplot.plot`.  This function 
    differs from :func:`~.showProjection` by accepting modes from two different 
    models.
    
    :arg ensemble: Ensemble for which deviations will be projected
    :type ensemble: :class:`~.Ensemble`
    :arg mode_x: Projection onto this mode will be shown along x-axis. 
    :type mode_x: :class:`~.Mode`
    :arg mode_y: Projection onto this mode will be shown along y-axis.
    :type mode_y: :class:`~.Mode`
    :arg scale: Scale width of the projection onto one of modes. 
                ``x`` and ``y`` are accepted.
    :type scale: str
    :arg scalar: Scalar factor for ``x`` or ``y``.  If ``scalar=None`` is 
        passed, best scaling factor will be calculated and printed on the
        console.
    :type scalar: float
    
    The projected values are by default converted to RMSD. 
    Pass ``rmsd=False`` to calculate raw projection values.
    :class:`~.Vector` instances are accepted as *ensemble* argument to allow
    for projecting a deformation vector onto normal modes.  
    
    By default ``marker='o', ls='None'`` is passed to the plotting function 
    to disable lines.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5.2,4))
       showCrossProjection(p38_ensemble, p38_pca[0], p38_anm[2])
    
    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    |example| See :ref:`pca-xray-plotting` for a more elaborate example."""

    import matplotlib.pyplot as plt
    if not isinstance(ensemble, (Ensemble, Conformation, Vector)):
        raise TypeError('ensemble must be Ensemble, Conformation, or Vector, '
                        'not {0:s}'.format(type(ensemble)))
    if not isinstance(mode_x, VectorBase):
        raise TypeError('mode_x must be a Mode instance, not {0:s}'
                        .format(type(mode_x)))
    if not mode_x.is3d():
        raise ValueError('mode_x must be 3-dimensional')
    if not isinstance(mode_y, VectorBase):
        raise TypeError('mode_y must be a Mode instance, not {0:s}'
                        .format(type(mode_y)))
    if not mode_y.is3d():
        raise ValueError('mode_y must be 3-dimensional')
    if scale is not None:
        assert isinstance(scale, str), 'scale must be a string'
        scale = scale.lower()
        assert scale in ('x', 'y'), 'scale must be x or y'
    if scalar is not None:
        assert isinstance(scalar, float), 'scalar must be a float'
    xcoords = calcProjection(ensemble, mode_x, kwargs.get('rmsd', True))
    ycoords = calcProjection(ensemble, mode_y, kwargs.pop('rmsd', True))
    if scale:
        if scalar is None:
            scalar = ((ycoords.max() - ycoords.min()) / 
                      (xcoords.max() - xcoords.min())) 
            scalar = scalar * np.sign(calcOverlap(mode_x, mode_y))
            if scale == 'x':
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode_x, scalar))
            else:
                scalar = 1 / scalar
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode_y, scalar))
        if scale == 'x':
            xcoords = xcoords * scalar  
        else:
            ycoords = ycoords * scalar
    if 'ls' not in kwargs:
        kwargs['ls'] = 'None'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'
    show = plt.plot(xcoords, ycoords, *args, **kwargs)
    plt.xlabel('{0:s} coordinate'.format(mode_x))
    plt.ylabel('{0:s} coordinate'.format(mode_y))
    return show

def showOverlapTable(rows, cols, *args, **kwargs):
    """Show overlap table using :func:`~matplotlib.pyplot.pcolor`.  *rows* and 
    *cols* are sets of normal modes, and correspond to rows and columns of the 
    displayed matrix.  Note that mode indices are increased by 1.  List of 
    modes should contain a set of contiguous modes from the same model. 
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showOverlapTable( p38_pca[:6], p38_anm[:6] )
       plt.title('p38 PCA vs ANM')

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    if not isinstance(rows, (NMA, ModeSet)):
        raise TypeError('rows must be an NMA model or a ModeSet, not {0:s}'
                        .format(type(rows)))
    if not isinstance(rows, (NMA, ModeSet)):
        raise TypeError('cols must be an NMA model or a ModeSet, not {0:s}'
                        .format(type(cols)))
    overlap = abs(calcOverlap(rows, cols))
    if isinstance(rows, NMA):
        rows = rows[:]
    if isinstance(cols, NMA):
        cols = cols[:]
    show = (plt.pcolor(overlap, cmap=plt.cm.jet, *args, **kwargs), 
            plt.colorbar())
    x_range = np.arange(1, len(cols)+1)
    plt.xticks(x_range-0.5, x_range)
    plt.xlabel(str(cols))
    y_range = np.arange(1, len(rows)+1)
    plt.yticks(y_range-0.5, y_range)
    plt.ylabel(str(rows))
    plt.axis([0, len(cols), 0, len(rows)])
    return show

def showCrossCorr(modes, *args, **kwargs):
    """Show cross-correlations for given modes using :func:`~matplotlib.pyplot.
    imshow`.  By default, *origin=lower* and *interpolation=bilinear* keyword 
    arguments are passed to imshow function. User can overwrite these 
    parameters.  See also :func:`~.calcCrossCorr`.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,5))
       # Show cross-correlations for ANM modes 1-3
       showCrossCorr( p38_anm[:3] )
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    arange = np.arange(modes.numAtoms())
    cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
    cross_correlations[arange[0]+1:, 
                       arange[0]+1:] = calcCrossCorr(modes)
    if not kwargs.has_key('interpolation'):
        kwargs['interpolation'] = 'bilinear'
    if not kwargs.has_key('origin'):
        kwargs['origin'] = 'lower'
    show = plt.imshow(cross_correlations, *args, **kwargs), plt.colorbar()
    plt.axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
    plt.title('Cross-correlations for {0:s}'.format(str(modes))) 
    plt.xlabel('Indices')
    plt.ylabel('Indices')
    return show

def showMode(mode, *args, **kwargs):
    """Show mode array using :func:`~matplotlib.pyplot.plot`.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,4))
       showMode( p38_anm[0] )
       plt.grid()
       plt.legend(loc='lower right', prop={'size': 10})
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, '
                        'not {0:s}'.format(type(mode)))
    if mode.is3d():
        a3d = mode.getArrayNx3()
        show = plt.plot(a3d[:, 0], *args, label='x-component', **kwargs)
        plt.plot(a3d[:, 1], *args, label='y-component', **kwargs)
        plt.plot(a3d[:, 2], *args, label='z-component', **kwargs)
    else:
        show = plt.plot(mode._getArray(), *args, **kwargs)
    plt.title(str(mode))
    plt.xlabel('Indices')
    return show

def showSqFlucts(modes, *args, **kwargs):
    """Show square fluctuations using :func:`~matplotlib.pyplot.plot`.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,4))
       showSqFlucts( p38_anm[0] )
       showSqFlucts( p38_anm[1] )
       plt.legend(prop={'size': 10})
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    sqf = calcSqFlucts(modes)
    if not 'label' in kwargs:
        kwargs['label'] = str(modes) 
    show = plt.plot(sqf, *args, **kwargs)
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    plt.title(str(modes))
    return show

def showScaledSqFlucts(modes, *args, **kwargs):
    """Show scaled square fluctuations using :func:`~matplotlib.pyplot.plot`.
    Modes or mode sets given as additional arguments will be scaled to have
    the same mean squared fluctuations as *modes*. 
    
    .. plot::
       :context:
       :include-source:
       
       plt.figure(figsize=(5,4))
       showScaledSqFlucts(p38_pca[0], p38_anm[2])
       plt.legend(prop={'size': 10})

    .. plot::
       :context:
       :nofigs:

       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    sqf = calcSqFlucts(modes)
    mean = sqf.mean()
    args = list(args)
    modesarg = []
    i = 0
    while i < len(args):
        if isinstance(args[i], (VectorBase, ModeSet, NMA)):
            modesarg.append(args.pop(i))
        else:
            i += 1
    show = [plt.plot(sqf, *args, label=str(modes), **kwargs)]
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        scalar = mean / sqf.mean()
        show.append(plt.plot(sqf * scalar, *args, 
                           label='{0:s} (x{1:.2f})'.format(str(modes), scalar), 
                           **kwargs))
    return show

def showNormedSqFlucts(modes, *args, **kwargs):
    """Show normalized square fluctuations using :func:`~matplotlib.pyplot.
    plot`.
    
    .. plot::
       :context:
       :include-source:
       
       plt.figure(figsize=(5,4))
       showNormedSqFlucts(p38_pca[0], p38_anm[2])
       plt.legend(prop={'size': 10})

    .. plot::
       :context:
       :nofigs:

       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    sqf = calcSqFlucts(modes)
    args = list(args)
    modesarg = []
    i = 0
    while i < len(args):
        if isinstance(args[i], (VectorBase, ModeSet, NMA)):
            modesarg.append(args.pop(i))
        else:
            i += 1
    show = [plt.plot(sqf/(sqf**2).sum()**0.5, *args, 
                        label='{0:s}'.format(str(modes)), **kwargs)]    
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        show.append(plt.plot(sqf/(sqf**2).sum()**0.5, *args, 
                    label='{0:s}'.format(str(modes)), **kwargs))
    return show

def showContactMap(enm, *args, **kwargs):
    """Show Kirchhoff matrix using :func:`~matplotlib.pyplot.spy`.
    
    .. plot::
       :context:
       :include-source:
        
       p38_gnm = GNM('p38')
       p38_gnm.buildKirchhoff( p38_structure )
       plt.figure(figsize=(4,4))
       showContactMap( p38_gnm )

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    if not isinstance(enm, GNMBase):
        raise TypeError('model argument must be an ENM instance')
    kirchhoff = enm.getKirchhoff()
    if kirchhoff is None:
        LOGGER.warning('kirchhoff matrix is not set')
        return None
    show = plt.spy(kirchhoff, *args, **kwargs)
    plt.title('{0:s} contact map'.format(enm.getTitle())) 
    plt.xlabel('Residue index')
    plt.ylabel('Residue index')
    return show

def showOverlap(mode, modes, *args, **kwargs):
    """Show overlap :func:`~matplotlib.pyplot.bar`.
    
    :arg mode: a single mode/vector
    :type mode: :class:`~.Mode`, :class:`~.Vector` 
    :arg modes: multiple modes
    :type modes: :class:`~.ModeSet`, :class:`~.ANM`, :class:`~.GNM`, 
        :class:`~.PCA` 
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(4,4))
       showOverlap( p38_pca[0], p38_anm[:6] )

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be Mode or Vector, not {0:s}'
                        .format(type(mode)))
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA or ModeSet, not {0:s}'
                        .format(type(modes)))
    overlap = abs(calcOverlap(mode, modes))
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+0.5)
    else:
        arange = modes.getIndices() + 0.5
    show = plt.bar(arange, overlap, *args, **kwargs)
    plt.title('Overlap with {0:s}'.format(str(mode)))
    plt.xlabel('{0:s} mode index'.format(modes))
    plt.ylabel('Overlap')
    return show

def showCumOverlap(mode, modes, *args, **kwargs):
    """Show cumulative overlap using :func:`~matplotlib.pyplot.plot`.
    
    :type mode: :class:`~.Mode`, :class:`~.Vector` 
    :arg modes: multiple modes
    :type modes: :class:`~.ModeSet`, :class:`~.ANM`, :class:`~.GNM`, 
        :class:`~.PCA` 
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showCumOverlap( p38_pca[0], p38_anm )
       # Let's also show the overlap
       showOverlap( p38_pca[0], p38_anm )

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')"""
    
    import matplotlib.pyplot as plt
    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be NMA, ModeSet, Mode or Vector, not {0:s}'
                        .format(type(mode)))
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    cumov = (calcOverlap(mode, modes) ** 2).cumsum() ** 0.5
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+0.5)
    else:
        arange = modes.getIndices() + 0.5
    show = plt.plot(arange, cumov, *args, **kwargs)
    plt.title('Cumulative overlap with {0:s}'.format(str(mode)))
    plt.xlabel('{0:s} mode index'.format(modes))
    plt.ylabel('Cumulative overlap')
    plt.axis((arange[0]-0.5, arange[-1]+0.5, 0, 1))
    return show
    
def resetTicks(x, y=None):
    """Reset X (and Y) axis ticks using values in given *array*.  Ticks in the 
    current figure should not be fractional values for this function to work as
    expected."""
    
    import matplotlib.pyplot as plt
    if x is not None:
        try:    
            xticks = plt.xticks()[0]
            xlist = list(xticks.astype(int))
            if xlist[-1] > len(x):
                xlist.pop()
            if xlist:
                xlist = list(x[xlist]) 
                plt.xticks(xticks, xlist + [''] * (len(xticks) - len(xlist)))
        except:
            LOGGER.warning('xticks could not be reset.')
    if y is not None:
        try:    
            yticks = plt.yticks()[0]
            ylist = list(yticks.astype(int))
            if ylist[-1] > len(y):
                ylist.pop()
            if ylist:
                ylist = list(y[ylist]) 
                plt.yticks(yticks, ylist + [''] * (len(yticks) - len(ylist)))
        except:
            LOGGER.warning('xticks could not be reset.')
    
