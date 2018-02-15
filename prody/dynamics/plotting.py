# -*- coding: utf-8 -*-
"""This module defines plotting functions for protein dynamics analysis.

Plotting functions are called by the name of the plotted data/property
and are prefixed with ``show``.  Function documentations refers to the
:mod:`matplotlib.pyplot` function utilized for actual plotting.  Arguments
and keyword arguments are passed to the Matplotlib functions."""

from collections import defaultdict

import numpy as np

from prody import LOGGER, SETTINGS
from prody.utilities import showFigure

from .nma import NMA
from .gnm import GNMBase
from .mode import Mode, VectorBase, Vector
from .modeset import ModeSet
from .analysis import calcSqFlucts, calcProjection
from .analysis import calcCrossCorr, calcPairDeformationDist
from .analysis import calcFractVariance, calcCrossProjection 
from .perturb import calcPerturbResponse, calcPerturbResponseProfiles
from .compare import calcOverlap
from prody.atomic import AtomGroup, Selection

__all__ = ['showContactMap', 'showCrossCorr',
           'showCumulOverlap', 'showFractVars',
           'showCumulFractVars', 'showMode',
           'showOverlap', 'showOverlapTable', 'showProjection',
           'showCrossProjection', 'showEllipsoid', 'showSqFlucts',
           'showScaledSqFlucts', 'showNormedSqFlucts', 'resetTicks',
           'showDiffMatrix','showMechStiff','showNormDistFunct',
           'showPairDeformationDist','showMeanMechStiff', 
           'showPerturbResponse', 'showPerturbResponseProfiles',
           'showMatrix', 'showPlot', 'showAtomicData', 'showTree', 
           'showTree_networkx']


def showEllipsoid(modes, onto=None, n_std=2, scale=1., *args, **kwargs):
    """Show an ellipsoid using  :meth:`~mpl_toolkits.mplot3d.Axes3D
    .plot_wireframe`.

    Ellipsoid volume gives an analytical view of the conformational space that
    given modes describe.

    :arg modes: 3 modes for which ellipsoid will be drawn.
    :type modes: :class:`.ModeSet`, :class:`.PCA`, :class:`.ANM`, :class:`.NMA`

    :arg onto: 3 modes onto which ellipsoid will be projected.
    :type modes: :class:`.ModeSet`, :class:`.PCA`, :class:`.ANM`, :class:`.NMA`

    :arg n_std: Number of standard deviations to scale the ellipsoid.
    :type n_std: float

    :arg scale: Used for scaling the volume of ellipsoid. This can be
        obtained from :func:`.sampleModes`.
    :type scale: float"""

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be a NMA or ModeSet instance, '
                        'not {0}'.format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be from a 3-dimensional model')
    if len(modes) != 3:
        raise ValueError('length of modes is not equal to 3')
    if onto is not None:
        if not isinstance(onto, (NMA, ModeSet)):
            raise TypeError('onto must be a NMA or ModeSet instance, '
                            'not {0}'.format(type(onto)))
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
        onto = list(onto)
        show.set_xlabel('Mode {0} coordinate'.format(int(onto[0])+1))
        show.set_ylabel('Mode {0} coordinate'.format(int(onto[1])+1))
        show.set_zlabel('Mode {0} coordinate'.format(int(onto[2])+1))
    else:
        modes = list(modes)
        show.set_xlabel('Mode {0} coordinate'.format(int(modes[0])+1))
        show.set_ylabel('Mode {0} coordinate'.format(int(modes[1])+1))
        show.set_zlabel('Mode {0} coordinate'.format(int(modes[2])+1))
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showFractVars(modes, *args, **kwargs):
    """Show fraction of variances using :func:`~matplotlib.pyplot.bar`.  Note
    that mode indices are incremented by 1."""

    import matplotlib.pyplot as plt
    if not isinstance(modes, (ModeSet, NMA)):
        raise TypeError('modes must be NMA, or ModeSet, not {0}'
                        .format(type(modes)))
    
    if kwargs.pop('new_fig', True):
        plt.figure()

    fracts = calcFractVariance(modes)
    fracts = [(int(mode), fract) for mode, fract in zip(modes, fracts)]
    fracts = np.array(fracts)
    show = plt.bar(fracts[:,0]+0.5, fracts[:,1], *args, **kwargs)
    axis = list(plt.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    plt.axis(axis)
    plt.xlabel('Mode index')
    plt.ylabel('Fraction of variance')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showCumulFractVars(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`~matplotlib.pyplot.
    plot`.  Note that mode indices are incremented by 1.  See also
    :func:`.showFractVars` function."""

    import matplotlib.pyplot as plt
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0}'.format(type(modes)))
    if isinstance(modes, Mode):
        indices = modes.getIndices() + 0.5
        modes = [modes]
    elif isinstance(modes, ModeSet):
        indices = modes.getIndices() + 0.5
    else:
        indices = np.arange(len(modes)) + 0.5
    
    if kwargs.pop('new_fig', True):
        plt.figure()

    fracts = calcFractVariance(modes).cumsum()
    show = plt.plot(indices, fracts, *args, **kwargs)
    axis = list(plt.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    plt.axis(axis)
    plt.xlabel('Mode index')
    plt.ylabel('Fraction of variance')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showProjection(ensemble, modes, *args, **kwargs):
    """Show a projection of conformational deviations onto up to three normal
    modes from the same model.

    :arg ensemble: an ensemble, trajectory or a conformation for which
        deviation(s) will be projected, or a deformation vector
    :type ensemble: :class:`.Ensemble`, :class:`.Conformation`,
        :class:`.Vector`, :class:`.Trajectory`

    :arg modes: up to three normal modes
    :type modes: :class:`.Mode`, :class:`.ModeSet`, :class:`.NMA`

    :arg color: a color name or a list of color names or values, 
        default is ``'blue'``
    :type color: str, list

    :arg label: label or a list of labels
    :type label: str, list

    :arg marker: a marker or a list of markers, default is ``'o'``
    :type marker: str, list

    :arg linestyle: line style, default is ``'None'``
    :type linestyle: str

    :arg text: list of text labels, one for each conformation
    :type text: list

    :arg fontsize: font size for text labels
    :type fontsize: int

    :arg new_fig: if ``True`` then a new figure will be created before plotting.
        default is True
    :type new_fig: bool

    The projected values are by default converted to RMSD.  Pass ``rmsd=False``
    to use projection itself.

    Matplotlib function used for plotting depends on the number of modes:

      * 1 mode: :func:`~matplotlib.pyplot.hist`
      * 2 modes: :func:`~matplotlib.pyplot.scatter`
      * 3 modes: :meth:`~mpl_toolkits.mplot3d.Axes3D.scatter`"""

    import matplotlib.pyplot as plt

    cmap = kwargs.pop('cmap', None)

    if kwargs.pop('new_fig', True):
        fig, ax = plt.subplots() 
    projection = calcProjection(ensemble, modes, kwargs.pop('rmsd', True), kwargs.pop('norm', True))

    if projection.ndim == 1 or projection.shape[1] == 1:
        show = plt.hist(projection.flatten(), *args, **kwargs)
        plt.xlabel('{0} coordinate'.format(str(modes)))
        plt.ylabel('Number of conformations')
        return show
    elif projection.shape[1] > 3:
        raise ValueError('Projection onto up to 3 modes can be shown. '
                         'You have given {0} mode.'.format(len(modes)))

    num = projection.shape[0]

    markers = kwargs.pop('marker', 'o')
    if isinstance(markers, str) or markers is None:
        markers = [markers] * num
    elif isinstance(markers, list):
        if len(markers) != num:
            raise ValueError('length of marker must be {0}'.format(num))
    else:
        raise TypeError('marker must be a string or a list')

    colors = kwargs.pop('color', 'blue')
    if isinstance(colors, str) or colors is None:
        colors = [colors] * num
    elif isinstance(colors, list):
        if len(colors) != num:
            raise ValueError('length of color must be {0}'.format(num))
    else:
        raise TypeError('color must be a string or a list')

    labels = kwargs.pop('label', None)
    if isinstance(labels, str) or labels is None:
        labels = [labels] * num
    elif isinstance(labels, list):
        if len(labels) != num:
            raise ValueError('length of label must be {0}'.format(num))
    elif labels is not None:
        raise TypeError('label must be a string or a list')

    kwargs['linestyle'] = kwargs.pop('linestyle', None) or kwargs.pop('ls', 'None')

    texts = kwargs.pop('text', None)
    if texts:
        if not isinstance(texts, list):
            raise TypeError('text must be a list')
        elif len(texts) != num:
            raise TypeError('length of text must be {0}'.format(num))
        size = kwargs.pop('fontsize', None) or kwargs.pop('size', None)

    indict = defaultdict(list)
    for i, opts in enumerate(zip(markers, colors, labels)):  # PY3K: OK
        indict[opts].append(i)

    modes = [m for m in modes]
    if len(modes) == 2: 
        plot = plt.plot
        show = plt.gcf()
        text = plt.text
    else: 
        from mpl_toolkits.mplot3d import Axes3D
        cf = plt.gcf()
        show = None
        for child in cf.get_children():
            if isinstance(child, Axes3D):
                show = child
                break
        if show is None:
            show = Axes3D(cf)
        plot = show.plot
        text = show.text

    args = list(args)
    for opts, indices in indict.items():  # PY3K: OK
        marker, color, label = opts
        kwargs['marker'] = marker
        if cmap is not None:
            color = cmap.colors[int(float(color))]
        kwargs['c'] = color

        if label:
            kwargs['label'] = label
        else:
            kwargs.pop('label', None)

        plot(*(list(projection[indices].T) + args), **kwargs)

    if texts:
        kwargs = {}
        if size:
            kwargs['size'] = size
        for args in zip(*(list(projection.T) + [texts])):
            text(*args, **kwargs)

    if len(modes) == 2:
        plt.xlabel('{0} coordinate'.format(int(modes[0])+1))
        plt.ylabel('{0} coordinate'.format(int(modes[1])+1))
    elif len(modes) == 3:
        show.set_xlabel('Mode {0} coordinate'.format(int(modes[0])+1))
        show.set_ylabel('Mode {0} coordinate'.format(int(modes[1])+1))
        show.set_zlabel('Mode {0} coordinate'.format(int(modes[2])+1))

    if SETTINGS['auto_show']:
        showFigure()
    return show


def showCrossProjection(ensemble, mode_x, mode_y, scale=None, *args, **kwargs):
    """Show a projection of conformational deviations onto modes from
    different models using :func:`~matplotlib.pyplot.plot`.  This function
    differs from :func:`.showProjection` by accepting modes from two different
    models.

    :arg ensemble: an ensemble or a conformation for which deviation(s) will be
        projected, or a deformation vector
    :type ensemble: :class:`.Ensemble`, :class:`.Conformation`,
        :class:`.Vector`, :class:`.Trajectory`
    :arg mode_x: projection onto this mode will be shown along x-axis
    :type mode_x: :class:`.Mode`, :class:`.Vector`
    :arg mode_y: projection onto this mode will be shown along y-axis
    :type mode_y: :class:`.Mode`, :class:`.Vector`
    :arg scale: scale width of the projection onto mode ``x`` or ``y``,
        best scaling factor will be calculated and printed on the console,
        absolute value of scalar makes the with of two projection same,
        sign of scalar makes the projections yield a positive correlation
    :type scale: str
    :arg scalar: scalar factor for projection onto selected mode
    :type scalar: float
    :arg color: a color name or a list of color name, default is ``'blue'``
    :type color: str, list
    :arg label: label or a list of labels
    :type label: str, list
    :arg marker: a marker or a list of markers, default is ``'o'``
    :type marker: str, list
    :arg linestyle: line style, default is ``'None'``
    :type linestyle: str
    :arg text: list of text labels, one for each conformation
    :type text: list
    :arg fontsize: font size for text labels
    :type fontsize: int


    The projected values are by default converted to RMSD.  Pass ``rmsd=False``
    to calculate raw projection values.  See :ref:`pca-xray-plotting` for a
    more elaborate example."""

    import matplotlib.pyplot as plt

    if kwargs.pop('new_fig', True):
        plt.figure()

    norm = kwargs.pop('norm', True)
    xcoords, ycoords = calcCrossProjection(ensemble, mode_x, mode_y,
                                           scale=scale, norm=norm, **kwargs)

    num = len(xcoords)

    markers = kwargs.pop('marker', 'o')
    if isinstance(markers, str) or markers is None:
        markers = [markers] * num
    elif isinstance(markers, list):
        if len(markers) != num:
            raise ValueError('length of marker must be {0}'.format(num))
    else:
        raise TypeError('marker must be a string or a list')

    colors = kwargs.pop('color', 'blue')
    if isinstance(colors, str) or colors is None:
        colors = [colors] * num
    elif isinstance(colors, list):
        if len(colors) != num:
            raise ValueError('length of color must be {0}'.format(num))
    else:
        raise TypeError('color must be a string or a list')

    labels = kwargs.pop('label', None)
    if isinstance(labels, str) or labels is None:
        labels = [labels] * num
    elif isinstance(labels, list):
        if len(labels) != num:
            raise ValueError('length of label must be {0}'.format(num))
    elif labels is not None:
        raise TypeError('label must be a string or a list')

    kwargs['ls'] = kwargs.pop('linestyle', None) or kwargs.pop('ls', 'None')

    texts = kwargs.pop('text', None)
    if texts:
        if not isinstance(texts, list):
            raise TypeError('text must be a list')
        elif len(texts) != num:
            raise TypeError('length of text must be {0}'.format(num))
        size = kwargs.pop('fontsize', None) or kwargs.pop('size', None)

    indict = defaultdict(list)
    for i, opts in enumerate(zip(markers, colors, labels)):  # PY3K: OK
        indict[opts].append(i)

    for opts, indices in indict.items():  # PY3K: OK
        marker, color, label = opts
        kwargs['marker'] = marker
        kwargs['color'] = color
        if label:
            kwargs['label'] = label
        else:
            kwargs.pop('label', None)
        show = plt.plot(xcoords[indices], ycoords[indices], *args, **kwargs)
    if texts:
        kwargs = {}
        if size:
            kwargs['size'] = size
        for x, y, t in zip(xcoords, ycoords, texts):
            plt.text(x, y, t, **kwargs)
    plt.xlabel('{0} coordinate'.format(mode_x))
    plt.ylabel('{0} coordinate'.format(mode_y))
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showOverlapTable(modes_x, modes_y, **kwargs):
    """Show overlap table using :func:`~matplotlib.pyplot.pcolor`.  *modes_x*
    and *modes_y* are sets of normal modes, and correspond to x and y axes of
    the plot.  Note that mode indices are incremented by 1.  List of modes
    is assumed to contain a set of contiguous modes from the same model.

    Default arguments for :func:`~matplotlib.pyplot.pcolor`:

      * ``cmap=plt.cm.jet``
      * ``norm=matplotlib.colors.Normalize(0, 1)``"""

    import matplotlib.pyplot as plt
    import matplotlib

    overlap = abs(calcOverlap(modes_y, modes_x))
    if overlap.ndim == 0:
        overlap = np.array([[overlap]])
    elif overlap.ndim == 1:
        overlap = overlap.reshape((modes_y.numModes(), modes_x.numModes()))

    cmap = kwargs.pop('cmap', plt.cm.jet)
    norm = kwargs.pop('norm', matplotlib.colors.Normalize(0, 1))

    if kwargs.pop('new_fig', True):
        plt.figure()
    show = (plt.pcolor(overlap, cmap=cmap, norm=norm, **kwargs),
            plt.colorbar())
    x_range = np.arange(1, modes_x.numModes() + 1)
    plt.xticks(x_range-0.5, x_range)
    plt.xlabel(str(modes_x))
    y_range = np.arange(1, modes_y.numModes() + 1)
    plt.yticks(y_range-0.5, y_range)
    plt.ylabel(str(modes_y))
    plt.axis([0, modes_x.numModes(), 0, modes_y.numModes()])
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showCrossCorr(modes, *args, **kwargs):
    """Show cross-correlations using :func:`~matplotlib.pyplot.imshow`.  By
    default, *origin=lower* and *interpolation=bilinear* keyword  arguments
    are passed to this function, but user can overwrite these parameters.
    See also :func:`.calcCrossCorr`."""

    import matplotlib.pyplot as plt
    if kwargs.pop('new_fig', True):
        plt.figure()

    arange = np.arange(modes.numAtoms())
    cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
    cross_correlations[arange[0]+1:,
                       arange[0]+1:] = calcCrossCorr(modes)
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'
    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
    show = plt.imshow(cross_correlations, *args, **kwargs), plt.colorbar()
    plt.axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
    plt.title('Cross-correlations for {0}'.format(str(modes)))
    plt.xlabel('Indices')
    plt.ylabel('Indices')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showMode(mode, *args, **kwargs):
    """Show mode array using :func:`~matplotlib.pyplot.plot`."""
    
    import matplotlib.pyplot as plt

    show_hinges = kwargs.pop('show_hinges', False)
    show_zero = kwargs.pop('show_zero', True)
    overlay_chains = kwargs.get('overlay_chains',False)
    atoms = kwargs.get('atoms',None)

    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be a Mode or Vector instance, '
                        'not {0}'.format(type(mode)))
    if mode.is3d():
        a3d = mode.getArrayNx3()
        show = []
        show.append(showPlot(a3d[:, 0], *args, label='x-component', **kwargs))
        show.append(showPlot(a3d[:, 1], *args, label='y-component', **kwargs))
        show.append(showPlot(a3d[:, 2], *args, label='z-component', **kwargs))
        if atoms is not None:
            show[0][0].set_title(str(mode))
        else:
            show[0].set_title(str(mode))
    else:
        a1d = mode._getArray()
        show = showPlot(a1d, *args, **kwargs)
        if show_hinges and isinstance(mode, Mode):
            hinges = mode.getHinges()
            if hinges is not None:
                if atoms is not None:
                    if overlay_chains:
                        n = 0
                        chain_colors = 'gcmyrwbk'
                        for i in atoms.getHierView().iterChains():
                            for hinge in hinges:
                                if i.getResindices()[0] < hinge < i.getResindices()[-1]:
                                    show[0].plot(hinge - i.getResindices()[0], \
                                                 a1d[hinge], '*', color=chain_colors[n], \
                                                 markeredgecolor='k', markersize=10)
                            n += 1
                    else:
                        show[0].plot(hinges, a1d[hinges], 'r*')
                else:
                    show.plot(hinges, a1d[hinges], 'r*')
        if atoms is not None:
            show[0].set_title(str(mode))
        else:
            show.set_title(str(mode))
    if show_zero:
        if not mode.is3d():
            if atoms is not None:
                show[0].plot(show[0].get_xlim(), (0,0), '--', color='grey')
            else:
                show.plot(show.get_xlim(), (0,0), '--', color='grey')
        else:
            if atoms is not None:
                show[0][0].plot(show[0][0].get_xlim(), (0,0), '--', color='grey')
                show[1][0].plot(show[1][0].get_xlim(), (0,0), '--', color='grey')
                show[2][0].plot(show[2][0].get_xlim(), (0,0), '--', color='grey')
            else:
                show[0].plot(show[0].get_xlim(), (0,0), '--', color='grey')
                show[1].plot(show[1].get_xlim(), (0,0), '--', color='grey')
                show[2].plot(show[2].get_xlim(), (0,0), '--', color='grey')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showSqFlucts(modes, *args, **kwargs):
    """Show square fluctuations using :func:`~matplotlib.pyplot.plot`.  See
    also :func:`.calcSqFlucts`."""

    import matplotlib.pyplot as plt
    show_hinge = kwargs.pop('hinge', False)
    sqf = calcSqFlucts(modes)
    if not 'label' in kwargs:
        kwargs['label'] = str(modes)
    show = showPlot(sqf, *args, **kwargs)
    atoms = kwargs.get('atoms',None)
    if atoms is not None:
        show[0].set_ylabel('Square fluctuations')
        show[0].set_title(str(modes))
    else:
        show.set_ylabel('Square fluctuations')
        show.set_title(str(modes))
    if show_hinge and not modes.is3d():
        hinges = modes.getHinges()
        if hinges is not None:
            show[0].plot(hinges, sqf[hinges], 'r*')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showScaledSqFlucts(modes, *args, **kwargs):
    """Show scaled square fluctuations using :func:`~matplotlib.pyplot.plot`.
    Modes or mode sets given as additional arguments will be scaled to have
    the same mean squared fluctuations as *modes*."""

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
    show = showPlot(sqf, *args, label=str(modes), **kwargs)
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        scalar = mean / sqf.mean()
        show.append(plt.plot(sqf * scalar, *args,
                             label='{0} (x{1:.2f})'.format(str(modes), scalar),
                             **kwargs))
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showNormedSqFlucts(modes, *args, **kwargs):
    """Show normalized square fluctuations via :func:`~matplotlib.pyplot.plot`.
    """

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
    show = showPlot(sqf/(sqf**2).sum()**0.5, *args,
                     label='{0}'.format(str(modes)), **kwargs)
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        show.append(plt.plot(sqf/(sqf**2).sum()**0.5, *args,
                    label='{0}'.format(str(modes)), **kwargs))
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showContactMap(enm, *args, **kwargs):
    """Show Kirchhoff matrix using :func:`~matplotlib.pyplot.spy`."""

    import matplotlib.pyplot as plt
    if kwargs.pop('new_fig', True):
        plt.figure()
        
    if not isinstance(enm, GNMBase):
        raise TypeError('model argument must be an ENM instance')
    kirchhoff = enm.getKirchhoff()
    if kirchhoff is None:
        LOGGER.warning('kirchhoff matrix is not set')
        return None
    show = plt.spy(kirchhoff, *args, **kwargs)
    plt.title('{0} contact map'.format(enm.getTitle()))
    plt.xlabel('Residue index')
    plt.ylabel('Residue index')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showOverlap(mode, modes, *args, **kwargs):
    """Show overlap :func:`~matplotlib.pyplot.bar`.

    :arg mode: a single mode/vector
    :type mode: :class:`.Mode`, :class:`.Vector`
    :arg modes: multiple modes
    :type modes: :class:`.ModeSet`, :class:`.ANM`, :class:`.GNM`, :class:`.PCA`
    """

    import matplotlib.pyplot as plt

    if kwargs.pop('new_fig', True):
        plt.figure()

    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be Mode or Vector, not {0}'
                        .format(type(mode)))
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA or ModeSet, not {0}'
                        .format(type(modes)))
    overlap = abs(calcOverlap(mode, modes))
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+0.5)
    else:
        arange = modes.getIndices() + 0.5
    show = plt.bar(arange, overlap, *args, **kwargs)
    plt.title('Overlap with {0}'.format(str(mode)))
    plt.xlabel('{0} mode index'.format(modes))
    plt.ylabel('Overlap')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showCumulOverlap(mode, modes, *args, **kwargs):
    """Show cumulative overlap using :func:`~matplotlib.pyplot.plot`.

    :type mode: :class:`.Mode`, :class:`.Vector`
    :arg modes: multiple modes
    :type modes: :class:`.ModeSet`, :class:`.ANM`, :class:`.GNM`, :class:`.PCA`
    """

    import matplotlib.pyplot as plt
    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be NMA, ModeSet, Mode or Vector, not {0}'
                        .format(type(mode)))
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0}'
                        .format(type(modes)))
    cumov = (calcOverlap(mode, modes) ** 2).cumsum() ** 0.5
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+0.5)
    else:
        arange = modes.getIndices() + 0.5
    
    if kwargs.pop('new_fig', True):
        plt.figure()
    show = plt.plot(arange, cumov, *args, **kwargs)
    plt.title('Cumulative overlap with {0}'.format(str(mode)))
    plt.xlabel('{0} mode index'.format(modes))
    plt.ylabel('Cumulative overlap')
    plt.axis((arange[0]-0.5, arange[-1]+0.5, 0, 1))
    if SETTINGS['auto_show']:
        showFigure()
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


def showDiffMatrix(matrix1, matrix2, *args, **kwargs):
    """Show the difference between two cross-correlation matrices from
    different models. For given *matrix1* and *matrix2* show the difference
    between them in the form of (matrix2 - matrix1) and plot the difference
    matrix using :func:`showMatrix`. When :class:`.NMA` models
    are passed instead of matrices, the functions could call
    :func:`.calcCrossCorr` function to calculate the matrices for given modes.

    To display the absolute values in the difference matrix, user could set
    *abs* keyword argument **True**.

    By default, *origin=lower* and *interpolation=bilinear* keyword arguments
    are passed to this function, but user can overwrite these parameters.
    """

    import matplotlib.pyplot as plt
    try:
        dim1, shape1 = matrix1.ndim, matrix1.shape
    except AttributeError:
        matrix1 = calcCrossCorr(matrix1)
        dim1, shape1 = matrix1.ndim, matrix1.shape
    try:
        dim2, shape2 = matrix2.ndim, matrix2.shape
    except AttributeError:
        matrix2 = calcCrossCorr(matrix2)
        dim2, shape2 = matrix2.ndim, matrix2.shape
    if (not ((dim1 == dim2 == 2) and (shape1 == shape2))):
        raise ValueError('Matrices must have same square shape.')
    if shape1[0] * shape1[1] == 0:
        raise ValueError('There are no data in matrices.')
    diff = matrix2 - matrix1
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'
    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
    if kwargs.pop('abs', False):
        diff = np.abs(diff)
    if kwargs.pop('new_fig', True):
        plt.figure()
    show = showMatrix(diff, *args, **kwargs)
    show.im3.axis([-.5, shape1[1] - .5, -.5, shape1[0] - .5])
    plt.title('Difference Matrix')
    if SETTINGS['auto_show']:
        showFigure()
    plt.xlabel('Indices')
    plt.ylabel('Indices')
    return show


def showMechStiff(model, coords, *args, **kwargs):
    """Show mechanical stiffness matrix using :func:`~matplotlib.pyplot.imshow`.
    By default, *origin=lower* keyword  arguments are passed to this function, 
    but user can overwrite these parameters."""

    import math
    import matplotlib
    import matplotlib.pyplot as plt
    arange = np.arange(model.numAtoms())
    model.buildMechStiff(coords)

    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
    if 'jet_r' in kwargs:
        import matplotlib.cm as plt
        kwargs['jet_r'] = 'cmap=cm.jet_r'
        
    MechStiff = model.getStiffness()
    matplotlib.rcParams['font.size'] = '14'

    if kwargs.pop('new_fig', True):
        fig = plt.figure(num=None, figsize=(10,8), dpi=100, facecolor='w')
    vmin = math.floor(np.min(MechStiff[np.nonzero(MechStiff)]))
    vmax = round(np.amax(MechStiff),1)
    show = showMatrix(MechStiff, vmin=vmin, vmax=vmax, *args, **kwargs)
    plt.title('Mechanical Stiffness Matrix')# for {0}'.format(str(model)))
    plt.xlabel('Indices', fontsize='16')
    plt.ylabel('Indices', fontsize='16')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showNormDistFunct(model, coords, *args, **kwargs):
    """Show normalized distance fluctuation matrix using 
    :func:`~matplotlib.pyplot.imshow`. By default, *origin=lower* 
    keyword  arguments are passed to this function, 
    but user can overwrite these parameters."""

    import math
    import matplotlib
    import matplotlib.pyplot as plt
    normdistfunct = model.getNormDistFluct(coords)

    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
        
    matplotlib.rcParams['font.size'] = '14'

    if kwargs.pop('new_fig', True):
        fig = plt.figure(num=None, figsize=(10,8), dpi=100, facecolor='w')
    vmin = math.floor(np.min(normdistfunct[np.nonzero(normdistfunct)]))
    vmax = round(np.amax(normdistfunct),1)
    show = showMatrix(normdistfunct, vmin=vmin, vmax=vmax, *args, **kwargs)
    #plt.clim(math.floor(np.min(normdistfunct[np.nonzero(normdistfunct)])), \
    #                                       round(np.amax(normdistfunct),1))
    plt.title('Normalized Distance Fluctution Matrix')
    plt.xlabel('Indices', fontsize='16')
    plt.ylabel('Indices', fontsize='16')
    if SETTINGS['auto_show']:
        showFigure()
    return show


def showPairDeformationDist(model, coords, ind1, ind2, *args, **kwargs):
    """Show distribution of deformations in distance contributed by each mode
    for selected pair of residues *ind1* *ind2*
    using :func:`~matplotlib.pyplot.plot`. """

    import matplotlib
    import matplotlib.pyplot as plt
    if not isinstance(model, NMA):
        raise TypeError('model must be a NMA instance, '
                        'not {0}'.format(type(model)))
    elif not model.is3d():
        raise TypeError('model must be a 3-dimensional NMA instance')
    elif len(model) == 0:
        raise ValueError('model must have normal modes calculated')
    elif model.getStiffness() is None:
        raise ValueError('model must have stiffness matrix calculated')

    d_pair = calcPairDeformationDist(model, coords, ind1, ind2)
    with plt.style.context('fivethirtyeight'):
        matplotlib.rcParams['font.size'] = '16'
        fig = plt.figure(num=None, figsize=(12,8), dpi=100, facecolor='w')
        #plt.title(str(model))
        plt.plot(d_pair[0], d_pair[1], 'k-', linewidth=1.5, *args, **kwargs)
        plt.xlabel('mode (k)', fontsize = '18')
        plt.ylabel('d$^k$' '($\AA$)', fontsize = '18')
    if SETTINGS['auto_show']:
        showFigure()
    return plt.show


def showMeanMechStiff(model, coords, header, chain='A', *args, **kwargs):
    """Show mean value of effective spring constant with secondary structure
    taken from MechStiff. Header is needed to obatin secondary structure range.
    Using ``'jet_r'`` as argument color map will be reverse (similar to VMD 
    program coding).
    """
    meanStiff = np.array([np.mean(model.getStiffness(), axis=0)])
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    fig=plt.figure(figsize=[18,6], facecolor='w', dpi=100)
    
    if 'jet_r' in kwargs:
       import matplotlib.cm as plt
       kwargs['jet_r'] = 'cmap=cm.jet_r'
    if 'nearest' in kwargs:
        kwargs['nearest'] = 'interpolation=nearest'

    with plt.style.context('fivethirtyeight'):
        ax = fig.add_subplot(111)
        matplotlib.rcParams['font.size'] = '24'
        plt.plot(np.arange(len(meanStiff[0]))+coords.getResnums()[0],meanStiff[0], 'k-', linewidth = 3)
        plt.xlim(coords.getResnums()[0], coords.getResnums()[-1])
        ax_top=round(np.max(meanStiff[0])+((np.max(meanStiff[0])-np.min(meanStiff[0]))/3))
        ax_bottom=np.floor(np.min(meanStiff[0]))
        LOGGER.info('The range of mean effective force constant is: {0} to {1}.'
                                           .format(min(meanStiff[0]), max(meanStiff[0])))
        plt.ylim(ax_bottom,ax_top)
        plt.xlabel('residue', fontsize = '22')
        plt.ylabel('mean $\kappa$ [a.u.]', fontsize = '22')

    ax = fig.add_subplot(411, aspect='equal')
    plt.imshow(meanStiff, *args, **kwargs)
    header_ss = header['sheet_range'] + header['helix_range']
    for i in range(len(header_ss)):
        if header_ss[i][1] == chain:
            beg = int(header_ss[i][-2])-coords.getResnums()[0]
            end = int(header_ss[i][-1])-coords.getResnums()[0]
            add_beg = end - beg
            if header_ss[i][0] == 'H':
                ax.add_patch(patches.Rectangle((beg+1,-0.7),add_beg,\
                1.4,fill=False, linestyle='solid',edgecolor='#b22683', linewidth=2))    
            elif header_ss[i][0] == 'E':
                if header_ss[i][2] == -1:    
                    ax.add_patch(patches.Arrow(beg+1,0,add_beg,0,width=4.65, \
                    fill=False, linestyle='solid',edgecolor='black', linewidth=2))
                else: 
                    ax.add_patch(patches.Arrow(end+1,0,add_beg*(-1),0,width=4.65, \
                    fill=False, linestyle='solid',edgecolor='black', linewidth=2))
    plt.axis('off')
    ax.set_ylim(-1.7,1.7)
    if SETTINGS['auto_show']:
        showFigure()
    return plt.show

def showPerturbResponse(**kwargs):
    """ Plot the PRS matrix with the profiles along the right and bottom.

    If no PRS matrix or profiles are provided, these will be calculated first
    using the provided options with a provided model (e.g. ANM, GNM or EDA).
    So as to obtain different sensors and effectors, normMatrix=True by default.

    If atoms are provided then residue numbers can be used from there.
    *model* and *atoms* must have the same number of atoms. *atoms* must be an
    :class:`.AtomGroup` instance.

    :arg prs_matrix: a perturbation response matrix
    :type prs_matrix: array

    :arg effectiveness: an effectiveness profile from a PRS matrix
    :type effectiveness: array

    :arg sensitivity: a sensitivity profile from a PRS matrix
    :type sensitivity: array

    :arg model: any object with a calcCovariance method
        e.g. :class:`.ANM` instance
    :type model: NMA

    :arg atoms: a :class: `AtomGroup` instance
    :type atoms: AtomGroup

    :arg returnData: whether to return data for further analysis
        default is False
    :type returnData: bool
	
    :arg percentile: percentile argument for showMatrix
    :type percentile: float

    Return values are prs_matrix, effectiveness, sensitivity, ax1, ax2, im, ax3, ax4
    The PRS matrix, effectiveness and sensitivity will not be returned if provided. 
    If returnData is False then only the last five objects are returned.
    """

    import matplotlib.pyplot as plt
    import matplotlib

    prs_matrix = kwargs.get('prs_matrix')
    effectiveness = kwargs.get('effectiveness')
    sensitivity = kwargs.get('sensitivity')
    model = kwargs.pop('model')
    atoms = kwargs.get('atoms')
    returnData = kwargs.pop('returnData',False)

    if atoms is None:

        if prs_matrix is None:
            if model is None:
                raise ValueError('Please provide a PRS matrix or model.')
            else:
                prs_matrix = calcPerturbResponse(model=model)

        if effectiveness is None or sensitivity is None:
            effectiveness, sensitivity = calcPerturbResponseProfiles(prs_matrix)
        else:
            returnData = False

        showMatrix_returns = showMatrix(prs_matrix, effectiveness, sensitivity, **kwargs)

    else:
        if not isinstance(atoms, AtomGroup) and not isinstance(atoms, Selection):
            raise TypeError('atoms must be an AtomGroup instance')
        elif model is not None and atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')

        if prs_matrix is None: 
            if model is None:
                raise ValueError('Please provide a PRS matrix or model.')
            atoms, prs_matrix = calcPerturbResponse(model=model,atoms=atoms)

        if effectiveness is None or sensitivity is None:
            atoms, effectiveness, sensitivity = calcPerturbResponseProfiles(prs_matrix,atoms)

        showMatrix_returns = showMatrix(prs_matrix, effectiveness, sensitivity, **kwargs)

    if not returnData:
        return showMatrix_returns
    elif kwargs.get('prs_matrix') is not None:
       if atoms is not None:
           return atoms, effectiveness, sensitivity, showMatrix_returns
       else:
           return effectiveness, sensitivity, showMatrix_returns
    else:
       if atoms is not None:
           return atoms, prs_matrix, effectiveness, sensitivity, showMatrix_returns
       else:
           return prs_matrix, effectiveness, sensitivity, showMatrix_returns

def showPerturbResponseProfiles(prs_matrix,atoms=None,**kwargs):
    """Plot as a line graph the average response to perturbation of
    a particular residue (a row of a perturbation response matrix)
    or the average effect of perturbation of a particular residue
    (a column of a normalized perturbation response matrix).

    If no PRS matrix or profiles are provided, these will be calculated first
    using the provided options with a provided model (e.g. ANM, GNM or EDA).
    So as to obtain different sensitivity and effectiveness, normMatrix=True by default.

    If no residue number is given then the effectiveness and sensitivity
    profiles will be plotted instead. These two profiles are also returned
    as arrays for further analysis if they aren't already provided.

    :arg prs_matrix: a perturbation response matrix
    :type prs_matrix: ndarray

    :arg atoms: a :class: `AtomGroup` instance for matching 
        residue numbers and chain IDs. 
    :type atoms: AtomGroup

    :arg effectiveness: an effectiveness profile from a PRS matrix
    :type effectiveness: list

    :arg sensitivity: a sensitivity profile from a PRS matrix
    :type sensitivity: list

    :arg model: any object with a calcCovariance method
        e.g. :class:`.ANM` instance
        *model* and *atoms* must have the same number of atoms.
    :type model: NMA

    :arg chain: chain identifier for the residue of interest
        default is to make a plot for each chain in the protein
    :type chain: str

    :arg resnum: residue number for the residue of interest
    :type resnum: int

    :arg direction: the direction you want to use to read data out
        of the PRS matrix for plotting: the options are 'effect' or 'response'.
        Default is 'effect'.
        A row gives the effect on each residue of peturbing the specified 
        residue.
        A column gives the response of the specified residue to perturbing 
        each residue.
        If no residue number is provided then this option will be ignored
    :type direction: str

    :arg returnData: whether to return profiles for further analysis
        default is False
    :type returnProfiles: bool
    """
    model = kwargs.get('model')
    if not type(prs_matrix) is np.ndarray:
        if prs_matrix is None:
            if model is None:
                raise ValueError('Please provide a PRS matrix or model.')
            else:
                if kwargs.get('normMatrix') is None:
                    kwargs.set('normMatrix',True)
                prs_matrix = calcPerturbResponse(**kwargs)
        else:
            raise TypeError('Please provide a valid PRS matrix (as array).')

    if atoms is None:
        raise ValueError('Please provide an AtomGroup object for matching ' \
                         'residue numbers and chain IDs.')
    else:
        if not isinstance(atoms, AtomGroup) and not isinstance(atoms, Selection):
            raise TypeError('atoms must be an AtomGroup instance')
        elif model is not None and atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')

    chain = kwargs.get('chain')
    hv = atoms.getHierView()
    chains = []
    for i in range(len(list(hv))):
        chainAg = list(hv)[i]
        chains.append(chainAg.getChids()[0])

    chains = np.array(chains)
    if chain is None:
        chain = ''.join(chains)

    resnum = kwargs.get('resnum', None)
    direction = kwargs.get('direction','effect')

    if resnum is not None: 
        timesNotFound = 0
        for n in range(len(chain)):
            if not chain[n] in chains:
                raise PRSMatrixParseError('Chain {0} was not found in {1}'.format(chain[n], pdbIn))

            chainNum = int(np.where(chains == chain[n])[0])
            chainAg = list(hv)[chainNum]
            if not resnum in chainAg.getResnums():
                LOGGER.info('A residue with number {0} was not found' \
                            ' in chain {1}. Continuing to next chain.' \
                            .format(resnum, chain[n]))
                timesNotFound += 1
                continue

        profiles = []
        for n in range(len(chain)):
            chainNum = int(np.where(chains == chain[n])[0])
            i = np.where(atoms.getResnums() == resnum)[0][chainNum-timesNotFound] 
            if direction is 'effect':
                profiles.append(prs_matrix[i,:])
            else:
                profiles.append(prs_matrix[:,i])

    else:
        effectiveness = kwargs.get('effectiveness')
        sensitivity = kwargs.get('sensitivity')
        if effectiveness is None or sensitivity is None:
            effectiveness, sensitivity = calcPerturbResponseProfiles(prs_matrix)
        profiles = [effectiveness, sensitivity]

    for profile in profiles:
        show = showPlot(profile,atoms=atoms,**kwargs)

    returnData = kwargs.get('returnData',False)
    if returnData:
        return show, profiles
    else:
        return show

def showMatrix(matrix=None, x_array=None, y_array=None, **kwargs):
    """Show a matrix using :meth:`~matplotlib.axes.Axes.imshow`. Curves on x- and y-axis can be added.
    The first return value is the :class:`~matplotlib.axes.Axes` object for the upper plot, and the second
    return value is equivalent object for the left plot. The third return value is 
    the :class:`~matplotlib.image.AxesImage` object for the matrix plot. The last return value is the 
    :class:`~matplotlib.axes.Axes` object for the color bar.

    :arg matrix: Matrix to be displayed.
    :type matrix: :class:`~numpy.ndarray`

    :arg x_array: Data to be plotted above the matrix.
    :type x_array: :class:`~numpy.ndarray`

    :arg y_array: Data to be plotted on the left side of the matrix.
    :type y_array: :class:`~numpy.ndarray`

    :arg percentile: A percentile threshold to remove outliers, i.e. only showing data within *p*-th 
                     to *100-p*-th percentile.
    :type percentile: float

    :arg vmin: Minimum value that can be used together with vmax 
               as an alternative way to remove outliers
    :type vmin: float

    :arg vmax: Maximum value that can be used together with vmin 
               as alternative way to remove outliers
    :type vmax: float

    :arg atoms: a :class: `AtomGroup` instance for matching 
        residue numbers and chain IDs. 
    :type atoms: :class: `AtomGroup`

    :arg num_div: the number of divisions for each chain
        default 2
    :type num_div: int

    :arg resnum_tick_labels: residue number labels in place of num_div.
         A list can be used to set the same labels on all chains or 
         a dictionary of lists to set different labels for each chain
    :type resnum_tick_labels: list or dictionary

    :arg add_last_resi: whether to add a label for the last residue
        default False
    :type add_last_resi: bool

    :arg label_size: size for resnum labels
        default is 6, which works well for 4 residues on 4 chains
    :type label_size: int
    """ 

    num_div = kwargs.pop('num_div',2)
    resnum_tick_labels = kwargs.pop('resnum_tick_labels',None)
    add_last_resi = kwargs.pop('add_last_resi',False)

    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.collections import LineCollection
    from matplotlib.pyplot import figure, imshow

    if matrix is None:
        raise TypeError('You need to provide a matrix.') 
    elif len(np.shape(matrix)) != 2:
        raise ValueError('The matrix must be a 2D array.')

    p = kwargs.pop('percentile', None)
    vmin = kwargs.pop('vmin', None)
    vmax = kwargs.pop('vmax', None)

    if vmin is None and vmax is None and p is not None:
        vmin = np.percentile(matrix, p)
        vmax = np.percentile(matrix, 100-p)
    
    W = 10
    H = 10
    aspect = 'auto'

    if x_array is not None and y_array is not None:
        nrow = 3; ncol = 5
        i = 1; j = 1
        width_ratios = [1, W, 0.2]
        height_ratios = [1, H, 0.2]
    elif x_array is not None and y_array is None:
        nrow = 3; ncol = 4
        i = 1; j = 0
        width_ratios = [W, 0.2]
        height_ratios = [1, H, 0.2]
    elif x_array is None and y_array is not None:
        nrow = 2; ncol = 5
        i = 0; j = 1
        width_ratios = [1, W, 0.2]
        height_ratios = [H, 0.2]
    else:
        nrow = 2; ncol = 4
        i = 0; j = 0
        width_ratios = [W, 0.2]
        height_ratios = [H, 0.2]

    main_index = (i,j)
    upper_index = (i-1,j)
    lower_index = (i+1,j)
    left_index = (i,j-1)
    right_index = (i,j+1)

    outer = GridSpec(1, 3, width_ratios = [sum(width_ratios), 1, 4], hspace=0., wspace=0.2) 

    gs = GridSpecFromSubplotSpec(nrow, ncol-2, subplot_spec = outer[0], width_ratios=width_ratios,
                                 height_ratios=height_ratios, hspace=0., wspace=0.)

    gs_bar = GridSpecFromSubplotSpec(nrow-1, 1, subplot_spec = outer[1], height_ratios=height_ratios[:-1], hspace=0., wspace=0.)
    gs_legend = GridSpecFromSubplotSpec(nrow-1, 1, subplot_spec = outer[2], height_ratios=height_ratios[:-1], hspace=0., wspace=0.)

    new_fig = kwargs.pop('new_fig', True)
    if new_fig:
        fig = plt.figure(figsize=[9.5,6]) 
    axes = []

    atoms = kwargs.pop('atoms', None)
    if atoms is not None:
        if not isinstance(atoms, AtomGroup) and not isinstance(atoms, Selection):
            raise TypeError('atoms must be an AtomGroup instance')
    
    cmap = kwargs.pop('cmap', 'jet')
    label_size = kwargs.pop('label_size', 6)
 
    ax1 = ax2 = ax3 = ax4 = ax5 = ax6 = ax7 = None
    if nrow > 2:
        y1 = x_array
        x1 = np.arange(len(y1))
        ax1 = plt.subplot(gs[upper_index])
        points1 = np.array([x1, y1]).T.reshape(-1, 1, 2)
        segments1 = np.concatenate([points1[:-1], points1[1:]], axis=1)
        lc1 = LineCollection(segments1, array=y1, linewidths=1, cmap=cmap)
        ax1.add_collection(lc1)

        ax1.set_xlim(x1.min(), x1.max())
        ax1.set_ylim(y1.min(), y1.max())
        ax1.axis('off')

    if ncol > 4:
        x2 = y_array
        y2 = np.arange(len(x2))
        ax2 = plt.subplot(gs[left_index])
        points2 = np.array([x2, y2]).T.reshape(-1, 1, 2)
        segments2 = np.concatenate([points2[:-1], points2[1:]], axis=1)
        lc2 = LineCollection(segments2, array=x2, linewidths=1, cmap=cmap)
        ax2.add_collection(lc2)

        ax2.set_xlim(x2.min(), x2.max())
        ax2.set_ylim(y2.min(), y2.max())
        ax2.axis('off')
        ax2.invert_xaxis()

    ax3 = plt.subplot(gs[main_index])
    im = imshow(matrix, aspect=aspect, vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)
    ax3.yaxis.tick_right()

    ax4 = plt.subplot(gs_bar[-1])
    plt.colorbar(cax=ax4)
 
    if atoms is not None:

        # Add bars along the bottom and right that are colored by chain and numbered with residue

        ax5 = plt.subplot(gs[lower_index])
        ax6 = plt.subplot(gs[right_index])

        n = 0
        resnum_tick_locs = []
        resnum_tick_labels_list = []

        if resnum_tick_labels is None:
            resnum_tick_labels = []
            user_set_labels = False
        elif type(resnum_tick_labels) is list:
            user_set_labels = list
        elif type(resnum_tick_labels) is dict:
            user_set_labels = dict
        else:
            raise TypeError('The resnum tick labels should be a list or dictionary of lists')

        chain_colors = 'gcmyrwbk'
        chain_handles = []
        for i in atoms.getHierView().iterChains():
            
            chain_handle, = ax5.plot([i.getResindices()[0], i.getResindices()[-1]], [0, 0], \
                                     '-', linewidth=3, color=chain_colors[n], label=str(i))
            chain_handles.append(chain_handle)

            ax6.plot([0,0], [np.flip(i.getResindices(),0)[0], np.flip(i.getResindices(),0)[-1]], \
                     '-', linewidth=3, color=chain_colors[n], label=str(i))

            if not user_set_labels:
                for j in range(num_div):
                    resnum_tick_locs.append(i.getResindices()[i.numAtoms()/num_div*j])
                    resnum_tick_labels.append(i.getResnums()[i.numAtoms()/num_div*j])
            elif user_set_labels is list:
                for j in resnum_tick_labels:
                    resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                    resnum_tick_labels_list.append(j)
            else:
                for k in resnum_tick_labels.keys():
                    if i.getChids()[0] == k:
                       for j in resnum_tick_labels[k]:
                           resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                           resnum_tick_labels_list.append(j)

            n += 1

        ax7 = plt.subplot(gs_legend[-1])
        plt.legend(handles=chain_handles, loc=2, bbox_to_anchor=(0.25, 1))
        ax7.axis('off')

        if add_last_resi:
            resnum_tick_locs.append(atoms.getResindices()[-1])
            resnum_tick_labels_list.append(atoms.getResnums()[-1])

        resnum_tick_locs = np.array(resnum_tick_locs)
        resnum_tick_labels = np.array(resnum_tick_labels_list)

        ax3.axis('off')

        ax5.set_xticks(resnum_tick_locs)
        ax5.set_xticklabels(resnum_tick_labels)
        ax5.tick_params(labelsize=label_size)
        ax5.set_yticks([])

        ax6.set_xticks([])
        ax6.yaxis.tick_right()
        ax6.set_yticks(resnum_tick_locs)
        ax6.set_yticklabels(resnum_tick_labels)
        ax6.tick_params(labelsize=label_size)

        ax5.set_xlim([-0.5, len(matrix)+0.5])
        ax6.set_ylim([-0.5, len(matrix.T)+0.5])

    if SETTINGS['auto_show']:
        showFigure()
 
    return ax1, ax2, im, ax3, ax4, ax5, ax6, ax7

def showAtomicData(y, atoms=None, **kwargs):
    """
    Show a plot with the option to include chain color bars using provided atoms.
    
    :arg atoms: a :class: `AtomGroup` instance for matching 
        residue numbers and chain IDs. 
    :type atoms: :class: `AtomGroup`

    :arg add_last_resi: whether to add a label for the last residue
        default False
    :type add_last_resi: bool

    :arg label_size: size for resnum labels
        default is 6, which works well for 4 residues on 4 chains
    :type label_size: int

    :arg overlay_chains: overlay the chains rather than having them one after another
        default False
    :type overlay_chains: bool

    :arg domain_bar: color the bar at the bottom by domains rather than chains
        default False
    :type domain_bar: bool
    """

    overlay_chains = kwargs.pop('overlay_chains', False)
    show_domains = kwargs.pop('show_domains', None)
    domain_bar = kwargs.pop('domain_bar', False)

    from prody.utilities import showData
    from matplotlib.pyplot import figure, imshow, ylim
    from matplotlib import ticker

    new_fig = kwargs.pop('new_fig', True)
    if new_fig:
        figure()
    
    ticklabels = None
    if atoms is not None:
        hv = atoms.getHierView()
        if hv.numChains() == 0:
            raise ValueError('atoms should contain at least one chain.')
        elif hv.numChains() == 1:
            if show_domains is None:
                show_domains = False
            ticklabels = atoms.getResnums()
        else:
            ticklabels = []
            chids = atoms.getChids()
            resnums = atoms.getResnums()
            ticklabels = ['%s:%d'%(c, n) for c, n in zip(chids, resnums)]
            
    ax = showData(y, ticklabels=ticklabels)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    if show_domains is None:
        show_domains = atoms is not None

    if show_domains and atoms is not None:
        yl = ylim()
        d_loc = yl[0]
        D = []
        chids = atoms.getChids()
        uni_chids = np.unique(chids)
        for chid in uni_chids:
            d = chids == chid
            D.append(d)
        D = np.vstack(D).T
        F = np.zeros(D.shape)
        F[~D] = np.nan
        F[D] = d_loc

        ax.plot(F, linewidth=5)
        ylim(yl)

    if SETTINGS['auto_show']:
        showFigure()
    return ax

def showPlot(y, **kwargs):

    """
    Show a plot with the option to include chain color bars using provided atoms.
    
    :arg atoms: a :class: `AtomGroup` instance for matching 
        residue numbers and chain IDs. 
    :type atoms: :class: `AtomGroup`
    
    :arg num_div: the number of divisions for each chain
        default 2
    :type num_div: int

    :arg resnum_tick_labels: residue number labels in place of num_div.
         A list can be used to set the same labels on all chains or 
         a dictionary of lists to set different labels for each chain
    :type resnum_tick_labels: list or dictionary

    :arg add_last_resi: whether to add a label for the last residue
        default False
    :type add_last_resi: bool

    :arg label_size: size for resnum labels
        default is 6, which works well for 4 residues on 4 chains
    :type label_size: int

    :arg overlay_chains: overlay the chains rather than having them one after another
        default False
    :type overlay_chains: bool

    :arg domain_bar: color the bar at the bottom by domains rather than chains
        default False
    :type domain_bar: bool
    """
    atoms = kwargs.pop('atoms',None)
    overlay_chains = kwargs.pop('overlay_chains',False)
    domain_bar = kwargs.pop('domain_bar',False)

    num_div = kwargs.pop('num_div',2)
    resnum_tick_labels = kwargs.pop('resnum_tick_labels',None)
    add_last_resi = kwargs.pop('add_last_resi',False)
    label_size = kwargs.pop('label_size',6)

    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.collections import LineCollection
    from matplotlib.pyplot import figure, imshow

    if y is None:
        raise TypeError('You need to provide data for the y-axis.')
    elif len(np.shape(y)) != 1:
        raise ValueError('The data must be a 1D array.')

    new_fig = kwargs.pop('new_fig', True)
    if new_fig:
        fig = plt.figure(figsize=[9.5,6])
    axes = [] 

    if atoms is not None:
        height_ratios = [15,0.2]
        nrows = 2
    else:
        height_ratios = None
        nrows = 1

    outer = GridSpec(1, 2, width_ratios = [16, 4], hspace=0., wspace=0.)

    gs = GridSpecFromSubplotSpec(nrows, 1, subplot_spec = outer[0], \
                                 height_ratios=height_ratios, hspace=0., wspace=0.)

    gs_legend = GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[1], hspace=0., wspace=0.)
    
    ax1 = plt.subplot(gs[0])

    chain_colors = 'gcmyrwbk'
    chain_handles = []

    if overlay_chains:
        n = 0
        for i in atoms.getHierView().iterChains():
            chain_handle, = ax1.plot(y[i.getResindices()[0]:i.getResindices()[-1]], color=chain_colors[n], label=str(i), **kwargs)
            chain_handles.append(chain_handle)
            n += 1
    else:
        ax1.plot(y, **kwargs)

    if nrows > 1:
        ax2 = plt.subplot(gs[1])

        resnum_tick_locs = []
        resnum_tick_labels_list = []

        if resnum_tick_labels is None:
            resnum_tick_labels = []
            user_set_labels = False
        elif type(resnum_tick_labels) is list:
            user_set_labels = list
        elif type(resnum_tick_labels) is dict:
            user_set_labels = dict
        else:
            raise TypeError('The resnum tick labels should be a list or dictionary of lists')

        n = 0
        for i in atoms.getHierView().iterChains():
            if not overlay_chains:
                chain_handle, = ax2.plot([i.getResindices()[0], i.getResindices()[-1]], [0, 0], \
                                         '-', linewidth=3, color=chain_colors[n], label=str(i))
                chain_handles.append(chain_handle)

            if not user_set_labels:
                for j in range(num_div):
                    resnum_tick_locs.append(i.getResindices()[i.numAtoms()/num_div*j])
                    resnum_tick_labels.append(i.getResnums()[i.numAtoms()/num_div*j])
            elif user_set_labels is list:
                for j in resnum_tick_labels:
                    resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                    resnum_tick_labels_list.append(j)
            else:
                for k in resnum_tick_labels.keys():
                    if i.getChids()[0] == k:
                       for j in resnum_tick_labels[k]: 
                           resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                           resnum_tick_labels_list.append(j)

            n += 1

        if domain_bar:
            try:
                atoms.getData('domain')[0]
            except:
                raise ValueError('A domain bar can only be generated if \
                                  there is domain data associated with \
                                  the atoms.')

            borders = {}
            for i in range(atoms.numAtoms()/atoms.getHierView().numChains()):
                if atoms.getData('domain')[i] != atoms.getData('domain')[i-1]:
                    if i != 0:
                        borders[atoms.getData('domain')[i-1]][-1].append(i-1)
                    if not atoms.getData('domain')[i] in borders.keys():
                        borders[atoms.getData('domain')[i]] = []
                    borders[atoms.getData('domain')[i]].append([])
                    borders[atoms.getData('domain')[i]][-1].append(i)

            hsv = plt.get_cmap('hsv')
            colors = hsv(np.linspace(0, 1.0, len(borders.keys())))

            for chain in atoms.getHierView().iterChains():
                domains_found = []
                for i in range(chain.numAtoms()):
                    if not atoms.getData('domain')[i] in domains_found and str(atoms.getData('domain')[i]) is not '':
                        n = 0
                        for j in borders[atoms.getData('domain')[i]]:
                            m = 0
                            if m == 0:
                                domain_handle, = ax2.plot([j[0], j[-1]], [0, 0], '-', linewidth=3, \
                                                          color=colors[n], label=str(atoms.getData('domain')[i]))
                                chain_handles.append(domain_handle)
                            else:
                                ax2.plot([j[0], j[-1]], [0, 0], '-', linewidth=3, color=colors[n])
                            m += 1
                        n += 1
 
        ax3 = plt.subplot(gs_legend[-1])
        plt.legend(handles=chain_handles, loc=2, bbox_to_anchor=(0.25, 1))
        ax3.axis('off')

        if not user_set_labels:
            resnum_tick_labels_list = resnum_tick_labels

        if add_last_resi:
            resnum_tick_locs.append(atoms.getResindices()[-1])
            resnum_tick_labels_list.append(atoms.getResnums()[-1])

        resnum_tick_locs = np.array(resnum_tick_locs)
        resnum_tick_labels = np.array(resnum_tick_labels_list)

        ax1.set_xticks([])

        if overlay_chains:
            ax1.set_xlim(-0.5,atoms.numAtoms()/atoms.getHierView().numChains()+0.5)

        ax2.set_xticks(resnum_tick_locs)
        ax2.set_xticklabels(resnum_tick_labels)
        ax2.tick_params(labelsize=label_size)
        ax2.set_yticks([])

        ax2.set_xlim(ax1.get_xlim())

    if atoms is not None:
        return ax1, ax2, ax3
    else:
        return ax1

def showTree(tree, format='ascii', **kwargs):
    """ Given a tree, creates visualization in different formats. 
    arg tree: Tree needs to be unrooted and should be generated by tree generator from Phylo in biopython. 
    type tree: Bio.Phylo.BaseTree.Tree
    arg format: Depending on the format, you will see different forms of trees. Acceptable formats are `plt`
    and `ascii`.
    type format: str
    arg font_size: Font size for branch labels
    type: float
    arg line_width: The line width for each branch
    type: float
    """
    try: 
        from Bio import Phylo
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')
    font_size = float(kwargs.get('font_size', 8.0))
    line_width = float(kwargs.get('line_width', 1.5))
    if format == 'ascii':
        obj = Phylo.draw_ascii(tree)
    elif format == 'pylab' or format == 'matplotlib': 
        try:
            import pylab
        except:
            raise ImportError("Pylab or matplotlib is not installed.")
        pylab.rcParams["font.size"]=font_size
        pylab.rcParams["lines.linewidth"]=line_width
        obj = Phylo.draw(tree, do_show=False)
        pylab.xlabel('distance')
        pylab.ylabel('proteins')
    elif format == 'networkx':
        node_size = kwargs.pop('node_size', 20)
        node_color = kwargs.pop('node_color', 'red')
        withlabels = kwargs.pop('withlabels', True)
        scale = kwargs.pop('scale', 1.)
        iterations = kwargs.pop('iterations', 500)
        obj = showTree_networkx(tree, node_size, node_color, 
                                withlabels, scale, iterations, **kwargs)

    return obj

def showTree_networkx(tree, node_size=20, node_color='red', withlabels=True, scale=1., iterations=500, **kwargs):
    from Bio import Phylo
    import matplotlib.pyplot as mpl

    try:
        import networkx
    except ImportError:
        raise ImportError('Please install networkx to use this function.')

    G = Phylo.to_networkx(tree)
    labels = {}
    colors = []
    sizes = []

    for node in G.nodes():
        lbl = node.name
        if lbl is None:
            lbl = ''
            colors.append('black')
            sizes.append(0)
        else:
            sizes.append(node_size)
            if isinstance(node_color, str):
                nc = node_color
            else:
                nc = node_color[lbl] if lbl in node_color else 'red'
            colors.append(nc)
        labels[node] = lbl

    if kwargs.pop('new_fig', True):
        mpl.figure()

    layout = networkx.spring_layout(G, scale=scale, iterations=iterations)
    networkx.draw(G, pos=layout, withlabels=False, node_size=sizes, node_color=colors)

    if withlabels:
        fontdict = kwargs.pop('fontdict', None)
        if fontdict is None:
            fontsize = kwargs.pop('font_size', 6)
            fontcolor = kwargs.pop('font_color', 'black')
            fontdict = {'size': fontsize, 'color': fontcolor}

        for node, pos in layout.iteritems():
            mpl.text(pos[0], pos[1], labels[node], fontdict=fontdict)

    if SETTINGS['auto_show']:
        showFigure()

    return mpl.gca()
