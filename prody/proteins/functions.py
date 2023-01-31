# -*- coding: utf-8 -*-
"""This module defines miscellaneous functions dealing with protein data."""

import numpy as np

from prody.atomic import Atomic, Atom, AtomGroup, Selection, HierView
from prody.utilities import openFile, showFigure, createStringIO, wrapModes
from prody import SETTINGS, PY3K

__all__ = ['view3D', 'showProtein']

def view3D(*alist, **kwargs):
    """Return a py3Dmol view instance for interactive visualization in
    Jupyter notebooks. Available arguments are: width, height (of
    the viewer), backgroundColor, zoomTo (a py3Dmol selection to center
    around), and style, which is a py3Dmol style object that will be
    applied to all atoms in the scene. More complex styling can be achieved
    by manipulating the view object directly.
    
    The default style is to show the protein in a rainbow cartoon and
    hetero atoms in sticks/spheres.
    
    GNM/ANM Coloring
    
    An array of fluctuation values can be provided with the *data* kwarg
    for visualization of GNM/ANM calculations.  The array is assumed to 
    correpond to a calpha selection of the provided protein.
    The default color will be set to a RWB color scheme on a per-residue
    basis.  If the fluctuation vector contains negative values, the
    midpoint (white) will be at zero.  Otherwise the midpoint is the mean.
    
    An array of displacement vectors can be provided with the *mode* kwarg.
    The animation of these motions can be controlled with frames (number
    of frames to animate over), amplitude (scaling factor), and animate
    (3Dmol.js animate options).  If animation isn't enabled, by default
    arrows are drawn.  Drawing of arrows is controlled by the boolean arrows
    option and the arrowcolor option.
    
    If multiple structures are provided with the data or mode arguments, these
    arguments must be provided as lists of arrays of the appropriate dimension.

    If a 3Dmol.js viewer as specified as the view argument, that viewer will be
    modified and returned.  After modification, update instead of show should
    be called on the viewer object if it is desired to update in-place
    instead of instantiating a new viewer.
    """

    try:
        import py3Dmol
    except:
        raise ImportError('py3Dmol needs to be installed to use view3D')
    
    from .pdbfile import writePDBStream
    
    width = kwargs.get('width', 400)
    height = kwargs.get('height', 400)
    data_list = kwargs.pop('data', None)
    modes = kwargs.pop('mode', None)
    style = kwargs.pop('style', [])
    zoomto = kwargs.pop('zoomto', {})
    bgcolor = kwargs.pop('backgroundcolor', 'white')
    bgcolor = kwargs.pop('backgroundColor', bgcolor)
    frames = kwargs.pop('frames', 30)
    interval = kwargs.pop('interval', 1)
    anim = kwargs.pop('anim', False)
    scale = kwargs.pop('scale', 100)
    arrows = kwargs.pop('arrows',True)
    arrowcolor = kwargs.pop('arrowcolor', 'darkgrey')
    arrowcolor = kwargs.pop('arrowColor', arrowcolor)

    if modes is None:
        n_modes = 0
    else:
        modes = wrapModes(modes)
        n_modes = len(modes)

    if data_list is None:
        n_data = 0
    else:
        data_list = wrapModes(data_list)
        n_data = len(data_list)

    view = kwargs.get('view',py3Dmol.view(width=width, height=height, js=kwargs.get('js','https://3dmol.csb.pitt.edu/build/3Dmol-min.js')))

    def _mapData(atoms, data):
        # construct map from residue to data property
        propmap = []
        for j, a in enumerate(atoms.calpha):
            propmap.append({'model': -1, 'chain': a.getChid(), 'resi': int(a.getResnum()), 
                            'props': {'data': float(data[j]) } })
        # set the atom property 
        # TODO: implement something more efficient on the 3Dmol.js side (this is O(n*m)!)
        view.mapAtomProperties(propmap)
        
        # color by property using gradient
        extreme = np.max(np.abs(data))
        lo = -extreme if np.min(data) < 0 else 0
        mid = np.mean(data) if np.min(data) >= 0 else 0
        view.setColorByProperty({'model': -1}, 'data', 'rwb', [extreme,lo,mid])
        view.setStyle({'model': -1},{'cartoon':{'style':'oval'}})    


    for i, atoms in enumerate(alist):
        pdb = createStringIO()
        writePDBStream(pdb, atoms)
        view.addAsOneMolecule(pdb.getvalue(), 'pdb')
        view.setStyle({'model': -1}, {'cartoon': {'color':'spectrum'}})
        view.setStyle({'model': -1, 'hetflag': True}, {'stick':{}})
        view.setStyle({'model': -1, 'bonds': 0}, {'sphere':{'radius': 0.5}})    

        if n_data:
            data = data_list[i]
            try:
                data = data.getArray()
            except AttributeError:
                pass

            if atoms.calpha.numAtoms() != len(data):
                raise RuntimeError("Atom count mismatch: {} vs {}. data styling assumes a calpha selection."
                                   .format(atoms.calpha.numAtoms(), len(data)))
            _mapData(atoms, data)
            
        if n_modes:
            mode = modes[i]
            try:
                arr = mode.getArray()
                is3d = mode.is3d()
            except AttributeError:
                arr = mode
                is3d = len(arr) == atoms.calpha.numAtoms()*3

            if is3d:
                if atoms.calpha.numAtoms()*3 != len(arr):
                    raise RuntimeError("Atom count mismatch: {} vs {}. mode animation assumes a calpha selection."
                                       .format(atoms.calpha.numAtoms(), len(arr)//3))

                if anim:
                    # construct map from residue to anm property and dx,dy,dz vectors
                    propmap = []
                    for j, a in enumerate(atoms.calpha):
                        propmap.append({'model': -1, 'chain': a.getChid(), 'resi': int(a.getResnum()),
                                        'props': {'dx': arr[3*j], 'dy': arr[3*j+1], 'dz': arr[3*j+2] } })
                    # set the atom property 
                    # TODO: implement something more efficient on the 3Dmol.js side (this is O(n*m)!)
                    view.mapAtomProperties(propmap)
                elif arrows:
                    for j, a in enumerate(atoms.calpha):
                        start = a._getCoords()
                        dcoords = arr[3*j:3*j+3]
                        end = start + dcoords * scale
                        view.addArrow({'start': {'x':start[0], 'y':start[1], 'z':start[2]},
                                       'end': {'x':end[0], 'y':end[1], 'z':end[2]},
                                       'radius': 0.3, 'color': arrowcolor})
            else:
                if atoms.calpha.numAtoms() != len(arr):
                    raise RuntimeError("Atom count mismatch: {} vs {}. mode styling assumes a calpha selection."
                                       .format(atoms.calpha.numAtoms(), len(arr)))

                _mapData(atoms, arr)

                
    # setting styles ...
    view.setBackgroundColor(bgcolor)

    if n_modes:
        #create vibrations
        view.vibrate(frames, scale)
        
        animate = kwargs.get('animate', {'loop':'rock', 'interval':interval})
        view.animate(animate)     
        
    if isinstance(style, dict):
        style = ({}, style)
    if isinstance(style, tuple):
        styles = [style]
    else:
        styles = style
    for sel, style in styles:
        view.setStyle(sel, style)
        
    view.zoomTo(zoomto)
            
    return view

def showProtein(*atoms, **kwargs):
    """Show protein representation using :meth:`~mpl_toolkits.mplot3d.Axes3D`.
    This function is designed for generating a quick view of the contents of a
    :class:`~.AtomGroup` or :class:`~.Selection`.    

    Protein atoms matching ``"calpha"`` selection are displayed using solid
    lines by picking a random and unique color per chain.  Line with can
    be adjusted using *lw* argument, e.g. ``lw=12``. Default width is 4.
    Chain colors can be overwritten using chain identifier as in ``A='green'``.

    Water molecule oxygen atoms are represented by red colored circles.  Color
    can be changed using *water* keyword argument, e.g. ``water='aqua'``.
    Water marker and size can be changed using *wmarker* and *wsize* keywords,
    defaults values are ``wmarker='.', wsize=6``.

    Hetero atoms matching ``"hetero and noh"`` selection are represented by
    circles and unique colors are picked at random on a per residue basis.
    Colors can be customized using residue name as in ``NAH='purple'``.  Note
    that this will color all distinct residues with the same name in the same
    color.  Hetero atom marker and size can be changed using *hmarker* and
    *hsize* keywords, default values are ``hmarker='o', hsize=6``.

    ProDy will set the size of axis so the representation is not distorted when
    the shape of figure window is close to a square.  Colors are picked at
    random, except for water oxygens which will always be colored red.
    
    *** Interactive 3D Rendering in Jupyter Notebook ***
    
    If py3Dmol has been imported then it will be used instead to display 
    an interactive viewer.  See :func:`view3D`
    
    
    """

    from prody.dynamics.mode import Mode

    method = kwargs.pop('draw', None)
    modes = kwargs.get('mode', None)
    scale = kwargs.get('scale', 100)

    # modes need to be specifically a list or a tuple (cannot be an array)
    if modes is None:
        n_modes = 0
    else:
        modes = wrapModes(modes)
        n_modes = len(modes)

    if method is None:
        import sys        
        if 'py3Dmol' in sys.modules: 
            method = 'py3Dmol'
        else:
            method = 'matplotlib'
    method = method.lower()
        
    alist = atoms
    for atoms in alist:
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance')
            
    if n_modes and n_modes != len(alist):
        raise RuntimeError('the number of proteins ({0}) does not match that of the modes ({1}).'
                            .format(len(alist), n_modes))

    if '3dmol' in method:
        mol = view3D(*alist, **kwargs)
        return mol
    else:
        kwargs.pop('mode', None)
        kwargs.pop('scale', 100)

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        cf = plt.gcf()
        show = None
        for child in cf.get_children():
            if isinstance(child, Axes3D):
                show = child
                break
        if show is None:
            show = Axes3D(cf)
        from matplotlib import colors
        cnames = dict(colors.cnames)
        wcolor = kwargs.get('water', 'red').lower()
        avoid = np.array(colors.hex2color(cnames.pop(wcolor, cnames.pop('red'))))
        for cn, val in cnames.copy().items():  # PY3K: OK
            clr = np.array(colors.hex2color(val))
            if clr.sum() > 2.4:
                cnames.pop(cn)
            elif np.abs(avoid - clr).sum() <= 0.6:
                cnames.pop(cn)
        cnames = list(cnames)
        import random
        random.shuffle(cnames)
        cnames_copy = list(cnames)
        min_ = list()
        max_ = list()
        for i, atoms in enumerate(alist):
            if isinstance(atoms, AtomGroup):
                title = atoms.getTitle()
            else:
                title = atoms.getAtomGroup().getTitle()
            calpha = atoms.select('calpha')
            if calpha:
                partition = False
                mode = modes[i] if n_modes else None
                if mode is not None:
                    is3d = False
                    try:
                        arr = mode.getArray()
                        is3d = mode.is3d()
                        n_nodes = mode.numAtoms()
                    except AttributeError:
                        arr = mode
                        is3d = len(arr) == len(calpha)*3
                        n_nodes = len(arr)//3 if is3d else len(arr)
                    if n_nodes != len(calpha):
                        raise RuntimeError('size mismatch between the protein ({0} residues) and the mode ({1} nodes).'
                                            .format(len(calpha), n_nodes))
                    partition = not is3d

                if partition:
                    xyz = calpha._getCoords()
                    chids = calpha.getChids()
                    rbody = []
                    last_sign = np.sign(arr[0])
                    last_chid = chids[0]
                    rcolor = ['red', 'red', 'blue']
                    n = 1
                    for i,a in enumerate(arr):
                        s = np.sign(a)
                        ch = chids[i]
                        if s == 0: s = last_sign
                        if last_sign != s or i == len(arr)-1 or last_chid != ch:
                            if last_chid == ch:
                                rbody.append(i)
                            show.plot(xyz[rbody, 0], xyz[rbody, 1], xyz[rbody, 2],
                                      label=title + '_regid%d'%n,
                                      color=rcolor[int(last_sign+1)],
                                      lw=kwargs.get('lw', 4))
                            rbody = []
                            n += 1
                            last_sign = s
                            last_chid = ch
                        rbody.append(i)
                else:
                    for ch in HierView(calpha, chain=True):
                        xyz = ch._getCoords()
                        chid = ch.getChid()
                        if len(cnames) == 0:
                            cnames = list(cnames_copy)
                        show.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                                label=title + '_' + chid,
                                color=kwargs.get(chid, cnames.pop()).lower(),
                                lw=kwargs.get('lw', 4))
                    
                    if mode is not None:
                        from prody.utilities.drawtools import drawArrow3D
                        XYZ = calpha._getCoords()
                        arr = arr.reshape((n_nodes, 3))
                        XYZ2 = XYZ + arr * scale
                        for i, xyz in enumerate(XYZ):
                            xyz2 = XYZ2[i]
                            mutation_scale = kwargs.pop('mutation_scale', 10)
                            drawArrow3D(xyz, xyz2, mutation_scale=mutation_scale, **kwargs)

            water = atoms.select('water and noh')
            if water:
                xyz = atoms.select('water')._getCoords()
                show.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], label=title + '_water',
                          color=wcolor,
                          ls='None', marker=kwargs.get('wmarker', '.'),
                          ms=kwargs.get('wsize', 6))
            hetero = atoms.select('not protein and not nucleic and not water and not dummy')
            if hetero:
                for res in HierView(hetero).iterResidues():
                    xyz = res._getCoords()
                    resname = res.getResname()
                    resnum = str(res.getResnum())
                    chid = res.getChid()
                    if len(cnames) == 0:
                        cnames = list(cnames_copy)
                    show.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], ls='None',
                              color=kwargs.get(resname, cnames.pop()).lower(),
                              label=title + '_' + chid + '_' + resname + resnum,
                              marker=kwargs.get('hmarker', 'o'),
                              ms=kwargs.get('hsize', 6))
            xyz = atoms._getCoords()
            min_.append(xyz.min(0))
            max_.append(xyz.max(0))

        show.set_xlabel('x')
        show.set_ylabel('y')
        show.set_zlabel('z')
        min_ = np.array(min_).min(0)
        max_ = np.array(max_).max(0)
        center = (max_ + min_) / 2
        half = (max_ - min_).max() / 2
        show.set_xlim3d(center[0]-half, center[0]+half)
        show.set_ylim3d(center[1]-half, center[1]+half)
        show.set_zlim3d(center[2]-half, center[2]+half)
        if kwargs.get('legend', False):
            show.legend(prop={'size': 10})
        if SETTINGS['auto_show']:
            showFigure()
        return show
