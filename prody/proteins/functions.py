# -*- coding: utf-8 -*-
"""This module defines miscellaneous functions dealing with protein data."""

import numpy as np

from prody.atomic import Atomic, Atom, AtomGroup, Selection, HierView
from prody.utilities import openFile, showFigure
from prody import SETTINGS

__all__ = ['showProtein', 'writePQR', ]


def writePQR(filename, atoms):
    """Write *atoms* in PQR format to a file with name *filename*.  Only
    current coordinate set is written.  Returns *filename* upon success.  If
    *filename* ends with :file:`.gz`, a compressed file will be written."""

    if not isinstance(atoms, Atomic):
        raise TypeError('atoms does not have a valid type')
    if isinstance(atoms, Atom):
        atoms = Selection(atoms.getAtomGroup(), [atoms.getIndex()],
                          atoms.getACSIndex(),
                          'index ' + str(atoms.getIndex()))
    stream = openFile(filename, 'w')
    n_atoms = atoms.numAtoms()
    atomnames = atoms.getNames()
    if atomnames is None:
        raise RuntimeError('atom names are not set')
    for i, an in enumerate(atomnames):
        lenan = len(an)
        if lenan < 4:
            atomnames[i] = ' ' + an
        elif lenan > 4:
            atomnames[i] = an[:4]

    s_or_u = np.array(['a']).dtype.char

    resnames = atoms._getResnames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms
    resnums = atoms._getResnums()
    if resnums is None:
        resnums = np.ones(n_atoms, int)
    chainids = atoms._getChids()
    if chainids is None:
        chainids = np.zeros(n_atoms, s_or_u + '1')
    charges = atoms._getCharges()
    if charges is None:
        charges = np.zeros(n_atoms, float)
    radii = atoms._getRadii()
    if radii is None:
        radii = np.zeros(n_atoms, float)
    icodes = atoms._getIcodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, s_or_u + '1')
    hetero = ['ATOM'] * n_atoms
    heteroflags = atoms._getFlags('hetatm')
    if heteroflags is None:
        heteroflags = atoms._getFlags('hetero')
    if heteroflags is not None:
        hetero = np.array(hetero, s_or_u + '6')
        hetero[heteroflags] = 'HETATM'
    altlocs = atoms._getAltlocs()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, s_or_u + '1')

    format = ('{0:6s}{1:5d} {2:4s}{3:1s}' +
              '{4:4s}{5:1s}{6:4d}{7:1s}   ' +
              '{8:8.3f}{9:8.3f}{10:8.3f}' +
              '{11:8.4f}{12:7.4f}\n').format
    coords = atoms._getCoords()
    write = stream.write
    for i, xyz in enumerate(coords):
        write(format(hetero[i], i+1, atomnames[i], altlocs[i],
                     resnames[i], chainids[i], int(resnums[i]),
                     icodes[i], xyz[0], xyz[1], xyz[2], charges[i], radii[i]))
    write('TER\nEND')
    stream.close()
    return filename


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
    an interactive viewer.  Available arguments are: width, height (of
    the viewer), backgroundColor, zoomTo (a py3Dmol selection to center
    around), and styles, which should be a list of (selection, style) 
    tuple objects in py3Dmol format.  A single style may also be applied
    to all atoms with the 'style' keyword.
    
    The default style is to show the protein in a rainbow cartoon and
    hetero atoms in sticks/spheres.
    
    GNM/ANM Coloring
    
    An array of fluctuation values can be provided with the flucts kwarg
    for visualization of GNM/ANM calculations.  The array is assumed to 
    correpond to a calpha selection of the provided protein.
    The default color will be set to a RWB color scheme on a per-residue
    basis.  
    
    An array of displacement vectors can be provided with the vecs kwarg.
    The animation of these motions can be controlled with frames (number
    of frames to animate over), amplitude (scaling factor), and animate
    (3Dmol.js animate options).
    
    
    """

    alist = atoms
    for atoms in alist:
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance')
    
    import sys        
    if 'py3Dmol' in sys.modules:
        
        import StringIO, py3Dmol
        from pdbfile import writePDBStream

        
        pdb = StringIO.StringIO()
        
        for mol in alist:
            writePDBStream(pdb, mol)
        
        width = kwargs.get('width',400)
        height = kwargs.get('height',400)
        view = py3Dmol.view(width=width,height=height)
        
        #case insensitive kwargs..
        bgcolor = kwargs['backgroundcolor'] if 'backgroundcolor' in kwargs else kwargs.get('backgroundColor','white')
        view.setBackgroundColor(bgcolor)

        view.addModels(pdb.getvalue(),'pdb')
        view.setStyle({'cartoon': {'color':'spectrum'}})
        view.setStyle({'hetflag': True}, {'stick':{}})
        view.setStyle({'bonds': 0}, {'sphere':{'radius': 0.5}})    
    
        if 'flucts' in kwargs:
            garr = kwargs['flucts']
            #note we are only getting info from last set of atoms..
            if atoms.calpha.numAtoms() != len(garr):
                print "Atom count mismatch: %d vs %d.  GNM styling assume a calpha selection." % (atoms.calpha.numAtoms(), len(garr))
            else:
                #construct map from residue to flucts property
                propmap = []
                for (i,a) in enumerate(atoms.calpha):
                    propmap.append({'chain': a.getChid(), 'resi':a.getResnum(), 'props': {'flucts': garr[i] } })
                #set the atom property 
                #TODO: implement something more efficient on the 3Dmol.js side (this is O(n*m)!)
                view.mapAtomProperties(propmap)
                
                #color by property using gradient
                extreme = np.abs(garr).max()
                lo = -extreme if garr.min() < 0 else 0
                view.setColorByProperty({}, 'flucts', 'rwb', [extreme,lo])
                view.setStyle({'cartoon':{'style':'trace'}})
                
        if 'vecs' in kwargs:
            aarr = kwargs['vecs']  #has xyz coordinates

            #note we are only getting info from last set of atoms..
            if atoms.calpha.numAtoms()*3 != len(aarr):
                print "Atom count mismatch: %d vs %d.  GNM styling assume a calpha selection." % (atoms.calpha.numAtoms(), len(aarr)/3)
            else:
                #construct map from residue to anm property and dx,dy,dz vectors
                propmap = []
                for (i,a) in enumerate(atoms.calpha):
                    propmap.append({'chain': a.getChid(), 'resi':a.getResnum(),
                        'props': {'dy': aarr[3*i+1], 'dz': aarr[3*i+2] } });
                #set the atom property 
                #TODO: implement something more efficient on the 3Dmol.js side (this is O(n*m)!)
                view.mapAtomProperties(propmap)
                
                #create vibrations
                frames = kwargs.get('frames',10)
                amplitude = kwargs.get('amplitude',100)
                view.vibrate(frames, amplitude)
                
                animate = kwargs.get('animate',{'loop':'rock'})
                view.animate(animate)                
    
        if 'style' in kwargs: # this is never a list
            view.setStyle({},kwargs['style'])
            
        if 'styles' in kwargs:
            #allow simpler forms - convert them into a list
            styles = kwargs['styles']
            if type(styles) == dict:
                styles = ({},styles)
            if type(styles) == tuple:
                styles = [styles]
            for (sel, style) in styles:
                view.setStyle(sel, style)
        
        zoomto = kwargs['zoomto'] if 'zoomto' in kwargs else kwargs.get('zoomTo',{})
        view.zoomTo(zoomto)
                
        return view.show()

    else:
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
        for cn, val in cnames.items():  # PY3K: OK
            clr = np.array(colors.hex2color(val))
            if clr.sum() > 2.4:
                cnames.pop(cn)
            elif np.abs(avoid - clr).sum() <= 0.6:
                cnames.pop(cn)
        cnames = list(cnames)
        import random
        random.shuffle(cnames)
        min_ = list()
        max_ = list()
        for atoms in alist:
            if isinstance(atoms, AtomGroup):
                title = atoms.getTitle()
            else:
                title = atoms.getAtomGroup().getTitle()
            calpha = atoms.select('calpha')
            if calpha:
                for ch in HierView(calpha, chain=True):
                    xyz = ch._getCoords()
                    chid = ch.getChid()
                    show.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                              label=title + '_' + chid,
                              color=kwargs.get(chid, cnames.pop()).lower(),
                              lw=kwargs.get('lw', 4))
            water = atoms.select('water and noh')
            if water:
                xyz = atoms.select('water')._getCoords()
                show.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], label=title + '_water',
                          color=wcolor,
                          ls='None', marker=kwargs.get('wmarker', '.'),
                          ms=kwargs.get('wsize', 6))
            hetero = atoms.select('not protein and not nucleic and not water')
            if hetero:
                for res in HierView(hetero).iterResidues():
                    xyz = res._getCoords()
                    resname = res.getResname()
                    resnum = str(res.getResnum())
                    chid = res.getChid()
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
