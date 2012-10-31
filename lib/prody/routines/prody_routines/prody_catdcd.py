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

"""Concatenate, slice, and/or select DCD files."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..actions import *

__all__ = ['prody_catdcd']

def prody_catdcd(*dcd, **kwargs):
    """Concatenate *dcd* files.
    
    :arg select: atom selection 
    
    :arg align: atom selection for aligning frames 
    
    :arg pdb: PDB file used in atom selections and as reference for alignment
        
    :arg psf: PSF file used in atom selections
    
    :arg output: output filename
    
    :arg first: index of the first output frame
    
    :arg last: index of the last output frame
    
    :arg stride: number of steps between output frames"""
    
    import prody
    LOGGER = prody.LOGGER
    if kwargs.get('numframes', False):
        for fn in dcd:
            print(prody.DCDFile(fn).numFrames())
        return

    from os.path import splitext
    
    ag = None
    psf = kwargs.get('psf', None)
    if psf:
        ag = prody.parsePSF(psf)
    
    pdb = kwargs.get('pdb', None)
    if pdb:
        if ag:
            ag = prody.parsePDB(pdb, ag=ag)
        else:
            ag = prody.parsePDB(pdb)

    align = kwargs.get('align', None)
    select = kwargs.get('select', None)
    if select == 'all':
        select = None
    if ag is None and (align or select):
        raise ValueError('one of PSF or PDB files must be provided for '
                         'align and select options to work')
    dcd = list(dcd)
    traj = prody.Trajectory(dcd.pop(0))
    while dcd:
        traj.addFile(dcd.pop(0))
    
    if ag:
        traj.link(ag)
        if ag.numCoordsets():
            traj.setCoords(ag.getCoords())
        if select:
            select = ag.select(select)
            if select is None:
                raise ValueError('{0:s} did not match any atoms'
                                 .format(repr(kwargs.get('select'))))
        else:
            select = ag
        LOGGER.info('{0:d} atoms are selected for writing output.'
                    .format(len(select)))
        if align:
            align = ag.select(align)
            traj.setAtoms(align)
            LOGGER.info('{0:d} atoms are selected for aligning frames.'
                        .format(len(align)))
    
    output = kwargs.get('output', 'trajectory.dcd')
    out = prody.DCDFile(output, 'w')
    count = 0
    stride = kwargs.get('stride', 1)
    goto = stride != 1
    slc = slice(kwargs.get('first', 0), kwargs.get('last', -1), 
                stride).indices(len(traj)+1)
    for i in range(*slc):
        if goto:
            traj.goto(i)
        frame = traj.next()
        if align:
            frame.superpose()
        if select:
            out.write(select._getCoords(), frame.getUnitcell())
        else:
            out.write(frame._getCoords(), frame.getUnitcell())
        count += 1
    traj.close()
    out.close()
    LOGGER.info("{0:d} frames are written into {1:s}."
                .format(count, output))

def addCommand(commands):
    
    subparser = commands.add_parser('catdcd', 
        help='concatenate dcd files')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Concatenate two DCD files and output all atmos: 
      
  $ prody catdcd mdm2.dcd mdm2sim2.dcd
  
Concatenate two DCD files and output backbone atoms: 
    
  $ prody catdcd mdm2.dcd mdm2sim2.dcd --pdb mdm2.pdb -s bb""")

    subparser.add_argument('-s', '--select', default='all', type=str, 
        dest='select', metavar='SEL', 
        help='atom selection (default: %(default)s)')

    subparser.add_argument('-o', '--output', type=str, metavar='FILE', 
        default='trajectory.dcd',
        help='output filename (default: %(default)s)')

    subparser.add_argument('-n', '--num', default=False, action='store_true',
        dest='numframes', 
        help='print the number of frames in each file and exit')

    subparser.add_argument('--psf', 
        help='PSF filename (must have same number of atoms as DCDs)')
    subparser.add_argument('--pdb', 
        help='PDB filename (must have same number of atoms as DCDs)')

    subparser.add_argument('--first', metavar='INT', type=int, default=0,
        help='index of the first output frame, default: %(default)s')
    subparser.add_argument('--last', metavar='INT', type=int, default=-1,
        help='index of the last output frame, default: %(default)s')
    subparser.add_argument('--stride', metavar='INT', type=int, default=1,
        help='number of steps between output frames, default: %(default)s')

    subparser.add_argument('--align', metavar='SEL', type=str,
        help='atom selection for aligning frames, a PSF or PDB file must be '
             'provided, if a PDB is provided frames will be superposed onto '
             'PDB coordinates')

    subparser.add_argument('dcd', nargs='+',
        help='DCD filename(s) (all must have same number of atoms)')

    subparser.set_defaults(subparser=subparser)
    subparser.set_defaults(func=lambda ns: prody_catdcd(
                                    *ns.__dict__.pop('dcd'), **ns.__dict__))
