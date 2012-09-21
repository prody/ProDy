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

import os.path

from actions import *

def prody_catdcd(opt):
    """Concatenate DCD files."""
    
    import prody
    LOGGER = prody.LOGGER
    if opt.num:
        num = [] 
        for dcd in opt.dcd:
            dcd = prody.DCDFile(dcd)
            num.append(dcd.numFrames())
        for n in num:
            print(n)
        print(sum(num))
        return
    align = opt.align
    ag = opt.psf or opt.pdb
    if ag:
        if os.path.splitext(ag)[1].lower() == '.psf':
            ag = prody.parsePSF(ag)
        else:
            ag = prody.parsePDB(ag)
    elif align:
        raise ValueError('one of PSF or PDB files must be provided for '
                         'align option to work')
    
    dcd = opt.dcd
    traj = prody.Trajectory(dcd.pop(0))
    while dcd:
        traj.addFile(dcd.pop(0))
    if ag:
        traj.setAtoms(ag)
        select = traj.select(opt.select)
        LOGGER.info('{0:d} atoms are selected for writing output.'
                    .format(len(select)))
        if align:
            _ = traj.select(align)
            LOGGER.info('{0:d} atoms are selected for aligning frames.'
                        .format(len(_)))

    out = prody.DCDFile(opt.output, 'w')
    count = 0
    goto = False
    if opt.stride != 1:
        goto = True
    slc = slice(opt.first, opt.last, opt.stride).indices(len(traj)+1)
    for i in range(*slc):
        if goto:
            traj.goto(i)
        frame = traj.next()
        if align:
            frame.superpose()
            out.write(select._getCoords(), frame.getUnitcell())
        else:
            out.write(frame._getCoords(), frame.getUnitcell())
        count += 1
    traj.close()
    out.close()
    LOGGER.info("{0:d} frames are written into '{1:s}'."
                .format(count, opt.output))

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
        dest='select', metavar='SELSTR', 
        help='atom selection (default: "%(default)s")')

    subparser.add_argument('-o', '--output', type=str, metavar='FILENAME', 
        default='trajectory.dcd',
        help='output filename (default: "%(default)s")')

    subparser.add_argument('-n', '--num', default=False, action='store_true',
        dest='num',
        help='print the number of frames in each file and exit')

    group = subparser.add_mutually_exclusive_group()
    group.add_argument('--psf', 
        help='PSF filename (must have same number of atoms as DCDs)')
    group.add_argument('--pdb', 
        help='PDB filename (must have same number of atoms as DCDs)')

    subparser.add_argument('--first', metavar='INT', type=int, default=0,
        help="the first frame to be written to the output file "
             "(default: %(default)s, first frame)")
    subparser.add_argument('--last', metavar='INT', type=int, default=-1,
        help="the last frame to be written to the output file "
             "(default: %(default)s, last frame)")
    subparser.add_argument('--stride', metavar='INT', type=int, default=1,
        help="number of frames to skip when writing "
             "(default: %(default)s, skip none)")

    subparser.add_argument('--align', metavar='SELSTR', type=str,
        help="atom selection for aligning frames, one of PSF or PDB files "
             "must be provided")

    subparser.add_argument('dcd', nargs='+',
        help='DCD filename(s) (all must have same number of atoms)')

    subparser.set_defaults(subparser=subparser)
    subparser.set_defaults(func=prody_catdcd)
