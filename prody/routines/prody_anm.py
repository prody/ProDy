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

"""This module defines :program:`prody anm` command."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

from actions import *
from nmaoptions import *

def prody_anm(opt):
    """Perform ANM calculations based on command line arguments."""
    
    outdir = opt.outdir
    if not os.path.isdir(outdir):
        opt.subparser.error('{0:s} is not a valid path'.format(outdir))
        
    import numpy as np
    import prody
    LOGGER = prody.LOGGER


    pdb = opt.pdb
    prefix = opt.prefix
    cutoff, gamma = opt.cutoff, opt.gamma, 
    nmodes, selstr, model = opt.nmodes, opt.select, opt.model
    
    pdb = prody.parsePDB(pdb, model=model)
    if prefix == '_anm':
        prefix = pdb.getTitle() + '_anm'

    select = pdb.select(selstr)
    if select is None:
        opt.subparser('Selection "{0:s}" do not match any atoms.'
                       .format(selstr))
    LOGGER.info('{0:d} atoms will be used for ANM calculations.'
                .format(len(select)))

    anm = prody.ANM(pdb.getTitle())
    anm.buildHessian(select, cutoff, gamma)
    anm.calcModes(nmodes)
    LOGGER.info('Writing numerical output.')
    if opt.npz:
        prody.saveModel(anm)
    prody.writeNMD(os.path.join(outdir, prefix + '.nmd'), anm, select)

    outall = opt.all
    delim, ext, format = opt.delim, opt.ext, opt.numformat

    if outall or opt.eigen:
        prody.writeArray(os.path.join(outdir, prefix + '_evectors'+ext), 
                         anm.getArray(), delimiter=delim, format=format)
        prody.writeArray(os.path.join(outdir, prefix + '_evalues'+ext), 
                         anm.getEigenvalues(), delimiter=delim, format=format)
    if outall or opt.beta:
        fout = prody.openFile(prefix + '_beta.txt', 'w', folder=outdir)
        fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4s} {0[3]:5s} {0[4]:5s}\n'
                       .format(['C', 'RES', '####', 'Exp.', 'The.']))
        for data in zip(select.getChids(), select.getResnames(), 
                        select.getResnums(), select.getBetas(), 
                        prody.calcTempFactors(anm, select)):
            fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4d} {0[3]:5.2f} {0[4]:5.2f}\n'
                       .format(data))
        fout.close()
    if outall or opt.covar:
        prody.writeArray(os.path.join(outdir, prefix + '_covariance'+ext), 
                         anm.getCovariance(), delimiter=delim, format=format)
    if outall or opt.ccorr:
        prody.writeArray(os.path.join(outdir, prefix + '_cross-correlations' 
                                                     + ext), 
                         prody.calcCrossCorr(anm), delimiter=delim, 
                         format=format)
    if outall or opt.hessian:
        prody.writeArray(os.path.join(outdir, prefix + '_hessian'+ext), 
                         anm.getHessian(), delimiter=delim, format=format)
    if outall or opt.kirchhoff:
        prody.writeArray(os.path.join(outdir, prefix + '_kirchhoff'+ext), 
                         anm.getKirchhoff(), delimiter=delim, format=format)
    if outall or opt.sqflucts:
        prody.writeArray(os.path.join(outdir, prefix + '_sqflucts'+ext), 
                         prody.calcSqFlucts(anm), delimiter=delim, 
                         format=format)
          
    figall, cc, sf, bf, cm = opt.figures, opt.cc, opt.sf, opt.bf, opt.cm

    if figall or cc or sf or bf or cm: 
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warning('Matplotlib could not be imported. '
                           'Figures are not saved.')
        else:
            LOGGER.info('Saving graphical output.')
            format, width, height, dpi = \
                opt.figformat, opt.width, opt.height, opt.dpi
            format = format.lower()
            if figall or cc:
                plt.figure(figsize=(width, height))
                prody.showCrossCorr(anm)
                plt.savefig(os.path.join(outdir, prefix + '_cc.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or cm:
                plt.figure(figsize=(width, height))
                prody.showContactMap(anm)
                plt.savefig(os.path.join(outdir, prefix + '_cm.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                prody.showSqFlucts(anm)
                plt.savefig(os.path.join(outdir, prefix + '_sf.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or bf:
                plt.figure(figsize=(width, height))
                bexp = select.getBetas()
                bcal = prody.calcTempFactors(anm, select)
                plt.plot(bexp, label='Experimental')
                plt.plot(bcal, label=('Theoretical (R={0:.2f})'
                                        .format(np.corrcoef(bcal, bexp)[0,1])))
                plt.legend(prop={'size': 10})
                plt.xlabel('Node index')
                plt.ylabel('Experimental B-factors')
                plt.title(pdb.getTitle() + ' B-factors')
                plt.savefig(os.path.join(outdir, prefix + '_bf.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
                
def addCommand(commands):

    subparser = commands.add_parser('anm', 
        help='perform anisotropic network model calculations')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')
    subparser.set_defaults(usage_example=
    """This command performs ANM calculations for given PDB structure and outputs \
    results in NMD format. If an identifier is passed, structure file will be \
    downloaded from the PDB FTP server.

    Fetch PDB 1p38, run ANM calculations using default parameters, and write \
    NMD file:
        
      $ prody anm 1p38
        
    Fetch PDB 1aar, run ANM calculations using default parameters for chain A \
    carbon alpha atoms with residue numbers less than 70, and save all of the \
    graphical output files:

      $ prody anm 1aar -s "calpha and chain A and resnum < 70" -A"""
    )

    group = addNMAParameters(subparser)
    group.add_argument('-c', '--cutoff', dest='cutoff', type=float, 
                       default=15.0, metavar='FLOAT', 
                       help='cutoff distance (A) (default: %(default)s)')
    group.add_argument('-g', '--gamma', dest='gamma', type=float, 
                       default=1.0, metavar='FLOAT', 
                       help='spring constant (default: %(default)s)')
    group.add_argument('-m', '--model', dest='model', type=int, 
                       default=1, metavar='INT', 
                       help=('model that will be used in the calculations')) 
                        
    group = addNMAOutput(subparser)
    group.add_argument('-b', '--beta-factors', dest='beta', action='store_true', 
                       default=False, help='write B-factors')
    group.add_argument('-l', '--hessian', dest='hessian', action='store_true', 
                       default=False, help='write Hessian matrix')
    group.add_argument('-k', '--kirchhoff', dest='kirchhoff', action='store_true', 
                       default=False, help='write Kirchhoff matrix')

    group = addNMAOutputOptions(subparser, '_anm')

    group = addNMAFigures(subparser)
    group.add_argument('-B', '--beta-factors-figure', dest='bf', 
                       action='store_true', default=False, 
                       help='save beta-factors')
    group.add_argument('-K', '--contact-map', dest='cm', action='store_true', 
                       default=False, 
                       help='save contact map (Kirchhoff matrix)')        

    group = addNMAFigureOptions(subparser)

    subparser.add_argument('pdb', help='PDB identifier or filename')

    subparser.set_defaults(func=prody_anm)
    subparser.set_defaults(subparser=subparser)
