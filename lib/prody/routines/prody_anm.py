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

"""Perform ANM calculations and output the results in plain text, NMD, and 
graphical formats."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

from actions import *
from nmaoptions import *

defaultVals = {'number-of-modes': 10, 'select':
    'protein and name CA or nucleic and name P C4\' C2', 'cutoff': 15.0, 
    'gamma': 1.0, 'model': 1, 'output-dir': '.', 'delimiter': ' ', 'extension':
    '.txt', 'number-format': '%12g', 'figure-format': 'pdf', 'resolution': 300,
    'height': 6.0, 'width': 8.0}

def prody_anm(pdbname, **kwargs):
    """Perform ANM calculations based on command line arguments.
    
    :arg pdb:  :term:`PDB` identifier or filename

    **NMA Parameters:**
    
    :arg number-of-modes: number of non-zero eigenvalues/vectors to
        calculate, default is 10
    :type number-of-modes: int
                        
    :arg select: atom selection string, default is :
        term:`protein and name CA or nucleic and name P C4' C2`
    :type select: str
    
    :arg cutoff: cutoff distance (A), default is 15.0
    :type cutoff: float
    
    :arg gamma: spring constant, default is 1.0
    :type gamma: float
    
    :arg model: model that will be used in the calculations, default is 1
    :type model: int

    **NMA Output options:**
    
    :arg all-output: write all outputs
    
    :arg output-dir: output directory, default is "."
    
    :arg eigenvs: write eigenvectors/values
    
    :arg cross-correlations: write cross-correlations
    
    :arg square-fluctuations: write square-fluctuations
    
    :arg covariance: write covariance matrix
    
    :arg npz: write compressed ProDy data file
    
    :arg extend:  output NMD files for the model extended to
        "backbone" ("bb") or "all" atoms of the residue, model
        must have one node per residue
    :type extend: str
    
    :arg beta-factors: write B-factors
    
    :arg hessian: write Hessian matrix
    
    :arg kirchhoff: write Kirchhoff matrix

    **Output options:**
    
    :arg file-prefix: prefix for output files, default is :file:`_anm`
    :type file-prefix: str
    
    :arg number-format: number format, default is `%12g`
    :type number-format: str
    
    :arg delimiter: delimiter, default is ` `
    :type delimiter: str
    
    :arg extension: file extension, default is `.txt`
    :type extension: str

    **Figures output options:**
    
    :arg all-figures: save all figures  
    
    :arg cross-correlations-figure: save cross-correlations
    
    :arg square-fluctuations-figure: save square-fluctuations
    
    :arg beta-factors-figure: save beta-factors
    
    :arg contact-map: save contact map (Kirchhoff matrix)

    **Figure options:**
    
    :arg figure-format: figure format, one of eps, pdf, png, ps, raw, rgba,
        svg, svgz, default is `pdf`
    :type figure-format: str
    
    :arg resolution: figure resolution (dpi), default is 300
    :type resolution: int
    
    :arg width: figure width (inch), default is 8.0
    :type width: float
    
    :arg height: figure height (inch), default is 6.0
    :type height: float    
    """
    
    outdir = kwargs.get('output-dir',defaultVals.get('outpur-dir'))
    print outdir
    if not os.path.isdir(outdir):   
        raise ValueError('{0:s} is not a valid path'.format(repr(outdir)))
    
        
    import numpy as np
    import prody
    LOGGER = prody.LOGGER


    prefix = kwargs.get('file-prefix', None)
    cutoff = kwargs.get('cutoff',defaultVals.get('cutoff'))
    gamma = kwargs.get('gamma',defaultVals.get('gamma')) 
    nmodes = kwargs.get('number-of-modes',defaultVals.get('number-of-modes'))
    model = kwargs.get('model',defaultVals.get('model'))
    selstr = kwargs.get('select',defaultVals.get('select'))     
    
    pdb = prody.parsePDB(pdbname, model=model)
    #if prefix == '_anm':
    #    prefix = pdb.getTitle() + '_anm'
    prefix = prefix or (pdb.getTitle() + '_anm')

    select = pdb.select(selstr)
    if select is None:
        LOGGER.warn('Selection {0:s} did not match any atoms.'
                        .format(repr(selstr)))
        return
    LOGGER.info('{0:d} atoms will be used for ANM calculations.'
                .format(len(select)))

    anm = prody.ANM(pdb.getTitle())
    anm.buildHessian(select, cutoff, gamma)
    anm.calcModes(nmodes)
    LOGGER.info('Writing numerical output.')
    if kwargs.get('npz'):
        prody.saveModel(anm)
    prody.writeNMD(os.path.join(outdir, prefix + '.nmd'), anm, select)
    
    extend = kwargs.get('extend')
    if extend:
        if extend == 'all':
            extended = prody.extendModel(anm, select, pdb)        
        else:
            extended = prody.extendModel(anm, select, select | pdb.bb)
        prody.writeNMD(os.path.join(outdir, prefix + '_extended_' + 
                       extend + '.nmd'), *extended)
        
    outall = kwargs.get('all-output')
    delim = kwargs.get('delimiter',defaultVals.get('delimiter')),
    ext = kwargs.get('extension',defaultVals.get('extension'))
    format = kwargs.get('number-format',defaultVals.get('%12g'))

    if outall or kwargs.get('eigenvs'):
        prody.writeArray(os.path.join(outdir, prefix + '_evectors'+ext), 
                         anm.getArray(), delimiter=delim, format=format)
        prody.writeArray(os.path.join(outdir, prefix + '_evalues'+ext), 
                         anm.getEigvals(), delimiter=delim, format=format)
    if outall or kwargs.get('beta-factors'):
        from prody.utilities import openFile
        fout = openFile(prefix + '_beta.txt', 'w', folder=outdir)
        fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4s} {0[3]:5s} {0[4]:5s}\n'
                       .format(['C', 'RES', '####', 'Exp.', 'The.']))
        for data in zip(select.getChids(), select.getResnames(), 
                        select.getResnums(), select.getBetas(), 
                        prody.calcTempFactors(anm, select)):
            fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4d} {0[3]:5.2f} {0[4]:5.2f}\n'
                       .format(data))
        fout.close()
    if outall or kwargs.get('covariance'):
        prody.writeArray(os.path.join(outdir, prefix + '_covariance'+ext), 
                         anm.getCovariance(), delimiter=delim, format=format)
    if outall or kwargs.get('cross-correlations'):
        prody.writeArray(os.path.join(outdir, prefix + '_cross-correlations' 
                                                     + ext), 
                         prody.calcCrossCorr(anm), delimiter=delim, 
                         format=format)
    if outall or kwargs.get('hessian'):
        prody.writeArray(os.path.join(outdir, prefix + '_hessian'+ext), 
                         anm.getHessian(), delimiter=delim, format=format)
    if outall or kwargs.get('kirchhoff'):
        prody.writeArray(os.path.join(outdir, prefix + '_kirchhoff'+ext), 
                         anm.getKirchhoff(), delimiter=delim, format=format)
    if outall or kwargs.get('square-fluctuations'):
        prody.writeArray(os.path.join(outdir, prefix + '_sqflucts'+ext), 
                         prody.calcSqFlucts(anm), delimiter=delim, 
                         format=format)
          
    figall, cm = kwargs.get('all-figures'), kwargs.get('cm')
    cc = kwargs.get('cc')
    sf = kwargs.get('sq')
    bf = kwargs.get('bf') 

    if figall or cc or sf or bf or cm: 
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warning('Matplotlib could not be imported. '
                           'Figures are not saved.')
        else:
            LOGGER.info('Saving graphical output.')
            format, width, height, dpi = \
                kwargs.get('figure-format',defaultVals.get('figure-format')),\
                kwargs.get('width',defaultVals.get('width')),\
                kwargs.get('height',defaultVals.get('height')),\
                kwargs.get('resolution',defaultVals.get('resolution'))
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
    """This command performs ANM calculations for given PDB structure and \
outputs results in NMD format. If an identifier is passed, structure file \
will be downloaded from the PDB FTP server.

Fetch PDB 1p38, run ANM calculations using default parameters, and write \
NMD file:
    
  $ prody anm 1p38
    
Fetch PDB 1aar, run ANM calculations using default parameters for chain A \
carbon alpha atoms with residue numbers less than 70, and save all of the \
graphical output files:

  $ prody anm 1aar -s "calpha and chain A and resnum < 70" -A"""
    )
    
    group = subparser.add_argument_group('parameters')
    #group = addNMAParameters(subparser)
    group.add_argument('-c', '--cutoff', dest='cutoff', type=float, 
                       default=defaultVals.get('cutoff'), metavar='FLOAT', 
                       help='cutoff distance (A) (default: %(default)s)')
    group.add_argument('-g', '--gamma', dest='gamma', type=float, 
                       default=defaultVals.get('gamma'), metavar='FLOAT', 
                       help='spring constant (default: %(default)s)')
    group.add_argument('-m', '--model', dest='model', type=int, 
                       default=defaultVals.get('model'), metavar='INT', 
                       help=('model that will be used in the calculations'))
    group.add_argument('-n', '--number-of-modes', dest='number-of-modes',
                       type=int, default=defaultVals.get('umber-of-modes'),
                       metavar='INT', help=('number of non-zero eigenvalues\
                                            /vectors to calculate'))
    group.add_argument('-s', '--select', dest='select', type=str, 
                        default=defaultVals.get('select'), metavar='SEL',
                        help='atom selection  (default: %(default)s)')
                        
    group = subparser.add_argument_group('output')
    group.add_argument('-a', '--all-output', dest='all-output',
                       action='store_true', default=False,
                       help='write all outputs')
    subparser.add_argument('-o', '--output-dir', dest='output-dir',
                        type=str, default=defaultVals.get('output-dir'),
                        metavar='PATH', help=('output directory for\
                                              saving file(s)'))
    group.add_argument('-b', '--beta-factors', dest='beta', action='store_true', 
                       default=False, help='write B-factors')
    group.add_argument('-l', '--hessian', dest='hessian', action='store_true', 
                       default=False, help='write Hessian matrix')
    group.add_argument('-k', '--kirchhoff', dest='kirchhoff', action='store_true', 
                       default=False, help='write Kirchhoff matrix')
    group.add_argument('-e', '--eigenvs', dest='eigenvs', action='store_true', 
                       default=False, help='write eigenvectors/values')
    group.add_argument('-r', '--cross-correlations', dest='cross-correlations',
                       action='store_true', default=False,
                       help='write cross-correlations')
    group.add_argument('-q', '--square-fluctuations', dest='square-fluctuations',
                       action='store_true', default=False,
                       help='write square-fluctuations')
    group.add_argument('-v', '--covariance', dest='covariance',
                       action='store_true', default=False,
                       help='write covariance')
    group.add_argument('-z', '--npz', dest='npz',action='store_true',
                       default=False, help='write compressed ProDy data file')
    group.add_argument('-t', '--extend', dest='extend', type=str, 
                        default='', metavar='STR', 
                        help=('output NMD files for the model extended to\
                        "backbone" ("bb") or "all" atoms of the residue,\
                        model must have one node per residue'))
    

    group = subparser.add_argument_group('output options')
    group.add_argument('-p', '--file-prefix', dest='file-prefix', type=str, 
                        default='', metavar='STR', 
                        help=('output filename prefix (default: pdb_anm)'))
    group.add_argument('-f', '--number-format', dest='number-format', type=str, 
                        default=defaultVals.get('number-format'), metavar='STR', 
                        help=('number_format (default: %(default)s)'))
    group.add_argument('-d', '--delimiter', dest='delimiter', type=str, 
                        default=defaultVals.get('delimiter'), metavar='STR', 
                        help=('delimiter (default: %(default)s)'))
    group.add_argument('-x', '--extension', dest='extension', type=str, 
                        default=defaultVals.get('extension'), metavar='STR', 
                        help=('file extension (default: %(default)s)'))
    

    #group = addNMAFigures(subparser)
    group = subparser.add_argument_group('figures')
    group.add_argument('-A', '--all-figures', dest='all-figures', 
                       action='store_true', default=False, 
                       help='save allfigures')
    group.add_argument('-B', '--beta-factors-figure', dest='bf', 
                       action='store_true', default=False, 
                       help='save beta-factors')
    group.add_argument('-K', '--contact-map', dest='cm', action='store_true', 
                       default=False, 
                       help='save contact map (Kirchhoff matrix)')
    group.add_argument('-Q', '--square-fluctuations-figure', dest='sq',
                       action='store_true', default=False, 
                       help='save square-fluctuations-figure')
    group.add_argument('-R', '--cross-correlations-figure', dest='cc',
                       action='store_true', default=False, 
                       help='save cross-correlations-figure')

    group = subparser.add_argument_group('figure options')
    group.add_argument('-F', '--figure-format', dest='figure-format', type=str, 
                        default=defaultVals.get('figure-format'), metavar='STR', 
                        help=('figure format, one of eps, pdf, png, ps, raw,\
                        rgba,svg, svgz (default: %(default)s)'))
    group.add_argument('-D', '--resolution', dest='resolution', type=int, 
                        default=defaultVals.get('resolution'), metavar='INT', 
                        help=('figure resolution (dpi) (default: %(default)s)'))
    group.add_argument('-W', '--width', dest='width', type=float, 
                        default=defaultVals.get('width'), metavar='FLOAT', 
                        help=('figure width (inch) (default: %(default)s)'))
    group.add_argument('-H', '--height', dest='height', type=float, 
                        default=defaultVals.get('height'), metavar='FLOAT', 
                        help=('figure height (inch) (default: %(default)s)'))


    subparser.add_argument('pdb', help='PDB identifier or filename')

    subparser.set_defaults(func=lambda ns: prody_anm(ns.pdb, **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
