# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan
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

"""This module contains functions which are used as command line programs."""

__author__ = 'Ahmet Bakan, Lidio Meireles'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan, Lidio Meireles'

import sys
import os.path

import numpy as np

from prody import *
from prody import ProDyLogger as LOGGER

__all__ = ['anm', 'gnm', 'pca', 'alignmodels', 'biomolecule', 'blastpdb',
           'fetchpdb', 'pdbselect']

PY3K = sys.version_info[0] > 2

if PY3K:
    raise NotImplemented('ProDy scripts are not yet implemented for Python 3')
else:
    from optparse import OptionParser, OptionGroup

DEFAULT_FIGURE_FORMAT = 'pdf'
DEFAULT_FILE_EXT = '.txt'
DEFAULT_NUMBER_FORMAT = '%12g'
DEFAULT_DELIMITER = ' '

def addOptions(parser):
    parser.add_option('', '--examples', dest='examples', action='store_true', 
                      default=False, 
                      help='show usage examples and exit')
    parser.add_option('', '--quiet', dest='silent', action='store_true', 
                      default=False, 
                      help='don\'t print status messages to stdout')

def addParameters(parser):
    parser.add_option('-n', '--number-of-modes', dest='nmodes', type='int', 
                      default=10, metavar='INT', 
                      help=('number of non-zero eigenvalues/vectors to '
                            'calculate, default is %default'))    
    parser.add_option('-s', '--select', dest='select', type='string', 
                      default='calpha', metavar='STRING',
                      help='selection string, default is "%default"')

def addOutput(parser):
    parser.add_option('-a', '--all-output', dest='all', action='store_true', 
                      default=False, help='write all outputs')
    parser.add_option('-e', '--eigenvs', dest='eigen', action='store_true', 
                      default=False, help='write eigenvectors/values')
    parser.add_option('-r', '--cross-correlations', dest='ccorr', 
                      action='store_true', 
                      default=False, help='write cross-correlations')
    parser.add_option('-q', '--square-fluctuations', dest='sqflucts', 
                      action='store_true', 
                      default=False, help='write square-fluctuations')
    parser.add_option('-v', '--covariance', dest='covar', action='store_true', 
                     default=False, help='write covariance matrix')

def addNumericalOptions(parser, prefix):
    parser.add_option('-p', '--file-prefix', dest='prefix', type='string', 
                      default=prefix, metavar='STRING', 
                      help=('prefix for output files, default is "PDB%default"'))
    parser.add_option('-f', '--number-format', dest='numformat', type='string', 
                     default=DEFAULT_NUMBER_FORMAT, metavar='STRING',
                     help=('delimiter, default is "%default"'))
    parser.add_option('-d', '--delimiter', dest='delim', type='string', 
                     default=DEFAULT_DELIMITER, metavar='STRING',
                     help=('delimiter, default is "%default"'))
    parser.add_option('-x', '--extension', dest='ext', type='string', 
                     default=DEFAULT_FILE_EXT, metavar='STRING',
                     help=('file extension, default is "%default"'))

def addFigures(parser):
    parser.add_option('-A', '--all-figures', dest='figures', 
                      action='store_true', 
                     default=False, help='save all figures')
    parser.add_option('-R', '--cross-correlations-figure', dest='cc', 
                      action='store_true', default=False, 
                      help='save cross-correlations')
    parser.add_option('-Q', '--square-fluctuations-figure', dest='sf', 
                      action='store_true', default=False, 
                      help='save square-fluctuations')

def addFigureOptions(parser):
    parser.add_option('-F', '--figure-format', dest='figformat', type='string', 
                      default=DEFAULT_FIGURE_FORMAT, metavar='STRING',
                      help=('figure format, one of eps, pdf, png, ps, raw, '
                           'rgba, svg, svgz, default is "%default"'))
    parser.add_option('-D', '--resolution', dest='dpi', type='int', 
                      default=300, metavar='INT', 
                      help='figure resolution (dpi), default is %default')
    parser.add_option('-W', '--width', dest='width', type='float', 
                      default=8.0, metavar='FLOAT', 
                      help='figure width (inch), default is %default')
    parser.add_option('-H', '--height', dest='height', type='float', 
                      default=6.0, metavar='FLOAT', 
                      help='figure height (inch), default is %default')
    

def anm():
    """Perform ANM calculations based on command line arguments."""
    
    usage = """%prog [options] PDB

ProDy v{0:s} - Anisotropic Network Model

Perform ANM calculations for given PDB structure and output results in NMD 
format. If an identifier is passed, file from the PDB FTP server will be 
downloaded.""".format(prody.__version__)
 
    parser = OptionParser(usage=usage)
    addOptions(parser)
    
    group = OptionGroup(parser, 'Parameters')
    group.add_option('-c', '--cutoff', dest='cutoff', type='float', 
                     default=15.0, metavar='FLOAT', 
                     help='cutoff distance, default is %default A')
    group.add_option('-g', '--gamma', dest='gamma', type='float', 
                     default=1.0, metavar='FLOAT', 
                     help='spring constant, default is %default')
    group.add_option('-m', '--model', dest='model', type='int', 
                     default=1, metavar='INT', 
                     help=('model that will be used in the calculations'))   
    addParameters(group)
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Output')
    addOutput(group)
    group.add_option('-b', '--beta-factors', dest='beta', action='store_true', 
                     default=False, help='write B-factors')
    group.add_option('-l', '--hessian', dest='hessian', action='store_true', 
                     default=False, help='write Hessian matrix')
    group.add_option('-k', '--kirchhoff', dest='kirchhoff', action='store_true', 
                     default=False, help='write Kirchhoff matrix')
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Output Options')
    addNumericalOptions(group, '_anm')
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Figures')
    addFigures(group)
    group.add_option('-B', '--beta-factors-figure', dest='bf', 
                     action='store_true', default=False, 
                     help='save beta-factors')
    group.add_option('-K', '--contact-map', dest='cm', action='store_true', 
                     default=False, 
                     help='save contact map (Kirchhoff matrix)')        
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Figure Options')
    addFigureOptions(group)
    parser.add_option_group(group)
    
    usage_examples="""
Fetch PDB 1p38, run ANM calculations using default parameters, and write 
NMD file:
    
  $ anm.py 1p38
    
Fetch PDB 1aar, run ANM calculations using default parameters for chain A 
carbon alpha atoms with residue numbers less than 70, and save all of the
graphical output files:

  $ anm.py 1aar -s "calpha and chain A and resnum < 70" -A
"""
    
    opt, args = parser.parse_args()
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) < 1:
        parser.print_help()
        print "\nError: PDB missing\n"
        sys.exit(-1)
    if opt.silent:
        changeVerbosity('warning')
                
    pdb = args[0]
    prefix = opt.prefix
    cutoff, gamma = opt.cutoff, opt.gamma, 
    nmodes, selstr, model = opt.nmodes, opt.select, opt.model
    
    pdb = parsePDB(pdb, model=model)
    if prefix == '_anm':
        prefix = pdb.getName() + '_anm'

    select = pdb.select(selstr)
    if select is None:
        LOGGER.warning('Selection "{0:s}" do not match any atoms.'
                       .format(selstr))
        sys.exit(-1)
    LOGGER.info('{0:d} atoms will be used for ANM calculations.'
                .format(len(select)))

    anm = ANM(pdb.getName())
    anm.buildHessian(select, cutoff, gamma)
    anm.calcModes(nmodes)
    LOGGER.info('Writing numerical output.')
    writeNMD(prefix + '.nmd', anm, select)

    outall = opt.all
    delim, ext, format = opt.delim, opt.ext, opt.numformat
    eigen, sqflucts, ccorr, covar = \
        opt.eigen, opt.sqflucts, opt.ccorr, opt.covar
    beta, hessian, kirchhoff = opt.beta, opt.hessian, opt.kirchhoff 

    if outall or eigen:
        writeArray(prefix + '_evectors'+ext, anm.getArray(), 
            delimiter=delim, format=format)
        writeArray(prefix + '_evalues'+ext, anm.getEigenvalues(), 
            delimiter=delim, format=format)
    
    if outall or beta:
        fout = open(prefix + '_beta.txt', 'w')
        fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4s} {0[3]:5s} {0[4]:5s}\n'
                       .format(['C', 'RES', '####', 'Exp.', 'The.']))
        for data in zip(select.getChainIdentifiers(),
                    select.getResidueNames(), select.getResidueNumbers(),
                    select.getTempFactors(), calcTempFactors(anm, select)):
            fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4d} {0[3]:5.2f} {0[4]:5.2f}\n'
                       .format(data))
        fout.close()
    if outall or covar:
        writeArray(prefix + '_covariance'+ext, anm.getCovariance(), 
            delimiter=delim, format=format)
    if outall or ccorr:
        writeArray(prefix + '_cross-correlations'+ext, 
            calcCrossCorrelations(anm), delimiter=delim, format=format)
    if outall or hessian:
        writeArray(prefix + '_hessian'+ext, anm.getHessian(), 
            delimiter=delim, format=format)
    if outall or kirchhoff:
        writeArray(prefix + '_kirchhoff'+ext, anm.getKirchhoff(), 
            delimiter=delim, format=format)
    if outall or sqflucts:
        writeArray(prefix + '_sqflucts'+ext, calcSqFlucts(anm), 
            delimiter=delim, format=format)
          
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
                showCrossCorrelations(anm)
                plt.savefig(prefix + '_cc.'+format, dpi=dpi, format=format)
                plt.close('all')
            if figall or cm:
                plt.figure(figsize=(width, height))
                showContactMap(anm)
                plt.savefig(prefix + '_cm.'+format, dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                showSqFlucts(anm)
                plt.savefig(prefix + '_sf.'+format, dpi=dpi, format=format)
            if figall or bf:
                plt.figure(figsize=(width, height))
                bexp = select.getTempFactors()
                bcal = calcTempFactors(anm, select)
                plt.plot(bexp, label='Experimental')
                plt.plot(bcal, label=('Theoretical (R={0:.2f})'
                                        .format(np.corrcoef(bcal, bexp)[0,1])))
                plt.legend()
                plt.xlabel('Node index')
                plt.ylabel('Experimental B-factors')
                plt.title(pdb.getName() + ' B-factors')
                plt.savefig(prefix + '_bf.'+format, dpi=dpi, format=format)

def gnm():
    """Perform GNM calculations based on command line arguments."""

    usage = """%prog [options] PDB

ProDy v{0:s} - Gaussian Network Model

Perform GNM calculations for given PDB structure and output eigenvalues and
eigenvectors. If an identifier is passed, file from the PDB FTP server will be 
downloaded.""".format(prody.__version__)
 
    parser = OptionParser(usage=usage)
    addOptions(parser)
    
    group = OptionGroup(parser, 'Parameters')
    group.add_option('-c', '--cutoff', dest='cutoff', type='float', 
                     default=10.0, metavar='FLOAT', 
                     help='cutoff distance, default is %default A')
    group.add_option('-g', '--gamma', dest='gamma', type='float', 
                     default=1.0, metavar='FLOAT', 
                     help='spring constant, default is %default')
    group.add_option('-m', '--model', dest='model', type='int', 
                     default=1, metavar='INT', 
                     help=('model that will be used in the calculations'))   
    addParameters(group)
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Output')
    addOutput(group)
    group.add_option('-b', '--beta-factors', dest='beta', action='store_true', 
                     default=False, help='write B-factors')
    group.add_option('-k', '--kirchhoff', dest='kirchhoff', action='store_true', 
                     default=False, help='write Kirchhoff matrix')
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Output Options')
    addNumericalOptions(group, '_gnm')
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Figures')
    addFigures(group)
    group.add_option('-B', '--beta-factors-figure', dest='bf', 
                     action='store_true', default=False, 
                     help='save beta-factors')
    group.add_option('-K', '--contact-map', dest='cm', action='store_true', 
                     default=False, 
                     help='save contact map (Kirchhoff matrix)')    
    group.add_option('-M', '--mode-figure', dest='modes', type='string', 
                     default='', metavar='STRING',
                     help=('save mode shape figures for specified modes, '
                           'e.g. "1-3 5" for modes 1, 2, 3 and 5'))
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Figure Options')
    addFigureOptions(group)
    parser.add_option_group(group)
    
    usage_examples="""
Fetch PDB 1p38, run GNM calculations using default parameters, and results:
    
  $ gnm.py 1p38
    
Fetch PDB 1aar, run GNM calculations with cutoff distance 7 angstrom for 
chain A carbon alpha atoms with residue numbers less than 70, and 
save all of the graphical output files:

  $ gnm.py 1aar -c 7 -s "calpha and chain A and resnum < 70" -A
"""
    opt, args = parser.parse_args()
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) < 1:
        parser.print_help()
        print "\nError: PDB missing\n"
        sys.exit(-1)
    if opt.silent:
        changeVerbosity('warning')
        
    pdb = args[0]
    prefix = opt.prefix

    cutoff, gamma = opt.cutoff, opt.gamma, 
    nmodes, selstr, model = opt.nmodes, opt.select, opt.model
    
    pdb = parsePDB(pdb, model=model)
    if prefix == '_gnm':
        prefix = pdb.getName() + '_gnm'

    select = pdb.select(selstr)
    if select is None:
        LOGGER.warning('Selection "{0:s}" do not match any atoms.'
                       .format(selstr))
        sys.exit(-1)
    LOGGER.info('{0:d} atoms will be used for GNM calculations.'
                .format(len(select)))

    gnm = GNM(pdb.getName())
    gnm.buildKirchhoff(select, cutoff, gamma)
    gnm.calcModes(nmodes)
    LOGGER.info('Writing numerical output.')
    outall = opt.all
    delim, ext, format = opt.delim, opt.ext, opt.numformat
    
    writeArray(prefix + '_evectors'+ext, gnm.getArray(), 
        delimiter=delim, format=format)
    writeArray(prefix + '_evalues'+ext, gnm.getEigenvalues(), 
        delimiter=delim, format=format)
    
    beta, sqflucts = opt.beta, opt.sqflucts
    ccorr, covar = opt.ccorr, opt.covar
    kirchhoff = opt.kirchhoff 

    
    if outall or beta:
        fout = open(prefix + '_beta.txt', 'w')
        fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4s} {0[3]:5s} {0[4]:5s}\n'
                       .format(['C', 'RES', '####', 'Exp.', 'The.']))
        for data in zip(select.getChainIdentifiers(),
                    select.getResidueNames(), select.getResidueNumbers(),
                    select.getTempFactors(), calcTempFactors(gnm, select)):
            fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4d} {0[3]:5.2f} {0[4]:5.2f}\n'
                       .format(data))
        fout.close()
    if outall or covar:
        writeArray(prefix + '_covariance'+ext, gnm.getCovariance(), 
            delimiter=delim, format=format)
    if outall or ccorr:
        writeArray(prefix + '_cross-correlations'+ext, 
                   calcCrossCorrelations(gnm), 
                   delimiter=delim, format=format)
    if outall or kirchhoff:
        writeArray(prefix + '_kirchhoff'+ext, gnm.getKirchhoff(), 
            delimiter=delim, format=format)
    if outall or sqflucts:
        writeArray(prefix + '_sqfluct'+ext, calcSqFlucts(gnm), 
            delimiter=delim, format=format)
          
    figall, cc, sf, bf, cm, modes = \
        opt.figures, opt.cc, opt.sf, opt.bf, opt.cm, opt.modes
    if figall or cc or sf or bf or cm or modes: 
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
                showCrossCorrelations(gnm)
                plt.savefig(prefix + '_cc.'+format, dpi=dpi, format=format)
                plt.close('all')
            if figall or cm:
                plt.figure(figsize=(width, height))
                showContactMap(gnm)
                plt.savefig(prefix + '_cm.'+format, dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                showSqFlucts(gnm)
                plt.savefig(prefix + '_sf.'+format, dpi=dpi, format=format)
            if figall or bf:
                plt.figure(figsize=(width, height))
                bexp = select.getTempFactors()
                bcal = calcTempFactors(gnm, select)
                plt.plot(bexp, label='Experimental')
                plt.plot(bcal, label=('Theoretical (corr coef = {0:.2f})'
                                        .format(np.corrcoef(bcal, bexp)[0,1])))
                plt.legend()
                plt.xlabel('Node index')
                plt.ylabel('Experimental B-factors')
                plt.title(pdb.getName() + ' B-factors')
                plt.savefig(prefix + '_bf.'+format, dpi=dpi, format=format)
            if modes: 
                modes = modes.split('-')
                try:
                    if len(modes) == 1:
                        modes = gnm(int(modes[0])-1)
                    else: 
                        modes = gnm[int(modes[0])-1:int(modes[1])]
                except:
                    LOGGER.warning('An error occured, and mode shapes were'
                                   'not plotted. See help: gnm.py -h')
                else:
                    for mode in modes:
                        plt.figure(figsize=(width, height))
                        showMode(mode)
                        plt.grid()
                        plt.savefig(prefix + '_mode_'+str(mode.getIndex()+1)+
                                    '.'+format, dpi=dpi, format=format)
                        


def pca():
    """Perform PCA calculations based on command line arguments."""
    
    usage = """%prog [options] PDB

ProDy v{0:s} - Principal Component Analysis

Perform PCA calculations for given PDB structure with multiple models and save 
results in NMD format. If an identifier is passed, file from the PDB FTP 
server will be downloaded.""".format(prody.__version__)
 
    parser = OptionParser(usage=usage)
    addOptions(parser)
    
    group = OptionGroup(parser, 'Parameters')
    addParameters(group)
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Output')
    addOutput(group)
    group.add_option('-j', '--projection', dest='proj', action='store_true', 
                     default=False, help='write projections onto PCs')
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Output Options')
    addNumericalOptions(group, '_pca')
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Figures')
    addFigures(group)
    group.add_option('-J', '--projection-figure', dest='sp', type='string', 
                     default='', metavar='STRING',
                     help=('save projections onto specified subspaces, e.g. '
                           '"1,2" for projections onto PCs 1 and 2; '
                           '"1,2 1,3" for projections onto PCs 1,2 and 1, 3; '
                           '"1 1,2,3" for projections onto PCs 1 and 1, 2, 3'))
    parser.add_option_group(group)
    
    group = OptionGroup(parser, 'Figure Options')
    addFigureOptions(group)
    parser.add_option_group(group)
    
    usage_examples="""
Fetch pdb 2k39, perform PCA calculations, and output NMD file:
    
    $ pca.py 2k39
    
Fetch pdb 2k39 and perform calculations for backbone of residues up to 71,
and save all output and figure files:

    $ pca.py 2k39 --select "backbone and resnum < 71" -a -A
"""
    
    opt, args = parser.parse_args()
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) != 1:
        parser.print_help()
        print "\nError: PDB missing\n"
        sys.exit(-1)
    if opt.silent:
        changeVerbosity('warning')
        
    pdb = args[0]
    prefix = opt.prefix
    nmodes, selstr = opt.nmodes, opt.select
    
    pdb = parsePDB(pdb)
    if prefix == '_pca':
        prefix = pdb.getName() + '_pca'
    select = pdb.select(selstr)
    if select is None:
        LOGGER.warning('Selection "{0:s}" do not match any atoms.'
                       .format(selstr))
        sys.exit(-1)
    LOGGER.info('{0:d} atoms will be used for PCA calculations.'
                .format(len(select)))
    ensemble = Ensemble(select)
    ensemble.iterpose()
    pca = PCA(pdb.getName())
    pca.buildCovariance(ensemble)
    pca.calcModes(nmodes)
    
    LOGGER.info('Writing numerical output.')
    writeNMD(prefix + '.nmd', pca, select)

    outall = opt.all
    delim, ext, format = opt.delim, opt.ext, opt.numformat
    eigen, sqflucts, proj = opt.eigen, opt.sqflucts, opt.proj
    ccorr, covar = opt.ccorr, opt.covar
    
    if outall or eigen:
        writeArray(prefix + '_evectors'+ext, pca.getArray(), 
            delimiter=delim, format=format)
        writeArray(prefix + '_evalues'+ext, pca.getEigenvalues(), 
            delimiter=delim, format=format)
    
    if outall or covar:
        writeArray(prefix + '_covariance'+ext, pca.getCovariance(), 
            delimiter=delim, format=format)
    if outall or ccorr:
        writeArray(prefix + '_cross-correlations'+ext, 
                   calcCrossCorrelations(pca), 
                   delimiter=delim, format=format)
    if outall or sqflucts:
        writeArray(prefix + '_sqfluct'+ext, calcSqFlucts(pca), 
            delimiter=delim, format=format)
    if outall or proj:
        writeArray(prefix + '_proj'+ext, calcProjection(ensemble, pca), 
                   delimiter=delim, format=format)
          
    figall, cc, sf, sp = opt.figures, opt.cc, opt.sf, opt.sp

    if figall or cc or sf: 
        format = format.lower()
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
                showCrossCorrelations(pca)
                plt.savefig(prefix + '_cc.'+format, dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                showSqFlucts(pca)
                plt.savefig(prefix + '_sf.'+format, dpi=dpi, format=format)
        
            if figall or sp:
                subsps = sp.split()
                for j, ss in enumerate(subsps):
                    try:
                        ss = ss.split(',')
                        modes = [int(i)-1 for i in ss]
                    except:
                        LOGGER.warning(ss + ' is not understood')
                    else:
                        if 1 <= len(modes) <= 3:
                            plt.figure(figsize=(width, height))
                            showProjection(ensemble, pca[modes])
                            plt.savefig(prefix + '_proj_' + '_'.join(ss) + 
                                '.' + format, dpi=dpi, format=format)
                            plt.close('all')
                    
    
def alignmodels():
    """Align models in a PDB file based on command line arguments."""
    
    usage = """Usage: {0:s} [options] PDB  

ProDy v{1:s} - Align Models

Align models in the PDB file using selected atoms and save aligned coordinate
sets.

""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    addOptions(parser)
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is PDB_aligned'))
    parser.add_option('-s', '--select', dest='select', type='string', 
                      default='calpha', metavar='STRING',
                      help='selection string, default is "%default"')
    parser.add_option('-m', '--model', dest='model', type='int', 
                      default=1, metavar='INT',
                      help=('model index onto which other models will be ' 
                            'superposed, default is %default'))
    usage_examples="""
Fetch pdb 2k39 and align models:
    
    $ alignmodels.py 2k39
    
Fetch pdb 2k39 and align models using backbone of residues with number smaller
than 71:

    $ alignmodels.py 2k39 --select "backbone and resnum < 71"
"""
        
    opt, args = parser.parse_args()
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) != 1:
        parser.print_help()
        print "\nError: PDB missing\n"
        sys.exit(-1)
    if opt.silent:
        changeVerbosity('warning')
        
    pdb = args[0]
    selstr, prefix, model = opt.select, opt.prefix, opt.model
    pdb = parsePDB(pdb)
    if prefix == '':
        prefix = pdb.getName() + '_aligned'
    pdbselect = pdb.select(selstr)
    if pdbselect is None:
        LOGGER.warning('Selection "{0:s}" do not match any atoms.'
                       .format(selstr))
        sys.exit(-1)
    prody.ProDyLogger.info('{0:d} atoms will be used for alignment.'
                           .format(len(pdbselect)))
    pdb.setActiveCoordsetIndex(model-1)
    alignCoordsets(pdb, selstr=selstr)
    rmsd = calcRMSD(pdb)
    LOGGER.info('Max RMSD: {0:0.2f} Mean RMSD: {1:0.2f}'
          .format(rmsd.max(), rmsd.mean()))
    outfn = prefix + '.pdb'
    LOGGER.info('Writing file: ' + outfn)
    writePDB(outfn, pdb)

def biomolecule():
    """Generate biomolecule coordinates based on command line arguments."""
    
    usage = """Usage: {0:s} [options] PDB  

ProDy v{1:s} - Generate Biomolecule Coordinates""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    addOptions(parser)
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is PDB_biomol_'))
    parser.add_option('-b', '--biomol', dest='biomol', type='int', 
                      default=None, metavar='INT',
                      help='index of the biomolecule, default is "%default"')
    usage_examples="""
Fetch pdb 2bfu and generate the biomolecular assembly:
    
  $ biomolecule.py 2bfu
"""
        
    opt, args = parser.parse_args()
    
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) != 1:
        parser.print_help()
        print "\nError: PDB missing\n"
        sys.exit(-1)
    if opt.silent:
        changeVerbosity('warning')

    pdb = args[0]
    prefix, biomol = opt.prefix, opt.biomol
        
    pdb, header = parsePDB(pdb, header=True)
    if prefix == '':
        prefix = pdb.getName()
        
    biomols = applyBiomolecularTransformations(header, pdb, biomol=biomol)
    if not isinstance(biomols, list):
        biomols = [biomols]
    
    for i, biomol in enumerate(biomols):
        if isinstance(biomol, Atomic):
            outfn = '{0:s}_biomol_{1:d}.pdb'.format(prefix, i+1)
            LOGGER.info('Writing {0:s}'.format(outfn))
            writePDB(outfn, biomol)
        elif isinstance(biomol, tuple):
            for j, part in enumerate(biomol):
                outfn = ('{0:s}_biomol_{1:d}_part_{2:d}.pdb'
                         .format(prefix, i+1, j+1))
                LOGGER.info('Writing {0:s}'.format(outfn))
                writePDB(outfn, part)

def readFirstSequenceFasta(seqfn):
    """Return first sequence from a file."""
    f = open(seqfn)
    lines = []
    n = 0
    for line in f.xreadlines():
        n += 1
        if n == 1 and line.startswith('>'): continue
        if n  > 1 and line.startswith('>'): break
        lines.append( line.strip() )
    f.close()
    return ''.join(lines)  

def blastpdb():
    """Blast search PDB based on command line arguments."""
    
    usage = """Usage: {0:s} [options] SEQUENCE
    
ProDy v{1:s} - Blast PDB

SEQUENCE can be a sequence string or a file in fasta format.  
""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    addOptions(parser)
    parser.add_option('-c', '--coverage', dest='coverage', type='float', 
                      default=90.0, metavar='FLOAT', 
                      help='percent coverage, default is %default%')
    parser.add_option('-d', '--dir', dest='folder', type='string',
                      default='', metavar='PATH', 
                      help=('if given, download PDB files to the folder'))
    parser.add_option('-i', '--identity', dest='identity', type='float', 
                      default=90.0, metavar='FLOAT', 
                      help='percent sequence identity, default is %default%')
    usage_examples="""
Blast search PDB for the first sequence in a fasta file:
    
  $ blastpdb.py seq.fasta -i 70

Blast search PDB for the sequence argument:

  $ blastpdb.py MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
"""
        
    opt, args = parser.parse_args()
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) < 1:
        parser.print_help()
        print "\nError: SEQUENCE missing\n"
        sys.exit(-1)
    if opt.silent:
        changeVerbosity('warning')
        
    seqfn = args[0]
    seq = seqfn
    if os.path.isfile(seq):
        seq = readFirstSequenceFasta(seqfn)
    if not seq.isalpha() or not seq.isupper():
        parser.print_help()
        print "\nError: {0:s} is not a valid sequence or a file\n".format(seq)
        sys.exit(-1)
        
    
    folder, identity, coverage, silent = (
        opt.folder, opt.identity, opt.coverage, opt.silent)
    assert 0 < identity < 100, 'identity must be between 0 and 100'
    assert 0 < coverage < 100, 'coverage must be between 0 and 100'
    if silent:
        changeVerbosity('warning')
    
    blast_results = blastPDB(seq)
    hits = blast_results.getHits(percent_identity=identity, 
                                 percent_coverage=coverage)
    
    #sort hits by decreasing percent identity
    hits2 = []
    for pdb in hits:
        hits2.append( (-hits[pdb]['percent_identity'], pdb) )
    hits2.sort()
    
    # print hits
    for identity,pdb in hits2:
        chain = hits[pdb]['chain_id']
        percent_identity = hits[pdb]['percent_identity']
        title = hits[pdb]['pdb_title']
        LOGGER.info(pdb + ' ' + chain + ' ' + ('%5.1f%%' % (percent_identity)) + ' ' + title)
    
    # download hits if --folder is given
    if opt.folder != '':
        LOGGER.info('Downloading hits to ' + opt.folder)
        pdblist = [ pdb for identity, pdb in hits2 ]
        pdblist2 = fetchPDB(pdblist, opt.folder)

def fetchpdb():
    """Fetch PDB files from PDB FTP server."""
    
    usage = """Usage: {0:s} [options] PDB PDB2 PDB3 ...  

ProDy v{1:s} - Fetch PDB

Download PDB files specified by their identifiers.
""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    addOptions(parser)
    parser.add_option('-d', '--dir', dest='folder', type='string',
                      default='.', metavar='PATH', 
                      help=('target directory saving downloaded PDB files'))
    parser.add_option('-f', '--file', dest="listfn", type='string', 
                      default='', metavar='FILE', 
                      help='file that contains PDB identifiers')
    usage_examples="""
Fetch PDB files for given identifiers:
    
  $ fetchpdb.py 1mkp 1p38
"""
    
    opt, args = parser.parse_args()
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) == 0 and opt.listfn == '':
        parser.print_help()
        print "\nError: PDB missing\n"
        sys.exit(-1)    
    if opt.silent:
        changeVerbosity('warning')
    
    folder, listfn = opt.folder, opt.listfn

    pdblist = []
    pdblist += args
    if opt.listfn != '':
        f = open(listfn)
        for line in f.xreadlines():
            line = line.strip()
            for s in line.split(','):
                for pdb in s.split():
                    if len(pdb) == 4: pdblist.append(pdb)
        f.close()
    
    pdblist2 = fetchPDB(pdblist, folder)
    
def pdbselect():
    """Write selected atoms from a PDB file in PDB format."""
    
    usage = """Usage: {0:s} [options] PDB SELECTION  

ProDy v{1:s} - PDB Select

Select atoms specified by SELECTION from PDB and write them in a file.
""".format(sys.argv[0], prody.__version__)
    parser = OptionParser(usage=usage)
    addOptions(parser)
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is "PDB_selected"'))
    usage_examples="""
Fetch PDB 1aar and write chain A carbon alpha atoms in a file:
        
  $ pdbselect.py 2bfu "backbone"
"""
    
    opt, args = parser.parse_args()
    if opt.examples:
        print 'Usage Examples:\n', usage_examples
        sys.exit(-1)
    if len(args) != 2:
        parser.print_help()
        print "\nError: PDB or SELECTION missing\n"
        sys.exit(-1)
    if opt.silent:
        changeVerbosity('warning')
        
    prefix = opt.prefix

    pdb = parsePDB(args[0])
    if prefix == '':
        prefix = pdb.getName() + '_selected'
    pdbselect = pdb.select(args[1])
    if pdbselect is None:
        LOGGER.warning('Selection "{0:s}" do not match any atoms.'
                       .format(args[1]))
        sys.exit(-1)
    LOGGER.info('Writing ' + prefix + '.pdb')
    writePDB(prefix + '.pdb', pdbselect)
    
