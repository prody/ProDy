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

"""Perform PCA/EDA calculations and output the results in plain text, NMD, and 
graphical formats.

Download example :download:`MDM2 trajectory files </doctest/mdm2.tar.gz>`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

from actions import *
from nmaoptions import *

def prody_pca(opt):
    """Perform PCA calculations based on command line arguments."""
    
    outdir = opt.outdir
    if not os.path.isdir(outdir):
        opt.subparser.error('{0:s} is not a valid path'.format(repr(outdir)))
        
    import prody
    LOGGER = prody.LOGGER
        
    coords = opt.coords
    prefix = opt.prefix
    nmodes, selstr = opt.nmodes, opt.select
    
    if os.path.splitext(coords)[1].lower() == '.dcd':     
        pdb = opt.psf or opt.pdb
        if pdb:
            if os.path.splitext(pdb)[1].lower() == '.psf':
                pdb = prody.parsePSF(pdb)
            else:
                pdb = prody.parsePDB(pdb)
        dcd = prody.DCDFile(opt.coords)
        if prefix == '_pca' or prefix == '_eda':
            prefix = dcd.getTitle() + prefix

        if len(dcd) < 2:
            opt.subparser.error("DCD file must contain multiple frames.")
        if pdb:
            if pdb.numAtoms() == dcd.numAtoms():
                dcd.setAtoms(pdb)
                select = dcd.select(selstr)
                LOGGER.info('{0:d} atoms are selected for calculations.'
                            .format(len(select)))
            else:
                select = pdb.select(selstr)
                if select.numAtoms() != dcd.numAtoms():
                    raise ValueError('number of selected atoms ({0:d}) does '
                                     'not match number of atoms in the DCD '
                                     'file ({1:d})'.format(select.numAtoms(),
                                                           dcd.numAtoms()))
                
        else:
            select = prody.AtomGroup()
            select.setCoords(dcd.getCoords())
        pca = prody.PCA(dcd.getTitle())
        if len(dcd) > 1000:
            pca.buildCovariance(dcd)
            pca.calcModes(nmodes)
        else:
            ens = dcd[:]
            ens.iterpose()
            pca.performSVD(ens)
    else:
        pdb = prody.parsePDB(opt.coords)
        if pdb.numCoordsets() < 2:
            opt.subparser.error("PDB file must contain multiple models.")
        if prefix == '_pca' or prefix == '_eda':
            prefix = pdb.getTitle() + prefix
        select = pdb.select(selstr)
        LOGGER.info('{0:d} atoms are selected for calculations.'
                    .format(len(select)))
        if select is None:
            opt.subparser.error('Selection {0:s} do not match any atoms.'
                                .format(repr(selstr)))
        LOGGER.info('{0:d} atoms will be used for PCA calculations.'
                    .format(len(select)))
        ensemble = prody.Ensemble(select)
        pca = prody.PCA(pdb.getTitle())
        ensemble.iterpose()
        pca.performSVD(ensemble)

    LOGGER.info('Writing numerical output.')
    if opt.npz:
        prody.saveModel(pca)
    prody.writeNMD(os.path.join(outdir, prefix + '.nmd'), pca[:nmodes], select)
    if opt.extend:
        if pdb:
            if opt.extend == 'all':
                extended = prody.extendModel(pca[:nmodes], select, pdb)        
            else:
                extended = prody.extendModel(pca[:nmodes], select, 
                                             select | pdb.bb)
            prody.writeNMD(os.path.join(outdir, prefix + '_extended_' + 
                           opt.extend + '.nmd'), *extended)
        else:
            prody.LOGGER.warn('Model could not be extended, provide a PDB or '
                              'PSF file.')
    outall = opt.all
    delim, ext, format = opt.delim, opt.ext, opt.numformat
    if outall or opt.eigen:
        prody.writeArray(os.path.join(outdir, prefix + '_evectors'+ext), 
                         pca.getArray(), delimiter=delim, format=format)
        prody.writeArray(os.path.join(outdir, prefix + '_evalues'+ext), 
                         pca.getEigvals(), delimiter=delim, format=format)
    if outall or opt.covar:
        prody.writeArray(os.path.join(outdir, prefix + '_covariance'+ext), 
                         pca.getCovariance(), delimiter=delim, format=format)
    if outall or opt.ccorr:
        prody.writeArray(os.path.join(outdir, prefix + '_cross-correlations' + 
                                              ext), prody.calcCrossCorr(pca), 
                         delimiter=delim, format=format)
    if outall or opt.sqflucts:
        prody.writeArray(os.path.join(outdir, prefix + '_sqfluct'+ext), 
                         prody.calcSqFlucts(pca), delimiter=delim, 
                         format=format)
    if outall or opt.proj:
        prody.writeArray(os.path.join(outdir, prefix + '_proj'+ext), 
                         prody.calcProjection(ensemble, pca), delimiter=delim, 
                         format=format)
          
    figall, cc, sf, sp = opt.figures, opt.cc, opt.sf, opt.sp

    if figall or cc or sf or sp: 
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
                prody.showCrossCorr(pca)
                plt.savefig(os.path.join(outdir, prefix + '_cc.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                prody.showSqFlucts(pca)
                plt.savefig(os.path.join(outdir, prefix + '_sf.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')                    
            if figall or sp:
                indices = []
                for item in sp.split():
                    try:
                        if '-' in item:
                            item = item.split('-')
                            if len(item) == 2:
                                indices.append(range(int(item[0])-1, 
                                               int(item[1])))
                        elif ',' in item:
                            indices.append([int(i)-1 for i in item.split(',')])
                        else:
                            indices.append(int(item)-1)
                    except:
                        pass
                for index in indices:
                        plt.figure(figsize=(width, height))
                        prody.showProjection(ensemble, pca[index])
                        if isinstance(index, int):
                            index = [index]
                        index = [str(i+1) for i in index]
                        plt.savefig(os.path.join(outdir, prefix + '_proj_' + 
                            '_'.join(index) + '.' + format),
                            dpi=dpi, format=format)
                        plt.close('all')
                        
def addCommand(commands):

    subparser = commands.add_parser('pca',
        help='perform principal component analysis calculations')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """This command performs PCA (or EDA) calculations for given multi-model \
PDB structure or DCD format trajectory file and outputs results in NMD \
format.  If a PDB identifier is given, structure file will be downloaded from \
the PDB FTP server.  DCD files may be accompanied with PDB or PSF files to \
enable atoms selections.

Fetch pdb 2k39, perform PCA calculations, and output NMD file:
    
  $ prody pca 2k39
    
Fetch pdb 2k39 and perform calculations for backbone of residues up to 71, \
and save all output and figure files:

  $ prody pca 2k39 --select "backbone and resnum < 71" -a -A
  
Perform EDA of MDM2 trajectory:
    
  $ prody eda mdm2.dcd
  
Perform EDA for backbone atoms:

  $ prody eda mdm2.dcd --pdb mdm2.pdb --select backbone"""
    )

    group = addNMAParameters(subparser)

    group = addNMAOutput(subparser)

    group.add_argument('-j', '--projection', dest='proj', action='store_true', 
        default=False, help='write projections onto PCs')

    group = addNMAOutputOptions(subparser, '_pca')

    group = addNMAFigures(subparser)

    group.add_argument('-J', '--projection-figure', dest='sp', type=str, 
        default='', metavar='STR',
        help=('save projections onto specified subspaces, e.g. '
              '"1,2" for projections onto PCs 1 and 2; '
              '"1,2 1,3" for projections onto PCs 1,2 and 1, 3; '
              '"1 1,2,3" for projections onto PCs 1 and 1, 2, 3'))

    group = addNMAFigureOptions(subparser)

    group = subparser.add_mutually_exclusive_group()
    group.add_argument('--psf', 
        help='PSF filename (must have same number of atoms as DCDs)')
    group.add_argument('--pdb', 
        help='PDB filename (must have same number of atoms as DCDs)')

    subparser.add_argument('coords', help='PDB or DCD filename')

    subparser.set_defaults(func=prody_pca)
    subparser.set_defaults(subparser=subparser)
    
    subparser = commands.add_parser('eda', parents=[subparser], 
        help='perform essential dynamics analysis calculations', add_help=False)

