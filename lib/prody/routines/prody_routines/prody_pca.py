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


from ..actions import *
from .nmaoptions import *
from . import nmaoptions

DEFAULTS = {}
HELPTEXT = {} 
for key, txt, val in [
    ('aligned', 'trajectory is already aligned', False),
    ('outproj', 'write projections onto PCs', False),
    ('figproj', 'save projections onto specified subspaces, e.g. '
                '"1,2" for projections onto PCs 1 and 2; '
                '"1,2 1,3" for projections onto PCs 1,2 and 1, 3; '
                '"1 1,2,3" for projections onto PCs 1 and 1, 2, 3', ''),]:
    
    DEFAULTS[key] = val
    HELPTEXT[key] = txt

DEFAULTS.update(nmaoptions.DEFAULTS)
HELPTEXT.update(nmaoptions.HELPTEXT)

DEFAULTS['prefix'] = '_pca'

__all__ = ['prody_pca']

def prody_pca(coords, **kwargs):
    """Perform PCA calculations for PDB or DCD format *coords* file.
    
    """
    
    for key in DEFAULTS:
        if not key in kwargs:
            kwargs[key] = DEFAULTS[key]
    
    from os.path import isdir, splitext, join
    outdir = kwargs['outdir']
    if not isdir(outdir):
        raise IOError('{0:s} is not a valid path'.format(repr(outdir)))
        
    import prody
    LOGGER = prody.LOGGER
        
    prefix = kwargs['prefix']
    nmodes = kwargs['nmodes']
    selstr = kwargs['select']
    
    ext = splitext(coords)[1].lower() 
    if ext == '.gz':
        ext = splitext(coords[:-3])[1].lower()
    
    if ext == '.dcd':     
        pdb = kwargs.get('psf') or kwargs.get('pdb')
        if pdb:
            if splitext(pdb)[1].lower() == '.psf':
                pdb = prody.parsePSF(pdb)
            else:
                pdb = prody.parsePDB(pdb)
        dcd = prody.DCDFile(coords)
        if prefix == '_pca' or prefix == '_eda':
            prefix = dcd.getTitle() + prefix

        if len(dcd) < 2:
            raise ValueError('DCD file must have multiple frames')
        if pdb:
            if pdb.numAtoms() == dcd.numAtoms():
                select = pdb.select(selstr)
                dcd.setAtoms(select)
                LOGGER.info('{0:d} atoms are selected for calculations.'
                            .format(len(select)))
            else:
                select = pdb.select(selstr)
                if select.numAtoms() != dcd.numAtoms():
                    raise ValueError('number of selected atoms ({0:d}) does '
                                     'not match number of atoms in the DCD '
                                     'file ({1:d})'.format(select.numAtoms(),
                                                           dcd.numAtoms()))
                if pdb.numCoordsets():
                    dcd.setCoords(select.getCoords())
                    
        else:
            select = prody.AtomGroup()
            select.setCoords(dcd.getCoords())
        pca = prody.PCA(dcd.getTitle())
        if len(dcd) > 1000:
            pca.buildCovariance(dcd, aligned=kwargs['aligned'])
            pca.calcModes(nmodes)
            ensemble = dcd
        else:
            ensemble = dcd[:]
            if not kwargs['aligned']:
                ensemble.iterpose()
            pca.performSVD(ensemble)
    
    else:
        pdb = prody.parsePDB(coords)
        if pdb.numCoordsets() < 2:
            raise ValueError('PDB file must contain multiple models')
            
        if prefix == '_pca' or prefix == '_eda':
            prefix = pdb.getTitle() + prefix
            
        select = pdb.select(selstr)
        LOGGER.info('{0:d} atoms are selected for calculations.'
                    .format(len(select)))
        if select is None:
            raise ValueError('selection {0:s} do not match any atoms'
                                .format(repr(selstr)))
        LOGGER.info('{0:d} atoms will be used for PCA calculations.'
                    .format(len(select)))
        ensemble = prody.Ensemble(select)
        pca = prody.PCA(pdb.getTitle())
        if not kwargs['aligned']:
            ensemble.iterpose()
        pca.performSVD(ensemble)
    
        
    LOGGER.info('Writing numerical output.')
    if kwargs['outnpz']:
        prody.saveModel(pca, join(outdir, prefix))
        
    prody.writeNMD(join(outdir, prefix + '.nmd'), pca[:nmodes], select)
    
    extend = kwargs['extend']
    if extend:
        if pdb:
            if extend == 'all':
                extended = prody.extendModel(pca[:nmodes], select, pdb)        
            else:
                extended = prody.extendModel(pca[:nmodes], select, 
                                             select | pdb.bb)
            prody.writeNMD(join(outdir, prefix + '_extended_' + 
                           extend + '.nmd'), *extended)
        else:
            prody.LOGGER.warn('Model could not be extended, provide a PDB or '
                              'PSF file.')
    outall = kwargs['outall']
    delim = kwargs['numdelim']
    ext = kwargs['numext']
    format = kwargs['numformat']
    
    if outall or kwargs['outeig']:
        prody.writeArray(join(outdir, prefix + '_evectors'+ext), 
                         pca.getArray(), delimiter=delim, format=format)
        prody.writeArray(join(outdir, prefix + '_evalues'+ext), 
                         pca.getEigvals(), delimiter=delim, format=format)
    if outall or kwargs['outcov']:
        prody.writeArray(join(outdir, prefix + '_covariance'+ext), 
                         pca.getCovariance(), delimiter=delim, format=format)
    if outall or kwargs['outcc']:
        prody.writeArray(join(outdir, prefix + '_cross-correlations' + 
                                              ext), prody.calcCrossCorr(pca), 
                         delimiter=delim, format=format)
    if outall or kwargs['outsf']:
        prody.writeArray(join(outdir, prefix + '_sqfluct'+ext), 
                         prody.calcSqFlucts(pca), delimiter=delim, 
                         format=format)
    if outall or kwargs['outproj']:
        prody.writeArray(join(outdir, prefix + '_proj'+ext), 
                         prody.calcProjection(ensemble, pca), delimiter=delim, 
                         format=format)
          
    figall = kwargs['figall']
    cc = kwargs['figcc']
    sf = kwargs['figsf']
    sp = kwargs['figproj']

    if figall or cc or sf or sp: 
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warning('Matplotlib could not be imported. '
                           'Figures are not saved.')
        else:
            LOGGER.info('Saving graphical output.')
            format = kwargs['figformat']
            width = kwargs['figwidth']
            height = kwargs['figheight']
            dpi = kwargs['figdpi']
            
            format = format.lower()
            if figall or cc:
                plt.figure(figsize=(width, height))
                prody.showCrossCorr(pca)
                plt.savefig(join(outdir, prefix + '_cc.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                prody.showSqFlucts(pca)
                plt.savefig(join(outdir, prefix + '_sf.'+format), 
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
                        plt.savefig(join(outdir, prefix + '_proj_' + 
                            '_'.join(index) + '.' + format),
                            dpi=dpi, format=format)
                        plt.close('all')
                        

_ = list(HELPTEXT)
_.sort()
for key in _:
    
    prody_pca.__doc__ += """
    :arg {0}: {1}, default is ``{2!r}``""".format(key, HELPTEXT[key], 
                                                  DEFAULTS[key])


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

  $ prody eda mdm2.dcd --pdb mdm2.pdb --select backbone""",
    test_examples=[0, 1])

    group = addNMAParameters(subparser)

    group = addNMAOutput(subparser)

    group.add_argument('-j', '--projection', dest='outproj', 
        action='store_true', 
        default=DEFAULTS['outproj'], help=HELPTEXT['outproj'])

    group = addNMAOutputOptions(subparser, '_pca')

    group = addNMAFigures(subparser)

    group.add_argument('-J', '--projection-figure', dest='figproj', type=str, 
        default=DEFAULTS['figproj'], metavar='STR',
        help=HELPTEXT['figproj'])

    group = addNMAFigureOptions(subparser)

    group = subparser.add_mutually_exclusive_group()
    group.add_argument('--psf', help='PSF filename')
    group.add_argument('--pdb', help='PDB filename')
    subparser.add_argument('--aligned', dest='aligned', action='store_true', 
        default=DEFAULTS['aligned'], help=HELPTEXT['aligned'])

    subparser.add_argument('dcd', help='file in DCD or PDB format')

    subparser.set_defaults(func=lambda ns: prody_pca(ns.dcd, **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
