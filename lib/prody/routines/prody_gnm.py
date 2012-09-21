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

"""Perform GNM calculations and output the results in plain text NMD, and 
graphical formats."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'


from .actions import *
from .nmaoptions import *
from . import nmaoptions

__all__ = ['prody_gnm']

DEFAULTS = {}
HELPTEXT = {} 
for key, txt, val in [
    ('model', 'index of model that will be used in the calculations', 1),
    ('cutoff', 'cutoff distance (A)', 10.),
    ('gamma', 'spring constant', 1.),
    
    ('outbeta', 'write beta-factors calculated from GNM modes', False),
    ('kirchhoff', 'write Kirchhoff matrix', False),
    ('figcmap', 'save contact map (Kirchhoff matrix) figure', False),
    ('figbeta', 'save beta-factors figure', False),
    ('figmode', 'save mode shape figures for specified modes, '
                 'e.g. "1-3 5" for modes 1, 2, 3 and 5', '')]:
    
    DEFAULTS[key] = val
    HELPTEXT[key] = txt

DEFAULTS.update(nmaoptions.DEFAULTS)
HELPTEXT.update(nmaoptions.HELPTEXT)

DEFAULTS['prefix'] = '_gnm'

def prody_gnm(pdb, **kwargs):
    """Perform GNM calculations for *pdb*.
    
    """
    
    for key in DEFAULTS:
        if not key in kwargs:
            kwargs[key] = DEFAULTS[key]
    
    from os.path import isdir, splitext, join
    outdir = kwargs['outdir']
    if not isdir(outdir):
        raise IOError('{0:s} is not a valid path'.format(repr(outdir)))
        
    import numpy as np
    import prody
    LOGGER = prody.LOGGER

    selstr = kwargs['select']
    prefix = kwargs['prefix']
    cutoff = kwargs['cutoff']
    gamma = kwargs['gamma'] 
    nmodes = kwargs['nmodes']
    selstr = kwargs['select']
    model = kwargs['model']
    
    pdb = prody.parsePDB(pdb, model=model)
    if prefix == '_gnm':
        prefix = pdb.getTitle() + '_gnm'

    select = pdb.select(selstr)
    if select is None:
        raise ValueError('selection {0:s} do not match any atoms'
                          .format(repr(selstr)))
    LOGGER.info('{0:d} atoms will be used for GNM calculations.'
                .format(len(select)))

    gnm = prody.GNM(pdb.getTitle())
    gnm.buildKirchhoff(select, cutoff, gamma)
    gnm.calcModes(nmodes)
    
    LOGGER.info('Writing numerical output.')
    
    if kwargs['outnpz']:
        prody.saveModel(gnm, join(outdir, prefix))
    
    prody.writeNMD(join(outdir, prefix + '.nmd'), gnm, select)
    
    extend = kwargs['extend']
    if extend:
        if extend == 'all':
            extended = prody.extendModel(gnm, select, pdb)        
        else:
            extended = prody.extendModel(gnm, select, select | pdb.bb)
        prody.writeNMD(join(outdir, prefix + '_extended_' + 
                       extend + '.nmd'), *extended)
    
    outall = kwargs['outall']
    delim = kwargs['numdelim']
    ext = kwargs['numext']
    format = kwargs['numformat']
    
    if outall or kwargs['outeig']:
        prody.writeArray(join(outdir, prefix + '_evectors'+ext), 
                         gnm.getArray(), delimiter=delim, format=format)
        prody.writeArray(join(outdir, prefix + '_evalues'+ext), 
                         gnm.getEigvals(), delimiter=delim, format=format)
    
    if outall or kwargs['outbeta']:
        from prody.utilities import openFile
        fout = openFile(prefix + '_beta.txt', 'w', folder=outdir)
        fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4s} {0[3]:5s} {0[4]:5s}\n'
                       .format(['C', 'RES', '####', 'Exp.', 'The.']))
        for data in zip(select.getChids(), select.getResnames(), 
                        select.getResnums(), select.getBetas(), 
                        prody.calcTempFactors(gnm, select)):
            fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4d} {0[3]:5.2f} {0[4]:5.2f}\n'
                       .format(data))
        fout.close()
        
    if outall or kwargs['outcov']:
        prody.writeArray(join(outdir, prefix + '_covariance'+ext), 
                         gnm.getCovariance(), delimiter=delim, format=format)
    
    if outall or kwargs['outcc']:
        prody.writeArray(join(outdir, prefix + '_cross-correlations' 
                                                     + ext), 
                         prody.calcCrossCorr(gnm), delimiter=delim, 
                         format=format)
    
    if outall or kwargs['kirchhoff']:
        prody.writeArray(join(outdir, prefix + '_kirchhoff'+ext), 
                         gnm.getKirchhoff(), delimiter=delim, format=format)
    
    if outall or kwargs['outsf']:
        prody.writeArray(join(outdir, prefix + '_sqfluct'+ext), 
                         prody.calcSqFlucts(gnm), delimiter=delim, 
                         format=format)
          
    figall = kwargs['figall']
    cc = kwargs['figcc']
    sf = kwargs['figsf']
    bf = kwargs['figbeta']
    cm = kwargs['figcmap']
    modes = kwargs['figmode']
    
    if figall or cc or sf or bf or cm or modes: 
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
                prody.showCrossCorr(gnm)
                plt.savefig(join(outdir, prefix + '_cc.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            
            if figall or cm:
                plt.figure(figsize=(width, height))
                prody.showContactMap(gnm)
                plt.savefig(join(outdir, prefix + '_cm.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            
            if figall or sf:
                plt.figure(figsize=(width, height))
                prody.showSqFlucts(gnm)
                plt.savefig(join(outdir, prefix + '_sf.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            
            if figall or bf:
                plt.figure(figsize=(width, height))
                bexp = select.getBetas()
                bcal = prody.calcTempFactors(gnm, select)
                plt.plot(bexp, label='Experimental')
                plt.plot(bcal, label=('Theoretical (corr coef = {0:.2f})'
                                        .format(np.corrcoef(bcal, bexp)[0,1])))
                plt.legend(prop={'size': 10})
                plt.xlabel('Node index')
                plt.ylabel('Experimental B-factors')
                plt.title(pdb.getTitle() + ' B-factors')
                plt.savefig(join(outdir, prefix + '_bf.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            
            if modes: 
                indices = []
                items = modes.split()
                items = sum([item.split(',') for item in items], [])
                for item in items:
                    try:
                        item = item.split('-')
                        if len(item) == 1:
                            indices.append(int(item[0])-1)
                        elif len(item) == 2:
                            indices.extend(range(int(item[0])-1, int(item[1])))
                    except:
                        pass
                for index in indices:
                    try:
                        mode = gnm[index]
                    except: 
                        pass
                    else:
                        plt.figure(figsize=(width, height))
                        prody.showMode(mode)
                        plt.grid()
                        plt.savefig(join(outdir, prefix + '_mode_' + 
                            str(mode.getIndex()+1) + '.' + format), 
                            dpi=dpi, format=format)
                        plt.close('all')
           

_ = list(HELPTEXT)
_.sort()      
for key in _:
    
    prody_gnm.__doc__ += """
    :arg {0}: {1}, default is ``{2!r}``""".format(key, HELPTEXT[key], 
                                                  DEFAULTS[key])
                                                  
                        
def addCommand(commands):

    subparser = commands.add_parser('gnm', 
        help='perform Gaussian network model calculations')
        
    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """This command performs GNM calculations for given PDB structure and \
outputs results in NMD format. If an identifier is passed, structure file \
will be downloaded from the PDB FTP server.

Fetch PDB 1p38, run GNM calculations using default parameters, and results:
    
  $ prody gnm 1p38
    
Fetch PDB 1aar, run GNM calculations with cutoff distance 7 angstrom for \
chain A carbon alpha atoms with residue numbers less than 70, and \
save all of the graphical output files:

  $ prody gnm 1aar -c 7 -s "calpha and chain A and resnum < 70" -A""",
    test_examples=[0, 1])

    group = addNMAParameters(subparser)

    group.add_argument('-c', '--cutoff', dest='cutoff', type=float, 
        default=DEFAULTS['cutoff'], metavar='FLOAT', 
        help=HELPTEXT['cutoff'] + ' (default: %(default)s)')
    
    group.add_argument('-g', '--gamma', dest='gamma', type=float, 
        default=DEFAULTS['gamma'], metavar='FLOAT', 
        help=HELPTEXT['gamma'] + ' (default: %(default)s)')
    
    group.add_argument('-m', '--model', dest='model', type=int, 
        metavar='INT', default=DEFAULTS['model'], help=HELPTEXT['model']) 
                        
    group = addNMAOutput(subparser)
    
    group.add_argument('-b', '--beta-factors', dest='outbeta', 
        action='store_true', default=DEFAULTS['outbeta'], 
        help=HELPTEXT['outbeta'])
    
    group.add_argument('-k', '--kirchhoff', dest='kirchhoff', 
        action='store_true', 
        default=DEFAULTS['kirchhoff'], help=HELPTEXT['kirchhoff'])

    group = addNMAOutputOptions(subparser, '_gnm')

    group = addNMAFigures(subparser)
    
    group.add_argument('-B', '--beta-factors-figure', dest='figbeta', 
        action='store_true', 
        default=DEFAULTS['figbeta'], help=HELPTEXT['figbeta'])
        
    group.add_argument('-K', '--contact-map', dest='figcmap', 
        action='store_true', 
        default=DEFAULTS['figcmap'], help=HELPTEXT['figcmap'])
        
    group.add_argument('-M', '--mode-shape-figure', dest='figmode', type=str, 
        default=DEFAULTS['figmode'], help=HELPTEXT['figmode'], metavar='STR')

    group = addNMAFigureOptions(subparser)

    subparser.add_argument('pdb', help='PDB identifier or filename')
        
    subparser.set_defaults(func=lambda ns: prody_gnm(ns.__dict__.pop('pdb'), 
                                                     **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
