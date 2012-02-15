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

"""This module contains functions which are used as command line programs."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import sys
import os.path

__all__ = ['prody_anm', 'prody_gnm', 'prody_pca', 'prody_align',
           'prody_blast', 'prody_biomol', 'prody_catdcd', 'prody_fetch',  
           'prody_select', ]

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


def prody_gnm(opt):
    """Perform GNM calculations based on command line arguments."""
    
    outdir = opt.outdir
    if not os.path.isdir(outdir):
        opt.subparser('{0:s} is not a valid path'.format(outdir))
        
    import numpy as np
    import prody
    LOGGER = prody.LOGGER

    pdb, prefix = opt.pdb, opt.prefix

    cutoff, gamma = opt.cutoff, opt.gamma, 
    nmodes, selstr, model = opt.nmodes, opt.select, opt.model
    
    pdb = prody.parsePDB(pdb, model=model)
    if prefix == '_gnm':
        prefix = pdb.getTitle() + '_gnm'

    select = pdb.select(selstr)
    if select is None:
        opt.subparser('Selection "{0:s}" do not match any atoms.'
                      .format(selstr))
    LOGGER.info('{0:d} atoms will be used for GNM calculations.'
                .format(len(select)))

    gnm = prody.GNM(pdb.getTitle())
    gnm.buildKirchhoff(select, cutoff, gamma)
    gnm.calcModes(nmodes)
    LOGGER.info('Writing numerical output.')
    if opt.npz:
        prody.saveModel(gnm)
    prody.writeNMD(os.path.join(outdir, prefix + '.nmd'), gnm, select)
    outall = opt.all
    delim, ext, format = opt.delim, opt.ext, opt.numformat
    
    if outall or opt.eigen:
        prody.writeArray(os.path.join(outdir, prefix + '_evectors'+ext), 
                         gnm.getArray(), delimiter=delim, format=format)
        prody.writeArray(os.path.join(outdir, prefix + '_evalues'+ext), 
                         gnm.getEigenvalues(), delimiter=delim, format=format)
    
    if outall or opt.beta:
        fout = prody.openFile(prefix + '_beta.txt', 'w', folder=outdir)
        fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4s} {0[3]:5s} {0[4]:5s}\n'
                       .format(['C', 'RES', '####', 'Exp.', 'The.']))
        for data in zip(select.getChids(), select.getResnames(), 
                        select.getResnums(), select.getBetas(), 
                        prody.calcTempFactors(gnm, select)):
            fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4d} {0[3]:5.2f} {0[4]:5.2f}\n'
                       .format(data))
        fout.close()
    if outall or opt.covar:
        prody.writeArray(os.path.join(outdir, prefix + '_covariance'+ext), 
                         gnm.getCovariance(), delimiter=delim, format=format)
    if outall or opt.ccorr:
        prody.writeArray(os.path.join(outdir, prefix + '_cross-correlations' 
                                                     + ext), 
                         prody.calcCrossCorr(gnm), delimiter=delim, 
                         format=format)
    if outall or opt.kirchhoff:
        prody.writeArray(os.path.join(outdir, prefix + '_kirchhoff'+ext), 
                         gnm.getKirchhoff(), delimiter=delim, format=format)
    if outall or opt.sqflucts:
        prody.writeArray(os.path.join(outdir, prefix + '_sqfluct'+ext), 
                         prody.calcSqFlucts(gnm), delimiter=delim, 
                         format=format)
          
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
                prody.showCrossCorr(gnm)
                plt.savefig(os.path.join(outdir, prefix + '_cc.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or cm:
                plt.figure(figsize=(width, height))
                prody.showContactMap(gnm)
                plt.savefig(os.path.join(outdir, prefix + '_cm.'+format), 
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                prody.showSqFlucts(gnm)
                plt.savefig(os.path.join(outdir, prefix + '_sf.'+format), 
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
                plt.savefig(os.path.join(outdir, prefix + '_bf.'+format), 
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
                        plt.savefig(os.path.join(outdir, prefix + '_mode_' + 
                            str(mode.getIndex()+1) + '.' + format), 
                            dpi=dpi, format=format)
                        plt.close('all')

def prody_pca(opt):
    """Perform PCA calculations based on command line arguments."""
    
    outdir = opt.outdir
    if not os.path.isdir(outdir):
        opt.subparser.error('{0:s} is not a valid path'.format(outdir))
        
    import prody
    LOGGER = prody.LOGGER
        
    coords = opt.coords
    prefix = opt.prefix
    nmodes, selstr = opt.nmodes, opt.select
    
    if os.path.splitext(coords)[1].lower() == '.dcd':     
        ag = opt.psf or opt.pdb
        if ag:
            if os.path.splitext(ag)[1].lower() == '.psf':
                ag = prody.parsePSF(ag)
            else:
                ag = prody.parsePDB(ag)
        dcd = prody.DCDFile(opt.coords)
        if len(dcd) < 2:
            opt.subparser("DCD file must contain multiple frames.")
        if ag:
            dcd.setAtomGroup(ag)
            select = dcd.select(selstr)
            LOGGER.info('{0:d} atoms are selected for calculations.'
                        .format(len(select)))
        else:
            select = prody.AtomGroup()
            select.setCoords(dcd.getCoords())
        pca = prody.PCA(dcd.getTitle())
        if len(dcd) > 1000:
            pca.buildCovariance(dcd)
            pca.calcModes(dcd)
        else:
            pca.performSVD(dcd[:])
    else:
        pdb = prody.parsePDB(opt.coords)
        if pdb.numCoordsets() < 2:
            opt.subparser("PDB file must contain multiple models.")
        if prefix == '_pca':
            prefix = pdb.getTitle() + '_pca'
        select = pdb.select(selstr)
        LOGGER.info('{0:d} atoms are selected for calculations.'
                    .format(len(select)))
        if select is None:
            opt.subparser('Selection "{0:s}" do not match any atoms.'
                          .format(selstr))
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

    outall = opt.all
    delim, ext, format = opt.delim, opt.ext, opt.numformat
    if outall or opt.eigen:
        prody.writeArray(os.path.join(outdir, prefix + '_evectors'+ext), 
                         pca.getArray(), delimiter=delim, format=format)
        prody.writeArray(os.path.join(outdir, prefix + '_evalues'+ext), 
                         pca.getEigenvalues(), delimiter=delim, format=format)
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
    ag = opt.psf or opt.pdb
    if ag:
        if os.path.splitext(ag)[1].lower() == '.psf':
            ag = prody.parsePSF(ag)
        else:
            ag = prody.parsePDB(ag)
    
    dcd = opt.dcd
    traj = prody.Trajectory(dcd.pop(0))
    while dcd:
        traj.addFile(dcd.pop(0))
    if ag:
        traj.setAtomGroup(ag)
        select = traj.select(opt.select)
        LOGGER.info('{0:d} atoms are selected for writing output.'
                    .format(len(select)))

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
        out.write(frame._getCoords(), frame.getUnitcell())
        count += 1
    traj.close()
    out.close()
    LOGGER.info("{0:d} frames are written into '{1:s}'."
                .format(count, opt.output))
    
    

def prody_align(opt):
    """Align models in a PDB file or a PDB file onto others."""
            
    import prody
    LOGGER = prody.LOGGER

    args = opt.pdb
    if len(args) == 1:
        pdb = args[0]
        LOGGER.info('Aligning multiple models in: ' + pdb)
        selstr, prefix, model = opt.select, opt.prefix, opt.model
        pdb = prody.parsePDB(pdb)
        pdbselect = pdb.select(selstr)
        if pdbselect is None:
            LOGGER.warning('Selection "{0:s}" do not match any atoms.'
                           .format(selstr))
            sys.exit(-1)
        LOGGER.info('{0:d} atoms will be used for alignment.'
                               .format(len(pdbselect)))
        pdb.setACSIndex(model-1)
        prody.alignCoordsets(pdb, selstr=selstr)
        rmsd = prody.calcRMSD(pdb)
        LOGGER.info('Max RMSD: {0:0.2f} Mean RMSD: {1:0.2f}'
              .format(rmsd.max(), rmsd.mean()))
        if prefix == '':
            prefix = pdb.getTitle() + '_aligned'
        outfn = prefix + '.pdb'
        LOGGER.info('Writing file: ' + outfn)
        prody.writePDB(outfn, pdb)
    else:
        reffn = args.pop(0)
        LOGGER.info('Aligning structures onto: ' + reffn)
        ref = prody.parsePDB(reffn)
        for arg in args:
            if arg == reffn:
                continue
            if '_aligned.pdb' in arg:
                continue
            pdb = prody.parsePDB(arg)
            if prody.matchAlign(pdb, ref):
                outfn = pdb.getTitle() + '_aligned.pdb'
                LOGGER.info('Writing file: ' + outfn)
                prody.writePDB(outfn, pdb)
            else:
                LOGGER.warning('Failed to align ' + arg)

def prody_biomol(opt):
    """Generate biomolecule coordinates based on command line arguments."""
        
    import prody
    LOGGER = prody.LOGGER
    prefix, biomol = opt.prefix, opt.biomol
    pdb, header = prody.parsePDB(opt.pdb, header=True)
    if not prefix:
        prefix = pdb.getTitle()
        
    biomols = prody.buildBiomolecules(header, pdb, biomol=biomol)
    if not isinstance(biomols, list):
        biomols = [biomols]
    
    for i, biomol in enumerate(biomols):
        if isinstance(biomol, prody.Atomic):
            outfn = '{0:s}_biomol_{1:d}.pdb'.format(prefix, i+1)
            LOGGER.info('Writing {0:s}'.format(outfn))
            prody.writePDB(outfn, biomol)
        elif isinstance(biomol, tuple):
            for j, part in enumerate(biomol):
                outfn = ('{0:s}_biomol_{1:d}_part_{2:d}.pdb'
                         .format(prefix, i+1, j+1))
                LOGGER.info('Writing {0:s}'.format(outfn))
                prody.writePDB(outfn, part)

def readFirstSequenceFasta(filename):
    """Return first sequence from a file."""
    
    fasta = open(filename)
    seq = []
    title = '' 
    first = True
    for line in fasta:
        if line[0] == '>': 
            if first:
                title = line[1:].strip()
                first = False
            else:    
                break
        else:
            seq.append( line.strip() )
    fasta.close()
    return title, ''.join(seq)

def prody_blast(opt):
    """Blast search PDB based on command line arguments."""
    
    import prody
    LOGGER = prody.LOGGER
    seq = opt.seq
    title = None
    if os.path.isfile(seq):
        title, seq = readFirstSequenceFasta(seq)
        LOGGER.info("First sequence ({0:s}) is parsed from '{1:s}'."
                    .format(title, seq))
    if not seq.isalpha() or not seq.isupper():
        opt.subparser("{0:s} is not a valid sequence or a file".format(seq))
        
    folder, identity, coverage = opt.folder, opt.identity, opt.coverage
    if not 0 < identity < 100: 
        opt.subparser.error('identity must be between 0 and 100')
    if not 0 < coverage < 100:
        opt.subparser.error('overlap must be between 0 and 100')
    
    blast_results = prody.blastPDB(seq)
    hits = blast_results.getHits(percent_identity=identity, 
                                 percent_coverage=coverage)
    
    #sort hits by decreasing percent identity
    hits2 = []
    for pdb in hits:
        hits2.append( (-hits[pdb]['percent_identity'], pdb) )
    hits2.sort()
    
    for identity, pdb in hits2:
        chain = hits[pdb]['chain_id']
        percent_identity = hits[pdb]['percent_identity']
        title = hits[pdb]['title']
        print(pdb + ' ' + chain + ' ' + ('%5.1f%%' % (percent_identity)) + 
              ' ' + title)
    
    # download hits if --folder is given
    if opt.folder:
        LOGGER.info('Downloading hits to ' + opt.folder)
        pdblist = [ pdb for identity, pdb in hits2 ]
        pdblist2 = prody.fetchPDB(pdblist, opt.folder, 
                                  compressed=opt.gzip, copy=True)

def prody_fetch(opt):
    """Fetch PDB files from PDB FTP server."""
    
    import prody
    pdblist = opt.pdb
    listfn = opt.listfn
    if listfn:
        if os.path.isfile(listfn):
            inp = prody.openFile(listfn)
            for line in inp:
                line = line.strip()
                for s in line.split(','):
                    for pdb in s.split():
                        if len(pdb) == 4: 
                            pdblist.append(pdb)
            inp.close()
        else:    
            opt.subparser.error("No such file: '{0:s}'".format(listfn))
    prody.fetchPDB(pdblist, opt.folder, compressed=opt.gzip, copy=True)

def prody_select(opt):
    """Write selected atoms from a PDB file in PDB format."""

    import prody
    LOGGER = prody.LOGGER
    pdb = prody.parsePDB(opt.pdb)
    prefix = opt.output
    if not prefix:
        prefix = pdb.getTitle() + '_selected'
    pdbselect = pdb.select(opt.selstr)
    if pdbselect is None:
        opt.subparser('Selection "{0:s}" do not match any atoms.'
                      .format(opt.selstr))
    LOGGER.info('Writing ' + prefix + '.pdb')
    prody.writePDB(prefix + '.pdb', pdbselect)


    
