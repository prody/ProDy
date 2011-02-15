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
from prody import *
from prody import ProDyLogger as LOGGER
import os.path

__all__ = ['anm', 'gnm', 'pca', 'alignmodels', 'biomolecule', 'blastpdb',
           'fetchpdb', 'pdbselect']

PY3K = sys.version_info[0] > 2

if PY3K:
    raise NotImplemented('ProDy scripts are not yet implemented for Python 3')
else:
    from optparse import OptionParser


def anm():
    """Perform ANM calculations based on command line arguments."""
    
    usage = """%prog [options] <pdb>

ProDy v{0:s} - Anisotropic Network Model""".format(prody.__version__) 
    parser = OptionParser(usage=usage)
    parser.add_option('-n', '--nmodes', dest='nmodes', type='int', 
                      default=20, metavar='INT', 
                      help=('number of non-zero eigenvalues/vectors to '
                            'calculate, default is %default'))    
    parser.add_option('-c', '--cutoff', dest='cutoff', type='float', 
                      default=15.0, metavar='FLOAT', 
                      help='cutoff distance, default is %default A')
    parser.add_option('-g', '--gamma', dest='gamma', type='float', 
                      default=1.0, metavar='FLOAT', 
                      help='spring constant, default is %default')
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is pdb_anm'))
    parser.add_option('-s', '--select', dest='select', type='string', 
                      default='calpha', metavar='STRING',
                      help='selection string, default is "%default"')
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
    parser.add_option('-e', '--examples', dest='examples', action='store_true', 
                      default=False, help='show usage examples')
    
    usage_examples="""
Fetch PDB 1p38 and run ANM using default parameters:
    
    $ anm.py 1p38
    
Fetch PDB 1p38 and run ANM using default parameters on all CA atoms of protein 
chain A:

    $ anm.py 1p38 --select "protein and name CA and chain A"
"""
    
    opt, args = parser.parse_args()
    
    if opt.examples:
        print usage_examples
        sys.exit(-1)
    
    if len(args) != 1:
        parser.print_help()
        print "\nError: <pdb> missing\n"
        sys.exit(-1)
        
    pdb = args[0]
    
    nmodes, select, cutoff, gamma, prefix, silent = (
        opt.nmodes, opt.select, opt.cutoff, opt.gamma, opt.prefix, opt.silent)

    if silent:
        changeVerbosity('warning')
    
    pdb = parsePDB(pdb)
    if prefix == '':
        pdb.getName() + '_anm'
    select = pdb.select(select)
    prody.ProDyLogger.info('{0:d} atoms will be used for ANM calculations.'
                           .format(len(select)))
    anm = ANM(pdb.getName())
    anm.buildHessian(select, cutoff, gamma)
    anm.calcModes(nmodes)
    
    writeArray(prefix + '_evalues.csv', anm.getArray(), delimiter=',')
    writeArray(prefix + '_evectors.csv', anm.getEigenvalues(), delimiter=',')
    writeNMD(prefix + '.nmd', anm, select)


def gnm():
    """Perform GNM calculations based on command line arguments."""

    usage = """%prog [options] <pdb>

ProDy v{0:s} - Gaussian Network Model""".format(prody.__version__)

    parser = OptionParser(usage=usage)
    parser.add_option('-n', '--nmodes', dest='nmodes', type='int', 
                      default=20, metavar='INT', 
                      help=('number of non-zero eigenvalues/vectors to '
                            'calculate, default is %default'))    
    parser.add_option('-c', '--cutoff', dest='cutoff', type='float', 
                      default=10.0, metavar='FLOAT', 
                      help='cutoff distance, default is %default A')
    parser.add_option('-g', '--gamma', dest='gamma', type='float', 
                      default=1.0, metavar='FLOAT', 
                      help='spring constant, default is %default')
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is pdb_gnm'))
    parser.add_option('-s', '--select', dest='select', type='string', 
                      default='calpha', metavar='STRING',
                      help='selection string, default is "%default"')
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
    parser.add_option('-e', '--examples', dest='examples', action='store_true', 
                      default=False, help='show usage examples')
    
    usage_examples="""
Fetch pdb 1p38 and run GNM using default parameters:
    
    $> gnm 1p38

Fetch pdb 1p38 and run GNM using default parameters on all CA atoms of protein 
chain A:
    
    $> gnm 1p38 --select "protein and name CA and chain A"
"""
    
    opt, args = parser.parse_args()
    
    if opt.examples:
        print usage_examples
        sys.exit(-1)
    
    if len(args) != 1:
        parser.print_help()
        print "\nError: <pdb> missing\n"
        sys.exit(-1)
        
    pdb = args[0]
    
    nmodes, select, cutoff, gamma, prefix, silent = (
        opt.nmodes, opt.select, opt.cutoff, opt.gamma, opt.prefix, opt.silent)
    if silent:
        changeVerbosity('warning')
    
    pdb = parsePDB(pdb)
    if prefix == '':
        pdb.getName() + '_gnm'
    select = pdb.select(select)
    prody.ProDyLogger.info('{0:d} atoms will be used for GNM calculations.'
                           .format(len(select)))
    gnm = GNM(pdb.getName())
    gnm.buildKirchhoff(select, cutoff, gamma)
    gnm.calcModes(nmodes)
    
    writeArray(prefix + '_evalues.csv', gnm.getArray(), delimiter=',')
    writeArray(prefix + '_evectors.csv', gnm.getEigenvalues(), delimiter=',')


def pca():
    """Perform PCA calculations based on command line arguments."""
    
    usage = """%prog [options] <pdb>

ProDy v{0:s} - Principal Component Analysis""".format(prody.__version__) 
    parser = OptionParser(usage=usage)
    parser.add_option('-n', '--nmodes', dest='nmodes', type='int', 
                      default=20, metavar='INT', 
                      help=('number of non-zero eigenvalues/vectors to '
                            'calculate, default is %default'))    
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is pdb_pca'))
    parser.add_option('-s', '--select', dest='select', type='string', 
                      default='calpha', metavar='STRING',
                      help='selection string, default is "%default"')
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
    parser.add_option('-e', '--examples', dest='examples', action='store_true', 
                      default=False, help='show usage examples')
    
    usage_examples="""
Fetch pdb 2k39 and run PCA for carbon alpha atoms of all residues:
    
    $ pca.py 2k39
    
Fetch pdb 2k39 and run PCA for residues with numbers less than 71:

    $ pca.py 2k39 --select "protein and name CA and resnum <= 71"
"""
    
    opt, args = parser.parse_args()
    
    if opt.examples:
        print usage_examples
        sys.exit(-1)
    
    if len(args) != 1:
        parser.print_help()
        print "\nError: <pdb> missing\n"
        sys.exit(-1)
        
    pdb = args[0]
    
    nmodes, select, prefix, silent = (
        opt.nmodes, opt.select, opt.prefix, opt.silent)

    if silent:
        changeVerbosity('warning')
    
    pdb = parsePDB(pdb)
    if prefix == '':
        pdb.getName() + '_pca'
    alignCoordsets(pdb, select)
    select = pdb.select(select)
    prody.ProDyLogger.info('{0:d} atoms will be used for PCA calculations.'
                           .format(len(select)))
    # gnm calculation
    pca = PCA(pdb.getName())
    pca.buildCovariance(select)
    pca.calcModes(nmodes)
    
    writeArray(prefix + '_evalues.csv', pca.getArray(), delimiter=',')
    writeArray(prefix + '_evectors.csv', pca.getEigenvalues(), delimiter=',')
    writeNMD(prefix + '.nmd', pca, select)
    
    
def alignmodels():
    """Align models in a PDB file based on command line arguments."""
    
    usage = """Usage: {0:s} [options] <pdb>  

ProDy v{1:s} - Align Models""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is pdb_aligned'))
    parser.add_option('-s', '--select', dest='select', type='string', 
                      default='calpha', metavar='STRING',
                      help='selection string, default is "%default"')
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
    parser.add_option('-e', '--examples', dest='examples', action='store_true', 
                      default=False, help='show usage examples')
    
    usage_examples="""
Fetch pdb 2k39 and align models:
    
    $ alignmodels.py 2k39
    
Fetch pdb 2k39 and align models using backbone of residues with number smaller
than 71:

    $ alignmodels.py 2k39 --select "backbone and resnum <= 71"
"""
        
    opt, args = parser.parse_args()
    
    if opt.examples:
        print usage_examples
        sys.exit(-1)
    
    if len(args) != 1:
        parser.print_help()
        print "\nError: <pdb> missing\n"
        sys.exit(-1)
        
    pdb = args[0]
    
    select, prefix, silent = (opt.select, opt.prefix, opt.silent)
    if silent:
        changeVerbosity('warning')
    pdb = parsePDB(pdb)
    if prefix == '':
        prefix = pdb.getName() + '_aligned'
    alignCoordsets(pdb, selstr=select)
    prody.ProDyLogger.info('{0:d} atoms will be used for alignment.'
                           .format(len(pdb.select(select))))
    rmsd = calcRMSD(pdb)
    LOGGER.info('Max RMSD: {0:0.2f} Mean RMSD: {1:0.2f}'
          .format(rmsd.max(), rmsd.mean()))
    outfn = prefix + '.pdb'
    LOGGER.info('Writing file: ' + outfn)
    writePDB(outfn, pdb)

def biomolecule():
    """Generate biomolecule coordinates based on command line arguments."""
    
    usage = """Usage: {0:s} [options] <pdb>  

ProDy v{1:s} - Generate Biomolecule Coordinates""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is pdb_biomol_'))
    parser.add_option('-b', '--biomol', dest='biomol', type='int', 
                      default=None, metavar='INT',
                      help='selection string, default is "%default"')
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
    parser.add_option('-e', '--examples', dest='examples', action='store_true', 
                      default=False, help='show usage examples')
    
    usage_examples="""
Fetch pdb 2bfu and generate the biomolecular assembly:
    
    $ biomolecule.py 2bfu
"""
        
    opt, args = parser.parse_args()
    
    if opt.examples:
        print usage_examples
        sys.exit(-1)
    
    if len(args) != 1:
        parser.print_help()
        print "\nError: <pdb> missing\n"
        sys.exit(-1)
        
    pdb = args[0]
    
    prefix, biomol, silent = (opt.prefix, opt.biomol, opt.silent)
    if silent:
        changeVerbosity('warning')
        
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
                outfn = '{0:s}_biomol_{1:d}_part_{2:d}.pdb'.format(prefix, i+1, j+1)
                LOGGER.info('Writing {0:s}'.format(outfn))
                writePDB(outfn, part)

def readFastaSequence(seqfn):
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
    
    usage = """Usage: {0:s} [options] <sequence.fasta>  

ProDy v{1:s} - Blast PDB""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    parser.add_option('-b', '--biomol', dest='biomol', type='int', 
                      default=None, metavar='INT',
                      help='selection string, default is "%default"')
    parser.add_option('-f', '--folder', dest='folder', type='string',
                      default='', metavar='PATH', 
                      help=('if given, download PDB files to the folder'))
    parser.add_option('-i', '--identity', dest='identity', type='float', 
                      default=90.0, metavar='FLOAT', 
                      help='percent sequence identity, default is %default%')
    parser.add_option('-c', '--coverage', dest='coverage', type='float', 
                      default=90.0, metavar='FLOAT', 
                      help='percent coverage, default is %default%')
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
        
    opt, args = parser.parse_args()
    
    if len(args) != 1:
        parser.print_help()
        print "\nError: <sequence.fasta> missing\n"
        sys.exit(-1)
        
    seqfn = args[0]
    
    folder, identity, coverage, silent = (
        opt.folder, opt.identity, opt.coverage, opt.silent)
    assert 0 < identity < 100, 'identity must be between 0 and 100'
    assert 0 < coverage < 100, 'coverage must be between 0 and 100'
    if silent:
        changeVerbosity('warning')

    seq = readFastaSequence(seqfn)
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
    
    usage = """Usage: {0:s} [options] <sequence.fasta>  

ProDy v{1:s} - Fetch PDB""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    parser.add_option('-f', '--folder', dest='folder', type='string',
                      default='.', metavar='PATH', 
                      help=('if given, download PDB files to the folder'))
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
    parser.add_option('', '--list', dest="listfn", type='string', 
                      default='', metavar='FILE', 
                      help='file that contains PDB identifiers')
    
    opt, args = parser.parse_args()
    
    folder, listfn, silent = (opt.folder, opt.listfn, opt.silent)
    if len(args) == 0 and opt.listfn == '':
        parser.print_help()
        sys.exit(-1)
    if silent:
        changeVerbosity('warning')
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
    
    usage = """Usage: {0:s} [options] <pdb> <selection>  

ProDy v{1:s} - PDB select""".format(sys.argv[0], prody.__version__)
        
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--prefix', dest='prefix', type='string', 
                      default='', metavar='STRING', 
                      help=('prefix for output files, default is pdb_selected'))
    parser.add_option('', '--silent', dest='silent', action='store_true', 
                      default=False, 
                      help='omit verbose information, default is %default')
    
    opt, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        sys.exit(-1)
        
    prefix, silent = (opt.prefix, opt.silent)
    if silent:
        changeVerbosity('warning')

    pdb = parsePDB(args[0])
    if prefix == '':
        prefix = pdb.getName() + '_selected'
    pdbselect = pdb.select(args[1])
    LOGGER.info('Writing ' + prefix + '.pdb')
    writePDB(prefix + '.pdb', pdbselect)
    
