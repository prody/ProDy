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
import textwrap

import argparse

from routines import *

__all__ = ['main']

DEFAULT_FIGURE_FORMAT = 'pdf'
DEFAULT_FILE_EXT = '.txt'
DEFAULT_NUMBER_FORMAT = '%12g'
DEFAULT_DELIMITER = ' '

class Quiet(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        import prody
        prody.setVerbosity('warning')


class UsageExample(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        tw = textwrap.TextWrapper()
        for line in namespace.usage_example.splitlines():
            print("\n".join(tw.wrap(line)))
        parser.exit()

class ProDyCitation(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print("Bakan A, Meireles LM, Bahar I "
              "ProDy: Protein Dynamics Inferred from Theory and Experiments "
              "Bioinformatics 2011 27(11):1575-1577.")
        parser.exit()

class ProDyVersion(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        import prody
        print("ProDy " + prody.__version__)
        parser.exit()


###############################################################################
# add common arguments to subparsers
###############################################################################

def addNMAParameters(parser):
    parser = parser.add_argument_group('parameters')
    parser.add_argument('-n', '--number-of-modes', dest='nmodes', type=int, 
        default=10, metavar='INT', 
        help="number of non-zero eigenvalues/vectors to "
              "calculate (default: %(default)s)")    
    parser.add_argument('-s', '--select', dest='select', type=str, 
        default="protein and name CA or nucleic and name P C4' C2", 
        metavar='SELSTR', help='atom selection (default: "%(default)s")')
    return parser

def addNMAOutput(parser):
    parser = parser.add_argument_group('output')
    parser.add_argument('-a', '--all-output', dest='all', action='store_true', 
                        default=False, help='write all outputs')
    parser.add_argument('-o', '--output-dir', dest='outdir', type=str, 
                        default='.', metavar='PATH', 
                        help='output directory (default: "%(default)s")')
    parser.add_argument('-e', '--eigenvs', dest='eigen', action='store_true', 
                        default=False, help='write eigenvectors/values')
    parser.add_argument('-r', '--cross-correlations', dest='ccorr', 
                        action='store_true', 
                        default=False, help='write cross-correlations')
    parser.add_argument('-q', '--square-fluctuations', dest='sqflucts', 
                        action='store_true', 
                        default=False, help='write square-fluctuations')
    parser.add_argument('-v', '--covariance', dest='covar', action='store_true', 
                        default=False, help='write covariance matrix')
    parser.add_argument('-z', '--npz', dest='npz', action='store_true', 
                        default=False, help='write compressed ProDy data file')
    return parser

def addNMAOutputOptions(parser, prefix):
    parser = parser.add_argument_group('output options')
    parser.add_argument('-p', '--file-prefix', dest='prefix', type=str, 
                        default=prefix, metavar='STR', 
                        help=('prefix for output files (default: '
                              'pdb%(default)s)'))
    parser.add_argument('-f', '--number-format', dest='numformat', type=str, 
                       default=DEFAULT_NUMBER_FORMAT, metavar='STR',
                       help=('delimiter (default: "%(default)s")'))
    parser.add_argument('-d', '--delimiter', dest='delim', type=str, 
                        default=DEFAULT_DELIMITER, metavar='STR',
                        help='delimiter (default: "%(default)s")')
    parser.add_argument('-x', '--extension', dest='ext', type=str, 
                       default=DEFAULT_FILE_EXT, metavar='STR',
                       help=("file extension (default: %(default)s)"))
    return parser

def addNMAFigures(parser):
    parser = parser.add_argument_group('figures')
    parser.add_argument('-A', '--all-figures', dest='figures', 
                        action='store_true', 
                        default=False, help='save all figures')
    parser.add_argument('-R', '--cross-correlations-figure', dest='cc', 
                        action='store_true', default=False, 
                        help='save cross-correlations')
    parser.add_argument('-Q', '--square-fluctuations-figure', dest='sf', 
                        action='store_true', default=False, 
                        help='save square-fluctuations')
    return parser

def addNMAFigureOptions(parser):
    parser = parser.add_argument_group('figure options')
    parser.add_argument('-F', '--figure-format', dest='figformat', type=str, 
                        default=DEFAULT_FIGURE_FORMAT, metavar='STR',
                        help="figure format, one of eps, pdf, png, ps, raw, "
                              "rgba, svg, svgz (default: %(default)s)",
                        choices="eps pdf png ps raw rgba svg svgz".split())
    parser.add_argument('-D', '--resolution', dest='dpi', type=int, 
                        default=300, metavar='INT', 
                        help="figure resolution (dpi) (default: %(default)s)")
    parser.add_argument('-W', '--width', dest='width', type=float, 
                        default=8.0, metavar='FLOAT', 
                        help="figure width (inch) (default: %(default)s)")
    parser.add_argument('-H', '--height', dest='height', type=float, 
                        default=6.0, metavar='FLOAT', 
                        help="figure height (inch) (default: %(default)s)")
    return parser

###############################################################################
# prody
###############################################################################

parser = argparse.ArgumentParser(
    description="ProDy: A Python Package for Protein Dynamics Analysis",
    epilog="See 'prody <command> -h' for more information on a specific "
           "command."
    )

parser.add_argument('-c', '--cite', help="print citation info and exit",
    action=ProDyCitation, nargs=0)

parser.add_argument('-v', '--version', help="print ProDy version and exit",
    action=ProDyVersion, nargs=0)

commands = parser.add_subparsers(
    title='subcommands')

###############################################################################
# anm
###############################################################################

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

###############################################################################
# gnm
###############################################################################

subparser = commands.add_parser('gnm', 
    help='perform Gaussian network model calculations')
    
subparser.add_argument('--quiet', help="suppress info messages to stderr",
    action=Quiet, nargs=0)

subparser.add_argument('--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')

subparser.set_defaults(usage_example=
"""This command performs GNM calculations for given PDB structure and outputs \
results in NMD format. If an identifier is passed, structure file will be \
downloaded from the PDB FTP server.

Fetch PDB 1p38, run GNM calculations using default parameters, and results:
    
  $ prody gnm 1p38
    
Fetch PDB 1aar, run GNM calculations with cutoff distance 7 angstrom for \
chain A carbon alpha atoms with residue numbers less than 70, and \
save all of the graphical output files:

  $ prody gnm 1aar -c 7 -s "calpha and chain A and resnum < 70" -A"""
)

group = addNMAParameters(subparser)

group.add_argument('-c', '--cutoff', dest='cutoff', type=float, 
                   default=10.0, metavar='FLOAT', 
                   help='cutoff distance (A) (default: "%(default)s")')
group.add_argument('-g', '--gamma', dest='gamma', type=float, 
                   default=1.0, metavar='FLOAT', 
                   help='spring constant (default: %(default)s)')
group.add_argument('-m', '--model', dest='model', type=int, 
                   default=1, metavar='INT', 
                   help=('model that will be used in the calculations')) 
                    
group = addNMAOutput(subparser)
group.add_argument('-b', '--beta-factors', dest='beta', action='store_true', 
                   default=False, help='write B-factors')
group.add_argument('-k', '--kirchhoff', dest='kirchhoff', action='store_true', 
                   default=False, help='write Kirchhoff matrix')

group = addNMAOutputOptions(subparser, '_gnm')

group = addNMAFigures(subparser)
group.add_argument('-B', '--beta-factors-figure', dest='bf', 
                   action='store_true', default=False, 
                   help='save beta-factors')
group.add_argument('-K', '--contact-map', dest='cm', action='store_true', 
                   default=False, 
                   help='save contact map (Kirchhoff matrix)')        
group.add_argument('-M', '--mode-shape-figure', dest='modes', type=str, 
                   default='', metavar='STR',
                   help=('save mode shape figures for specified modes, '
                         'e.g. "1-3 5" for modes 1, 2, 3 and 5'))

group = addNMAFigureOptions(subparser)

subparser.add_argument('pdb', help='PDB identifier or filename')
    
subparser.set_defaults(func=prody_gnm)
subparser.set_defaults(subparser=subparser)

###############################################################################
# pca
###############################################################################

subparser = commands.add_parser('pca',
    help='perform principal component analysis calculations')

subparser.add_argument('--quiet', help="suppress info messages to stderr",
    action=Quiet, nargs=0)

subparser.add_argument('--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')

subparser.set_defaults(usage_example=
"""This command performs PCA (or EDA) calculations for given multi-model PDB \
structure or DCD format trajectory file and outputs results in NMD format.  \
If a PDB identifier is given, structure file will be downloaded from the PDB  \
FTP server.  DCD files may be accompanied with PDB or PSF files to enable \
atoms selections.

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

###############################################################################
# eda
###############################################################################

subparser = commands.add_parser('eda', parents=[subparser], 
    help='perform essential dynamics analysis calculations', add_help=False)


###############################################################################
# align
###############################################################################

subparser = commands.add_parser('align', 
    help='align models or structures')

subparser.add_argument('--quiet', help="suppress info messages to stderr",
    action=Quiet, nargs=0)

subparser.add_argument('--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')

subparser.set_defaults(usage_example=
"""Align models in PDB structure or multiple PDB structures and save aligned \
coordinate sets.  When multiple structures are aligned, ProDy will match \
chains and use best match for aligning the structures.  Note that options \
are not used when aligning multiple structures.

Fetch PDB structure 2k39 and align models:
    
    $ prody align 2k39
    
Fetch PDB structure 2k39 and align models using backbone of residues with \
number smaller than 71:

    $ prody align 2k39 --select "backbone and resnum < 71" 
    
Fetch PDB structures 1p38, 1r39 and 1zz2 and superpose 1r39 and 1zz2 onto 1p38:

    $ prody align 1p38 1r39 1zz2"""
)
    
subparser.add_argument('-p', '--prefix', dest='prefix', type=str, 
    default='', metavar='STR', 
    help=('prefix for output files, default is PDB_aligned'))
subparser.add_argument('-s', '--select', dest='select', type=str, 
    default='calpha', metavar='SELSTR',
    help='selection string (default: "%(default)s")')
subparser.add_argument('-m', '--model', dest='model', type=int, 
    default=1, metavar='INT',
    help=('for NMR files, reference model index (default: %(default)s)'))

subparser.add_argument('pdb', nargs='+', 
    help='PDB identifier(s) or filename(s)')
        
subparser.set_defaults(func=prody_align)
subparser.set_defaults(subparser=subparser)

###############################################################################
# biomol
###############################################################################

subparser = commands.add_parser('biomol', 
    help='build biomolecules')

subparser.add_argument('--quiet', help="suppress info messages to stderr",
    action=Quiet, nargs=0)

subparser.add_argument('--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')

subparser.set_defaults(usage_example=
"""Generate biomolecule coordinates:
    
  $ prody biomol 2bfu"""
)

subparser.add_argument('-p', '--prefix', dest='prefix', type=str, 
    default=None, metavar='STR', 
    help=('prefix for output files (default: pdb_biomol_)'))
subparser.add_argument('-b', '--biomol', dest='biomol', type=int, 
    default=None, metavar='INT',
    help='index of the biomolecule, by default all are generated')

subparser.add_argument('pdb', help='PDB identifier or filename')

subparser.set_defaults(func=prody_biomol)
subparser.set_defaults(subparser=subparser)

###############################################################################
# blast
###############################################################################

subparser = commands.add_parser('blast', 
    help='blast search Protein Data Bank')

subparser.add_argument('--quiet', help="suppress info messages to stderr",
    action=Quiet, nargs=0)

subparser.add_argument('--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')

subparser.set_defaults(usage_example=
"""Blast search PDB for the first sequence in a fasta file:
    
  $ prody blast seq.fasta -i 70

Blast search PDB for the sequence argument:

  $ prody blast MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ\
KESTLHLVLRLRGG

Blast search PDB for avidin structures, download files, and align all files \
onto the 2avi structure:
    
  $ prody blast -d . ARKCSLTGKWTNDLGSNMTIGAVNSRGEFTGTYITAVTATSNEIKESPLHGTQNTIN\
KRTQPTFGFTVNWKFSESTTVFT
  $ prody align 2avi.pdb *pdb """)

subparser.add_argument('-i', '--identity', dest='identity', type=float, 
    default=90.0, metavar='FLOAT', 
    help='percent sequence identity (default: %(default)s)')
subparser.add_argument('-o', '--overlap', dest='coverage', type=float, 
    default=90.0, metavar='FLOAT', 
    help='percent sequence overlap (default: %(default)s)')
subparser.add_argument('-d', '--dir', dest='folder', type=str,
    default=None, metavar='PATH', 
    help=('download uncompressed PDB files to given path'))

subparser.add_argument('-z', '--gzip', dest='gzip', action='store_true', 
                 default=False, help='write compressed PDB file')

subparser.add_argument('seq', type=str,  
    help=('sequence or file in fasta format'))

subparser.set_defaults(func=prody_blast)
subparser.set_defaults(subparser=subparser)

###############################################################################
# catdcd
###############################################################################

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
    
  $ prody catdcd mdm2.dcd mdm2sim2.dcd --pdb mdm2.pdb -s bb"""
)

subparser.add_argument('-s', '--select', default='all', type=str, 
    dest='select', metavar='SELSTR', 
    help='atom selection (default: "%(default)s")')

subparser.add_argument('-o', '--output', type=str, metavar='FILENAME', 
    default='trajectory.dcd',
    help='output filename (default: %(default)s)')

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

subparser.add_argument('dcd', nargs='+',
    help='DCD filename(s) (all must have same number of atoms)')

subparser.set_defaults(subparser=subparser)
subparser.set_defaults(func=prody_catdcd)

###############################################################################
# fetch
###############################################################################

subparser = commands.add_parser('fetch', 
    help='fetch a PDB file')
    
subparser.add_argument('--quiet', help="suppress info messages to stderr",
    action=Quiet, nargs=0)

subparser.add_argument('--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')

subparser.set_defaults(usage_example=
"""Download PDB file(s) by specifying identifiers:
    
  $ prody fetch 1mkp 1p38"""
)

subparser.add_argument('-d', '--dir', dest='folder', type=str,
                  default='.', metavar='PATH', 
                  help=('target directory for saving PDB files'))
subparser.add_argument('-f', '--file', dest="listfn", type=str, 
                  default='', metavar='FILE', 
                  help='file that contains PDB identifiers')
subparser.add_argument('-z', '--gzip', dest='gzip', action='store_true', 
                 default=False, help='write compressed PDB file')

subparser.add_argument('pdb', nargs='+', help='PDB identifier(s)')

subparser.set_defaults(func=prody_fetch)
subparser.set_defaults(subparser=subparser)

###############################################################################
# select
###############################################################################

subparser = commands.add_parser('select', 
    help='select atoms and write a PDB file')

subparser.add_argument('--quiet', help="suppress info messages to stderr",
    action=Quiet, nargs=0)

subparser.add_argument('--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')

subparser.set_defaults(usage_example=
"""This command selects specified atoms and writes them in a PDB file.

Fetch PDB 2bfu and write backbone atoms in a file:
        
  $ prody select 2bfu "backbone" """
)


subparser.add_argument('-o', '--output', dest='output', type=str, 
    metavar='STR', help="output filanem (default: 'pdb_selected.pdb')")
    
subparser.add_argument('pdb', help='PDB identifier or filename')
subparser.add_argument('selstr', help='atom selection string')

subparser.set_defaults(func=prody_select)
subparser.set_defaults(subparser=subparser)

def main():
    
    if len(sys.argv) == 1:    
        parser.print_help()
    else:
        args = parser.parse_args()
        args.func(args)

    
if __name__ == '__main__':
    main()
