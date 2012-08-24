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

"""This module defines options common to anm/gnm/pca commands."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

DEFAULT_FIGURE_FORMAT = 'pdf'
DEFAULT_FILE_EXT = '.txt'
DEFAULT_NUMBER_FORMAT = '%12g'
DEFAULT_DELIMITER = ' '

__all__ = ['addNMAParameters', 'addNMAOutput', 'addNMAOutputOptions',
           'addNMAFigures', 'addNMAFigureOptions']

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
    parser.add_argument('-t', '--extend', dest='extend', type=str,
                        metavar='STR', choices=set(['bb', 'all', 'backbone']),
                        help='output NMD files for the model extended to '
                             '"backbone" ("bb") or "all" atoms of the residue,'
                             ' model must have one node per residue')
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
