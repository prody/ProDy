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

DEFAULTS = {}
HELPTEXT = {} 

FIGFORMATS = 'eps pdf png ps raw rgba svg svgz'.split() 

for key, txt, val in [
    ('outdir', 'output directory', '.'),
    ('prefix', 'output file prefix', ''),
    ('select', 'atom selection',
         "protein and name CA or nucleic and name P C4' C2"),
    ('nmodes', 'number of non-zero eigenvectors (modes) to calculate', 10),

    ('outall', 'write all outputs', False),
    ('outeig', 'write eigenvalues/vectors', False),
    ('outcov', 'write covariance matrix', False),
    ('outnpz', 'write compressed ProDy data file', False),
    ('outcc', 'write cross-correlations', False),
    ('outhm', 'write cross-correlations heatmap file', False),
    ('outsf', 'write square-fluctuations', False),
    ('extend', 'write NMD file for the model extended to '
               '"backbone" ("bb") or "all" atoms of the residue,'
               ' model must have one node per residue', ''),

    ('numformat', 'number output format', '%12g'),
    ('numext', 'numeric file extension', '.txt'),
    ('numdelim', 'number delimiter', ' '),
    
    
    ('figformat', 'figure file format', 'pdf'),
    ('figwidth', 'figure width (inch)', 8.),
    ('figheight', 'figure height (inch)', 6.),
    ('figdpi', 'figure resolution (dpi)', 300),
    ('figall', 'save all figures', False),
    ('figcc', 'save cross-correlations figure', False),
    ('figsf', 'save square-fluctuations figure', False),]:
    
    DEFAULTS[key] = val
    HELPTEXT[key] = txt
    

__all__ = ['addNMAParameters', 'addNMAOutput', 'addNMAOutputOptions',
           'addNMAFigures', 'addNMAFigureOptions']


def addNMAParameters(parser):
    
    parser = parser.add_argument_group('parameters')

    parser.add_argument('-n', '--number-of-modes', dest='nmodes', type=int, 
        default=DEFAULTS['nmodes'], metavar='INT', 
        help=HELPTEXT['nmodes'] + ' (default: %(default)s)')    

    parser.add_argument('-s', '--select', dest='select', type=str, 
        default=DEFAULTS['select'], metavar='SEL', 
        help=HELPTEXT['select'] + ' (default: "%(default)s")')

    return parser


def addNMAOutput(parser):
    
    parser = parser.add_argument_group('output')
    
    parser.add_argument('-a', '--all-output', 
        dest='outall', action='store_true', 
        default=DEFAULTS['outall'], help=HELPTEXT['outall'])
    
    parser.add_argument('-o', '--output-dir', dest='outdir', type=str, 
        default=DEFAULTS['outdir'], metavar='PATH', 
        help=HELPTEXT['outdir'] + ' (default: %(default)s)')
    
    parser.add_argument('-e', '--eigenvs', dest='outeig', action='store_true', 
        default=DEFAULTS['outeig'], help=HELPTEXT['outeig'])
        
    parser.add_argument('-r', '--cross-correlations', dest='outcc', 
        action='store_true', 
        default=DEFAULTS['outcc'], help=HELPTEXT['outcc'])
        
    parser.add_argument('-u', '--heatmap', dest='outhm', 
        action='store_true', 
        default=DEFAULTS['outhm'], help=HELPTEXT['outhm'])

    parser.add_argument('-q', '--square-fluctuations', dest='outsf', 
        action='store_true', 
        default=DEFAULTS['outsf'], help=HELPTEXT['outsf'])
        
    parser.add_argument('-v', '--covariance', dest='outcov', 
        action='store_true', 
        default=DEFAULTS['outcov'], help=HELPTEXT['outcov'])
        
    parser.add_argument('-z', '--npz', dest='outnpz', action='store_true', 
        default=DEFAULTS['outnpz'], help=HELPTEXT['outnpz'])
        
    parser.add_argument('-t', '--extend', dest='extend', type=str,
        metavar='STR', choices=set(['bb', 'all', 'backbone']),
        help=HELPTEXT['extend'])
    
    return parser


def addNMAOutputOptions(parser, prefix):
    
    parser = parser.add_argument_group('output options')
    
    parser.add_argument('-p', '--file-prefix', dest='prefix', type=str, 
        default=prefix, metavar='STR', 
        help=HELPTEXT['prefix'] + ' (default: pdb%(default)s)')
        
    parser.add_argument('-f', '--number-format', dest='numformat', type=str, 
        default=DEFAULTS['numformat'], metavar='STR',
        help=HELPTEXT['numformat'] + ' (default: %(default)s)')
        
    parser.add_argument('-d', '--delimiter', dest='numdelim', type=str, 
        default=DEFAULTS['numdelim'], metavar='STR',
        help=HELPTEXT['numdelim'] + ' (default: "%(default)s")')
        
    parser.add_argument('-x', '--extension', dest='numext', type=str, 
        default=DEFAULTS['numext'], metavar='STR',
        help=HELPTEXT['numext'] + ' (default: %(default)s)')
       
    return parser


def addNMAFigures(parser):
    
    parser = parser.add_argument_group('figures')
    
    parser.add_argument('-A', '--all-figures', dest='figall', 
        action='store_true', 
        default=DEFAULTS['figall'], help=HELPTEXT['figall'])
        
    parser.add_argument('-R', '--cross-correlations-figure', dest='figcc', 
        action='store_true', default=DEFAULTS['figcc'], 
        help=HELPTEXT['figcc'])
        
    parser.add_argument('-Q', '--square-fluctuations-figure', dest='figsf', 
        action='store_true', default=DEFAULTS['figsf'], 
        help=HELPTEXT['figsf'])
        
    return parser


def addNMAFigureOptions(parser):
    
    parser = parser.add_argument_group('figure options')
    
    parser.add_argument('-F', '--figure-format', dest='figformat', type=str, 
        default=DEFAULTS['figformat'], metavar='STR',
        help=DEFAULTS['figformat'] + ' (default: %(default)s)', 
        choices=set(FIGFORMATS))
             
    parser.add_argument('-D', '--dpi', dest='figdpi', type=int, 
        default=DEFAULTS['figdpi'], metavar='INT', 
        help=HELPTEXT['figdpi'] + ' (default: %(default)s)')
        
    parser.add_argument('-W', '--width', dest='figwidth', type=float, 
        default=DEFAULTS['figwidth'], metavar='FLOAT', 
        help=HELPTEXT['figwidth'] + ' (default: %(default)s)')
        
    parser.add_argument('-H', '--height', dest='figheight', type=float, 
        default=DEFAULTS['figheight'], metavar='FLOAT', 
        help=HELPTEXT['figheight'] + ' (default: %(default)s)')
        
    return parser
