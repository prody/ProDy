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

"""Refine MSA application."""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..apptools import DevelApp
import prody
import numpy as np

__all__ = ['evol_rankorder']

APP = DevelApp('rankorder', 'identify highly coevolving pairs of residues')

APP.setExample(
"""This application identifies that top ranking pairs of residues that \
coevolve based on their mutual information.  By default coevolution is \
reported for pairs that are at least 3 residues apart in sequence. A z-score \
normalization can be applied to the mutinfo matrix to identify coevolving \
pairs.  The following examples show how to use with default as well as \
additional options:

    $ evol rankorder piwi_refined_mutinfo.txt -z
    
    $ evol rankorder piwi_refined_mutinfo.txt --msa piwi_refined.slx \
--label AGO6_ARATH""", [])


APP.addArgument('mutinfo', 
    help='mutual information matrix')

APP.addGroup('input', 'input options')
APP.addArgument('-z', '--zscore', 
    dest='zscore', 
    action='store_true', 
    help='apply zscore for identifying top ranked coevolving pairs',
    group='input'
    )
APP.addArgument('-d', '--delimiter', 
    dest='delimiter',  
    help='delimiter used in mutual information matrix file',
    type=str,
    metavar='STR',
    default=None,
    group='input'
    )
APP.addArgument('-p', '--pdb',
    dest='pdb',
    help='PDB file that contains same number of residues as the mutual '
    'information matrix, output residue numbers will be based on PDB file',
    default=None,
    type=str,
    metavar='STR',
    group='input'
    )
APP.addArgument('-m', '--msa',
    dest='msa',
    help='MSA file used for building the mutual info matrix, '
    'output residue numbers will be based on the most complete sequence '
    'in MSA if a PDB file or sequence label is not specified',
    default=None,
    type=str,
    metavar='STR',
    group='input'
    )
APP.addArgument('-l', '--label',
    dest='label',
    help='label in MSA file for output residue numbers',
    default=None,
    type=str,
    metavar='STR',
    group='input'
    )
APP.addGroup('output', 'output options')
APP.addArgument('-n', '--num-pairs',
    dest='numpairs',
    help='number of top ranking residue pairs to list',
    default=100,
    type=int,
    metavar='INT',
    group='output'
    )
APP.addArgument('-q', '--seq-sep',
    dest='seqsep',
    help='report coevolution for residue pairs that are sequentially '
    'separated by input value',
    default=3,
    type=int,
    metavar='INT',
    group='output'
    )
APP.addArgument('-t', '--min-dist',
    dest='dist',
    help='report coevolution for residue pairs whose CA atoms are spatially '
    'separated by at least the input value, used when a PDB file is given '
    'and --use-dist is true',
    default=10.0,
    type=float,
    metavar='FLOAT',
    group='output'
    )
APP.addArgument('-u', '--use-dist',
    dest='usedist',
    action='store_true', 
    help='use structural separation to report coevolving pairs',
    group='output'
    )
APP.addArgument('-o', '--outname',
    dest='outname',
    help='output filename, default is mutinfo_rankorder.txt',
    type=str,
    metavar='STR',
    group='output'
    )

def calcAllDist(coordset):
    
    from prody import calcDistance
    shape = coordset.shape
    distance = np.zeros((shape[1], shape[1]))
    
    for i in range(shape[1]):
        temp = np.tile(coordset[0,i,:], (1, shape[1], 1))
        distance[:,i] = calcDistance(coordset, temp)
    
    return distance

def evol_rankorder(mutinfo, **kwargs):
    from prody import parseMSA, LOGGER, parsePDB, calcMSAOccupancy
    from prody.utilities import openFile
    from os.path import splitext
    
    delimiter = kwargs.get('delimiter')
    mi = np.loadtxt(str(mutinfo), delimiter=delimiter)
    
    ndim, shape = mi.ndim, mi.shape
    if ndim != 2 or shape[0] != shape[1]:
        raise ValueError('mutinfo must contain a square matrix')
    
    msa, label = kwargs.get('msa'), kwargs.get('label')
    
    pdb, pdbflag = kwargs.get('pdb'), False
    
    resnum = None
    
    if pdb is not None:
        from prody import parsePDB
        try:
            pdb = parsePDB(pdb)
        except:
            LOGGER.info('Could not parse PDB, ignoring PDB input')
        else:
            chains = list(pdb.iterChains())
            for chain in chains:
                sel = chain.select('protein and name CA')
                if sel.numAtoms() == shape[0]:
                    resnum = sel.getResnums()
                    coordset = sel.getCoordsets()
                    distance = calcAllDist(coordset)
                    pdbflag = True
                    label = pdb.getTitle()
                    LOGGER.info('Residue numbers will be based on pdb: '
                                '{0}'.format(pdb.getTitle()))
                else:
                    LOGGER.info('Number of residues in PDB does not match '
                                'mutinfo matrix, ignoring PDB input')
    
    if not pdbflag:
        if msa is not None:
            msa = parseMSA(msa)
            if msa.numResidues() != shape[0]:
                LOGGER.info('Input MSA and mutinfo do not have similar no '
                            'of residues, ignoring MSA')
            else:
                index = msa.getIndex(label)   
                if index is None:
                    if label is not None:
                        LOGGER.info('Could not find given label in MSA, '
                                    'using complete sequence from MSA')
                    occ = calcMSAOccupancy(msa._msa, 'row')
                    index = np.where(occ == occ.max())[0][0]
                    label, seq, start, end = msa[index]
                else:
                    label, seq, start, end = msa[index]
                if (start and end is not None) and (start < end):
                    resnum = np.arange(start, end+1)
                    if len(resnum) != shape[0]:
                        LOGGER.info('Label: {0}/{1}-{2} and mutinfo do '
                                    'not have similar no of residues, using '
                                    'serial indexing'.format(label, start, end))
                        label = 'Serial Index'
                        resnum = np.arange(1, shape[0]+1)
                    else:
                        LOGGER.info('Residue numbers will be based on label: '
                                    '{0}'.format(label))
                else:
                    LOGGER.info('Could not identify residue indexes from MSA'
                                    ' using serial indexing')
                    label = 'Serial Index'
                    resnum = np.arange(1, shape[0]+1)
        else:
            LOGGER.info('MSA or PDB not given or does not match mutinfo, '
                        'using serial indexing')
            resnum = np.arange(1, shape[0]+1)
    
    LOGGER.info('Residue numbers start and end with {0}-{1}'.
                format(str(resnum[0]), str(resnum[-1])))
    
    outname = kwargs.get('outname')
    if outname is None:
        outname, ext = splitext(str(mutinfo))
        if ext.lower() == '.gz': 
            outname, _ = splitext(str(mutinfo))
    else:
        outname, ext = splitext(str(outname))
        if ext is None:
            ext = '.txt'
    
    outname += '_rankorder' + ext
    zscore = kwargs.get('zscore')
    if zscore:
        LOGGER.info('zscore normalization applied such that each column '
                    'has 0 mean and standard deviation 1')
        header = 'Serial\tRow\tColumn\tZscore'
        mi = (mi - mi.mean(0)) / mi.std(0)
    else:
        header = 'Serial\tRow\tColumn\tMI'
    
    mi_ind_start, mi_ind_end = np.tril_indices(shape[0], k=-1)
    mi_matrix = mi[mi_ind_start, mi_ind_end]
    sorted_index = mi_matrix.argsort(axis=None)[::-1]
    row = mi_ind_start[sorted_index]
    column = mi_ind_end[sorted_index]
    count = 1
    i = 0
    
    f = openFile(outname, 'wb')
    if label is None:
        label = 'Serial Index'
    f.write(('Label: '+ label + '\t' + 'Residue Numbers: ' +
             str(resnum[0]) + '-' + str(resnum[-1]) + '\n'))
    
    numpairs = kwargs.get('numpairs')
    size = len(row)
    seqsep = kwargs.get('seqsep')
    if not kwargs.get('usedist') or not pdbflag:
        if kwargs.get('usedist'):
            LOGGER.info('use-struct-sep set to true, but PDB not given or '
                        'incorrect residue number. Using sequence separation')
        else:
            if pdbflag:
                LOGGER.info('use-struct-sep not set, using sequence separation'
                            ' to report coevolving pairs')
        f.write((header + '\n'))
        while count <=numpairs  and i < size:        
            if row[i] > (column[i] + seqsep):
                f.write('{0}\t{1}\t{2}\t{3:.3f}\n'.
                        format(count, resnum[row[i]], resnum[column[i]],
                               mi[row[i], column[i]]))
                count += 1
            i += 1
    else:
        structsep = kwargs.get('dist')
        f.write((header + '\tDistance Cutoff: ' + str(structsep) + '\n'))        
        while count <=numpairs  and i < size:        
            if distance[row[i], column[i]] > structsep:
                f.write('{0}\t{1}\t{2}\t{3:.3f}\t{4:.2f}\n'.
                        format(count, resnum[row[i]], resnum[column[i]],
                               mi[row[i], column[i]],
                               distance[row[i], column[i]]))
                count += 1                
            i += 1
    f.close()
    
APP.setFunction(evol_rankorder)
