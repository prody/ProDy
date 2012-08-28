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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import array, invert, ones, tile 

from prody import Ensemble, PDBEnsemble, LOGGER
from prody.tests.test_datafiles import parseDatafile, DATA_FILES

LOGGER.verbosity = 'none'

ATOL = 1e-5
RTOL = 0

ATOMS = parseDatafile('multi_model_truncated', subset='ca')
ALLATOMS = parseDatafile('multi_model_truncated')
DCD = parseDatafile('dcd')
COORDS = ATOMS.getCoords()
COORDSETS = ATOMS.getCoordsets()
ENSEMBLE = Ensemble(ATOMS)
CONF = ENSEMBLE[0]
DATA = DATA_FILES['multi_model_truncated']
ENSEMBLE_RMSD = DATA['rmsd_ca']
ENSEMBLE_SUPERPOSE = DATA['rmsd_ca_aligned']

ENSEMBLEW = Ensemble(ATOMS)
ENSEMBLEW.setWeights(ones(len(ATOMS), dtype=float))
CONFW = ENSEMBLEW[0]
        
PDBENSEMBLE = PDBEnsemble('PDBEnsemble')
PDBENSEMBLE.setCoords(COORDS)
WEIGHTS = []
for i, xyz in enumerate(ATOMS.iterCoordsets()):
    weights = ones((len(ATOMS), 1), dtype=float)
    if i > 0:
        weights[i] = 0
        weights[-i] = 0 
        PDBENSEMBLE.addCoordset(xyz, weights=weights)
    else:
        PDBENSEMBLE.addCoordset(xyz)
    WEIGHTS.append(weights)
PDBCONF = PDBENSEMBLE[0]
WEIGHTS = array(WEIGHTS)
WEIGHTS_BOOL = tile(WEIGHTS.astype(bool), (1,1,3))  
WEIGHTS_BOOL_INVERSE = invert(WEIGHTS_BOOL)
