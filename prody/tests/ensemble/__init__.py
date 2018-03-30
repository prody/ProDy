from numpy import array, invert, ones, tile

from prody import Ensemble, PDBEnsemble, LOGGER
from prody.tests.datafiles import parseDatafile, DATA_FILES

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

PDBENSEMBLEA = PDBEnsemble('PDBEnsemble 2')
PDBENSEMBLEA.setAtoms(ATOMS)
PDBENSEMBLEA.setCoords(COORDS)
for i in range(ATOMS.numCoordsets()):
    weights = ones((len(ATOMS), 1), dtype=float)
    ATOMS.setACSIndex(i)
    if i > 0:
        weights[i] = 0
        weights[-i] = 0
        PDBENSEMBLEA.addCoordset(ATOMS, weights=weights, degeneracy=True)
    else:
        PDBENSEMBLEA.addCoordset(ATOMS, degeneracy=True)

ATOMS.setACSIndex(0)
