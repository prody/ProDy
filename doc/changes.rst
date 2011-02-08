.. _changes:

.. currentmodule:: prody

*******************************************************************************
Changes
*******************************************************************************

Release 0.6 (in development)
===============================================================================

**Improvements**:

  * A few bugs in ensemble and dynamics classes are fixed.
  * Improvements in :class:`~proteins.RCSB_PDBFetcher` allow it not to miss a 
    PDB file if it exists in the target folder.
  * :func:`~dynamics.writeNMD` is fixed to output B-factors (Thanks to 
    Dan Holloway for pointing it out).

**New Features**:

  * Membership, equality, and non-equality test operation are defined for all 
    :mod:`~prody.atomic` classes. See :ref:`selection-operations`.
  * Two functions are implemented for dealing with anisotropic temperature 
    factors: :func:`~measure.calcADPAxes` and :func:`~measure.buildADPMatrix`.
  * :meth:`~dynamics.NMA.setEigens` method is implemented for analyzing
    normal modes calculated using external software.
  * :func:`~dynamics.parseArray` is implementing for reading numeric data,
    particularly normal mode data calculated using other software for analysis
    using ProDy. 
  * The method in [BH02]_ to calculate overlap between covariance matrices is 
    implemented as :func:`~dynamics.calcCovarianceOverlap` function.
  * A method to trim :class:`~ensemble.Ensemble` instances is implemented: 
    :func:`~ensemble.trimEnsemble`

**Changes**:

  * Change in default behavior of :func:`~proteins.parsePDB` function. When
    alternate locations exist, those indicated by A are parsed. For parsing
    all alternate locations user needs to pass ``altloc=True`` argument.    
  * :func:`getSumOfWeights` is renamed as :func:`~ensemble.calcSumOfWeights`.
  * :func:`mapAtomsToChain` is renamed as :func:`~compare.mapOntoChain`.


Release 0.5.2 (Jan 12, 2011)
===============================================================================

**Bugfixes**:
  
  * An important fix in :func:`~dynamics.sampleModes` function was made
    (Thanks to Alberto Perez for finding the bug and suggesting a solution). 
    
**Improvements**:
  
  * Improvements in :meth:`dynamics.ANM.calcModes`, :meth:`dynamics.GNM.calcModes`, 
    and :meth:`dynamics.PCA.calcModes` methods prevent Numpy/Scipy throwing an
    exception when more than available modes are requested by the user.
  * Improvements in :func:`proteins.blastPDB` enable ProDy throw an 
    exception when no internet connection is found, and warn user when
    downloads fail due to restriction in network regulations (Thanks to Serkan
    Apaydin for helping identify these improvements).
  * New example :ref:`writepdb`.

Release 0.5.1 (Dec 31, 2010)
===============================================================================

**Changes in dependencies**:

* Scipy (linear algebra module) is not required package anymore. 
  When available it replaces Numpy (linear algebra module) for greater flexibility
  and efficiency. A warning message is printed when Scipy is not found.
* Biopython KDTree module is not required for ENM calculations (specifically
  for building Hessian (ANM) or Kirchoff (GNM) matrices). When available it 
  is used to increase the performance. A warning message is printed when 
  KDTree is not found.

Release 0.5 (Dec 21, 2010)
===============================================================================

**New Features**:

  * :class:`~atomic.AtomPointer` base class for classes pointing to
    atoms in an :class:`~atomic.AtomGroup`.
  * :class:`~atomic.AtomPointer` instances (Selection, Residue, etc.)
    can be added. See :ref:`selection-operations` for examples.
  * :meth:`~select.Select.getIndices` and :meth:`~select.Select.getBoolArray` 
    methods to expand the usage of :class:`~select.Select`.
  * :meth:`~dynamics.sliceVector` and :meth:`~dynamics.sliceMode` functions.
  * :meth:`~dynamics.saveModel` and :meth:`~dynamics.loadModel` functions
    for saving and loading NMA data.
  * :func:`~proteins.parsePDBStream` can now parse specific chains or
    alternate locations from a PDB file.
  * :func:`~measure.alignCoordsets` is implemented to superimpose
    coordinate sets of an :class:`~atomic.AtomGroup` instance.

**Bugfixes**:

  * A bug in :func:`proteins.parsePDBStream` that caused unidentified errors 
    when a model in a multiple model file did not have the same number of 
    atoms is fixed.

**Changes**:

  * Iterating over a :class:`~atomic.Chain` instance yields :class:`atomic.Residue` 
    instances.
  * :class:`~dynamics.Vector` instantiation requires an *array* only. *name* 
    is an optional argument.
  * Functions starting with ``get`` and performing a calculations are renamed to
    start with ``calc``, e.g. :func:`getRMSD` is renamed to
    :func:`~prody.measure.calcRMSD`

Release 0.2 (Nov 16, 2010)
===============================================================================

**Important Changes**:


  * Single word keywords *not* followed by "and" logical operator are not 
    accepted, e.g. "protein within 5 of water" will raise an SelectionError, 
    use "protein and within 5 of water" instead.
  * :func:`findMatchingChains` is renamed to 
    :func:`~compare.matchChains`.
  * :func:`showOverlapMatrix` is renamed to 
    :func:`~dynamics.showOverlapTable`.
  * Modules are reorganized.

**New Features**:

  * :class:`~atomic.Atomic` for easy type checking.
  * :class:`~select.Contacts` for faster intermolecular contact 
    identification.
  * :class:`~select.Select` can identify intermolecular contacts. See
    :ref:`contacts` for an examples and details.
  * :func:`~dynamics.getCumulativeOverlapArray`.
  * :func:`~dynamics.sampleModes` implemented for sampling conformations 
    along normal modes.

**Improvements**:

  * :mod:`~prody.compare` functions are improved. Now they perform sequence
    alignment if simple residue number/identity based matchin does not work,
    or if user passes ``pwalign=True`` argument. This impacts the speed 
    of X-ray ensemble analysis.
  * :class:`~select.Select` can cache data optionally. This results
    in speeds up from 2 to 50 folds depending on number of atoms and selection
    operations.
  * Implementation of :func:`~dynamics.showProjection` is completed.

Release 0.1.2 (Nov 9, 2010)
===============================================================================

* Important bugfixes and improvements in NMA helper and plotting functions.
* Documentation updates and improvements.


Release 0.1.1 (Nov 8, 2010)
===============================================================================

* Important bugfixes and improvements in chain comparison functions.
* Bugfixes.
* Source clean up.
* Documentation improvements.

Release 0.1 (Nov 7, 2010)
===============================================================================

* First release.
