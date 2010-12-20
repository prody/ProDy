.. currentmodule:: prody

*******************************************************************************
Changes
*******************************************************************************

Release 0.5
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

  * Iterating over a :class:`atomic.Chain` instance yields :class:`atomic.Residue` 
    instances.
  * :class:`dynamics.Vector` instantiation requires an *array* only. *name* 
    is an optional argument.
  * Functions starting with ``get`` and performing a calculations are renamed to
    start with ``calc``, e.g. :func:`measure.getRMSD` is renamed to
    :func:`~prody.measure.calcRMSD`

Release 0.2
===============================================================================

**Important Changes**:


  * Single word keywords *not* followed by "and" logical operator are not 
    accepted, e.g. "protein within 5 of water" will raise an SelectionError, 
    use "protein and within 5 of water" instead.
  * :func:`compare.findMatchingChains` is renamed to 
    :func:`compare.matchChains`.
  * :func:`dynamics.showOverlapMatrix` is renamed to 
    :func:`dynamics.showOverlapTable`.
  * Modules are reorganized.

**New Features**:

  * :class:`atomic.Atomic` for easy type checking.
  * :class:`select.Contacts` for faster intermolecular contact 
    identification.
  * :class:`select.Select` can identify intermolecular contacts. See
    :ref:`contacts` for an examples and details.
  * :func:`dynamics.getCumulativeOverlapArray`.
  * :func:`dynamics.sampleModes` implemented for sampling conformations 
    along normal modes.

**Improvements**:

  * :mod:`prody.compare` functions are improved. Now they perform sequence
    alignment if simple residue number/identity based matchin does not work,
    or if user passes ``pwalign=True`` argument. This impacts the speed 
    of X-ray ensemble analysis.
  * :class:`select.Select` can cache data optionally. This results
    in speeds up from 2 to 50 folds depending on number of atoms and selection
    operations.
  * Implementation of :func:`dynamics.showProjection` is completed.

Release 0.1.2
===============================================================================

* Important bugfixes and improvements in NMA helper and plotting functions.
* Documentation updates and improvements.


Release 0.1.1
===============================================================================

* Important bugfixes and improvements in chain comparison functions.
* Bugfixes.
* Source clean up.
* Documentation improvements.
