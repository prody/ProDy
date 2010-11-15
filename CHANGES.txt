*******************************************************************************
Changes
*******************************************************************************

Release 0.2.0a1
===============================================================================

**Important Changes**:


  * Modules are reorganized.
  * :func:`prody.compare.findMatchingChains` is renamed to 
    :func:`prody.compare.matchChains`.
  * :func:`prody.dynamics.showOverlapMatrix` is renamed to 
    :func:`prody.dynamics.showOverlapTable`.
  * Single word keywords *not* followed by "and" logical operator are not 
    accepted, e.g. "protein within 5 of water" will raise an SelectionError, 
    use "protein and within 5 of water" instead.

**New Features**:

  * :class:`prody.select.Select` can identify intermolecular contacts. See
    :ref:`contacts` for an examples and details.
  * :func:`prody.dynamics.sampleModes` implemented for sampling conformations 
    along normal modes.

**Improvements**:

  * :mod:`prody.compare` functions are improved. Now they perform sequence
    alignment if simple residue number/identity based matchin does not work,
    or if user passes ``pwalign=True`` argument. This impacts the speed 
    of X-ray ensemble analysis.
  * :class:`prody.select.Select` can cache data optionally. This results
    in speeds up from 2 to 50 folds depending on number of atoms and selection
    operations.

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
