*******************************************************************************
Changes
*******************************************************************************

Release 0.2.0
===============================================================================

**Important Changes**:


  * Modules are reorganized.
  * :func:`prody.compare.findMatchingChains` is renamed to 
    :func:`prody.compare.matchChains`.
  * :func:`prody.dynamics.showOverlapMatrix` is renamed to 
    :func:`prody.dynamics.showOverlapTable`.

**New Features**:

  * :class:`prody.select.Select` can not identify intermolecular contacts. See
    :ref:`contacts` for an example

**Improvements**:

  * :mod:`prody.compare` functions are improved. Now they perform sequence
    alignment if simple residue number/identity based matchin does not work,
    or if user passes ``pwalign=True`` argument. This impacts the speed 
    of X-ray ensemble analysis.
  * :class:`prody.select.Select` can cache data optionally. This results
    in speeds up from 2 to 40 folds depending on number of atoms and selection
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
