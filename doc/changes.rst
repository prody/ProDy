*******************************************************************************
Changes
*******************************************************************************

Release 0.2.0
===============================================================================

**Important Changes**:


  * Modules are reorganized.
  * :func:`prody.compare.findMatchingChains` is renamed as 
    :func:`prody.compare.matchChains`.

**New Features**:

  * :class:`prody.select.Select` can not identify intermolecular contacts. See
    :ref:`contacts` for an example

**Improvements**:

  * :mod:`prody.compare` functions are improved. Now they perform sequence
    alignment if simple residue number/identity based matchin does not work,
    or if user passes ``pwalign=True`` argument. This impacts the speed 
    of X-ray ensemble analysis.

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
