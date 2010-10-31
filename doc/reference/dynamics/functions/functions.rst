.. module:: prody.dynamics

.. _dyfunctions:

*******************************************************************************
Function Library
*******************************************************************************

Many of the functions documented in this page accepts a *modes* argument (may also have different names).
This argument may be:

  * an NMA model, which may be an instance of one of :class:`ANM`, :class:`GNM`, :class:`PCA`.
  * a :class:`Mode` instance obtained by indexing an NMA model, e.g. ``nma[0]``
  * a list of :class:`Mode` instances obtained by slicing an NMA model, e.g. ``nma[:10]``

Some of these functions may also accept :class:`Vector` instances as *modes* argument.

**List of Functions**:

  * Short-hand functions:

    * :func:`getANM`
    * :func:`getGNM`

  * Analysis:

    * :func:`getCrossCorrelations`
    * :func:`getSqFlucts`  

  * Write data:

    * :func:`writeArray`
    * :func:`writeModes`
    * :func:`writeNMD`
        
  * Visualization:

    * :func:`viewNMDinVMD`
    * :func:`setVMDpath`
    * :func:`getVMDpath`
    
  * Compare NMA models:

    * :func:`getOverlap`
    * :func:`getCumulativeOverlap`
    * :func:`getSubspaceOverlap`
    * :func:`printOverlapTable`
  
  * Other:

    * :func:`getProjection`
    * :func:`getSumOfWeights`
    * :func:`reduceModel`
    
    
.. autofunction:: getANM

.. autofunction:: getGNM

.. autofunction:: getCrossCorrelations

.. autofunction:: getSqFlucts

.. autofunction:: writeArray

.. autofunction:: writeModes

.. autofunction:: writeNMD
        
.. autofunction:: viewNMDinVMD

.. autofunction:: setVMDpath

.. autofunction:: getVMDpath
    
.. autofunction:: getOverlap
    
.. autofunction:: getCumulativeOverlap

.. autofunction:: getSubspaceOverlap

.. autofunction:: printOverlapTable

.. autofunction:: getProjection

.. autofunction:: getSumOfWeights

.. autofunction:: reduceModel
