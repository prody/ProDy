.. currentmodule:: prody

.. _features:

*******************************************************************************
Feature Index
*******************************************************************************

**Protein Structural Data:**
  
  * :ref:`blastpdb` (:func:`~proteins.blastPDB`)
  * :ref:`fetchpdb` (:func:`~proteins.fetchPDB`)
  * :ref:`parsepdb` (:func:`~proteins.parsePDB`)

**Atom selections and contacts:**

  * :ref:`selections` (:class:`~select.Select`)
  * :ref:`contacts` (:class:`~select.Contacts`)
  * :ref:`selection-operations`

**Elastic Network Models:**

  * Anisotropic Network Model (:class:`~dynamics.ANM`)
 
    * :ref:`anm`
    * :ref:`p38-xray-calculations`
    * :ref:`pca-nmr`
    * :ref:`scripts-anm`
    
  * Gaussian Network Model (:class:`~dynamics.GNM`)
 
    * :ref:`gnm`
    * :ref:`scripts-gnm`
      
**Principal Component Analysis:** (:class:`~dynamics.PCA`)
   
  * :ref:`pca-xray`
  * :ref:`pca-nmr`
  * :ref:`eda-md`
   
**Generate conformations along normal modes:**

  * Sample random conformations along a set of modes (:func:`~dynamics.sampleModes`)
  * Generate evenly spaced conformations along a single mode (:func:`~dynamics.traverseMode`)

**Compare experimental and theoretical models:**
  
  * Calculate overlap (:func:`~dynamics.calcOverlap`)
  * Calculate cumulative overlap (:func:`~dynamics.calcCumulativeOverlap`)
  * Calculate subspace overlap (:func:`~dynamics.calcSubspaceOverlap`)
  * Print overlap table (:func:`~dynamics.printOverlapTable`)

**Write output files:**
  
  * Write PDB file (:func:`~proteins.writePDB`)
  * Write normal mode data in NMD format (:func:`~dynamics.writeNMD`)
  * Write Numpy arrays in file (:func:`~dynamics.writeArray`)
  * Write normal mode eigenvectors (:func:`~dynamics.writeModes`)
  * Write overlap table (:func:`~dynamics.writeOverlapTable`)
  
**Visualization**:
  
  * Get path to VMD executable (:func:`~dynamics.getVMDpath`)
  * Set path to VMD executable (:func:`~dynamics.setVMDpath`)
  * Visualize normal modes in VMD (:func:`~dynamics.viewNMDinVMD`)
