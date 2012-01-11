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
    * :ref:`pca-xray-calculations`
    * :ref:`pca-nmr`
    * :ref:`command-anm`
    
  * Gaussian Network Model (:class:`~dynamics.GNM`)
 
    * :ref:`gnm`
    * :ref:`command-gnm`
      
**Principal Component Analysis:** (:class:`~dynamics.PCA`)
   
  * :ref:`pca-xray`
  * :ref:`pca-nmr`
  * :ref:`pca-blast`
  * :ref:`eda`
  * :ref:`command-pca`

**Trajectory Analysis:**

  * :ref:`eda`
  * :ref:`trajectory`
  
**Execute/parse secondary structure analysis programs/data:**
   
  * :func:`~proteins.execDSSP`
  * :func:`~proteins.parseDSSP`
  * :func:`~proteins.performDSSP`
  * :func:`~proteins.execSTRIDE`
  * :func:`~proteins.parseSTRIDE`
  * :func:`~proteins.performSTRIDE`

**Generate conformations along normal modes:**

  * Sample random conformations along a set of modes (:func:`~dynamics.sampleModes`)
  * Generate evenly spaced conformations along a single mode (:func:`~dynamics.traverseMode`)

**Analyze dynamics models:**
  * Show ANM/GNM contact map (:func:`~dynamics.showContactMap`)
  * Show cross-correlations (:func:`~dynamics.showCrossCorrelations`)
  * Show fraction of variances (:func:`~dynamics.showFractOfVariances`)
  * Show cumulative fraction of variances (:func:`~dynamics.showCumFractOfVariances`)
  * Show projection onto normal modes (:func:`~dynamics.showProjection`)
  * Show square-fluctuations (:func:`~dynamics.showSqFlucts`)
  * Show mode shape (:func:`~dynamics.showMode`)

**Compare experimental and theoretical models:**
  
  * Calculate overlap (:func:`~dynamics.calcOverlap`)
  * Calculate cumulative overlap (:func:`~dynamics.calcCumulativeOverlap`)
  * Calculate subspace overlap (:func:`~dynamics.calcSubspaceOverlap`)
  * Print overlap table (:func:`~dynamics.printOverlapTable`)
  * Show overlap (:func:`~dynamics.showOverlap`)
  * Show cumulative overlap (:func:`~dynamics.showCumulativeOverlap`)
  * Show overlap table (:func:`~dynamics.showOverlapTable`)
  * Project conformations onto modes from different models (:func:`~dynamics.showCrossProjection`)
  * Show square fluctuations from two different models (:func:`~dynamics.showScaledSqFlucts`)
  * Show square fluctuations from two different models (:func:`~dynamics.showNormedSqFlucts`)

**Model reduction**:
    
  * Reduce an ANM, PCA, or GNM model to a subset of atoms (:func:`~dynamics.reduceModel`)
  * Get a slice of a mode matching a selection of atoms (:func:`~dynamics.sliceMode`)
  * Get a slice of a vector matching a selection of atoms (:func:`~dynamics.sliceVector`)

**Write output files:**
  
  * Write PDB file (:func:`~proteins.writePDB`)
  * Write normal mode data in NMD format (:func:`~dynamics.writeNMD`)
  * Write Numpy arrays in file (:func:`~dynamics.writeArray`)
  * Write normal mode eigenvectors (:func:`~dynamics.writeModes`)
  * Write overlap table (:func:`~dynamics.writeOverlapTable`)
  
**Read input files:**
  
  * Parse PDB files (:func:`~proteins.parsePDB`)
  * Parse NMD files (:func:`~dynamics.parseNMD`)
  * Parse arrays from file (:func:`~dynamics.parseArray`)
  * Parse normal mode data (:func:`~dynamics.parseModes`)
  * Parse DCD files (:func:`~ensemble.parseDCD`)
  * Parse DSSP output (:func:`~proteins.parseDSSP`)
  * Parse STRIDE output (:func:`~proteins.parseSTRIDE`)
  
**Visualization**:
  
  * Get path to VMD executable (:func:`~dynamics.getVMDpath`)
  * Set path to VMD executable (:func:`~dynamics.setVMDpath`)
  * Visualize NMD file in VMD (:func:`~dynamics.viewNMDinVMD`)

**Save/load ProDy session data**:
    
  * Save/load ANM/GNM/PCA data (:func:`~dynamics.saveModel`, :func:`~dynamics.loadModel`)
  * Save/load arbitrary vectors (:func:`~dynamics.saveVector`, :func:`~dynamics.loadVector`)
  * Save/load ensembles (:func:`~ensemble.saveEnsemble`, :func:`~ensemble.loadEnsemble`)
  * Save/load atomic data (:func:`~atomic.saveAtoms`, :func:`~atomic.loadAtoms`)
