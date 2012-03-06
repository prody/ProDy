.. _changes:

.. currentmodule:: prody

*******************************************************************************
Changes
*******************************************************************************

Release 1.0 (in development)
===============================================================================

**Improvements**:

  * Added ``bonded to ...`` selection method that expands a selection to 
    immediately bound atoms.  See :ref:`selections` for its description.
  
  * :func:`~.fetchPDBLigand` parses bond data from the XML file.

  * :func:`~.fetchPDBLigand` can optionally save compressed XML files into
    ProDy package folder so that frequent access to same files will be more
    rapid. See :func:`~.confProDy` function for setting this option.

  * :class:`~.Select` class is revised. All exceptions are handled delicately
    to increase the stability of the class.

  * Distance based atom selection is 10 to 15% faster for atom groups with
    more than 5K atoms.
    
  * Added uncompressed file saving option to :ref:`prody-blast` command.

**Changes**:

  * Deprecated method and functions are removed.
  
  * The relation between :class:`~.AtomGroup`, :class:`~.Trajectory`, and
    :class:`~.Frame` instances have changed. See :ref:`trajectory2` and
    :ref:`outputtraj`, and :ref:`atomsframes` usage examples.
  
  * :meth:`Mode.getCollectivity` method is removed, use 
    :func:`~.calcCollectivity` function instead.
  
  * :meth:`Mode.getFractOfVariance` method is removed, use the new 
    :func:`~.calcFractVariance` function instead.
  
  * :meth:`Mode.getSqFlucts` method is removed, use :func:`~.calcSqFlucts` 
    function instead.

  * Renamed :func:`showFractOfVar` function as :func:`~.showFractVariance` 
    function instead.
    
  * Removed :func:`calcCumOverlapArray`, use :func:`~.calcCumulOverlap`
    with ``array=True`` argument instead.
    
  * :class:`~.AtomGroup` cannot be deformed by direct addition with a vector 
    instance.  

  * Unmapped atoms in :class:`~.AtomMap` instances are called dummies.  
    :meth:`AtomMap.numUnmapped` method, for example, is renamed as 
    :meth:`.AtomMap.numDummies`.

  * Renamed :func:`extrapolateModel` as :func:`~.extendModel`.

  * :func:`~.fetchPDBLigand` accepts only *filename* (instead of *save* and 
    *folder*) argument to save an XML file.

**Bugfix**:

  * A problem in distance based atom selection which would could cause problems
    when a distance based selection is made from a selection is fixed.
    
  * Changed :ref:`prody-blast` so that when a path for downloading files
    are given files are not save to local PDB folder. 
    

Release 0.9.4 (Feb 4, 2012)
===============================================================================

**Changes**:

  * :meth:`setAtomGroup` and :meth:`getAtomGroup` methods are renamed as 
    :meth:`.Ensemble.setAtoms` and  :meth:`.Ensemble.getAtoms`.
    
  * :class:`~.AtomGroup` class trajectory methods, i.e.
    :meth:`~.AtomGroup.setTrajectory`, 
    :meth:`~.AtomGroup.getTrajectory`, 
    :meth:`~.AtomGroup.nextFrame`,
    :meth:`~.AtomGroup.nextFrame`, and 
    :meth:`~.AtomGroup.gotoFrame` 
    methods are deprecated. Version 1.0 will feature a better integration
    of :class:`~.AtomGroup` and :class:`~.Trajectory` classes.


**Bugfix**:

  * Bugfixes in :meth:`~.Bond.setACSIndex`, :func:`~.saveAtoms`,
    and :meth:`~.HierView.getSegment`.
    
  * Bugfixes in :class:`~.GammaVariableCutoff` and 
    :class:`~.GammaStructureBased` classes.

  * Bugfix in :func:`~.calcCrossCorr` function.

  * Bugfixes in :meth:`~.EnsembleBase.getWeights`, 
    :func:`~.showOccupancies`, :meth:`~.DCDFile.flush`.
    
  * Bugfixes in ProDy commands :ref:`prody-blast`, :ref:`prody-fetch`, and
    :ref:`prody-pca`.

  * Bugfix in :func:`~.calcCenter` function.


Release 0.9.3 (Feb 1, 2012)
===============================================================================

**New Features**:

  * :class:`~.DBRef` class is implemented for storing references 
    to sequence databases parsed from PDB header records.
    
  * Methods for storing coordinate set labels in :class:`~.AtomGroup` 
    instances are implemented: :meth:`~.AtomGroup.getACSLabel`, and 
    :meth:`~.AtomGroup.getACSLabel`.

  * :func:`~.calcCenter` and :func:`~.moveAtoms` functions 
    are implemented for dealing with coordinate translation.

  * Hierarchical view, :class:`~.HierView`, is completely redesigned.  
    PDB files that contain non-empty segment name column (or when such 
    information is parsed from a PSF file), new design delicately handles this 
    information to identify distinct chains and residues.  This prevents 
    merging distinct chains in different segments but with same identifiers 
    and residues in those with same numbers.  New design is also using ordered 
    dictionaries :class:`collections.OrderedDict` and lists so that chain and 
    residue iterations yield them in the order they are parsed from file.  
    These improvements also bring modest improvements in speed.

  * :class:`~.Segment` class is implemented for handling segments
    of atoms defined in molecular dynamics simulations setup, using 
    :program:`psfgen` for example. 

  * Context manager methods are added to trajectory classes.  A trajectory
    file can be opened as follows::
       
      with Trajectory('mdm2.dcd') as traj:
          for frame in traj:
              calcGyradius(frame)

  * :class:`~.Chain` slicing is implemented::
    
      p38 = parsePDB('1p38')
      chA = p38['A']
      res_4to10 = chA[4:11]
      res_100toLAST = chA[100:]
  
  * Some support for bonds is implemented to :class:`~.AtomGroup` class.
    Bonds can be set using :meth:`~.AtomGroup.setBonds` method.  All 
    bonds must be set at once.  :meth:`~.AtomGroup.iterBonds` or
    :meth:`~.Atom.iterBonds` methods can be used to iterate over bonds
    in an AtomGroup or an Atom.
      
  * :func:`~.parsePSF` parses bond information and sets to the
    atom group. 

  * :meth:`.Selection.update` method is implemented, which may be useful to 
    update a distance based selection after coordinate changes.  
    
  * :func:`~.buildKDTree` and :func:`~.iterNeighbors` methods
    are implemented for facilitating identification of pairs of atoms that 
    are proximal.  
    
  * :meth:`~.AtomGroup.iterAtoms` method is implemented to all 
    :mod:`~prody.atomic` classes to provide uniformity for atom iterations.
    
  * :meth:`~.calcAngle`, :meth:`~.calcDihedral`, 
    :meth:`~.calcPhi`, :meth:`~.calcPsi`, and 
    :meth:`~.calcOmega` methods are implemented.

**Improvements**:

  * :meth:`~.Chain.getSelstr` methods of :class:`~.Chain` and 
    :class:`~.Residue` classes are improved to include the selection
    string of a :class:`~.Selection` when they are built using one.

**Changes**:

  * :class:`~.Residue` methods :meth:`~.Residue.getNumber`, 
    :meth:`~.Residue.setNumber`, :meth:`~.Residue.getName`,
    :meth:`~.Residue.setName` methods are deprecated and will be 
    removed in v1.0.

  * :class:`~.Chain` methods :meth:`~.Chain.getIdentifier` and 
    :meth:`~.Chain.setIdentifier` methods are deprecated and will be 
    removed in v1.0.

  * :class:`~.Polymer` attribute :attr:`~.Polymer.identifier`
    is renamed as :attr:`~.Polymer.chid`.
  * :class:`~.Chemical` attribute :attr:`~.Chemical.identifier`
    is renamed as :attr:`~.Chemical.resname`.
    
  * :meth:`getACSI` and :meth:`setACSI` are renamed as 
    :meth:`~.AtomGroup.getACSIndex` and  
    :meth:`~.AtomGroup.setACSIndex`, respectively.

  * :func:`calcRadiusOfGyration` is deprecated and will be removed
    in v1.0.  Use :func:`~.calcGyradius` instead.


**Bugfix**:

  * Fixed a problem in :func:`~.parsePDB` that caused loosing existing
    coordinate sets in an :class:`~.AtomGroup` when passed as *ag* 
    argument.
    
  * Fixed a problem with ``"same ... as ..."`` argument of :class:`~.Select`
    that selected atoms when followed by an incorrect atom selection.
    
  * Fixed another problem with ``"same ... as ..."`` which result in selecting
    multiple chains when same chain identifier is found in multiple segments
    or multiple residues when same residue number is found in multiple 
    segments.
    
  * Improved handling of negative integers in indexing :class:`~.AtomGroup` 
    instances.


Release 0.9.2 (Jan 11, 2012)
===============================================================================

**New Features**:

  * :program:`prody catdcd` command is implemented for concatenating and/or 
    slicing :file:`.dcd` files.  See :ref:`prody-catdcd` for usage examples.

  * :class:`~.DCDFile` can be opened in write or append mode, and 
    coordinate sets can be added using :meth:`~.DCDFile.write` method. 
  
  * :class:`~select.getReservedWords` can be used to get a list of words
    that cannot be used to label user data.

  * :class:`~prody.confProDy` function is added for configuring ProDy.
  
  * ProDy can optionally backup existing files with :file:`.BAK` (or another) 
    extension instead of overwriting them.  This behavior can be activated 
    using :class:`~prody.confProDy` function. 
    
**Improvements**:

  * :func:`~.writeDCD` file accepts :class:`~.AtomGroup` or 
    other :class:`~.Atomic` instances as *trajectory* argument.

  * :program:`prody align` command can be used to align multiple PDB structures.
  
  * :program:`prody pca` command allows atom selections for DCD files that are 
    accompanied with a PDB or PSF file.

**Changes**:

  * :class:`~.DCDFile` instances, when closed, raise exception, similar
    to behavior of :class:`file` objects in Python. 
    
  * Title of :class:`~.AtomGroup` instances resulting from copying an 
    :class:`~.Atomic` instances does not start with 'Copy of'.
    
  * :func:`changeVerbosity` and :func:`getVerbosityLevel` are renamed as 
    :func:`~.setVerbosity` and :func:`~.getVerbosity`, respectively. 
    Old names will be removed in v1.0.
    
  * ProDy routines (commands) module is rewritten to use new :mod:`argparse`
    module. See :ref:`commands` for details of changes.
    
  * :mod:`argparse` module is added to the package for Python versions 2.6
    and older.
    

**Bugfix**:

  * Fixed problems in :func:`~.loadAtoms` and :func:`~.saveAtoms`
    functions.

  * Bugfixes in :func:`~.parseDCD` and :func:`~.writeDCD`
    functions for Windows compatability.


Release 0.9.1 (Nov 9, 2011)
===============================================================================

**Bug Fixes**:

  * Fixed problems with reading and writing configuration files.
  * Fixed problem with importing nose for testing.

Release 0.9 (Nov 8, 2011)
===============================================================================

**New Features**:

  * `PDBML <http://pdbml.pdb.org/>`_ and `mmCIF <http://mmcif.pdb.org/>`_ files
    can be retrieved using :func:`~.fetchPDB` function.

  * :func:`~.getPDBLocalFolder` and :func:`~.setPDBLocalFolder`
    functions are implemented for local PDB folder management.
    
  * :func:`~.parsePDBHeader` is implemented for convenient parsing of
    header data from :file:`.pdb` files.
    
  * :func:`~.showProtein` is implemented to allow taking a quick look
    at protein structure.

  * :class:`~.Chemical` and :class:`~.Polymer` classes are 
    implemented for storing chemical and polymer component data parsed from
    PDB header records.


**Changes**:

  .. warning::  This release introduces numerous changes in method and function
     names all aiming to improve the interactive usage experience.  All changes
     are listed below.  Currently these functions and methods are present in
     both old and new names, so code using ProDy must not be affected.  Old
     function names will be removed from version 1.0, which is expected to 
     happen late in the first quarter of 2012.
      
     Old function names are marked as deprecated, but ProDy will not issue any
     warnings until the end of 2011.  In 2012, ProDy will will automatically 
     start issuing :class:`DeprecationWarning` upon calls using old names to 
     remind the user of the name change.  
     
     For deprecated methods that are present in multiple classes, only the 
     affected modules are listed for brevity.
     
  .. note::  When modifying code using ProDy to adjust the name changes, 
     turning on deprecation warnings may help locating all use cases of the 
     deprecated names.  See :meth:`~.turnonDepracationWarnings` for this
     purpose.

  **Functions**:
  
  The following function name changes are mainly to reduce the length of the 
  name in order to make them more suitable for interactive sessions:

  ========================================  =====================================
  Old name                                  New name
  ========================================  =====================================
  :func:`applyBiomolecularTransformations`  :func:`~.buildBiomolecules` 
  :func:`assignSecondaryStructure`          :func:`~.assignSecstr`
  :func:`scanPerturbationResponse`          :func:`~.calcPerturbResponse`
  :func:`calcCrossCorrelations`             :func:`~.calcCrossCorr`
  :func:`calcCumulativeOverlap`             :func:`~.calcCumulOverlap`
  :func:`calcCovarianceOverlap`             :func:`~.calcCovOverlap`
  :func:`showFractOfVariances`              :func:`~.showFractVars`
  :func:`showCumFractOfVariances`           :func:`~.showCumulFractVars`
  :func:`showCrossCorrelations`             :func:`~.showCrossCorr`
  :func:`showCumulativeOverlap`             :func:`~.showCumulOverlap`
  :func:`deform`                            :func:`~.deformAtoms`
  :func:`calcSumOfWeights`                  :func:`~.calcOccupancies`
  :func:`showSumOfWeights`                  :func:`~.showOccupancies`
  :func:`trimEnsemble`                      :func:`~.trimPDBEnsemble`
  :func:`getKeywordResidueNames`            :func:`~.getKeywordResnames`
  :func:`setKeywordResidueNames`            :func:`~.setKeywordResnames`
  :func:`getPairwiseAlignmentMethod`        :func:`~.getAlignmentMethod`
  :func:`setPairwiseAlignmentMethod`        :func:`~.setAlignmentMethod`
  :func:`getPairwiseMatchScore`             :func:`~.getMatchScore`
  :func:`setPairwiseMatchScore`             :func:`~.setMatchScore`
  :func:`getPairwiseMismatchScore`          :func:`~.getMismatchScore`
  :func:`setPairwiseMismatchScore`          :func:`~.setMismatchScore`
  :func:`getPairwiseGapOpeningPenalty`      :func:`~.getGapPenalty`
  :func:`setPairwiseGapOpeningPenalty`      :func:`~.setGapPenalty`
  :func:`getPairwiseGapExtensionPenalty`    :func:`~.getGapExtPenalty`
  :func:`setPairwiseGapExtensionPenalty`    :func:`~.setGapExtPenalty`
  ========================================  =====================================

  **Coordinate methods**:
  
  All :meth:`getCoordinates` and :meth:`setCoordinates` methods in 
  :mod:`~prody.atomic` and :mod:`~prody.ensemble` classes are renamed as 
  :meth:`getCoords` and :meth:`setCoords`, respectively.
  
  ``getNumOf`` **methods**:
  
  All method names starting with ``getNumOf`` now start with ``num``.  This
  change brings two advantages: method names (i) are considerably shorter, 
  and (ii) do not suggest that there might also be corresponding ``set``
  methods. 
  
  ============================  ====================  =========================
  Old name                      New name              Affected modules
  ============================  ====================  =========================
  :meth:`getNumOfAtoms`         :meth:`numAtoms`      :mod:`~prody.atomic`, 
                                                      :mod:`~prody.ensemble`, 
                                                      :mod:`~prody.dynamics`
  :meth:`getNumOfChains`        :meth:`numChains`     :mod:`~prody.atomic`
  :meth:`getNumOfConfs`         :meth:`numConfs`      :mod:`~prody.ensemble`
  :meth:`getNumOfCoordsets`     :meth:`numCoordsets`  :mod:`~prody.atomic`, 
                                                      :mod:`~prody.ensemble`
  :meth:`getNumOfDegOfFreedom`  :meth:`numDOF`        :mod:`~prody.dynamics`
  :meth:`getNumOfFixed`         :meth:`numFixed`      :mod:`~prody.ensemble`
  :meth:`getNumOfFrames`        :meth:`numFrames`     :mod:`~prody.ensemble`
  :meth:`getNumOfResidues`      :meth:`numResidues`   :mod:`~prody.atomic`
  :meth:`getNumOfMapped`        :meth:`numMapped`     :mod:`~prody.atomic`
  :meth:`getNumOfModes`         :meth:`numModes`      :mod:`~prody.dynamics`
  :meth:`getNumOfSelected`      :meth:`numSelected`   :mod:`~prody.ensemble`
  :meth:`getNumOfUnmapped`      :meth:`numUnmapped`   :mod:`~prody.atomic`
  ============================  ====================  =========================
    
  ``getName`` **method**:
  
  :meth:`getName` methods are renamed as :meth:`getTitle` to avoid confusions 
  that might arise from changes in :mod:`~prody.atomic` method names listed 
  below.  All classes in :mod:`~prody.atomic`, :mod:`~prody.ensemble`, and 
  :mod:`~prody.dynamics` are affected from this change. 
  
  In line with this change, :func:`~.parsePDB` and 
  :func:`~.parsePQR` *name* arguments are changed to *title*, but 
  *name* argument will also work until release 1.0.
  
  This name change conflicted with :meth:`getTitle` method of :class:`~
  .DCDFile`.  The conflict is resolved in favor of the general :meth:`getTitle`
  method.  An alternative method will be implemented to handle title strings 
  in :file:`DCD` files.   
  
  ``get/set`` **methods of atomic classes**:
  
  Names of ``get`` and ``set`` methods allowing access to atomic data are all
  shortened as follows:
  
  ===========================  =======================
  Old name                     New name
  ===========================  =======================
  :meth:`getAtomNames`	       :meth:`getNames`
  :meth:`getAtomTypes`	       :meth:`getTypes`
  :meth:`getAltLocIndicators`  :meth:`getAltlocs`
  :meth:`getAnisoTempFactors`  :meth:`getAnisos`
  :meth:`getAnisoStdDevs`      :meth:`getAnistds`
  :meth:`getChainIdentifiers`  :meth:`getChains`
  :meth:`getElementSymbols`    :meth:`getElements`
  :meth:`getHeteroFlags`       :meth:`getHeteros`
  :meth:`getInsertionCodes`    :meth:`getIcodes`
  :meth:`getResidueNames`      :meth:`getResnames`
  :meth:`getResidueNumbers`    :meth:`getResnums`
  :meth:`getSecondaryStrs`     :meth:`getSecstrs`            
  :meth:`getSegmentNames`      :meth:`getSegnames`
  :meth:`getSerialNumbers`     :meth:`getSerials`
  :meth:`getTempFactors`	   :meth:`getBetas`
  ===========================  =======================
                       
  This change affects all :mod:`~prody.atomic` classes, 
  :class:`~.AtomGroup`, :class:`~.Atom`, :class:`~.Chain`, 
  :class:`~.Residue`, :class:`~.Selection` and 
  :class:`~.AtomMap`.     


  **Other changes in atomic methods**:
  
  * :meth:`getSelectionString` renamed as :meth:`getSelstr`

  Methods handling user data (which was previously called attribute) are 
  renamed as follows: 

  ====================  =======================
  Old name              New name
  ====================  =======================
  :meth:`getAttribute`  :meth:`getData`
  :meth:`getAttrNames`  :meth:`getDataLabels`
  :meth:`getAttrType`   :meth:`getDataType`
  :meth:`delAttribute`  :meth:`delData`
  :meth:`isAttribute`	:meth:`isData`
  :meth:`setAttribute`  :meth:`setData`
  ====================  =======================
  
  **To be removed**:
  
  Finally, the following methods will be removed, but other suitable methods
  are overloaded to perform their action:
  
  * removed :meth:`AtomGroup.getBySerialRange`, overloaded 
    :meth:`.AtomGroup.getBySerial`
  * removed :func:`~.getProteinResidueNames`, overloaded
    :func:`~select.getKeywordResnames` 
  * removed :func:`~.setProteinResidueNames`, overloaded
    :func:`~select.setKeywordResnames` 


**Scripts**:
  
  The way ProDy scripts work has changed. See :ref:`commands` for details.
  Using older scripts will start issuing deprecation warnings in 2012.

**Bug Fixes**:

  * Bugs in :func:`~.execDSSP` and :func:`~.execSTRIDE`
    functions that caused exceptions when compressed files were passed 
    is fixed.
  
  * A problem in scripts for PCA of DCD files is fixed. 
    
  
:ref:`nmwiz`
-------------------------------------------------------------------------------

Development of NMWiz is finalized and it will not be distributed in the ProDy
installation package anymore.  See :ref:`nmwiz` pages for instructions on
installing it. 


Release 0.8.3 (Oct 16, 2011)
===============================================================================

**New Features**:

  * Functions to read and write PQR files: :func:`~.parsePQR` and
    :func:`~.writePQR`.
    
  * Added :meth:`~.PDBEnsemble.getIdentifiers` method that returns
    identifiers of all conformations in the ensemble. 

  * ProDy tests are incorporated to the package installer.  If you are using 
    Python version 2.7, you can run the tests by calling ``prody.test()``.


**Improvements**:

  * :func:`~.blastPDB` function and :class:`~.PDBBlastRecord`
    class are rewritten to use faster and more compact code.
    
  * New :class:`PackageLogger` function is implemented to unify logging
    and reporting task progression.
  
  * Improvements in PDB ensemble support functions, e.g. 
    :func:`~.trimPDBEnsemble`, are made.
   
  * Improvements in ensemble concatenations are made.

**Bug Fixes**:

  * Bugfixes in :func:`~.PDBEnsemble` slicing operation.  This may
    have affected users when slicing a PDB ensemble for plotting projections
    in color for different forms of the protein.
    
Release 0.8.2 (Oct 14, 2011)
===============================================================================

**New Features**:

  * :func:`~.fetchPDBClusters`, :func:`~.loadPDBClusters`, and
    :func:`~.getPDBCluster` functions are implemented for handling
    PDB sequence cluster data. These functions can be used instead of 
    :func:`~.blastPDB` function for fast access to structures of 
    the same protein (at 95% sequence identity level) or similar proteins.
    
  * Perturbation response scanning method described in [CA09]_ is implemented 
    as :func:`~.scanPerturbationResponse` based on the code provided 
    by Ying Liu.
  

**Changes**:

  * :func:`~.fetchLigandData` returns the URL of the XML file in the
    ligand data dictionary.
    
  * Name of the ProDy configuration file in user :file:`home` directory 
    is renamed as :file:`.prodyrc` (used to be :file:`.prody`).
    
  * :func:`~.applyBiomolecularTransformations` and
    :func:`~.assignSecondaryStructure` functions raise 
    :class:`ValueError` when the function fails to perform its action
    due to missing data in header dictionary.
    
  * :func:`~.fetchPDB` decompresses PDB files found in the working
    directory when user asks for decompressed files.

  * :func:`~.parsePDB` appends *chain* and *subset* arguments to
    :func:`~.AtomGroup` name.

  * *chain* argument is added to :meth:`~.PDBBlastRecord.getHits`.

**Improvements**:

  * Atom selection class :class:`~select.Select` is completely redesigned
    to prevent breaking of the parser when evaluating invalid selection
    strings.
    
  * Improved type checking in :func:`~.parsePDB` function.
    
**Bug Fixes**:

  * Bugfixes in :func:`~.parseDSSP`: one emerged problems in lines
    indicating chain breaks, another did not parse bridge-partners correctly.
    Both fixes are contributed by Kian Ho.
    
  * Bugfix in :func:`~.parsePDB` function. When only header is desired
    (``header=True, model=0``), would return a tuple containing an empty 
    atom group and the header.

**Developmental**:

  * Unit tests for :mod:`~prody.proteins` and :mod:`~prody.select` modules are 
    developed. 

Release 0.8.1 (Sep 16, 2011)
===============================================================================

**New Features**:

  * :func:`~.fetchLigandData` is implemented for fetching ligand data
    from Ligand Expo.

  * :func:`~.parsePSF` function is implemented for parsing X-PLOR 
    format PSF files.
    
**Changes**:

  * __slots__ is used in :class:`~.AtomGroup` and :class:`~.Atomic`
    classes. This change prevents user from assigning new variables to 
    instances of all classes derived from the base :class:`~.Atomic`.
    
  * :mod:`pyparsing` is updated to version 1.5.6.

**Bug Fixes**:

  * A bug in :meth:`~.AtomGroup.copy` method is fixed. When AtomGroup
    instance itself is copied, deep copies of data arrays were not made.

  * A bug in :class:`~select.Select` class raising exceptions when negative
    residue number values are present is fixed.

  * Another bug in :class:`~select.Select` class misinterpreting
    ``same residue as ...`` statement when specific chains are involved is
    fixed.

  * A bug in :meth:`~.AtomGroup.addCoordset` method duplicating 
    coordinates when no coordsets are present in the instance is fixed.

:ref:`nmwiz`
-------------------------------------------------------------------------------

**Changes**:

  * Version number in main window is iterated.
  
  * Mode graphics material is stored for individual modes.

  * Mode scaling factor is printed when active mode or RMSD is changed.

  * All selections are deleted to avoid memory leaks.


Release 0.8 (Aug 24, 2011)
===============================================================================

.. note::
   After installing v0.8, you may need to make a small change in your 
   existing scripts. If you are using :class:`~.Ensemble` class 
   for analyzing PDB structures, rename it as :class:`~.PDBEnsemble`.
   See the other changes that may affect your work below and the class 
   documentation for more information.
    

**New Features**:

  * :class:`~.DCDFile` is implemented for handling DCD files.
    
  * :class:`~.Trajectory` is implemented for handling multiple 
    trajectory files.
  
  * :func:`~.writeDCD` is implemented for writing DCD files.
    
  * :ref:`trajectory` example to illustrate usage of new classes for handling
    DCD files. :ref:`eda` example is updated to use new ProDy classes.
  
  * :class:`~.PCA` supports :class:`~.Trajectory` and
    :class:`~.DCDFile` instances. 
  
  * :class:`~.Ensemble` and :class:`~.PDBEnsemble` classes
    can be associated with :class:`~.AtomGroup` instances. This allows
    selecting and evaluating coordinates of subset of atoms. See
    :meth:`~.EnsembleBase.setAtomGroup`, 
    :meth:`~.EnsembleBase.select`,
    :meth:`~.EnsembleBase.getAtomGroup`, and 
    :meth:`~.EnsembleBase.getSelection` methods.
  
  * :func:`~.execDSSP`, :func:`~.parseDSSP`, and 
    :func:`~.performDSSP` functions are implemented for 
    executing and parsing DSSP calculations.

  * :func:`~.execSTRIDE`, :func:`~.parseSTRIDE`, and 
    :func:`~.performSTRIDE` functions are implemented for 
    executing and parsing DSSP calculations.

  * :func:`~.parsePDB` function parses atom serial numbers. Atoms
    can be retrieved from an :class:`~.AtomGroup` instance by their
    serial numbers using :meth:`~.AtomGroup.getBySerial` and
    :meth:`~.AtomGroup.getBySerialRange` methods. 
    
  * :func:`~.calcADPs` function can be used to calculate anisotropic
    displacement parameters for atoms with anisotropic temperature factor
    data.
    
  * :meth:`~.Ensemble.getRMSFs` is implemented for calculating
    root mean square fluctuations.
    
  * :class:`~.AtomGroup` and :class:`~.Mode` or 
    :class:`~.Vector` additions are supported. This adds a new
    coordinate set to the :class:`~.AtomGroup` instance.
    
  * :meth:`~.AtomGroup.getAttrNames` is implemented for listing
    user set attribute names.
    

**Improvements**:

  * :func:`~.calcProjection`, :func:`~.showProjection`, 
    and :func:`~.showCrossProjection` functions can optionally 
    calculate/display RMSD along the normal mode. 
  
  * ANM, GNM, and PCA routines can optionally write compressed ProDy data files.
  
  * :func:`~.fetchPDB` function can optionally write decompressed 
    files and force copying a file from local mirror to target folder.

  * :meth:`~.PCA.buildCovariance` and :meth:`~.PCA.performSVD` 
    methods accept Numpy arrays as coordinate sets.

  * Performance of :meth:`~.PCA.buildCovariance` method is optimized 
    for evaluation of PDB ensembles.
  
  * :func:`~.calcRMSD` and :func:`~.superpose` functions are
    optimed for speed and memory usage.
  
  * :meth:`~.Ensemble.getMSFs` is optimized for speed and memory usage.

  * Improvements in memory operations in :mod:`~prody.atomic`, 
    :mod:`~prody.ensemble`, and :mod:`~prody.dynamics` modules for 
    faster data (PDB/NMD) output.
    
  * Optimizations in :class:`~select.Select` and :class:`~select.Contacts`
    classes.

**Changes**:

  * :class:`~.Ensemble` does not store conformation names. Instead, 
    newly implemented :class:`~.PDBEnsemble` class stores identifiers 
    for individual conformations (PDB IDs). This class should be used in cases 
    where source of individual conformations is important.
    
  * :func:`~.calcProjection`, :func:`~.showProjection`, and
    :func:`~.showCrossProjection` function calculate/display 
    root mean square deviations, by default. 
    
  * Oxidized cysteine residue abbreviation ``CSO`` is added to the definition
    of ``protein`` keyword.
    
  * :meth:`getMSF` method is renamed as :meth:`~.Ensemble.getMSFs`.
  
  * :func:`~.parseDCD` function returns :class:`~.Ensemble`
    instances.

**Bug Fixes**:
  
  * A bug in :mod:`~prody.select` module causing exceptions when regular 
    expressions are used is fixed.
  
  * Another bug in :mod:`select` module raising exception when "(not ..," is
    passed is fixed.

  * Various bugfixes in :mod:`~prody.ensemble` module.
  
  * Problem in :file:`fetchpdb.py` that occurred when a file is found in a
    local mirror is fixed.
    
  * Bugfix in :meth:`~.AtomPointer.copy` method.

:ref:`nmwiz`
-------------------------------------------------------------------------------

**New Features**:

  * NMWiz can be used to compare two structures by calculating and depicting
    structural changes.

  * Arrow graphics is scaled based on a user specified RMSD value.

**Improvements**:

  * NMWiz writes DCD format trajectories for PCA using ProDy. This provides
    significant speed up in cases where IO rate is the bottleneck.

**Changes**:

  * Help is provided in a text window to provide a cleaner GUI.


Release 0.7.2 (Jun 21, 2011)
===============================================================================

**New Features**:

  * :func:`~.parseDCD` is implemented for parsing coordinate sets 
    from DCD files.


**Improvements**:

  * :func:`~.parsePDB` parses ``SEQRES`` records in header sections.

**Changes**:

  * Major classes can be instantiated without passing a name argument.
  * Default selection in NMWiz ProDy interface is changed to ensure selection
    only protein Cα atoms.


**Bug Fixes**:

  * A bug in :func:`~.writeNMD` function causing problems when writing
    a single mode is fixeed.
  * Other bugfixes in :mod:`~prody.dynamics` module functions.


Release 0.7.1 (Apr 28, 2011)
===============================================================================

**Highlights**:

  * :class:`~.Atomic` :meth:`__getattribute__` is overloaded to interpret
    atomic selections following the dot operator. For example,
    ``atoms.calpha`` is interpreted as ``atoms.select('calpha')``. See
    :ref:`` for more details.
  * :class:`~.AtomGroup` class is integrated with 
    :class:`~.HierView` class. Atom group instances now can be indexed
    to get chains or residues and number of chains/residues can be retrieved.
    A hierarchical view is generated and updated when needed. See
    :ref:`` for more details.
      
     

**New Features**:

  * :func:`~.matchAlign` is implemented for quick alignment of protein
    structures. See :ref:`extract-ligands` usage example.
  * :meth:`~.AtomGroup.setAttribute`, 
    :meth:`~.AtomGroup.getAttribute`,
    :meth:`~.AtomGroup.delAttribute`, and
    :meth:`~.AtomGroup.isAttribute` functions are implemented for
    :class:`~.AtomGroup` class to facilitate storing user provided 
    atomic data. See :ref:`attributes` example.
  * :func:`~.saveAtoms` and :func:`~.loadAtoms` functions 
    are implemented to allow for saving atomic data and loading it 
    This saves custom atomic attributes and much faster than parsing
    data from PDB files. 
  * :func:`~.calcCollectivity` function is implemented to allow
    for calculating collectivity of deformation vectors.
    
**Improvements**:

  * :func:`~.parsePDB` can optionally return biomolecule when 
    ``biomol=True`` keyword argument is passed.
  * :func:`~.parsePDB` can optionally make secondary structure
    assignments when ``secondary=True`` keyword argument is passed.
  * :func:`~.calcSqFlucts` function is changed to accept 
    :class:`~.Vector` instances, e.g. deformation vectors.

**Changes**:

  * Changes were made in :func:`~.calcADPAxes` function to follow 
    the conventions in analysis ADPs. See its documentation.

**Bug Fixes**:

  * A in :class:`~.Ensemble` slicing operations is fixed. Weights are
    now copied to the new instances obtained by slicing.
  * Bug fixes in :mod:`~prody.dynamics` plotting functions 
    :func:`~.showScaledSqFlucts`, :func:`~.showNormedSqFlucts`, 

Release 0.7 (Apr 4, 2011)
===============================================================================

**New Features**:

  * Regular expressions can be used in atom selections. See 
    :mod:`~prody.select` module for details.
  
  * User can define selection macros using :func:`~select.defSelectionMacro`
    function. Macros are saved in ProDy configuration and loaded in later
    sessions. See :mod:`~prody.select` module for other related functions.
  
  * :func:`~.parseSparseMatrix` function is implemented for parsing
    matrices in sparse format. See the usage example in :ref:`external-matrix`. 

  * :func:`~.deform` function is implemented for deforming coordinate 
    sets along a normal mode or linear combination of multiple modes. 

  * :func:`~.sliceModel` function is implemented for slicing normal
    mode data to be used with functions calculating atomic properties using
    normal modes. 

**Improvements**:

  * Atom selections using bare keyword arguments is optimized. New keyword
    definitions are added. See :mod:`~prody.select` module for the complete 
    list.
  
  * A new keyword argument for :func:`~.calcADPAxes` allows for
    comparing largest axis to the second largest one.

**Changes**:

  * There are changes in function used to alter definitions of selection
    keywords. See :mod:`~prody.select` for details.
    
  * :func:`~.assignSecondaryStructure` function assigns SS identifiers
    to all atoms in a residue. Residues with no SS information specified is
    assigned coil conformation. 
  
  * When :func:`~.Ensemble` and :func:`~.NMA` classes are 
    instantiated with an empty string, instances are called "Unnamed".

  * :func:`~.sliceMode`, :func:`~.sliceVector` and
    :func:`~.reduceModel` functions return the atom selection 
    in addition to the sliced vector/mode/model instance.

**Bug Fixes**:

  * Default selection for :func:`~.calcGNM` function is set to
    "calpha".

:ref:`nmwiz`
-------------------------------------------------------------------------------

**New Features**:

  * NMWiz supports GNM data and can use ProDy for GNM calculations.
  * NMWiz can gather normal mode data from molecules loaded into VMD.
    This allows NMWiz to support all formats supported by VMD.
  * User can write data loaded into NMWiz in NMD format.
  * An Arrow Graphics option allows the user to draw arrows in both directions.
  * User can select Licorice representation for the protein if model is an
    all atom mode.
  * User can select Custom as the representation of the protein to prevent
    NMWiz from chancing a user set representation.
  * Trace is added as a protein backbone representation option.

**Improvements**:

  * NMWiz remembers all adjustments on arrow graphics for all modes.
  * Plotting :guilabel:`Clear` button clears only atom labels that are 
    associated with the dataset.
  * Removing a dataset removes all associated molecule objects.
  * Selected atom representations are turned on based on atom index. 
  * Padding around interface button has been standardized to provide a uniform
    experience between different platforms.

Release 0.6.2 (Mar 16, 2011)
===============================================================================

**New Features**:

  * :meth:`~.PCA.performSVD` function is implemented for faster
    and more memory efficient principal compoment analysis.
  * :func:`~.extrapolateModel` function is implemented for 
    extrapolating a coarse-grained model to an all atom model. See the 
    usage example :ref:`extendmodel`.
  * :func:`plog` is implemented for enabling users to make log entries.

**Improvements**:

  * :mod:`~prody.compare` functions are improved to handle insertion codes.
  * :class:`~.HierView` allows for indexing using chain identifier
    and residue numbers. See usage example :ref:`hierview`
  * :class:`~.Chain` allows for indexing using residue number and
    insertion code. See usage example :ref:`hierview`.
  * :meth:`~.AtomGroup.addCoordset` function accepts 
    :class:`~.Atomic` and :class:`~.Ensemble` instances
    as *coords* argument.
  * New method :meth:`~.HierView.getAtoms` is implemented.
  * :class:`~.AtomGroup` set functions check the correctness of 
    dimension of data arrays to prevent runtime problems.
  * :file:`pca.py` script is updated to use the faster PCA method
    that uses SVD.

**Changes**:

  * "backbone" definition now includes the backbone hydrogen atom 
    (Thanks to Nahren Mascarenhas for pointing to this discrepancy in the
    keyword definition). 

**Bug Fixes**:

  * A bug in :class:`~.PCA` allowed calculating covariance matrix
    for less than 3 coordinate sets is fixed.
  * A bug in :func:`~.mapOntoChain` function that caused problems
    when mapping all atoms is fixed.
    
    

Release 0.6.1 (Mar 2, 2011)
===============================================================================

**New Features**:

  * :func:`~.setWWPDBFTPServer` and :func:`~.getWWPDBFTPServer` 
    functions allow user to change or learn the WWPDB FTP server that ProDy 
    uses to download PDB files. Default server is RCSB PDB in USA. 
    User can change the default server to one in Europe or Japan.
  * :func:`~.setPDBMirrorPath` and :func:`~.getPDBMirrorPath` 
    functions allow user to specify or learn the path to a local PDB mirror.
    When specified, a local PDB mirror is preferred for accessing PDB files,
    over downloading them from FTP servers.
  * :func:`~.mapOntoChain` function is improved to map backbone or 
    all atoms.

**Improvements**:

  * :class:`~.WWPDB_PDBFetcher` can download PDB files from different
    WWPDB FTP servers.
  * :class:`~.WWPDB_PDBFetcher` can also use local PDB mirrors for
    accessing PDB files.

**Changes**:

  * :class:`RCSB_PDBFetcher` is renamed as :class:`~.WWPDB_PDBFetcher`.
  * :func:`~.mapOntoChain` and :func:`~.matchChains` functions
    accept ``"ca"`` and ``"bb"`` as *subset* arguments.
  * Definition of selection keyword "protein" is updated to include
    some non-standard amino acid abbreviations. 

**Bug Fixes**:

  * A bug in :class:`~.WWPDB_PDBFetcher` causing exceptions when
    non-string items passed in a list is fixed.
  * An important bug in :func:`~.parsePDB` is fixed. When parsing
    backbone or Cα atoms, residue names were not checked and this caused
    parsing water atoms with name ``"O"`` or calcium ions with name ``"CA"``.
    

Release 0.6 (Feb 22, 2011)
===============================================================================

**New Features**:

  * Biopython module pairwise2 and packages KDTree and Blast are incorporated
    in ProDy package to make installation easier. Only NumPy needs to be 
    installed before ProDy can be used. For plotting, Matplotlib is still 
    required.
  * :ref:`nmwiz` is distributed with ProDy source. On Linux, if VMD is 
    installed, ProDy installer locates VMD plugins folder and installs NMWiz.
    On Windows, user needs to follow a separate set of instructions (see
    :ref:`nmwiz`).
  * :class:`~.Gamma` class is implemented for facilitating use of  
    force constants based on atom type, residue type, or property. An
    example derived classes are :class:`~.GammaStructureBased` and 
    :class:`~.GammaVariableCutoff`.
  * :func:`~.calcTempFactors` function is implemented to 
    calculate theoretical temperature factors.
  * 5 new :ref:`commands` are implemented, and existing scripts are improved to
    output figures.
  * :meth:`~.NMABase.getModel` method is implemented to 
    make function development easier.
  * :func:`~.resetTicks` function is implemented to change X and/or Y
    axis ticks in plots when there are discontinuities in the plotted data. 

**Improvements**:

  * :meth:`~.ANM.buildHessian` and :meth:`~.GNM.buildKirchhoff`
    classes are improved to accept :class:`~.Gamma` instances
    or other custom function as *gamma* argument. See also :ref:`gamma`.
  * :class:`~select.Select` class is changed to treat single word keywords
    differently, e.g. "backbone" or "protein". 
    They are interpreted 10 times faster and in use achieve much higher 
    speed-ups when compared to composite selections. For example, using the 
    keyword "calpha" instead of the ``name CA and protein``,
    which returns the same selection, works >20 times faster. 
  * Optimizations in :class:`~select.Select` class to increase 
    performance (Thanks to Paul McGuire for providing several Pythonic tips
    and Pyparsing specific advice).
  * :func:`~.applyBiomolecularTransformations` function is improved
    to handle large biomolecular assemblies.
  * Performance optimizations in :func:`~.parsePDB` and other 
    functions.
  * :class:`~.Ensemble` class accepts :class:`.Atomic` 
    instances and automatically adds coordinate sets to the ensemble. 
  
**Changes**:
 
  * :class:`PDBlastRecord` is renamed as :class:`~.PDBBlastRecord`. 
  * :class:`~.NMABase` instances can be index using a list or tuple of 
    integers, e.g. ``anm[1,3,5]``.
  * "ca", "bb", and "sc" keywords are defined as short-hands 
    for "calpha", "backbone", and "sidechain", respectively.
  * Behavior of :func:`~.calcANM` and :func:`~.calcGNM` 
    functions have changed. They return the atoms used for calculation as well.

**Bug Fixes**:
    
  * A bug in :func:`~.assignSecondaryStructure` function is fixed.
  * Bug fixes in :ref:`prody-anm` and :ref:`prody-gnm`.
  * Bug fixes in :func:`~.showSqFlucts` and 
    :func:`~.showProjection` functions.
    
:ref:`nmwiz`
-------------------------------------------------------------------------------

  * NMWiz can be used as a graphical interface to ProDy. ANM or PCA 
    calculations can be performed for molecules that are loaded in VMD.
  * User can set default color for arrow graphics and paths to
    ANM and PCA scripts.
  * Optionally, NMWiz can preserve the current view in VMD display window when
    loading a new dataset. Check the box in the NMWiz GUI main window.
  * A bug that prevented selecting residues from plot window is fixed.

Release 0.5.3 (Feb 11, 2011)
===============================================================================

**New Features**:

  * Membership, equality, and non-equality test operation are defined for all 
    :mod:`~prody.atomic` classes. See :ref:`selection-operations`.
  * Two functions are implemented for dealing with anisotropic temperature 
    factors: :func:`~.calcADPAxes` and :func:`~.buildADPMatrix`.
  * :meth:`~.NMA.setEigens` and :meth:`~.NMA.addEigenpair` 
    methods are implemented to assist analysis of normal modes calculated using 
    external software.
  * :func:`~.parseNMD` is implemented for parsing NMD files.
  * :func:`~.parseModes` is implemented for parsing normal mode data.
  * :func:`~.parseArray` is implementing for reading numeric data,
    particularly normal mode data calculated using other software for analysis
    using ProDy. 
  * The method in [BH02]_ to calculate overlap between covariance matrices is 
    implemented as :func:`~.calcCovOverlap` function.
  * A method to trim :class:`~.Ensemble` instances is implemented: 
    :func:`~.trimEnsemble`
  * :func:`checkUpdates` to check for ProDy updates is implemented.

**Changes**:

  * Change in default behavior of :func:`~.parsePDB` function. When
    alternate locations exist, those indicated by A are parsed. For parsing
    all alternate locations user needs to pass ``altloc=True`` argument.    
  * :func:`getSumOfWeights` is renamed as :func:`~.calcSumOfWeights`.
  * :func:`mapAtomsToChain` is renamed as :func:`~.mapOntoChain`.
  * :func:`ProDyStartLogFile` is renamed as :func:`startLogfile`.
  * :func:`ProDyCloseLogFile` is renamed as :func:`closeLogfile`.
  * :func:`ProDySetVerbosity` is renamed as :func:`changeVerbosity`.

**Improvements**:

  * A few bugs in ensemble and dynamics classes are fixed.
  * Improvements in :class:`~.RCSB_PDBFetcher` allow it not to miss a 
    PDB file if it exists in the target folder.
  * :func:`~.writeNMD` is fixed to output B-factors (Thanks to 
    Dan Holloway for pointing it out).

Release 0.5.2 (Jan 12, 2011)
===============================================================================

**Bug Fixes**:
  
  * An important fix in :func:`~.sampleModes` function was made
    (Thanks to Alberto Perez for finding the bug and suggesting a solution). 
    
**Improvements**:
  
  * Improvements in :meth:`.ANM.calcModes`, :meth:`.GNM.calcModes`, 
    and :meth:`.PCA.calcModes` methods prevent Numpy/Scipy throwing an
    exception when more than available modes are requested by the user.
  * Improvements in :func:`.blastPDB` enable ProDy throw an 
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

  * :class:`~.AtomPointer` base class for classes pointing to
    atoms in an :class:`~.AtomGroup`.
  * :class:`~.AtomPointer` instances (Selection, Residue, etc.)
    can be added. See :ref:`selection-operations` for examples.
  * :meth:`.Select.getIndices` and :meth:`.Select.getBoolArray` 
    methods to expand the usage of :class:`~select.Select`.
  * :func:`~.sliceVector` and :func:`~.sliceMode` functions.
  * :func:`~.saveModel` and :func:`~.loadModel` functions
    for saving and loading NMA data.
  * :func:`~.parsePDBStream` can now parse specific chains or
    alternate locations from a PDB file.
  * :func:`~.alignCoordsets` is implemented to superimpose
    coordinate sets of an :class:`~.AtomGroup` instance.

**Bug Fixes**:

  * A bug in :func:`.parsePDBStream` that caused unidentified errors 
    when a model in a multiple model file did not have the same number of 
    atoms is fixed.

**Changes**:

  * Iterating over a :class:`~.Chain` instance yields :class:`~.Residue` 
    instances.
  * :class:`~.Vector` instantiation requires an *array* only. *name* 
    is an optional argument.
  * Functions starting with ``get`` and performing a calculations are renamed to
    start with ``calc``, e.g. :func:`getRMSD` is renamed to
    :func:`~.calcRMSD`

Release 0.2 (Nov 16, 2010)
===============================================================================

**Important Changes**:


  * Single word keywords *not* followed by "and" logical operator are not 
    accepted, e.g. "protein within 5 of water" will raise an SelectionError, 
    use "protein and within 5 of water" instead.
  * :func:`findMatchingChains` is renamed to 
    :func:`~.matchChains`.
  * :func:`showOverlapMatrix` is renamed to 
    :func:`~.showOverlapTable`.
  * Modules are reorganized.

**New Features**:

  * :class:`~.Atomic` for easy type checking.
  * :class:`~.Contacts` for faster intermolecular contact 
    identification.
  * :class:`~.Select` can identify intermolecular contacts. See
    :ref:`contacts` for an examples and details.
  * :func:`~.sampleModes` implemented for sampling conformations 
    along normal modes.

**Improvements**:

  * :mod:`~.proteins.compare` functions are improved. Now they perform sequence
    alignment if simple residue number/identity based matchin does not work,
    or if user passes ``pwalign=True`` argument. This impacts the speed 
    of X-ray ensemble analysis.
  * :class:`~.Select` can cache data optionally. This results
    in speeds up from 2 to 50 folds depending on number of atoms and selection
    operations.
  * Implementation of :func:`~.showProjection` is completed.
  
:ref:`nmwiz`
-------------------------------------------------------------------------------

**Release 0.2.3**

  * For each mode a molecule for drawing arrows and a molecule for showing 
    animation is formed in VMD on demand. NMWiz remembers a color associated 
    with a mode.
  * Deselecting a residue by clicking on a plot is possible.
  * A bug causing incorrect parsing of NMD files from ANM server is fixed.


**Release 0.2.2**

  * Selection string option allows user to show a subset of arrows matching 
    a VMD selection string. Optionally, this selection string may affect
    protein and animation representations.
  * A bug that caused problems when over plotting modes is removed.
  * A bug affecting line width changes in plots is removed.
  * Selected residue representations are colored according to the color of the
    plot.

**Release 0.2.1**

  * Usability improvements.
  * Loading the same data file more than once is prevented.
  * If a GUI window for a dataset is closed, it can be reloaded from the main
    window.
  * A dataset and GUI can be deleted from the VMD session via the main window.

**Release 0.2**

  * Instant documentation is improved.
  * Problem with clearing selections is fixed.
  * Plotting options frame is populated.
  * Multiple modes can be plotted on the same canvas.

Release 0.1.2 (Nov 9, 2010)
===============================================================================

* Important bug fixes and improvements in NMA helper and plotting functions.
* Documentation updates and improvements.


Release 0.1.1 (Nov 8, 2010)
===============================================================================

* Important bug fixes and improvements in chain comparison functions.
* Bug fixes.
* Source clean up.
* Documentation improvements.

Release 0.1 (Nov 7, 2010)
===============================================================================

* First release.
