.. _changes:

.. currentmodule:: prody

*******************************************************************************
Changes
*******************************************************************************

Release 0.8.2 (in development)
===============================================================================

**New Features**:

  * :func:`~proteins.fetchPDBClusters`, :func:`~proteins.loadPDBClusters`, and
    :func:`~proteins.getPDBCluster` functions are implemented for handling
    PDB sequence cluster data. These functions can be used instead of 
    :func:`~proteins.blastPDB` function for fast access to structures of 
    the same protein (at 95% sequence identity level) or similar proteins.
    
  * Perturbation response scanning method described in [CA09]_ is implemented 
    as :func:`~dynamics.scanPerturbationResponse` based on the code provided 
    by Ying Liu.
  

**Changes**:

  * :func:`~proteins.fetchLigandData` returns the URL of the XML file in the
    ligand data dictionary.
    
  * Name of the ProDy configuration file in user :file:`home` directory 
    is renamed as :file:`.prodyrc` (used to be :file:`.prody`).
    
  * :func:`~proteins.applyBiomolecularTransformations` and
    :func:`~proteins.assignSecondaryStructure` functions raise 
    :class:`ValueError` when the function fails to perform it action.
    
  * :func:`~proteins.fetchPDB` decompresses PDB files found in the working
    directory when user asks for decompressed files.

  * :func:`~proteins.parsePDB` appends *chain* and *subset* arguments to
    :func:`~atomic.AtomGroup` name.

  * *chain* argument is added to :meth:`~proteins.PDBBlastRecord.getHits`.

**Improvements**:

  * Atom selection class :class:`~select.Select` is completely redesigned
    to prevent breaking of the parser when evaluating invalid selection
    strings.
    
  * Improved type checking in :func:`~proteins.parsePDB` function.
    
**Bug Fixes**:

  * Bugfixes in :func:`~proteins.parseDSSP`: one emerged problems in lines
    indicating chain breaks, another did not parse bridge-partners correctly.
    Both fixes are contributed by Kian Ho.
    
  * Bugfix in :func:`~proteins.parsePDB` function. When only header is desired
    (``header=True, model=0``), would return a tuple containing an empty 
    atom group and the header.

**Developmental**:

  * Unit tests for :mod:`proteins` and :mod:`select` modules are developed. 

Release 0.8.1 (Sep 16, 2011)
===============================================================================

**New Features**:

  * :func:`~proteins.fetchLigandData` is implemented for fetching ligand data
    from Ligand Expo.

  * :func:`~proteins.parsePSF` function is implemented for parsing X-PLOR 
    format PSF files.
    
**Changes**:

  * __slots__ is used in :class:`~atomic.AtomGroup` and :class:`~atomic.Atomic`
    classes. This change prevents user from assigning new variables to 
    instances of all classes derived from the base :class:`~atomic.Atomic`.
    
  * :mod:`pyparsing` is updated to version 1.5.6.

**Bug Fixes**:

  * A bug in :meth:`~atomic.AtomGroup.copy` method is fixed. When AtomGroup
    instance itself is copied, deep copies of data arrays were not made.

  * A bug in :class:`~select.Select` class raising exceptions when negative
    residue number values are present is fixed.

  * Another bug in :class:`~select.Select` class misinterpreting
    ``same residue as ...`` statement when specific chains are involved is
    fixed.

  * A bug in :meth:`~atomic.AtomGroup.addCoordset` method duplicating 
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
   existing scripts. If you are using :class:`~ensemble.Ensemble` class 
   for analyzing PDB structures, rename it as :class:`~ensemble.PDBEnsemble`.
   See the other changes that may affect your work below and the class 
   documentation for more information.
    

**New Features**:

  * :class:`~ensemble.DCDFile` is implemented for handling DCD files.
    
  * :class:`~ensemble.Trajectory` is implemented for handling multiple 
    trajectory files.
  
  * :func:`~ensemble.writeDCD` is implemented for writing DCD files.
    
  * :ref:`trajectory` example to illustrate usage of new classes for handling
    DCD files. :ref:`eda` example is updated to use new ProDy classes.
  
  * :class:`~dynamics.PCA` supports :class:`~ensemble.Trajectory` and
    :class:`~ensemble.DCDFile` instances. 
  
  * :class:`~ensemble.Ensemble` and :class:`~ensemble.PDBEnsemble` classes
    can be associated with :class:`~atomic.AtomGroup` instances. This allows
    selecting and evaluating coordinates of subset of atoms. See
    :meth:`~ensemble.EnsembleBase.setAtomGroup`, 
    :meth:`~ensemble.EnsembleBase.select`,
    :meth:`~ensemble.EnsembleBase.getAtomGroup`, and 
    :meth:`~ensemble.EnsembleBase.getSelection` methods.
  
  * :func:`~proteins.execDSSP`, :func:`~proteins.parseDSSP`, and 
    :func:`~proteins.performDSSP` functions are implemented for 
    executing and parsing DSSP calculations.

  * :func:`~proteins.execSTRIDE`, :func:`~proteins.parseSTRIDE`, and 
    :func:`~proteins.performSTRIDE` functions are implemented for 
    executing and parsing DSSP calculations.

  * :func:`~proteins.parsePDB` function parses atom serial numbers. Atoms
    can be retrieved from an :class:`~atomic.AtomGroup` instance by their
    serial numbers using :meth:`~atomic.AtomGroup.getBySerial` and
    :meth:`~atomic.AtomGroup.getBySerialRange` methods. 
    
  * :func:`~measure.calcADPs` function can be used to calculate anisotropic
    displacement parameters for atoms with anisotropic temperature factor
    data.
    
  * :meth:`~ensemble.Ensemble.getRMSFs` is implemented for calculating
    root mean square fluctuations.
    
  * :class:`~atomic.AtomGroup` and :class:`~dynamics.Mode` or 
    :class:`~dynamics.Vector` additions are supported. This adds a new
    coordinate set to the :class:`~atomic.AtomGroup` instance.
    
  * :meth:`~atomic.AtomGroup.getAttrNames` is implemented for listing
    user set attribute names.
    

**Improvements**:

  * :func:`~dynamics.calcProjection`, :func:`~dynamics.showProjection`, 
    and :func:`~dynamics.showCrossProjection` functions can optionally 
    calculate/display RMSD along the normal mode. 
  
  * ANM, GNM, and PCA routines can optionally write compressed ProDy data files.
  
  * :func:`~proteins.fetchPDB` function can optionally write decompressed 
    files and force copying a file from local mirror to target folder.

  * :meth:`~dynamics.PCA.buildCovariance` and :meth:`~dynamics.PCA.performSVD` 
    methods accept Numpy arrays as coordinate sets.

  * Performance of :meth:`~dynamics.PCA.buildCovariance` method is optimized 
    for evaluation of PDB ensembles.
  
  * :func:`~measure.calcRMSD` and :func:`~measure.superpose` functions are
    optimed for speed and memory usage.
  
  * :meth:`~ensemble.Ensemble.getMSFs` is optimized for speed and memory usage.

  * Improvements in memory operations in :mod:`~prody.atomic`, 
    :mod:`~prody.ensemble`, and :mod:`~prody.dynamics` modules for 
    faster data (PDB/NMD) output.
    
  * Optimizations in :class:`~select.Select` and :class:`~select.Contacts`
    classes.

**Changes**:

  * :class:`~ensemble.Ensemble` does not store conformation names. Instead, 
    newly implemented :class:`~ensemble.PDBEnsemble` class stores identifiers 
    for individual conformations (PDB IDs). This class should be used in cases 
    where source of individual conformations is important.
    
  * :func:`~dynamics.calcProjection`, :func:`~dynamics.showProjection`, and
    :func:`~dynamics.showCrossProjection` function calculate/display 
    root mean square deviations, by default. 
    
  * Oxidized cysteine residue abbreviation ``CSO`` is added to the definition
    of ``protein`` keyword.
    
  * :meth:`getMSF` method is renamed as :meth:`~ensemble.Ensemble.getMSFs`.
  
  * :func:`~ensemble.parseDCD` function returns :class:`~ensemble.Ensemble`
    instances.

**Bug Fixes**:
  
  * A bug in :mod:`select` module causing exceptions when regular expressions
    are used is fixed.
  
  * Another bug in :mod:`select` module raising exception when "(not ..," is
    passed is fixed.

  * Various bugfixes in :mod:`ensemble` module.
  
  * Problem in :file:`fetchpdb.py` that occurred when a file is found in a
    local mirror is fixed.
    
  * Bugfix in :meth:`~atomic.AtomPointer.copy` method.

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

  * :func:`~ensemble.parseDCD` is implemented for parsing coordinate sets 
    from DCD files.


**Improvements**:

  * :func:`~proteins.parsePDB` parses ``SEQRES`` records in header sections.

**Changes**:

  * Major classes can be instantiated without passing a name argument.
  * Default selection in NMWiz ProDy interface is changed to ensure selection
    only protein Cα atoms.


**Bug Fixes**:

  * A bug in :func:`~dynamics.writeNMD` function causing problems when writing
    a single mode is fixeed.
  * Other bugfixes in :mod:`~dynamics` module functions.


Release 0.7.1 (Apr 28, 2011)
===============================================================================

**Highlights**:

  * :class:`~atomic.Atomic` :meth:`__getattribute__` is overloaded to interpret
    atomic selections following the dot operator. For example,
    ``atoms.calpha`` is interpreted as ``atoms.select('calpha')``. See
    :ref:`atomic` for more details.
  * :class:`~atomic.AtomGroup` class is integrated with 
    :class:`~atomic.HierView` class. Atom group instances now can be indexed
    to get chains or residues and number of chains/residues can be retrieved.
    A hierarchical view is generated and updated when needed. See
    :ref:`atomic` for more details.
      
     

**New Features**:

  * :func:`~compare.matchAlign` is implemented for quick alignment of protein
    structures. See :ref:`extract-ligands` usage example.
  * :meth:`~atomic.AtomGroup.setAttribute`, 
    :meth:`~atomic.AtomGroup.getAttribute`,
    :meth:`~atomic.AtomGroup.delAttribute`, and
    :meth:`~atomic.AtomGroup.isAttribute` functions are implemented for
    :class:`~atomic.AtomGroup` class to facilitate storing user provided 
    atomic data. See :ref:`attributes` example.
  * :func:`~atomic.saveAtoms` and :func:`~atomic.loadAtoms` functions 
    are implemented to allow for saving atomic data and loading it 
    This saves custom atomic attributes and much faster than parsing
    data from PDB files. 
  * :func:`~dynamics.calcCollectivity` function is implemented to allow
    for calculating collectivity of deformation vectors.
    
**Improvements**:

  * :func:`~proteins.parsePDB` can optionally return biomolecule when 
    ``biomol=True`` keyword argument is passed.
  * :func:`~proteins.parsePDB` can optionally make secondary structure
    assignments when ``secondary=True`` keyword argument is passed.
  * :func:`~dynamics.calcSqFlucts` function is changed to accept 
    :class:`~dynamics.Vector` instances, e.g. deformation vectors.

**Changes**:

  * Changes were made in :func:`~measure.calcADPAxes` function to follow 
    the conventions in analysis ADPs. See its documentation.

**Bug Fixes**:

  * A in :class:`~ensemble.Ensemble` slicing operations is fixed. Weights are
    now copied to the new instances obtained by slicing.
  * Bug fixes in :mod:`dynamics` plotting functions 
    :func:`~dynamics.showScaledSqFlucts`, :func:`~dynamics.showNormedSqFlucts`, 

Release 0.7 (Apr 4, 2011)
===============================================================================

**New Features**:

  * Regular expressions can be used in atom selections. See 
    :mod:`~prody.select` module for details.
  
  * User can define selection macros using :func:`~select.defSelectionMacro`
    function. Macros are saved in ProDy configuration and loaded in later
    sessions. See :mod:`~prody.select` module for other related functions.
  
  * :func:`~dynamics.parseSparseMatrix` function is implemented for parsing
    matrices in sparse format. See the usage example in :ref:`external-matrix`. 

  * :func:`~dynamics.deform` function is implemented for deforming coordinate 
    sets along a normal mode or linear combination of multiple modes. 

  * :func:`~dynamics.sliceModel` function is implemented for slicing normal
    mode data to be used with functions calculating atomic properties using
    normal modes. 

**Improvements**:

  * Atom selections using bare keyword arguments is optimized. New keyword
    definitions are added. See :mod:`~prody.select` module for the complete 
    list.
  
  * A new keyword argument for :func:`~measure.calcADPAxes` allows for
    comparing largest axis to the second largest one.

**Changes**:

  * There are changes in function used to alter definitions of selection
    keywords. See :mod:`~prody.select` for details.
    
  * :func:`~proteins.assignSecondaryStructure` function assigns SS identifiers
    to all atoms in a residue. Residues with no SS information specified is
    assigned coil conformation. 
  
  * When :func:`~ensemble.Ensemble` and :func:`~dynamics.NMA` classes are 
    instantiated with an empty string, instances are called "Unnamed".

  * :func:`~dynamics.sliceMode`, :func:`~dynamics.sliceVector` and
    :func:`~dynamics.reduceModel` functions return the atom selection 
    in addition to the sliced vector/mode/model instance.

**Bug Fixes**:

  * Default selection for :func:`~dynamics.calcGNM` function is set to
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

  * :meth:`~dynamics.PCA.performSVD` function is implemented for faster
    and more memory efficient principal compoment analysis.
  * :func:`~dynamics.extrapolateModel` function is implemented for 
    extrapolating a coarse-grained model to an all atom model. See the 
    usage example :ref:`extrapolate`.
  * :func:`plog` is implemented for enabling users to make log entries.

**Improvements**:

  * :mod:`~prody.compare` functions are improved to handle insertion codes.
  * :class:`~atomic.HierView` allows for indexing using chain identifier
    and residue numbers. See usage example :ref:`hierview`
  * :class:`~atomic.Chain` allows for indexing using residue number and
    insertion code. See usage example :ref:`hierview`.
  * :meth:`~atomic.AtomGroup.addCoordset` function accepts 
    :class:`~atomic.Atomic` and :class:`~ensemble.Ensemble` instances
    as *coords* argument.
  * New method :meth:`~atomic.HierView.getAtoms` is implemented.
  * :class:`~atomic.AtomGroup` set functions check the correctness of 
    dimension of data arrays to prevent runtime problems.
  * :file:`pca.py` script is updated to use the faster PCA method
    that uses SVD.

**Changes**:

  * "backbone" definition now includes the backbone hydrogen atom 
    (Thanks to Nahren Mascarenhas for pointing to this discrepancy in the
    keyword definition). 

**Bug Fixes**:

  * A bug in :class:`~dynamics.PCA` allowed calculating covariance matrix
    for less than 3 coordinate sets is fixed.
  * A bug in :func:`~compare.mapOntoChain` function that caused problems
    when mapping all atoms is fixed.
    
    

Release 0.6.1 (Mar 2, 2011)
===============================================================================

**New Features**:

  * :func:`~proteins.setWWPDBFTPServer` and :func:`~proteins.getWWPDBFTPServer` 
    functions allow user to change or learn the WWPDB FTP server that ProDy 
    uses to download PDB files. Default server is RCSB PDB in USA. 
    User can change the default server to one in Europe or Japan.
  * :func:`~proteins.setPDBMirrorPath` and :func:`~proteins.getPDBMirrorPath` 
    functions allow user to specify or learn the path to a local PDB mirror.
    When specified, a local PDB mirror is preferred for accessing PDB files,
    over downloading them from FTP servers.
  * :func:`~compare.mapOntoChain` function is improved to map backbone or 
    all atoms.

**Improvements**:

  * :class:`~proteins.WWPDB_PDBFetcher` can download PDB files from different
    WWPDB FTP servers.
  * :class:`~proteins.WWPDB_PDBFetcher` can also use local PDB mirrors for
    accessing PDB files.

**Changes**:

  * :class:`RCSB_PDBFetcher` is renamed as :class:`~proteins.WWPDB_PDBFetcher`.
  * :func:`~compare.mapOntoChain` and :func:`~compare.matchChains` functions
    accept ``"ca"`` and ``"bb"`` as *subset* arguments.
  * Definition of selection keyword "protein" is updated to include
    some non-standard amino acid abbreviations. 

**Bug Fixes**:

  * A bug in :class:`~proteins.WWPDB_PDBFetcher` causing exceptions when
    non-string items passed in a list is fixed.
  * An important bug in :func:`~proteins.parsePDB` is fixed. When parsing
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
    installed, ProDy installer locates VMD plug-ins folder and installs NMWiz.
    On Windows, user needs to follow a separate set of instructions (see
    :ref:`getnmwiz`).
  * :class:`~dynamics.Gamma` class is implemented for facilitating use of  
    force constants based on atom type, residue type, or property. An
    example derived classes are :class:`~dynamics.GammaStructureBased` and 
    :class:`~dynamics.GammaVariableCutoff`.
  * :func:`~dynamics.calcTempFactors` function is implemented to 
    calculate theoretical temperature factors.
  * 5 new :ref:`scripts` are implemented, and existing scripts are improved to
    output figures.
  * :meth:`~dynamics.NMABase.getModel` method is implemented to 
    make function development easier.
  * :func:`~dynamics.resetTicks` function is implemented to change X and/or Y
    axis ticks in plots when there are discontinuities in the plotted data. 

**Improvements**:

  * :meth:`~dynamics.ANM.buildHessian` and :meth:`~dynamics.GNM.buildKirchhoff`
    classes are improved to accept :class:`~dynamics.Gamma` instances
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
  * :func:`~proteins.applyBiomolecularTransformations` function is improved
    to handle large biomolecular assemblies.
  * Performance optimizations in :func:`~proteins.parsePDB` and other 
    functions.
  * :class:`~ensemble.Ensemble` class accepts :class:`atomic.Atomic` 
    instances and automatically adds coordinate sets to the ensemble. 
  
**Changes**:
 
  * :class:`PDBlastRecord` is renamed as :class:`~proteins.PDBBlastRecord`. 
  * :class:`~dynamics.NMABase` instances can be index using a list or tuple of 
    integers, e.g. ``anm[1,3,5]``.
  * "ca", "bb", and "sc" keywords are defined as short-hands 
    for "calpha", "backbone", and "sidechain", respectively.
  * Behavior of :func:`~dynamics.calcANM` and :func:`~dynamics.calcGNM` 
    functions have changed. They return the atoms used for calculation as well.

**Bug Fixes**:
    
  * A bug in :func:`~proteins.assignSecondaryStructure` function is fixed.
  * Bug fixes in :ref:`scripts-anm` and :ref:`scripts-gnm`.
  * Bug fixes in :func:`~dynamics.showSqFlucts` and 
    :func:`~dynamics.showProjection` functions.
    
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
    factors: :func:`~measure.calcADPAxes` and :func:`~measure.buildADPMatrix`.
  * :meth:`~dynamics.NMA.setEigens` and :meth:`~dynamics.NMA.addEigenpair` 
    methods are implemented to assist analysis of normal modes calculated using 
    external software.
  * :func:`~dynamics.parseNMD` is implemented for parsing NMD files.
  * :func:`~dynamics.parseModes` is implemented for parsing normal mode data.
  * :func:`~dynamics.parseArray` is implementing for reading numeric data,
    particularly normal mode data calculated using other software for analysis
    using ProDy. 
  * The method in [BH02]_ to calculate overlap between covariance matrices is 
    implemented as :func:`~dynamics.calcCovarianceOverlap` function.
  * A method to trim :class:`~ensemble.Ensemble` instances is implemented: 
    :func:`~ensemble.trimEnsemble`
  * :func:`checkUpdates` to check for ProDy updates is implemented.

**Changes**:

  * Change in default behavior of :func:`~proteins.parsePDB` function. When
    alternate locations exist, those indicated by A are parsed. For parsing
    all alternate locations user needs to pass ``altloc=True`` argument.    
  * :func:`getSumOfWeights` is renamed as :func:`~ensemble.calcSumOfWeights`.
  * :func:`mapAtomsToChain` is renamed as :func:`~compare.mapOntoChain`.
  * :func:`ProDyStartLogFile` is renamed as :func:`startLogfile`.
  * :func:`ProDyCloseLogFile` is renamed as :func:`closeLogfile`.
  * :func:`ProDySetVerbosity` is renamed as :func:`changeVerbosity`.

**Improvements**:

  * A few bugs in ensemble and dynamics classes are fixed.
  * Improvements in :class:`~proteins.RCSB_PDBFetcher` allow it not to miss a 
    PDB file if it exists in the target folder.
  * :func:`~dynamics.writeNMD` is fixed to output B-factors (Thanks to 
    Dan Holloway for pointing it out).

Release 0.5.2 (Jan 12, 2011)
===============================================================================

**Bug Fixes**:
  
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

**Bug Fixes**:

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
