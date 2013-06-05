Using NMWiz
===============================================================================

This chapter provides an introduction to NMWiz and its capabilities.

Main Window
-------------------------------------------------------------------------------

NMWiz is located in :menuselection:`VMD Main --> Extensions --> Analysis` menu
and labeled as :guilabel:`Normal Mode Wizard`, as shown below:

.. figure:: images/vmd_main_nmwiz.png
   :scale: 100 %

   Location of NMWiz plugin in VMD menus.


When you first open NMWiz, you will see the :guilabel:`NMWiz 1.2 - Main`
window shown below:

.. figure:: images/nmwiz_main.png
   :scale: 100 %

   NMWiz Main window.

:guilabel:`NMWiz 1.2 - Main` window offers options to load precalculated
data or to generate data for molecules loaded in VMD:

  * :guilabel:`Load NMD File` let's you load precalculated data (:file:`.nmd`
    files) for analysis.

  * :guilabel:`From Molecule` let's you load precalculated data stored in
    molecular data formats recognized by VMD.

  * :guilabel:`ProDy interface` allows for performing ANM, GNM, and PCA/EDA
    calculations.  This functionality is illustrated in the following chapters.

  * :guilabel:`Structure Comparison` allows for calculating, depicting, and
    animating deformation between two frames of the same molecule or two
    separate structures.


NMD Format
-------------------------------------------------------------------------------

NMD files (extension :file:`.nmd`) are plain text files that contain at least
normal mode and coordinate data. In addition to PCA, EDA, NMA, ANM, or GNM
data, arbitrary vectors can be stored in NMD files.  Following data fields
are recognized:

  * **coordinates**: Coordinates must be provided in one line as a list of
    decimal numbers. Number of atoms in the system is deduced from size of
    this data line.
  * **mode**: Normal mode array. Each normal mode array must be provided in
    one line as a list of decimal numbers. Mode array may be preceded by mode
    index and mode length (square root of variance or inverse frequency).
  * *title*: A title for the dataset.
  * *names*: Atom names. Default is "CA" for all atoms.
  * *resnames*: Residue names. Default value is "GLY".
  * *chainids*: Chain identifiers. Default value is "A".
  * *resids*: Residue numbers. If this data line if not found, residue numbers
    are started from 1 and incremented by one for each atom.
  * *betas*: Beta factors. Default value is 0 (zero). B-factors are used to
    color the protein representation.

.. note::
   All data types must be listed in a single line. The size of data lines must
   match the number of atoms in the system, i.e. the size of coordinates line.


For example, :download:`nmwiz_tutorial_files/1dlo_anm.nmd` contains
following lines::

   nmwiz_load 1dlo_anm.nmd
   name 1dlo
   atomnames CA CA CA CA CA CA CA CA CA CA ...
   resnames PRO ILE SER PRO ILE GLU THR VAL PRO VAL ...
   resids 1 2 3 4 5 6 7 8 9 10 ...
   chainids A A A A A A A A A A ...
   bfactors 71.47 63.68 46.35 42.28 36.52 31.66 20.85 18.49 24.94 29.01 ...
   coordinates 172.541 -33.818 72.141 172.591 -29.981 71.748 171.498 -29.271 ...
   mode 1 5.08 0.010 -0.014 -0.030 0.011 -0.013 -0.029 0.006 -0.012 -0.027 ...
   mode 2 3.99 -0.003 -0.014 0.007 0.001 -0.016 0.002 0.001 -0.019 0.000 ...
   mode 3 2.82 0.010 0.027 0.000 0.006 0.028 0.003 0.005 0.028 0.004 0.002 ...
   ...

The first line in this file contains ``nmwiz_load 1dlo_anm.nmd`` command.
With this, you can make :program:`vmd` load this file first thing when it
starts::

  $ vmd -e 1dlo_anm.nmd



.. note::
   Any data field name that is not recognized by NMWiz, such as
   ``nmwiz_load 1dlo_anm.nmd`` will be omitted.



NMWiz Window
-------------------------------------------------------------------------------

Now, let's load :file:`1dlo_anm.nmd` file to introduced the :guilabel:`NMWiz`
window:

.. figure:: images/nmwiz_initial.png
   :scale: 80 %

   NMWiz GUI


Changing and Editing Modes
^^^^^^^^^^^^^^^^^^^^^^^^^^

Active mode and color
"""""""""""""""""""""

Active mode can be changed using the drop down menu or :guilabel:`<=` and
:guilabel:`=>` buttons. Changing the active mode will hide previously
generated depiction.  Newly depicted mode will have a distinct color.

Mode direction
""""""""""""""

NMWiz depicts normal modes as they are read from the input file, but
:guilabel:`+/-`


Mode options
""""""""""""

.. figure:: images/nmwiz_mode_options.png
   :scale: 80 %

   Mode Graphics Options panel


:guilabel:`Mode Graphics Options` panel allows for customizing the look of
arrows. You can change the size, length, material, and resolution of arrows,
and also select to use porcupines instead of arrows.

.. note::
   Mode representation is not redrawn after manual changes in variables.
   When this is the case click one of :guilabel:`Draw` or :guilabel:`Redraw`
   buttons.


Scaling and RMSD
""""""""""""""""

Remember that normal modes are unit vectors, so when they are drawn as they are
arrows will not be visible, especially for large structures. Hence, modes are
automatically scaled to depict a given RMSD.  RMSD value and the corresponding
scaling factor is displayed in the top panel.

Atom selections
"""""""""""""""

Arrows are drawn for selected atoms. You can enter any valid atom selection
string to display arrows for a specific part of the molecule.

Molecule Representation
^^^^^^^^^^^^^^^^^^^^^^^

You can change molecule representations through NMWiz.
:guilabel:`Molecule` number in :guilabel:`Actions` panel shows the index
of the molecule that the changes will affect.  If you would like NMWiz
not to interfere with your customized visualization on this molecule,
you can turn this off in molecule :guilabel:`Options` panel.  Alternatively,
you can :guilabel:`Hide` this molecule, and make custom representations for
a molecule that you load into VMD manually.

Making Animations
^^^^^^^^^^^^^^^^^

Clicking :menuselection:`Animation --> Make` will start an animation along
the active mode.  The animated motion involves changes in the structure
along the mode in both directions as much as the RMSD displayed in top panel.

You can adjust animation settings, such as number of frames,
using the :menuselection:`Animation --> Options` panel.

Plotting and Heatmaps
^^^^^^^^^^^^^^^^^^^^^

You can generate plots and heatmaps that are linked to the protein display.
Click :guilabel:`Plot Mobility`, for example. You will see a plot


Main Controls
^^^^^^^^^^^^^

Use the following controls for:

  * :guilabel:`Main` - open NMWiz main window, useful when it's closed
  * :guilabel:`Save` - save current dataset in NMWiz format
  * :guilabel:`Remove` - remove current dataset and free memory
  * :guilabel:`Help` - get some help


Efficiency Tips
^^^^^^^^^^^^^^^

Displaying large number of primitive graphics for large molecules may make
your system unresponsive.  You can try some of the following tricks to
speed up drawing arrows and manipulating the view.

Show some of modes
""""""""""""""""""

  * Show subset of arrows evenly distributed in the molecular system
    by making a selection using ``%``:

    * ``index % 2 == 0`` - shows every other mode

    * ``resid % 4 == 0`` - shows every mode for every forth residue

    You need to enter these strings to :guilabel:`Selection` box.

  * Draw modes that are longer than a threshold length using the options
    in :guilabel:`Mode Graphics Options` in Figure 2.4.

Depending on the molecule one or both of these tricks can provide a more
responsive experience and a better (less crowded) depiction of motions.



