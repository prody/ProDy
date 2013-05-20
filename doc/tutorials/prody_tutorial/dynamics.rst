.. _dynamics-tutorial:


Dynamics Analysis
===============================================================================

In this section, we will show how to perform quick PCA and ANM analysis
using a solution structure of ubiquitin.  If you started a new Python session,
import ProDy contents:


.. ipython:: python
   :suppress:

   from prody import *; from pylab import *; ion()

.. ipython:: python

   structure = parsePDB('1p38')


PCA Calculations
-------------------------------------------------------------------------------

We perform principal component analysis (:class:`.PCA`) of NMR models
in PDB file 2k39 as follows:

.. ipython:: python

   ubi = parsePDB('2k39', subset='calpha')
   ubi_selection = ubi.select('resnum < 71')
   ubi_ensemble = Ensemble(ubi_selection)
   ubi_ensemble.iterpose()

Perform the PCA:

.. ipython:: python

   pca = PCA('Ubiquitin')
   pca.buildCovariance(ubi_ensemble)
   pca.calcModes()

Print the fraction of variance for top raking 4 PCs:

.. ipython:: python

   for mode in pca[:4]:
       print(calcFractVariance(mode).round(2))


PCA data can be saved on disck using :func:`.saveModel` function:

.. ipython:: python

   saveModel(pca)

This functions writes data in binary format, so is an efficient way of
storing data permanently.  In a later session, this data can be loaded using
:func:`.loadModel` function.

ANM Calculations
-------------------------------------------------------------------------------

Anisotropic network model (:class:`.ANM`) analysis can be
performed in two ways:

The shorter way, which may be suitable for interactive sessions:

.. ipython:: python

   anm, atoms = calcANM(ubi_selection, selstr='calpha')

The longer and more controlled way:

.. ipython:: python

   anm = ANM('ubi') # instantiate ANM object
   anm.buildHessian(ubi_selection) # build Hessian matrix for selected atoms
   anm.calcModes() # calculate normal modes
   saveModel(anm)


:ref:`anm` provides a more detailed discussion of ANM calculations.
The above longer way gives more control to the user. For example, instead of
building the Hessian matrix using uniform force constant and cutoff distance,
customized force constant functions (see :ref:`gamma`) or a pre-calculated matrix
(see :meth:`.ANM.setHessian`) may be used.

Individual :class:`.Mode` instances can be accessed by
indexing the :class:`.ANM` instance:

.. ipython:: python

   slowest_mode = anm[0]
   print( slowest_mode )
   print( slowest_mode.getEigval().round(3) )

Note that indices in Python start from zero (0).  0th mode is the 1st non-zero
mode in this case.

The :func:`.writeNMD` function writes ANM results in NMD format.
NMD files can be viewed using the :ref:`nmwiz` VMD plugin.

.. ipython:: python

   writeNMD('p38_anm.nmd', anm[:6], ubi_selection)

For more information on elastic network model calculations see
:ref:`enm-analysis` section.

Comparative Analysis
-------------------------------------------------------------------------------

ProDy comes with many built-in functions to facilitate a comparative analysis
of experimental and theoretical data. For example, using
:func:`.printOverlapTable` function you can see the agreement between
experimental (PCA) modes and theoretical (ANM) modes calculated above:

.. ipython:: python

   printOverlapTable(pca[:4], anm[:4])

Output above shows that PCA mode 2 and ANM mode 2 for ubiquitin show the
highest overlap (cosine-correlation).

.. ipython:: python

   @savefig prody_tutorial_dynamics_overlap.png width=4in
   showOverlapTable(pca[:4], anm[:4]);

This was a short example for a simple case. :ref:`pca` section contains more
comprehensive examples for heterogeneous datasets. :ref:`pca-xray-analysis`
shows more analysis function usage examples and :ref:`dynamics` module
documentation lists all of the analysis functions.

Output Data Files
-------------------------------------------------------------------------------

The :func:`.writeNMD` function writes PCA results in NMD format.
NMD files can be viewed using the :ref:`nmwiz` VMD plugin.

.. ipython:: python

   writeNMD('ubi_pca.nmd', pca[:3], ubi_selection)


Additionally, results can be written in plain text files for analysis with
other programs using the :func:`.writeArray` function:

.. ipython:: python

   writeArray('ubi_pca_modes.txt', pca.getArray(), format='%8.3f')


External Data
-------------------------------------------------------------------------------

Normal mode data from other NMA, EDA, or PCA programs can be parsed using
:func:`.parseModes` function for ProDy analysis.

In this case, we will parse ANM modes for p38 MAP Kinase calculated using
`ANM server <http://ignmtest.ccbb.pitt.edu/cgi-bin/anm/anm1.cgi>`_  as the
external software.  We use :download:`oanm.eigvals <prody_tutorial_files/oanm_eigvals.txt>`
and :download:`oanm.slwevs <prody_tutorial_files/oanm_slwevs.txt>` files from the ANM
server.

You can either download these files to your current working directory from here
or obtain them for another protein from the ANM server.

.. ipython:: python

   nma = parseModes(normalmodes='prody_tutorial_files/oanm_slwevs.txt',
    eigenvalues='prody_tutorial_files/oanm_eigvals.txt',
    nm_usecols=range(1,21), ev_usecols=[1], ev_usevalues=range(6,26))
.. ipython:: python

   nma
   nma.setTitle('1p38 ANM')
   slowmode = nma[0]
   print(slowmode.getEigval().round(2))


Plotting Data
-------------------------------------------------------------------------------

If you have `Matplotlib <http://matplotlib.sourceforge.net>`_, you can use
ProDy functions whose name start with ``show`` to plot data:

.. ipython:: python

   @savefig prody_tutorial_dynamics_sqflucts.png width=4in
   showSqFlucts(slowmode);


:ref:`pca-xray-plotting` shows more plotting examples and
:ref:`dynamics` module documentation lists all of the plotting functions.
