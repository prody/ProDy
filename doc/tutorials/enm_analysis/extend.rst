.. _extendmodel:


Extend a coarse-grained model
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example shows how to extend normal modes calculated for a
coarse-grained model to a larger set of atoms. Extended model can be
used to generate alternate conformers that can be save in PDB format.

Normal mode analysis


We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()

Conformers can be generated along any set of normal modes. In this example,
we will calculate normal modes for unbound structure of p38 MAP kinase and
generate backbone trace conformations.

.. ipython:: python

   p38 = parsePDB('1p38')
   p38_ca = p38.select('calpha')
   anm = ANM('1p38')
   anm.buildHessian(p38_ca)
   anm.calcModes()

Extrapolation
-------------------------------------------------------------------------------

ANM modes are extended using the :func:`.extendModel` function:

.. ipython:: python

   bb_anm, bb_atoms = extendModel(anm, p38_ca, p38.select('backbone'))
   bb_anm
   bb_atoms


Note that :class:`.GNM`, :class:`.PCA`, and :class:`.NMA` instances can also
be used as input to this function.

Write NMD file
-------------------------------------------------------------------------------

Extended modes can be visualized in VMD using :ref:`nmwiz` using
an NMD file:

.. ipython:: python

   writeNMD('p38_anm_backbone.nmd', bb_anm, bb_atoms)


Sample conformers
-------------------------------------------------------------------------------

We can use the extended model to sample backbone conformers:

.. ipython:: python

   ensemble = sampleModes(bb_anm[:3], bb_atoms, n_confs=40, rmsd=0.8)
   ensemble


Note that we made used of ANM modes beyond their theoretical limitations.


Write PDB file
-------------------------------------------------------------------------------

Generated conformers can be written in PDB format as follows:

.. ipython:: python

   backbone = bb_atoms.copy()
   backbone.addCoordset(ensemble)
   writePDB('p38_backbone_ensemble.pdb', backbone)