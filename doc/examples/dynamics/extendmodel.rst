.. _extendmodel:

*******************************************************************************
Extend a coarse-grained model
*******************************************************************************

Synopsis
=============================================================================

This example shows how to extend normal modes calculated for a 
coarse-grained model to a larger set of atoms. Extended model can be
used to generate alternate conformers. 

Input
-------------------------------------------------------------------------------

A PDB structure. 

Output
-------------------------------------------------------------------------------

Output is an extended model and conformers along selected normal modes. 
Conformers can be save in PDB format.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Normal mode analysis
-------------------------------------------------------------------------------

Conformers can be generated along any set of normal modes. In this example,
we will calculate normal modes for unbound structure of p38 MAP kinase and
generate backbone trace conformations. 

>>> p38 = parsePDB('1p38')
>>> p38_ca = p38.select('calpha')
>>> anm = ANM('1p38')
>>> anm.buildHessian(p38_ca)
>>> anm.calcModes()

Extrapolation
-------------------------------------------------------------------------------

ANM modes are extended using the :func:`~.extendModel` function: 

>>> bb_anm, bb_atoms = extendModel(anm, p38_ca, p38.select('backbone'))
>>> bb_anm
<NMA: Extended ANM 1p38 (20 modes, 1404 atoms)>
>>> bb_atoms
<AtomMap: Selection "backbone" from 1p38 from 1p38 (1404 atoms, 1404 mapped, 0 dummy)>

Note that :class:`~.GNM`, :class:`~.PCA`, and :class:`~.NMA` instances can also
be used as input to this function.

Write NMD file
-------------------------------------------------------------------------------

Extended modes can be visualized in VMD using :ref:`nmwiz` using 
an NMD file:

>>> writeNMD('p38_anm_backbone.nmd', bb_anm, bb_atoms)
'p38_anm_backbone.nmd'

Sample conformers
-------------------------------------------------------------------------------

We can use the extended model to sample backbone conformers:

>>> ensemble = sampleModes(bb_anm[:3], bb_atoms, n_confs=40, rmsd=0.8)
>>> ensemble
<Ensemble: Conformations along 3 modes from NMA Extended ANM 1p38 (40 conformations, 1404 atoms, 1404 selected)>

Note that we made used of ANM modes beyond their theoretical limitations.


Write PDB file
-------------------------------------------------------------------------------

Generated conformers can be written in PDB format as follows: 

>>> backbone = p38.copy( bb_atoms )
>>> backbone.addCoordset(ensemble)
>>> writePDB('p38_backbone_ensemble.pdb', backbone)
'p38_backbone_ensemble.pdb'

|questions|

|suggestions|
