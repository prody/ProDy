.. currentmodule:: prody.dynamics

.. _p38-xray-visualization:

*******************************************************************************
p38 X-ray Ensemble - Part IV: Visualization
*******************************************************************************

|more| Continued from :ref:`p38-xray-plotting`.

If you have left your session in the previous part, you will need to
load saved data::

    from prody import *

    pca = loadModel('p38_xray.pca.npz')
    anm = loadModel('1p38.anm.npz')
    ref_chain = parsePDB('p38_ref_chain.pdb')

It is possible to make a comparative visual analysis of PCA and ANM mode
that were calculated in the previous parts. We need |nmwiz| for this purpose.
NMWiz is a VMD plugin designed to complement ProDy. 

We will save PCA and ANM data in NMD format. 
NMWiz can read and visualize multiple NMD files at once. Interested
user is referred to NMWiz documentation for more information. NMD files
are saved as follows using :func:`writeNMD` functions::

    writeNMD('p38_pca.nmd', pca[:3], ref_chain)
    writeNMD('p38_anm.nmd', anm[:3], ref_chain)
   

It is also possible to load VMD to visualize normal mode data 
from within an interactive Python session. For this to work, you need
VMD and NMWiz plugin installed::

    # Check if VMD path is correct
    getVMDpath()
    
    viewNMDinVMD( '1p38_pca.nmd' )
