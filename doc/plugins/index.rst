.. _nmwiz:

*******************************************************************************
Normal Mode Wizard
*******************************************************************************

Normal Mode Wizard (NMWiz) is a `VMD <www.ks.uiuc.edu/Research/vmd/>`_ 
plugin designed for visual comparative analysis of normal mode data, 
i.e. modes may come from principal component, essential dynamics, normal 
mode analysis or may be any vector describing a molecular motion. 

NMWiz can be used for:

  * drawing normal modes arrows
  * making animations (conformations along a normal mode)
  * plotting square-fluctuations (labeling and highlighting residues)
  * comparing two structures and drawing deformation arrows
  
NMWiz is available with `VMD 1.9.1`_  or later.

.. _VMD 1.9.1: http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD
  
Following molecular representations are prepared using NMWiz:
  
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
|                                                  | Example figures                                |                                                |
+==================================================+================================================+================================================+
| .. image:: /_static/gallery/p38_modes_123_sm.png | .. image:: /_static/gallery/p38_anm_pca_sm.png | .. image:: /_static/gallery/p38_network_sm.png |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+
| ANM modes 1-3 for p38 MAPK                       | ANM and PCA modes for p38                      | p38 network model                              |
+--------------------------------------------------+------------------------------------------------+------------------------------------------------+

NMWiz can also be used to generate trajectories on the fly.  The movie shows 
normal mode representation and animation generated using NMWiz.  Anisotropic 
network model modes were calculated using ProDy.  Movie was generated using 
`VMD Movie Plugin <http://www.ks.uiuc.edu/Research/vmd/plugins/vmdmovie/>`_.

.. only:: html

   .. youtube:: 1OUzdzm68YY
      :width: 400

See `NMWiz documentation`_ or :ref:`nmwiz-tutorial` for usage details.  NMWiz 
recognizes :ref:`nmd-format`.

.. _NMWiz documentation: http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/

