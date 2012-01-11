.. _credits:

*******************************************************************************
Acknowledgments
*******************************************************************************

Contributors
===============================================================================

`Lidio Meireles <http://www.linkedin.com/in/lidio>`_ provided insightful 
comments on the design of ProDy classes and modules, and contributed to the 
development of early versions of :ref:`commands` for automated tasks.

`Ying Liu <http://www.linkedin.com/pub/ying-liu/15/48b/5a9>`_ provided the 
code for Perturbation Response Scanning method.   

`Kian Ho <https://github.com/kianho>`_ contributed with bug fixes and unit 
tests initiating the unit test development phase for ProDy that helped making
numerous improvements in atom selections and PDB parsing.


Funding
===============================================================================

Support from NIH grant 1R01GM086238-01 for development and maintenance of ProDy 
is acknowledged by `Ivet Bahar <http://www.ccbb.pitt.edu/Faculty/bahar/>`_.

The project described was supported by Grant Number UL1 RR024153 from the 
`National Center for Research Resources <http://www.ncrr.nih.gov/>`_ (NCRR), 
a component of the National 
Institutes of Health (NIH) and NIH Roadmap for Medical Research, and its 
contents are solely the responsibility of the authors and do not necessarily 
represent the official view of NCRR or NIH.  
Information on Re-engineering the Clinical Research Enterprise can be obtained 
from `here 
<http://nihroadmap.nih.gov/clinicalresearch/overview-translational.asp>`_.

Software
===============================================================================

Last but not least, ProDy makes use of the following great software:

* `Pyparsing <http://pyparsing.wikispaces.com/>`_, which is 
  distributed with ProDy, is used to define the sophisticated atom selection 
  grammar. This makes every user a power user by enabling fast access to and 
  easy handling of atomic data via simple selection statements.    

* `Biopython <http://biopython.org/>`_ KDTree package and pairwise2 module, 
  which are distributed ProDy, significantly enrich and improve the ProDy 
  user experience.  KDtree package allows for fast distance based selections
  making atom selections suitable for contact identification.  pairwise2 
  module enables performing sequence alignment for protein structure
  comparison and ensemble analysis.
     
* ProDy requires `Numpy <http://numpy.scipy.org/>`_ for almost all major 
  functionality including, but not limited to, storing atomic data and 
  performing normal mode calculations.  The power and speed of Numpy makes
  ProDy suitable for interactive and high-throughput structural analysis.
  
* Finally, ProDy can benefit from `Scipy <http://www.scipy.org/SciPy>`_ and
  `Matplotlib <http://matplotlib.sourceforge.net/>`_ packages.  Scipy
  makes ProDy normal calculations more flexible and on low memory machines 
  possible.  Matplotlib allows greatly enriches user experience by allowing
  plotting protein dynamics data calculated using ProDy. 
   
  
  
