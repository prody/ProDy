.. _pdbparser-performance-2:

PDB Parser Performance 2
===============================================================================

*Date: 31 Oct 2011*

ProDy v0.9 comes with improvements in :func:`~prody.proteins.parsePDB` and 
:class:`~prody.atomic.HierView`.  The improvements in speed is reported here.  


Dataset
-------------------------------------------------------------------------------
Three PDB files that contains different number of models and sized proteins 
were used in benchmarking:


.. csv-table:: Dataset
   :header: "PDB", "2k39", "1dlo", "3mrz"
   
   "# atoms", 1231, 7691, 90973
   "# models", 116, 1, 1


Results and Discussion
-------------------------------------------------------------------------------

The analysis was carried out using a desktop machine with Intel(R) Xeon(TM) CPU 
at 3.20GHz. ProDy was timed for (*i*) parsing all atoms, (*ii*) parsing Cα’s, 
and (*iii*) building a hierarchical view.

|

.. csv-table:: Parsing all atoms and models
   :header: "PDB", "2k39", "1dlo", "3mrz"
   
   "*v0.8.3*", 0.8919, 0.1629, 1.9691
   "*v0.9*", 0.4806, 0.0985, 1.2034
   
On average PDB parser is 1.7 times faster.

|

.. csv-table:: Parsing Cα atoms
   :header: "PDB", "2k39", "1dlo", "3mrz"
   
   "*v0.8.3*", 0.2834, 0.0337, 0.2157
   "*v0.9*", 0.2537, 0.0251, 0.1902

Performance of parsing atom subsets was not affected as much, it is only 1.2
times faster.

|

.. csv-table:: Building a hierarchical view
   :header: "PDB", "2k39", "1dlo", "3mrz"
   
   "*v0.8.3*", 0.0114, 0.1331, 1.0926
   "*v0.9*", 0.0041, 0.0412, 0.6422
   
Generating a hierarchical is 2.6 times faster on average.
 
Python code
-------------------------------------------------------------------------------

The following code was used for evaluation::

  from prody import *
  from time import time
  import numpy as np

  def timeit(func, *args, **kwargs):
      start = time()
      pdb = func(*args, **kwargs)
      return time() - start

  changeVerbosity('warning')

  N = 10 # repeats
  PDB = ['2k39.pdb', '1dlo.pdb', '3mrz.pdb']
  SUBSET = [None, 'ca']
  for pdb in PDB:
      for subset in SUBSET:
          tm = np.mean([timeit(parsePDB, pdb, subset=subset) for i in xrange(N)])
          print('{0:s} {1:s} {2:.4f}'.format(pdb[:4], str(subset), tm))
          if subset is None:
              ag = parsePDB(pdb)
              tm = np.mean([timeit(HierView, ag) for i in xrange(N)])
              print('{0:s} hv {1:.4f}'.format(pdb[:4], tm))
