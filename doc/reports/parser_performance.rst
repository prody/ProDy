.. _pdbparser-performance-2:

PDB Parser Performance 2
===============================================================================

*Date: 23 Nov 2010*

Performance of ProDy pdb parser :func:`prody.proteins.parsePDB` is compared to 
the :class:`Bio.PDB.PDBParser.PDBParser` using 4701 PDB structures. List of 
identifiers for a non-redundant set of PDB structures is obtained from 
http://bioinfo.tg.fh-giessen.de/pdbselect (:download:`pdb_select25`).


.. csv-table:: Results from parsing PDB select set of 4701 proteins.
   :header: "Parser", "Total time", "Per file time", "Failures"

   "*ProDy*", 636.40 s, 0.135 s, 0
   "*Bio.PDB*", 3125.89 s, 0.665 s, 0
   
ProDy PDB parser seems to be 4.9 times faster than Bio.PDB parser on average. 
   
Comparison was made using a desktop machine with Intel(R) Xeon(TM) CPU at 3.20GHz. 

The following code was used for evaluation::

  from time import time
  from glob import glob
  from prody import *
  from Bio.PDB import *

  def timeProDy():
      pdbfiles = glob('pdb_select/*pdb')
      fout = open('timer_failures_ProDy.txt', 'w')
      ProDySetVerbosity('critical')
      failures = 0
      st = time()
      for pdb in pdbfiles:
          try:
              pdb = parsePDB(pdb)
          except:
              failures += 1
              fout.write(pdb + '\n')
      print 'ProDy - Time (s): {0:.2f}'.format(time() - st)
      print 'ProDy - Failures: {0:d}'.format(failures)
      fout.close()
      
  def timeBioPDB():
      pdbfiles = glob('pdb_select/*pdb')
      fout = open('timer_failures_ProDy.txt', 'w')
      parser = PDBParser()
      failures = 0
      st = time()
      for pdb in pdbfiles:
          f = open(pdb)
          try:
              pdb = parser.get_structure('', f)
          except:
              failures += 1
              fout.write(pdb + '\n')
          f.close()
      print 'Bio.PDB - Time (s): {0:.2f}'.format(time() - st)
      print 'Bio.PDB - Failures: {0:d}'.format(failures)
      fout.close()


  if __name__ == '__main__':
      #timeProDy()
      timeBioPDB()
      
Output was::

  ProDy - Time (s): 636.40
  ProDy - Failures: 0
  Bio.PDB - Time (s): 3125.89
  Bio.PDB - Failures: 0
