.. _pdbparser-performance:

PDB Parser Performance
===============================================================================

*Date: 20 Dec 2010*

We tested the performance of the ProDy PDB parser :func:`~prody.proteins.parsePDB` 
by comparing it with the PDB parsers of Biopython 
(:class:`Bio.PDB.PDBParser.PDBParser`) and MMTK 
(:class:`MMTK.PDB.PDBConfiguration`). 

Dataset
-------------------------------------------------------------------------------
A non-redundant set of PDB structures was used. A list of PDB identifiers was 
obtained from http://bioinfo.tg.fh-giessen.de/pdbselect 
(:download:`pdb_select25.txt`). The dataset contained 4701 uncompressed files. 

Results
-------------------------------------------------------------------------------

.. csv-table:: Results from parsing PDB select set of 4701 proteins.
   :header: "", "ProDy HV", "ProDy All", "ProDy Ca", "ProDy m1", "Bio.PDB", "MMTK"

   "*Total*", 18.13 m, 11.16 m, 3.27 m, 2.23 m, 52.28 m, 155.6 m
   "*Per file*", 0.232 s, 0.142 s, 0.042 s, 0.028 s, 0.667 s, 1.986 s
   
The analysis was carried out using a desktop machine with Intel(R) Xeon(TM) CPU 
at 3.20GHz. ProDy was timed for (*i*) parsing all atoms and generating a 
hierarchical view (HV), parsing (*ii*) all atoms (All), (*iii*) Ca’s (Ca), and 
(*iv*) Ca’s from model 1 (m1). Note that by default Bio.PDB parser evaluates 
all models, and MMTK parser evaluates only the first model.
 
Discussion
-------------------------------------------------------------------------------

The ProDy PDB parser was 2.8 to 24 times faster than Bio.PDB parser on average. 
Note that Biopython and MMTK parsers perform additonal tasks when coordinates
are parsed, i.e. building a hierarchical view containing chains and residues.
ProDy parser evaluates coordinate lines and generates a plain view of atoms
to increase the speed of parsing action. A hiearchical view is generated
only when needed.  

Python code
-------------------------------------------------------------------------------

The following code was used for evaluation::

  from time import time
  from glob import glob
  from prody import *
  from Bio.PDB import *
  from MMTK.PDB import PDBConfiguration
  from MMTK.Proteins import Protein
  import numpy as np

  def getCAcoords_ProDy(pdb):
      return pdb.select('name CA').getCoordinates()

  def timeProDy(subset=None, model=None, hv=False):
      pdbfiles = glob('pdb_select/*pdb')
      fout = open('timer_failures_ProDy.txt', 'w')
      ProDySetVerbosity('critical')
      failures = 0
      st = time()
      for pdb in pdbfiles:
          try:
              structure = parsePDB(pdb, model=model, subset=subset)
              if hv:
                temp = structure.getHierView()
              caxyz = getCAcoords_ProDy(structure)
          except:
              failures += 1
              fout.write(pdb + '\n')
      print 'ProDy - Time (s): {0:.2f}'.format(time() - st)
      print 'ProDy - Failures: {0:d}'.format(failures)
      fout.close()

  def getCAcoords_BioPDB(pdb, model=0):
      """Return Cα coordinates from indicated model.
      
      Note that this function does not check whether a protein with name CA
      is from an amino acid residue.
      """
      CAxyz = []
      for chain in pdb[model]:
          for residue in chain:
              ca = residue.child_dict.get('CA', None)
              if ca is not None:
                  CAxyz.append(ca.coord)
      return np.array(CAxyz)
      
  def timeBioPDB():
      pdbfiles = glob('pdb_select/*pdb')
      fout = open('timer_failures_BioPDB.txt', 'w')
      parser = PDBParser()
      failures = 0
      st = time()
      for pdb in pdbfiles:
          f = open(pdb)
          try:
              structure = parser.get_structure('', f)
              caxyz = getCAcoords_BioPDB(structure)
          except:
              failures += 1
              fout.write(pdb + '\n')
          f.close()
      print 'Bio.PDB - Time (s): {0:.2f}'.format(time() - st)
      print 'Bio.PDB - Failures: {0:d}'.format(failures)
      fout.close()

  def getCAcoords_MMTK(filename):
      """Return Cα coordinates.
      
      Note that this function does not check whether a protein with name CA
      is from an amino acid residue.
      """
      pdb = PDBConfiguration(filename)
      CAxyz = []
      for res in pdb.residues:
          try:
              ca = res['CA']
              CAxyz.append(ca.position)
          except:
              pass
      return np.array(CAxyz)

  def getCAcoords_MMTK_2(filename):
      """Return Cα coordinates.
      This method was found to be slower, so is not reported."""
      protein = Protein(filename, model='calpha')
      return np.array([atom.position() for atom in protein.atoms])

  def timeMMTK():
      pdbfiles = glob('pdb_select/*pdb')
      fout = open('timer_failures_MMTK.txt', 'w')
      failures = 0
      st = time()
      for pdb in pdbfiles:
          try:
              caxyz = getCAcoords_MMTK(pdb)
              #caxyz = getCAcoords_MMTK_2(pdb)
          except:
              failures += 1
              fout.write(pdb + '\n')
      print 'MMTK - Time (s): {0:.2f}'.format(time() - st)
      print 'MMTK - Failures: {0:d}'.format(failures)
      fout.close()


  if __name__ == '__main__':
      #timeProDy()
      #timeBioPDB()
      timeMMTK()

