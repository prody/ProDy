Optimize Conformations
===============================================================================

In this part we will optimize the geometries of conformations generated
in the previous step using NAMD.


Configuration
-------------------------------------------------------------------------------

Let's find the location of NAMD executable:

.. ipython:: python

   from prody.utilities import which
   namd2 = which('namd2')
   namd2

We will need a force field file for energy minimization. VMD ships with
CHARMM force field files. We can find their location as follows:


.. ipython:: python

   with open('where_is_charmmpar.tcl', 'w') as inp:
       inp.write('''global env;
   puts $env(CHARMMPARDIR);
   exit;''')

When you run the following command, you will see more output than the
following, but the line that you need will be at the end:

.. ipython::
   :verbatim:

   In [1]: !vmd -e where_is_charmmpar.tcl
   /home/abakan/Programs/vmd-1.9.1/plugins/noarch/tcl/readcharmmpar1.2
   Info) VMD for LINUXAMD64, version 1.9.1 (February 1, 2012)
   Info) Exiting normally.

.. ipython:: python

   import os
   par = os.path.join('/home/abakan/Programs/vmd-1.9.1/'
                      'plugins/noarch/tcl/readcharmmpar1.2',
                      'par_all27_prot_lipid_na.inp')


Let's make a folder for writing optimization input and output files:

.. ipython::

   In [1]: mkdir -p p38_optimize

We will write an NAMD configuration file for each conformation based
on :download:`min.conf <conformational_sampling_files/min.conf>` file:

.. ipython:: python

   import glob
   conf = open('conformational_sampling_files/min.conf').read()
   for pdb in glob.glob(os.path.join('p38_ensemble', '*.pdb')):
       fn = os.path.splitext(os.path.split(pdb)[1])[0]
       pdb = os.path.join('..', pdb)
       out = open(os.path.join('p38_optimize', fn + '.conf'), 'w')
       out.write(conf.format(
           out=fn, pdb=pdb,
           par=par))
       out.close()


Optimization
-------------------------------------------------------------------------------

Now we will run NAMD to optimize each of these conformations. We make a list
of commands that we want to execute:

.. ipython:: python

   os.chdir('p38_optimize')  # we will run commands in this folder
   cmds = []
   for conf in glob.glob('*.conf'):
       fn = os.path.splitext(conf)[0]
       cmds.append('namd2 ' + conf + ' > ' + fn + '.log')

   cmds[:2]

We will run these commands using :mod:`multiprocessing` module.  We will
allocate 3 processors for the job:

.. ipython:: python

   from multiprocessing import Pool
   pool = Pool(3) # number of CPUs to use
   signals = pool.map(os.system, cmds)

``signals`` will collect the output from execution of NAMD. If everything goes
right, we should have only 0s.

.. ipython:: python

   set(signals)

All NAMD output should be in :file:`p38_optimize` folder.  We go back to
origional folder as follows:

.. ipython:: python


   os.chdir('..')
