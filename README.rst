.. image:: https://img.shields.io/github/actions/workflow/status/prody/prody/main.yml
   :target: https://github.com/prody/ProDy/actions/workflows/main.yml

.. image:: https://img.shields.io/pypi/v/ProDy.svg
   :target: https://pypi.org/project/ProDy/

.. image:: https://img.shields.io/github/commit-activity/m/prody/ProDy.svg
   :target: https://github.com/prody/ProDy/commits/master

.. image:: https://img.shields.io/pypi/dm/ProDy.svg
   :target: http://www.bahargroup.org/prody/downloads/

SYNOPSIS
--------

ProDy is a free and open-source Python package for protein structure, dynamics,
and sequence analysis.  It allows for comparative analysis and modeling of
protein structural dynamics and sequence co-evolution.  Fast and flexible ProDy
API is for interactive usage as well as application development.  ProDy also
comes with several analysis applications and a graphical user interface for
visual analysis.

Further details are described in the ProDy papers:

  | Bakan A, Meireles LM, Bahar I.
  | *ProDy*: Protein Dynamics Inferred from Theory and Experiments.
  | *Bioinformatics* **2011** 27(11):1575-1577.

  | Bakan A, Dutta A, Mao W, Liu Y, Chennubhotla C, Lezon TR, Bahar I.
  | *Evol* and *ProDy* for Bridging Protein Sequence Evolution and Structural Dynamics.
  | *Bioinformatics* **2014** 30(18):2681-2683.

  | Zhang S, Krieger JM, Zhang Y, Kaya C, Kaynak B, Mikulska-Ruminska K, Doruker P, Li H, Bahar I.
  | *ProDy* 2.0: Increased Scale and Scope after 10 Years of Protein Dynamics Modelling with Python.
  | *Bioinformatics* **2021** 37(20):3657-3659.

INSTALLING PRODY
________
ProDy is under active development, so we recommend installing it from source from GitHub to ensure everything works properly. 

We recommend downloading and installing the Anaconda package manager to handle dependencies in controlled environments. 

Then you should download the ProDy code either as a zipped folder to extract or using git as directed at the big green button says code. For example, if you have git installed then you can do the following.

git clone https://github.com/prody/ProDy.git

Then change directory into the ProDy directory tree and install it by running the following commands.

cd ProDy 
python setup.py build_ext --inplace --force
pip install -Ue .

More instructions can be found at http://bahargroup.org/prody/downloads/ and http://bahargroup.org/prody/manual/devel/develop.html

MODULES
--------

ProDy has a modular structure with modules inside various subpackages.

The main ones are:

- :mod:`~prody.atomic` handles all :class:`~prody.atomic.Atomic` objects including :class:`~prody.atomic.AtomGroup` and :class:`.Selection`

- :mod:`~prody.database` interfaces with databases such as CATH, DALI, UniProt and Pfam

- :mod:`~prody.dynamics` provides all the modules for normal mode analysis with various elastic network models, 
as well as PCA, SignDy (:mod:`~prody.dynamics.signature`), hybrid methods such as :mod:`~prody.dynamics.clustenm`, 
allosteric signal propagation methods :mod:`~prody.dynamics.perturb` (PRS) and :func:`~prody.dynamics.analysis.calcHitTime` (Markovian hitting time),
and various analysis and plotting functions.

- :mod:`~prody.ensemble` enables construction of heterogeneous structural ensembles for exploring dynamics from experiments and simulations

- :mod:`~prody.proteins` provides various modules for parsing different kinds of protein structure files including PDB, mmCIF, MMTF and maps,
as well as tools to align and compare structures, and analysis of :mod:`~prody.proteins.interactions` within and between proteins (InSty) and 
find :mod:`~prody.proteins.waterbridges` (WatFinder).

- :mod:`~prody.sequence` has all the sequence alignment and evolutionary analysis tools of Evol


Smaller ones include:

- :mod:`~prody.chromatin` specific to chromatin dynamics (ChromDy) including :mod:`~prody.chromatin.hic` and :mod:`~prody.chromatin.cluster`

- :mod:`~prody.compounds` for parsing small molecule compounds/ligands from the PDB and related databases

- :mod:`~prody.domain_decomposition` for Spectrus dynamical domain decomposition 

- :mod:`~prody.trajectory` for trajectories in DCD format

- :mod:`~prody.utilities`


GETTING PRODY
-------------

You can run ProDy on all major platforms.  For download and installation
instructions see:

* http://www.bahargroup.org/prody/downloads/


DOCUMENTATION
-------------

* Homepage: http://www.bahargroup.org/prody/

* Tutorials: http://www.bahargroup.org/prody/tutorials

* Reference: http://www.bahargroup.org/prody/manual

* Applications: http://www.bahargroup.org/prody/manual/apps

* NMWiz GUI: http://www.bahargroup.org/prody/nmwiz

* Changes: http://www.bahargroup.org/prody/manual/release

See also https://github.com/prody/ProDy-website for latest versions.


SOURCE CODE
-----------

* Source code: https://github.com/prody/ProDy

* Issue tracker: https://github.com/prody/ProDy/issues


LICENSE
-------

ProDy is available under MIT License. See LICENSE.txt for more details.

Biopython (http://biopython.org/) KDTree and TreeConstruction modules are distributed
with ProDy. Biopython is developed by The Biopython Consortium and is available
under the Biopython license (http://www.biopython.org/DIST/LICENSE).

Pyparsing (https://github.com/pyparsing/pyparsing) module is distributed with ProDy.
Pyparsing is developed by Paul T. McGuire and is available under the MIT
license (http://www.opensource.org/licenses/mit-license.php).

CEalign module (https://pymolwiki.org/index.php/Cealign_plugin) is distributed 
with ProDy. The original CE method was developed by Ilya Shindyalov and Philip 
Bourne. The Python version which is used by ProDy is developed by Jason Vertrees 
and available under the New BSD license. 

Hbp module: The calculation of hydrophobic interactions, solvent accessible surface 
area (SASA) and volume for each residue is using geometric methods based on the 
information of the atoms in the molecule. The methods have been programmed in C++ 
and can be compiled as a python module “hpb.so” which is then used by ProDy.
Files for compilation are stored at prody/proteins/hpbmodule folder and
required C++ and Fortran compiler. After compilation hpb.so file can be
stored in prody/proteins folder in ProDy or in the local directory which
is used to perform calulations. The precompiled versions for Python 2.7,
3.8, 3.9, and 3.10 are availabe in prody/proteins/hpbmodule. The user can
choose the correct version of hpb.so and copy to the prody/proteins or
local directory.

C++ code of hpb.so was developed by Xin Cao and Fortran code by Xin Cao, 
Michelle H. Hummel, Bihua Yu, and Evangelos A. Coutsias (License in 
prody/proteins/hpbmodule folder). Details of the method can be found 
in the Supplementary Material of InSty manuscript (https://doi.org/10.1016/j.jmb.2025.169009). 
