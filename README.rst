.. image:: https://img.shields.io/travis/prody/ProDy.svg
   :target: http://travis-ci.org/#!/prody/ProDy

.. image:: https://img.shields.io/pypi/v/ProDy.svg
   :target: https://pypi.org/project/ProDy/

.. image:: https://img.shields.io/github/commit-activity/m/prody/ProDy.svg
   :target: https://github.com/prody/ProDy/commits/master

.. image:: https://img.shields.io/pypi/dm/ProDy.svg
   :target: http://prody.csb.pitt.edu/downloads/

SYNOPSIS
--------

ProDy is a free and open-source Python package for protein structure, dynamics,
and sequence analysis.  It allows for comparative analysis and modeling of
protein structural dynamics and sequence co-evolution.  Fast and flexible ProDy
API is for interactive usage as well as application development.  ProDy also
comes with several analysis applications and a graphical user interface for
visual analysis.


GETTING PRODY
-------------

You can run ProDy on all major platforms.  For download and installation
instructions see:

* http://prody.csb.pitt.edu/downloads


DOCUMENTATION
-------------

* Homepage: http://prody.csb.pitt.edu/

* Tutorials: http://prody.csb.pitt.edu/tutorials

* Reference: http://prody.csb.pitt.edu/manual

* Applications: http://prody.csb.pitt.edu/manual/apps

* NMWiz GUI: http://prody.csb.pitt.edu/nmwiz

* Changes: http://prody.csb.pitt.edu/manual/release


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
in the Supplementary Material of InSty manuscript 
(soon will be submitted for publication). 