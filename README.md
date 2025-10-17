[![ProDy workflow status](https://img.shields.io/github/actions/workflow/status/prody/prody/main.yml)](https://github.com/prody/ProDy/actions/workflows/main.yml)

[![ProDy conda-forge page](https://anaconda.org/conda-forge/prody/badges/version.svg)](https://anaconda.org/conda-forge/prody)

[![ProDy pypi page](https://img.shields.io/pypi/v/ProDy.svg)](https://pypi.org/project/ProDy/)


## **SYNOPSIS**

ProDy is a free and open-source Python package for protein structure, dynamics, and sequence analysis. It allows for comparative analysis and modeling of protein structural dynamics and sequence co-evolution. The fast and flexible ProDy API is designed for interactive usage as well as application development. ProDy also comes with several analysis applications and a graphical user interface for visual analysis.

Further details are described in the ProDy papers:

Bakan A, Meireles LM, Bahar I.

ProDy: Protein Dynamics Inferred from Theory and Experiments.

Bioinformatics 2011 27(11):1575-1577.

Bakan A, Dutta A, Mao W, Liu Y, Chennubhotla C, Lezon TR, Bahar I.

Evol and ProDy for Bridging Protein Sequence Evolution and Structural Dynamics.

Bioinformatics 2014 30(18):2681-2683.

Zhang S, Krieger JM, Zhang Y, Kaya C, Kaynak B, Mikulska-Ruminska K, Doruker P, Li H, Bahar I.

ProDy 2.0: Increased Scale and Scope after 10 Years of Protein Dynamics Modelling with Python.

Bioinformatics 2021 37(20):3657-3659.

## **INSTALLING PRODY**

We recommend downloading and installing the Anaconda package manager to handle dependencies in controlled environments. ProDy can be installed with the following command:

conda install \-c conda-forge ProDy

ProDy is under active development. Install it from the source on GitHub if you want the most recent fixes.

First, download the ProDy code either as a zipped folder or using Git. For example, if you have Git installed, you can do the following:

git clone https://github.com/prody/ProDy.git

Then, change into the ProDy directory and install it by running the following commands:

cd ProDy  
python setup.py build\_ext \--inplace \--force  
pip install \-Ue .

More instructions can be found at [http://bahargroup.org/prody/downloads/](http://bahargroup.org/prody/downloads/) and [http://bahargroup.org/prody/manual/devel/develop.html](http://bahargroup.org/prody/manual/devel/develop.html).

## **MODULES**

ProDy has a modular structure with modules inside various subpackages.

The main ones are:

* **prody.atomic**: handles all Atomic objects including AtomGroup and Selection  
* **prody.database**: interfaces with databases such as CATH, DALI, UniProt and Pfam  
* **prody.dynamics**: provides all the modules for normal mode analysis with various elastic network models, as well as PCA, SignDy (prody.dynamics.signature), hybrid methods such as prody.dynamics.clustenm, allosteric signal propagation methods prody.dynamics.perturb (PRS) and calcHitTime (Markovian hitting time), and various analysis and plotting functions.  
* **prody.ensemble**: enables construction of heterogeneous structural ensembles for exploring dynamics from experiments and simulations  
* **prody.proteins**: provides various modules for parsing different kinds of protein structure files including PDB, mmCIF, MMTF and maps, as well as tools to align and compare structures, and analysis of prody.proteins.interactions within and between proteins (InSty) and find prody.proteins.waterbridges (WatFinder).  
* **prody.sequence**: has all the sequence alignment and evolutionary analysis tools of Evol

Smaller ones include:

* **prody.chromatin**: specific to chromatin dynamics (ChromDy) including prody.chromatin.hic and prody.chromatin.cluster  
* **prody.compounds**: for parsing small molecule compounds/ligands from the PDB and related databases  
* **prody.domain\_decomposition**: for Spectrus dynamical domain decomposition  
* **prody.trajectory**: for trajectories in DCD format  
* **prody.utilities**

## **GETTING PRODY**

You can run ProDy on all major platforms. For download and installation instructions see: [http://www.bahargroup.org/prody/downloads/](http://www.bahargroup.org/prody/downloads/)

## **DOCUMENTATION**

* **Homepage**: [http://www.bahargroup.org/prody/](http://www.bahargroup.org/prody/)  
* **Tutorials**: [http://www.bahargroup.org/prody/tutorials](http://www.bahargroup.org/prody/tutorials)  
* **Reference**: [http://www.bahargroup.org/prody/manual](http://www.bahargroup.org/prody/manual)  
* **Applications**: [http://www.bahargroup.org/prody/manual/apps](http://www.bahargroup.org/prody/manual/apps)  
* **NMWiz GUI**: [http://www.bahargroup.org/prody/nmwiz](http://www.bahargroup.org/prody/nmwiz)  
* **Changes**: [http://www.bahargroup.org/prody/manual/release](http://www.bahargroup.org/prody/manual/release)  
* See also [https://github.com/prody/ProDy-website](https://github.com/prody/ProDy-website) for latest versions.

## **SOURCE CODE**

* **Source code**: [https://github.com/prody/ProDy](https://github.com/prody/ProDy)  
* **Issue tracker**: [https://github.com/prody/ProDy/issues](https://github.com/prody/ProDy/issues)

## **LICENSE**

ProDy is available under the MIT License. See LICENSE.txt for more details.

* **Biopython**: ([http://biopython.org/](http://biopython.org/)) KDTree and TreeConstruction modules are distributed with ProDy. Biopython is developed by The Biopython Consortium and is available under the Biopython license ([http://www.biopython.org/DIST/LICENSE](http://www.biopython.org/DIST/LICENSE)).  
* **Pyparsing**: ([https://github.com/pyparsing/pyparsing](https://github.com/pyparsing/pyparsing)) module is distributed with ProDy. Pyparsing is developed by Paul T. McGuire and is available under the MIT license ([http://www.opensource.org/licenses/mit-license.php](http://www.opensource.org/licenses/mit-license.php)).  
* **CEalign module**: ([https://pymolwiki.org/index.php/Cealign\_plugin](https://pymolwiki.org/index.php/Cealign_plugin)) is distributed with ProDy. The original CE method was developed by Ilya Shindyalov and Philip Bourne. The Python version which is used by ProDy is developed by Jason Vertrees and available under the New BSD license.  
* **Hbp module**: The calculation of hydrophobic interactions, solvent accessible surface area (SASA) and volume for each residue is using geometric methods based on the information of the atoms in the molecule. The methods have been programmed in C++ and can be compiled as a python module “hpb.so” which is then used by ProDy. Files for compilation are stored at prody/proteins/hpbmodule folder and required C++ and Fortran compiler. After compilation hpb.so file can be stored in prody/proteins folder in ProDy or in the local directory which is used to perform calulations. The precompiled versions for Python 2.7, 3.8, 3.9, and 3.10 are availabe in prody/proteins/hpbmodule. The user can choose the correct version of hpb.so and copy to the prody/proteins or local directory. C++ code of hpb.so was developed by Xin Cao and Fortran code by Xin Cao, Michelle H. Hummel, Bihua Yu, and Evangelos A. Coutsias (License in prody/proteins/hpbmodule folder). Details of the method can be found in the Supplementary Material of InSty manuscript ([https://doi.org/10.1016/j.jmb.2025.169009](https://doi.org/10.1016/j.jmb.2025.169009)).
