.. _foreword:

*******************************************************************************
Foreword
*******************************************************************************

ProDy is designed for structure-based analysis of protein dynamics. It is 
suitable for a variety of uses from everyday parsing and analysis of protein 
structures to developing or prototyping software for structure-based analysis.  

Dynamics from experiments
===============================================================================

ProDy allows for quantitative analysis of *large* heterogeneous experimental 
structural datasets.  It takes, for example, under 6 seconds to parse, align, 
and analyze 75 p38 MAP kinase structures on a personal laptop computer 
(see :ref:`pca-xray-calculations`).  Such experimental datasets may contain 
sequence homologs, orthologs, mutant, or ligand bound forms of a protein.  
Dominant patterns in structural variability are extracted by principal 
component analysis (PCA) of the ensemble, and can be compared with 
theoretically predicted conformational dynamics using ProDy.

Dynamics from theory
===============================================================================

On the theoretical side, protein dynamics is predicted using normal mode 
analysis (NMA) based on elastic network models (ENMs) or is extracted from 
molecular dynamics (MD) trajectories using essential dynamics analysis (EDA).  
Numerous helper functions enable comparative analysis of thus obtained 
experimental and theoretical data, and visualize the principal changes 
in conformations that are accessible in different functional states.

Input for ProDy
===============================================================================

The input is for calculations is the set of atomic coordinates of the query 
protein in PDB file format, or simply the PDB id or single letter amino acid 
sequence of the protein.  Fast and flexible ProDy parsers are used to retrieve 
data from the `wwPDB <http://www.wwpdb.org/>`_ FTP servers and extract the
coordinate data and other relevant information. 

Additionally, PCA/NMA/EDA results calculated by other software can be parsed
and analyzed using rich analysis and plotting functions of ProDy.  User just 
needs to provide normal mode data in plain text files.  See the tutorial for 
an example.

Visualization
===============================================================================

Analysis and plotting capabilities of ProDy are complemented VMD plugin 
:ref:`NMWiz`.  NMWiz can be used to visualize and animate normal mode data 
calculated using ProDy or other software.  See :ref:`NMWiz` documentation for 
details. 
