.. currentmodule:: prody.select

.. _selections:

*******************************************************************************
Atom Selections
*******************************************************************************

ProDy offers a powerful atom selector which works very much the same way |vmd| 
does. There are only small differences which should not affect most practical 
uses of atom selections. This section describes the keywords and selection 
syntax.

.. _selkeys:

Selection Keywords
===============================================================================

Single words
-------------------------------------------------------------------------------

* :term:`all`, :term:`none`
* :term:`protein`, :term:`nucleic`, :term:`water`, :term:`hetero`
* :term:`calpha`, :term:`backbone`, :term:`sidechain`
* :term:`hydrogen`, :term:`noh`
* :term:`acidic`, :term:`basic`, :term:`charged`, :term:`polar`
* :term:`neutral`, :term:`aliphatic`, :term:`hydrophobic`
* :term:`aromatic`, :term:`cyclic`, :term:`acyclic`
* :term:`small`, :term:`medium`, :term:`large`

**Examples**:
  * "noh protein" selects non-hydrogen protein atoms
  * "charged and cyclic" selects histidine residues 

.. note::
   Definitions of these keywords can be obtained and changed using 
   corresponding get and set functions. These functions are mentioned in the
   definition of the keyword in this page and are also listed in :ref:`keywords`.

Keywords followed by characters or words
-------------------------------------------------------------------------------

* :term:`name`, :term:`element`, :term:`type`
* :term:`resname`
* :term:`chain`, :term:`segment`

**Examples**:
  * "name CA" selects atoms with name CA
  * "protein name CA and chain A" selects alpha carbons of chain A 


Keywords followed by integers and/or number ranges
-------------------------------------------------------------------------------

* :term:`index`, :term:`serial`
* :term:`resnum`, :term:`resid`

**Examples**:
  * "index 10" selects 11th atom
  * "serial 10" selects 10th atom
  * "resnum 1 to 10" selects residues with residue numbers from 1 to 10 

Keywords followed by floating point numbers and/or number ranges
-------------------------------------------------------------------------------

* :term:`x`, :term:`y`, :term:`z`
* :term:`beta`, :term:`occupancy`
* :term:`mass`, :term:`radius`, :term:`charge`

**Examples**:

* "x 0 to 20" selects atoms with x coordinates greater or equal to 0 and lesser or equal to 20
* "occupancy 1" selects atoms with occupancy values equal to 1 

Operations, Functions and Comparisons
-------------------------------------------------------------------------------

Operations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. csv-table::
   :header: "Operation", "Description"

   x ** y or x ^ y, "x to the power y"
   x * y, "x times y"
   x / y, "x divided by y"
   x // y, "x divided by y (floor devision)"
   x % y, "x modulo y"
   x + y, "x plus y" 
   x - y, "x minus y"

Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: "Function", "Description"

   abs(x), "absolute value of x" 
   acos(x), "arccos of x"
   asin(x), "arcsin of x"
   atan(x), "arctan of x"
   ceil(x), "smallest integer not less than x"
   cos(x), "cosine of x"
   cosh(x), "hyperbolic cosine of x"
   floor(x), "largest integer not greater than x" 
   exp(x), "e to the power x"
   log(x), "natural logarithm of x"
   log10(x), "base 10 logarithm of x"
   sin(x), "sine of x"
   sinh(x), "hyperbolic sine of x"
   sq(x), "square of x"
   sqrt(x), "square-root of x"
   tan(x), "tangent of x"
   tanh(x), "hyperbolic tangent of x"

Comparisons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. csv-table:: 
   :header: "Comparison", "Description"

   <, "less than"
   >, "greater than"
   <=, "less than or equal"
   >=, "greater than or equal"
   == or =, "equal"
   !=, "not equal"

**Examples**
  * "sqrt(x**2 + y**2 + z**2) < 10" selects atoms within 10 angstrom of the origin
  * "resnum <= 100" selects atoms with residue numbers less than or equal to 100  

.. _keydefs:

Keyword Definitions
===============================================================================

.. glossary:: 
    
    acidic
        Amino acid residues with acidic sidechains. Residue names include ASP and GLU.
        This definition can be changed using :func:`setAcidicResidueNames` 
        method.
    
    acyclic
        Non-:term:`cyclic` :term:`protein` residues.
    
    aliphatic
        Amino acid residues with aliphatic sidechains. Residue names include ALA, GLY, ILE, LEU, and VAL.
        This definition can be changed using :func:`setAliphaticResidueNames` 
        method.

    all
        All of the atoms in the molecule.
        
    altloc
        Alternative location identifier.

    aromatic
        Amino acid residues with aromatic sidechains. Residue names include HIS, PHE, TRP, TYR.
        This definition can be changed using :func:`setAromaticResidueNames` 
        method.

    backbone
        Group of :term:`protein` atoms whose names match one of CA, N, C, or O.
        Note that this definition contains only non-hydrogen atoms, but may
        be changed using :func:`setBackboneAtomNames` 
        method.

    basic
        Amino acid residues with basic sidechains. Residue names include ASP and GLU.
        This definition can be changed using :func:`setBasicResidueNames` 
        method.
    
    beta
        Atomic temperature (B/beta) factors.
    
    calpha
        Alpha carbon atoms of :term:`protein` residues.
    
    chain
        Poly-peptide/nucleotide/etc. chain identifier. "_" means atoms no chain
        identifier or a whitespace.
        
        e.g. "chain A B _" selects atoms whose chain identifiers are A, B, or a whitespace 
        
    charge
        Atomic partial charges.
    
    charged
        :term:`Acidic` and :term:`basic` residues.
    
    cyclic
        Amino acid residues with cyclic sidechians. 
        Residue names include HIS, PHE, PRO, TRP, TYR.
        This definition can be changed using :func:`setCyclicResidueNames` 
        method.
    
    element
        Chemical element symbols.
    
    hetero
        Not :term:`protein` or :term:`nucleic`.
    
    hydrogen
        Atoms with name matching the regular expression "[0-9]?H.*".
        This regular expression may be changed using 
        :func:`setHydrogenRegex`. See :mod:`re` module
        for more details on regular expressions.
        
    hydrophobic
        Not :term:`charged` or :term:`polar`.
    
    index 
        Atom numbers starting at 0.
        
    large
        :term:`Protein` residues that are not :term:`small` or :term:`medium`.

    mass
        Atomic mass.
    
    medium
        Amino acid residues with medium size sidechains. 
        Residue names include VAL, THR, ASP, ASN, PRO, CYS.
        This definition can be changed using :func:`setMediumResidueNames` 
        method.
        
    name
        Atom name.
    
    neutral
        Non-:term:`charged` :term:`protein` residues.
    
    noh
        Non-:term:`hydrogen` atoms.
    
    none
        None of the atoms in the molecule.
    
    nucleic
        Group of atoms whose residue names match one of GUA, ADE, CYT, THY, URA,
        DA, DC, DG, or DT.
        List of residue names in this definition can be chagned using
        :func:`setNucleicResidueNames` method.

    occupancy
        Atomic occupancy values.
    
    polar
        Amino acid residues with polar sidechains.
        
    protein
        Group of atoms whose residue names match 3-letter standard and 
        some non-standard amino acid abbreviations. List of residue names
        in the default definition is:
        ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, 
        LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, 
        or HSD, HSE, HSP.
        
        Note that this list of residue names can be changed using
        :func:`setProteinResidueNames` method.
        
    radius
        Atomic radius.
    
    resid
        Same as :term:`resnum`.

    residue
        A set pf atoms with the same residue number and chain identifier.

    resname
        Residue name abbreviation.
        
        e.g. resname ALA ARG ASN
    
    resnum
        Residue number. If there are multiple residues with same number but 
        distinguished with insertion codes, insertion code can be appended
        to the residue number. "_" stands for empty insertion code.
        
        Examples:
            
            * "resnum 5" selects residue 5 (all insertion codes)
            * "resnum 5A" selects residue 5 with insertion code A
            * "resnum 5\_" selects residue 5 with no insertion code
            * "resnum 5 10 to 15" selects residues 5, 10, 11, 12, 13, 14, and 15
            * "resnum 5 10:15" selects residues 5, 10, 11, 12, 13, and 14 (: works as it does in Python slicing operations)
            * "resnum 1:10:2" selects residues 1, 3, 5, 7, and 9
            
            
    
    sidechain
        Non-:term:`backbone` :term:`protein` atoms. Note that this
        defition includes backbone amide hydrogen.

    segment
        Group of atoms with same segment identifiers (segids).

    serial
        Atom numbers starting at 1.

    small
        Amino acid residues with small sidechains. Residue names include ALA, GLY, SER.
        This definition can be changed using :func:`setSmallResidueNames` 
        method.
        
    type
        Atom type (e.g. force field type).
    
    water
        Group of atoms whose residue names match one of HOH, WAT, TIP3, or H2O. This
        list may be expanded using :func:`setWaterResidueNames` 
        method.

    x
        x component of coordinates.
    
    y
        y component of coordinates.

    z
        z component of coordinates.

