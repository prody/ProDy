.. currentmodule:: prody.proteins

.. _water:

*******************************************************************************
Constructing AtomGroups
*******************************************************************************

This example shows how to construct an :class:`AtomGroup` 
instance from scratch. It is particularly useful for those who intend to design 
a molecular data file parser that returns parsed data in an :class:`AtomGroup`
instance.

.. literalinclude:: water.py
   
   
..    
    *******************************************************************************
    Constructing AtomGroups 2
    *******************************************************************************


    **Water AtomGroup Example**

    This example is useful if the user intends to design a molecular data file 
    parser using :class:`AtomGroup`.


    Instantiate an AtomGroup
    ===============================================================================

    >>> from prody import *
    >>> import numpy as np
    >>> 
    >>> # Instantiate an atom group with a name
    >>> wtr1 = AtomGroup('Water')
    >>> wtr1
    <AtomGroup: Water (0 atoms; 0 coordinate sets, active set index: 0)>


    Set/get atom attributes
    ===============================================================================

     
    >>> # Set coordinates
    >>> wtr1.setCoordinates( np.array( [ [1, 0, 0], [0, 0, 0], [0, 0, 1] ] ) )
    >>> wtr1
    <AtomGroup: Water (3 atoms; 1 coordinate sets, active set index: 0)>

    .. note::
       Atom numbers



    >>> # set atom names
    >>> wtr1.setAtomNames( ['H', 'O', 'H'] )
    >>> # data must always be passes in a list or an array
    >>> wtr1.setResidueNumbers( [1, 1, 1] )
    >>> wtr1.setResidueNames( ['WAT', 'WAT', 'WAT'] )
    >>> 
    >>> # Accessing data will return a copy of the data
    >>> print wtr1.getAtomNames()
    ['H', 'O', 'H']

    Individual atoms
    ===============================================================================

    Individual atoms are represented by instance of :class:`prody.proteins.atom.Atom`.

    **Iteration**

    >>> # Atoms in an atom group can be iterated over
    >>> for a in wtr1: print a
    ... 
    Atom H (index 0) from AtomGroup Water
    Atom O (index 1) from AtomGroup Water
    Atom H (index 2) from AtomGroup Water

    **Indexing**

    >>> # Atoms in an atom group can be accessed via indexing
    >>> a = wtr1[0]
    >>> # This returns an AtomGroup instance, which provides access to arbitrary 
    >>> # groups of atoms in an atom group
    >>> a
    <Atom: H from AtomGroup Water (index 0; 1 coordinate sets, active set index: 0)>
    >>> print a.getCoordinates()
    [1 0 0]


    Coordinate sets
    ===============================================================================

    >>> # Let's add a coordinate set to the atom group
    >>> wtr1.addCoordset( np.array( [ [0, 1, 0], [0, 0, 0], [0, 0, 1] ] ) )
    >>> # Note that number of coordinate sets is now 2, but active coordinate set index
    >>> # is still 0
    >>> wtr1
    <AtomGroup: Water (3 atoms; 2 coordinate sets, active set index: 0)>
    >>> # Active coordinate set incex can be changed for AtomGroups
    >>> a.setActiveCoordsetIndex(1)
    >>> a
    <AtomGroup from Water: Atom H (index 0) (1 atoms; 2 coordinate sets, active set index: 1)>
    >>> # Changing active coordinate set for an atom group, does not affect the 
    >>> # active coordinate set of the atom group
    >>> wtr1
    <AtomGroup: Water (3 atoms; 2 coordinate sets, active set index: 0)>
    >>> # Coordinates for the atom group will be returned from the active coordinate set
    >>> a.getCoordinates()
    array([0, 1, 0])

    **Iterations**

    >>> # coordinate sets can also be iterated over for atoms and atom groups
    >>> for xyz in a.iterCoordsets(): xyz
    ... 
    array([1, 0, 0])
    array([0, 1, 0])

    Clone atom groups
    ===============================================================================

    >>> # Now let's make another copy of this water
    >>> wtr2 = wtr1.copy()
    >>> wtr2
    <AtomGroup: Copy of Water (3 atoms; 2 coordinate sets, active set index: 0)>

    **Translate clone**

    >>> # let's translate the coordinates of wtr2 so that it does not overlap with wtr1
    >>> wtr2.setCoordinates( wtr2.getCoordinates() + 2 )
    >>> wtr2.getCoordinates()
    array([[3, 2, 2],
           [2, 2, 2],
           [2, 2, 3]])
    >>> # above operation only translated the coordinate set at index 0
    >>> wtr2.setActiveCoordsetIndex(1)
    >>> wtr2.getCoordinates()
    array([[0, 1, 0],
           [0, 0, 0],
           [0, 0, 1]])
    >>> wtr2.setCoordinates( wtr2.getCoordinates() + 2 ) # translate the second coordinate set as well

    **Change clone attributes**

    >>> # before we merge wtr1 and wtr2, let's change resid's of wtr2
    >>> wtr2.setResidueNumbers( [2, 2, 2] )
    >>> wtr2.getResidueNumbers()
    array([2, 2, 2])
    >>> # we can do this in an alternate way too
    >>> wtr2.select('all').setResidueNumbers(2)
    >>> wtr2.getResidueNumbers()
    array([2, 2, 2])
    >>> # note that the following won't work
    >>> wtr2.setResidueNumbers(2)
    Traceback (most recent call last):
      File "/usr/lib/python2.6/doctest.py", line 1248, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__[29]>", line 1, in <module>
        wtr2.resids = 2
      File "....", line 424, in set_resnums
        if len(resids) != self._n_atoms:
    TypeError: object of type 'int' has no len()

    Merge atom groups
    ===============================================================================

    >>> # let's merge two water atom groups
    >>> wtrs = wtr1 + wtr2
    >>> wtrs
    <AtomGroup: Water + Copy of Water (6 atoms; 2 coordinate sets, active set index: 0)>
    >>> wtrs.getCoordinates()
    array([[1, 0, 0],
           [0, 0, 0],
           [0, 0, 1],
           [3, 2, 2],
           [2, 2, 2],
           [2, 2, 3]])
    >>> print wtrs.getAtomNames()
    ['H' 'O' 'H' 'H' 'O' 'H']
    >>> wtrs.getResidueNumbers()
    [1 1 1 2 2 2]

    .. note::
       This hints to why :class:`AtomGroup` instead of Molecule is used. The entire 
       content of a PDB file is not a molecule in strict sense, even when so it is 
       not complete always. If it was, we don't store bond information anyhow. 
       We merely store coordinate, etc. data on some atoms in an AtomGroup.

    Hierarchical view
    ===============================================================================

    Hierarchical views of atom groups are represented by 
    :class:`prody.proteins.hierview.HierView`.

    >>> # Residues (and also chains) in an atom group can also be iterated over
    >>> for res in wtrs.getHierView().iterResidues(): print res
    <Residue: Residue WAT 1 (Water + Copy of Water; 3 atoms; 2 coordinate sets, active set index: 0)>
    <Residue: Residue WAT 2 (Water + Copy of Water; 3 atoms; 2 coordinate sets, active set index: 0)>
    >>>
    >>> # It's is also possible to change the name of AtomGroup "Water + Copy of Water"
    >>> wtrs.setName('2Waters')
    >>> wtrs
    <AtomGroup: 2Waters (6 atoms; 2 coordinate sets, active set index: 0)>
