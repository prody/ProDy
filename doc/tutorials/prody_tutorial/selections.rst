Atom Selections
===============================================================================

:class:`.AtomGroup` instances have a plain view of atoms for efficiency,
but they are coupled with a powerful atom selection engine.  You can get well
defined atom subsets by passing simple keywords or make rather sophisticated
selections using composite statements.  Selection keywords and grammar is very
much similar to those found in `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_.
Some examples are shown here:

Keyword selections
-------------------------------------------------------------------------------

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

.. ipython:: python

   structure = parsePDB('1p38')
   protein = structure.select('protein')
   protein

Using the "protein" keyword we selected 2833 atoms out of 2962 atoms.
:meth:`~.Atomic.select` method returned a :class:`.Selection` instance.
Note that all ``get`` and ``set`` methods defined for the :class:`.AtomGroup`
class are also defined for :class:`.Selection` class. For example:

.. ipython:: python

    protein.getResnames()


Select by name/type
-------------------------------------------------------------------------------

We select backbone atoms by passing atom names following "name" keyword:

.. ipython:: python

   backbone = structure.select('protein and name N CA C O')
   backbone


We can also use "backbone" to make the same selection.

We select acidic and basic residues by using residue names with
"resname" keyword:

.. ipython:: python

   charged = structure.select('resname ARG LYS HIS ASP GLU')
   charged

Alternatively, we can use predefined keywords "acidic" and "basic".

.. ipython:: python

   charged = structure.select('acidic or basic')
   charged
   set(charged.getResnames())

Composite selections
-------------------------------------------------------------------------------

Let's try a more sophisticated selection.  We first calculate the geometric
center of the protein atoms using :func:`.calcCenter` function.  Then, we
select the Cα and Cβ atoms of residues that have at least one atom within
10 Å away from the geometric center.

.. ipython:: python

   center = calcCenter(protein).round(3)
   center
   sel = structure.select('protein and name CA CB and same residue as '
                          '((x-1)**2 + (y-17.5)**2 + (z-40.0)**2)**0.5 < 10')
   sel

Alternatively, this selection could be done as follows:

.. ipython:: python

   sel = structure.select('protein and name CA CB and same residue as '
                          'within 10 of center', center=center)
   sel

Selection operations
-------------------------------------------------------------------------------

:class:`.Selection` instances can used with bitwise operators:

.. ipython:: python

   ca = structure.select('name CA')
   cb = structure.select('name CB')
   ca_or_cb = ca | cb
   ca_or_cb
   ca & cb # returns None, since there are no common atoms between the two

Selections simplified
-------------------------------------------------------------------------------

In interactive sessions, an alternative to typing in ``.select('protein')``
or ``.select('backbone')`` is using dot operator:

.. ipython:: python

   protein = structure.protein
   protein

You can use dot operator multiple times:

.. ipython:: python

   bb = structure.protein.backbone
   bb


This may go on and on:

.. ipython:: python

   ala_ca = structure.protein.backbone.resname_ALA.calpha
   ala_ca


More examples
-------------------------------------------------------------------------------

There is much more to what you can do with this flexible and fast atom
selection engine, without the need for writing nested loops with comparisons
or changing the source code.  See the following pages:

  * :ref:`selections` for description of all selection keywords
  * :ref:`selection-operations` for handy features of :class:`.Selection`
  * :ref:`contacts` for selecting interacting atoms