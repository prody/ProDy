Atom Selections
===============================================================================

This part gives more information on properties of :class:`.AtomGroup` objects.
We start with making necessary imports. Note that every documentation page
contains them so that the code within the can be executed independently.
You can skip them if you have already done them in a Python session.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Atom Selections
-------------------------------------------------------------------------------

:class:`.AtomGroup` instances have a plain view of atoms for efficiency,
but they are coupled with a powerful atom selection engine.  You can get well
defined atom subsets by passing simple keywords or make rather sophisticated
selections using composite statements.  Selection keywords and grammar are very
much similar to those found in `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_.
Some examples are shown here:

Keyword selections
^^^^^^^^^^^^^^^^^^

Now, we parse a structure. This could be any structure, one that you know
well from your research, for example.

.. ipython:: python

   structure = parsePDB('1p38')
   protein = structure.select('protein')
   protein

Using the ``"protein"`` keyword we selected 2833 atoms out of 2962 atoms.
:meth:`.Atomic.select` method returned a :class:`.Selection` instance.
Note that all ``get`` and ``set`` methods defined for the :class:`.AtomGroup`
objects are also defined for :class:`.Selection` objects. For example:

.. ipython:: python

    protein.getResnames()


Select by name/type
^^^^^^^^^^^^^^^^^^^

We can select backbone atoms by passing atom names following ``"name"`` keyword:

.. ipython:: python

   backbone = structure.select('protein and name N CA C O')
   backbone


Alternatively, we can use ``"backbone"`` to make the same selection:

.. ipython:: python

   backbone = structure.select('backbone')

We select acidic and basic residues by using residue names with
``"resname"`` keyword:

.. ipython:: python

   charged = structure.select('resname ARG LYS HIS ASP GLU')
   charged

Alternatively, we can use predefined keywords "acidic" and "basic".

.. ipython:: python

   charged = structure.select('acidic or basic')
   charged
   set(charged.getResnames())

Composite selections
^^^^^^^^^^^^^^^^^^^^

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

Selections simplified
^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^

There is much more to what you can do with this flexible and fast atom
selection engine, without the need for writing nested loops with comparisons
or changing the source code.  See the following pages:

  * :ref:`selections` for description of all selection keywords
  * :ref:`contacts` for selecting interacting atoms


.. _selection-operations:

Operations on Selections
-------------------------------------------------------------------------------

:class:`.Selection` objects can used with bitwise operators:

Union
^^^^^

Let's select β-carbon atoms for non-GLY amino acid residues, and
α-carbons for GLYs in two steps:

.. ipython:: python

   betas = structure.select('name CB and protein')
   len(betas)
   gly_alphas = structure.select('name CA and resname GLY')
   len(gly_alphas)

The above shows that the p38 structure contains 15 GLY residues.

These two selections can be combined as follows:

.. ipython:: python

   betas_gly_alphas = betas | gly_alphas
   betas_gly_alphas
   len(betas_gly_alphas)

The selection string for the union of selections becomes:

.. ipython:: python

   betas_gly_alphas.getSelstr()

Note that it is also possible to yield the same selection using selection
string ``(name CB and protein) or (name CA and resname GLY)``.


Intersection
^^^^^^^^^^^^

It is as easy to get the intersection of two selections. Let's find
charged and medium size residues in a protein:

.. ipython:: python

   charged = structure.select('charged')
   charged
   medium = structure.select('medium')
   medium

.. ipython:: python

   medium_charged = medium & charged
   medium_charged
   medium_charged.getSelstr()

Let's see which amino acids are considered charged and medium:

.. ipython:: python

   set(medium_charged.getResnames())

What about amino acids that are medium or charged:

.. ipython:: python

   set((medium | charged).getResnames())


Inversion
^^^^^^^^^

It is also possible to invert a selection:

.. ipython:: python

   only_protein = structure.select('protein')
   only_protein
   only_non_protein = ~only_protein
   only_non_protein
   water = structure.select('water')
   water

The above shows that 1p38 does not contain any non-water
hetero atoms.

Addition
^^^^^^^^

Another operation defined on the :class:`.Select` object is addition
(also on other :class:`.AtomPointer` derived classes).

This may be useful if you want to yield atoms in an :class:`.AtomGroup` in a
specific order.
Let's think of a simple case, where we want to output atoms in 1p38 in a
specific order:

.. ipython:: python

   protein = structure.select('protein')
   water = structure.select('water')
   water_protein = water + protein
   writePDB('1p38_water_protein.pdb', water_protein)

In the resulting file, the water atoms will precedes the
protein atoms.


Membership
^^^^^^^^^^

Selections also allows membership test operations:

.. ipython:: python

   backbone = structure.select('protein')
   calpha = structure.select('calpha')

Is :term:`calpha` a subset of :term:`backbone`?

.. ipython:: python

   calpha in backbone

Or, is water in protein selection?

.. ipython:: python

   water in protein

Other tests include:

.. ipython:: python

   protein in structure
   backbone in structure
   structure in structure
   calpha in calpha


Equality
^^^^^^^^

You can also check the equality of selections. Comparison will return
``True`` if both selections refer to the same atoms.

.. ipython:: python

   calpha = structure.select('protein and name CA')
   calpha2 = structure.select('calpha')
   calpha == calpha2

.. _attributes:

