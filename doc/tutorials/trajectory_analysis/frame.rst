.. _frame:

Frames and Atom Groups
===============================================================================

This part shows how to use :class:`.AtomGroup` in place of :class:`.Frame`.


Trajectory Frames
-------------------------------------------------------------------------------

:class:`.Frame` instances store only coordinate and some frame related data.
For each frame data, you will a different frame instance:

.. ipython:: python

   from prody import *
   dcd = Trajectory('trajectory_analysis_files/mdm2.dcd')
   dcd
   frame0 = dcd.next()
   frame0
   frame1 = dcd.next()
   frame1

These :class:`.Frame` instances are different objects:

.. ipython:: python

   frame0 is frame1

When you are not referring to any of these frames anymore in your code,
Python garbage collector will free or reuse the memory space that was used by
those frames.


Linking Atom Groups
-------------------------------------------------------------------------------

When trajectory is not linked to an :class:`.AtomGroup` (using
:meth:`~Trajectory.link`), :class:`.Frame` and :class:`.AtomGroup` objects
share the same coordinate data. Let's see how this works:

When an :class:`.AtomGroup` is linked to the trajectory as follows, things
work differently:

.. ipython:: python

   pdb = parsePDB('trajectory_analysis_files/mdm2.pdb')
   pdb
   dcd.link(pdb)
   dcd.reset()

We get :class:`.Frame` instances in the same way:

.. ipython:: python

   frame0 = dcd.next()
   frame0
   pdb.getACSLabel()

Note that the active coordinate set of the :class:`.AtomGroup` and its label
will change when we get the next frame:

.. ipython:: python

   frame1 = dcd.next()
   frame1
   pdb.getACSLabel()

Now the key difference is that the :class:`Frame` instances are the same
objects in this case:

.. ipython:: python

   frame0 is frame1

As you see, a new frame was not instantiated.  The same frame is reused and
it always points to the coordinates stored in the :class:`.AtomGroup`.
You can also make :class:`.Selection` instances that will point to the same
coordinate set.  This will allow making a more elaborate analysis of frames.
For an example see :ref:`trajectory2`.