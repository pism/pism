===============================
PISM Python Overview
===============================

PISM has two Python components: a set of bindings of a subset of 
the PISM C++ libraries, and a library of Python code with an emphasis
on script writing and on inversion algorithms. Both components are 
accessed by importing the ``PISM`` Python package.

For the most part, the use of PISM from Python is intuitive and familiar 
to a C++ PISM coder.  However, there are some additional features that
a PISM Python user should be aware of to simplify working with PISM.
The following is a minimal script that initializes PISM, creates an :cpp:class:`IceGrid`, and creates an :cpp:class:`IceModelVec2S` on the grid::

  import PISM
  
  context = PISM.Context()
  
  grid = context.newgrid()
  
  Lx = 10000.; Ly=10000;    # 10 km half-width grid
  Mx = 20; My = 20;         # 20x20 grid
  PISM.model.initShallowGrid(grid,Lx,Ly,Mx,My,PISM.XY_PERIODIC)
  
  vec = PISM.IceModelVec2S()
  stencil_width = 1
  vec.create(grid, "tauc", PISM.kHasGhosts, stencil_width)

It's worth walking through this example to see exactly what happens in the background.

  ::
  
    import PISM

  When the ``PISM`` package is imported, there is a check done to see if
  ``PETSc`` has already been initialized by the ``petsc4py`` package.  
  If not, then ``PETSc`` is initialized with the 
  system command line flags used to generate its options database.  You
  can override this behavior by initializing ``PETSc`` manually 
  before importing ``PISM``.

  ::
  
    context = PISM.Context()

  The PISM Python libraries introduce a new :class:`~PISM.Context` class that
  manages the following global information as public attributes:
  
    * An MPI Communicator ``com``, always ``PETSc.COMM_WORLD``.
    * MPI ``size`` and ``rank`` variables concerning the number of
      processors and our processor's number within the collection.
    * A :cpp:class:`NCConfigVariable`, ``config``, containing 
      global configuration parameters.
    
  There is only every one :class:`~PISM.Context` instance, and it can always
  be accessed by ``PISM.Context()``.  
  
  The :class:`~PISM.Context`\ 's ``config`` attribute is not created
  until it is first accessed.  By default, on first access it reads in 
  and merges a base and override config file specified by the flags 
  ``-config`` and ``-config_override`` by using the PISM C++
  :cpp:member:``init_config`` function, and then calls
  :cpp:member:``set_config_from_options`` to further overwrite some parameters
  from the command line.  This behavior can be modified somewhat by calling
  :meth:`PISM.Context.init_config` before first accessing the ``config``  
  attribute.
  
  ::
  
    grid = context.newgrid()
  
  A :class:`~PISM.Context` contains all the data needed to construct new 
  :cpp:class:`IceGrid` objects, so it contains a convenience function to do
  this.  This line is equivalent to ``grid = PISM.IceGrid(context.com,context.size,context.rank,context.config)``.
  
  ::
  
    Lx = 10000.; Ly=10000;    # 10 km half-width grid
    Mx = 20; My = 20;         # 20x20 grid
    PISM.model.initShallowGrid(grid,Lx,Ly,Mx,My,PISM.XY_PERIODIC)
  
  A number of explicit steps need to occur to fully initialize the grid.  
  The grid dimensions need to be specified, grid ownership needs to be
  partitioned between the processors, grid spacing needs to be 
  decided, and final allocation needs to occur.  To help with these
  steps, see: 
  
    * :func:`PISM.model.initShallowGrid`,
    * :func:`PISM.model.initGrid`,
    * :func:`PISM.model.initGridFromFile`.
  
  ::
  
    vec = PISM.IceModelVec2S()
    stencil_width = 1
    vec.create(grid, "tauc", PISM.kHasGhosts, stencil_width)
  
  Little needs to be said about this last snippet; it behaves exactly as the
  analogous C++ code
  
  .. code-block:: cpp
  
    IceModelVec2S vec;
    int stencil_width = 1;
    vec.create(grid,"tauc",kHasGhosts,stencil_width);



Options Database
----------------

The PISM Python library allows access to the PETSc options database 
through functions in the :mod:`PISM.options` module.
There is one key difference between C++ and Python code accessing options.
In C++ code, the ``PetscOptionsBegin`` and ``PetscOptionsEnd`` macros are
used to delimit a block of code that access a set of related command-line options that should be grouped and labeled nicely when using the ``-help``
command-line flag.  These macros are surprisingly technical (they implement a loop!), and there is a special idiom for their Python counterpart::

  com = PISM.Context().com
  for o in PISM.OptionsBegin(com,title="PISMI Options"):
    input_filename = PISM.optionsString("-i","input file")
    is_regional    = PISM.optionsFlag("-regional",
                       "Compute SIA/SSA using regional model semantics",
                       default=False)
    design_var     = PISM.optionsList(com,"-inv_design",
                       "design variable name", ['tauc','hardav']
                       default = 'tauc')

Note that the variable ``o`` is just a placeholder and can be named anything.
The :class:`~PISM.OptionsBegin` class is a Python iterator that ensures the
the loop repeats the requisite number of times, and Python syntax requires
the ``o``.

PETSc Error Codes
-----------------

Numerous PISM functions and methods return :cpp:type:`PetscErrorCode` values
to indicate error conditions.  These appear on the Python side as exceptions of type :class:`petsc4py.PETSc.Error`, so no explicit checks are needed.

Note that PISM is not exception safe.  In most circumstances, where Python is calling PISM code, this is not an issue; exceptions on the Python side (including :class:`petsc4py.PETSc.Error`\ s) never
cause a change in control flow on the C++ side.  However, if Python code 
is called from C++ via a callback, it is not safe to propagate a 
Python exception into a C++ exception. Callbacks in the :mod:`PISM.invert`
library are wrapped to to catch any exceptions that might be raised and print 
an error message instead.

:cpp:class:`IceModelVec` Access and Communication
-------------------------------------------------

Access to the contents of :cpp:class:`IceModelVec`\ s needs to be 
managed in parallel runs through
:cpp:member:`IceModelVec::begin_access`/:cpp:member:`IceModelVec::end_access` pairs.  Ghost communication is done with :cpp:member:`update_ghosts`.  This
leads to code blocks of the form

.. code-block:: cpp

  err=vec.begin_access(); CHKERRQ(err);
  
  /******************************************/
  /*  Do stuff with vec here                */
  /******************************************/
  
  err=vec.end_access(); CHKERRQ(err);
  err=vec.update_ghosts(); CHKERRQ(err); // If needed
  
Omitting either the :cpp:member:`end_access` or :cpp:member:`update_ghosts`
can lead to bugs, and the distance of these calls from the initial 
:cpp:member:`begin_access` can make it easy to forget these calls.

The Python :class:`PISM.vec.Access` class allows specification at
the start of a code block of the desired actions to take at the
end of the block. Assuming that ``vg1`` and ``vg2`` are vectors with
ghosts and ``w`` is unghosted::

  grid = w.get_grid()
  with PISM.vec.Access(comm=vg2,nocomm=[vg1,w]):
    for (i,j) in grid.points():
      vg2(i,j) = 2*vg1(i,j) - w(i,j)

On entry into the ``with`` block, there is a call to
:cpp:member:`begin_access` for all three vectors.
On exit, there is a 
call to :cpp:member:`end_access` for all three vectors,
and a call to :cpp:member:`update_ghosts` for ``vg2``.

Note the call to :meth:`grid.points`, which returns an iterator
over all grid points owned by the current processor.  We can 
iterate over all points including ghosted nodes with
:meth:`grid.points_with_ghosts`.

Of course, loops like this in Python will be much slower than in C++.
For many scripts, however, this is only a mild nuisance.
For future reference, Cython might have been a superior option to SWIG for
building the Python bindings, and would have given a good mechanism for
making fast loops in end-user code.


Other Bits
-------------------------------------------------

* The :mod:`PISM.model` module contains a number of methods for creating
  :cpp:class:`IceModelVec`\ s with metadata setup for standard PISM variables.
  See, e.g., :func:`PISM.model.createIceSurfaceVec` and its cousins.

* The :class:`PISM.vec.ToProcZero` class facilitates communication of 
  :cpp:class:`IceModelVec`\ s to processor zero and conversion into ``numpy``
  arrays.
  
* There is a basic message logging system in :mod:`PISM.logging` that 
  allows messages to be captured and stored in :file:`.nc` files.



