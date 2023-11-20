Changes from 0.6 (February 2014) to 0.7
=======================================

Most of the work on PISM since the release of v0.6 was aimed at making
PISM easier to maintain, test, extend, and couple to other models. These
changes consist of approximately 1200 individual Git commits. They
touched over 650 files, more than half of PISM's C++ code (44,000 out of
84,000 lines).

All changes listed below are included in v0.7 (May 12, 2015) and later.
Some changes listed below were included in v0.6.1 or v0.6.2.

User-visible changes
~~~~~~~~~~~~~~~~~~~~

Installing PISM
^^^^^^^^^^^^^^^

-  Support PETSc 3.5.
-  Support PETSc ``--with-64-bit-indices=1``
-  Require PETSc >= 3.5 for TAO-based code.
-  Require FFTW >= 3.1 to limit the amount of time spent choosing DFT
   algorithms.
-  Allow building PISM with GSL <= 1.15
   (`#304 <https://github.com/pism/pism/issues/304>`__).
-  Updated installation instructions for Cray systems
   (`#316 <https://github.com/pism/pism/issues/316>`__).
-  Put quick installation instructions in one spot
   (`#314 <https://github.com/pism/pism/issues/314>`__).
-  Allow building documentation on systems without full PISM
   prerequisites (`#251 <https://github.com/pism/pism/issues/251>`__).
   Give better warnings about missing tools needed to build the source
   code browser (`#137 <https://github.com/pism/pism/issues/137>`__).

New features
^^^^^^^^^^^^

-  Implement ``KirchhoffEnthalpyConverter``, an enthalpy converter which
   approximates L, the latent heat of fusion of ice, using Kirchhoff's
   law of thermochemistry.
-  Implement ``-atmosphere weather_station``. Reads scalar time-series
   of near-surface air temperature and precipitation and applies them to
   the whole domain. Use with lapse rate corrections to add spatial
   variability.
-  Re-implement and document the ocean model which uses the 3-equation
   sub-shelf melting parameterization (Hellmer and Olbers,
   ``-ocean th``).
-  The PDD model supports a spatially-variable standard deviation
   (`#179 <https://github.com/pism/pism/issues/179>`__) field used to
   model temperature variability and a parameterization of this standard
   deviation (`#265 <https://github.com/pism/pism/issues/265>`__).
-  Add ``-atmosphere ...,frac_P``, an atmosphere "modifier" that scales
   precipitation using a time-dependent factor read from a file.
   (`#271 <https://github.com/pism/pism/issues/271>`__).
-  Add a PETSc-based parallel version of the ``fill_missing`` script
   which can be used to fill gaps in high-resolution datasets.

New diagnostics
^^^^^^^^^^^^^^^

-  Add SIA-type shear stresses (``tauxz``, ``tauyz``) and hydrostatic
   pressure (``pressure``)
   (`#280 <https://github.com/pism/pism/issues/280>`__).
-  Add the bed parallel basal shear stress (``taub``) and its magnitude
   (``taub_mag``) (`#266 <https://github.com/pism/pism/issues/266>`__).
-  New names of vector diagnostic quantities:

   -  ``cbar`` was renamed to ``velbar_mag``,
   -  ``cbase`` was renamed to ``velbase_mag``,
   -  ``csurf`` was renamed to ``velsurf_mag``,
   -  ``cflx`` was renamed to ``flux_mag``.

-  Mass-conserving hydrology models add conservation-related scalar
   diagnostics (`#256 <https://github.com/pism/pism/issues/256>`__).
-  Add ``flux_divergence``
   (`#165 <https://github.com/pism/pism/issues/165>`__).
-  Add ``uflux``, ``vflux``, 3D horizontal ice fluxes in the X and Y
   direction.
-  Add hydrology diagnostics and CTS to ``-o_size big``
   (`#264 <https://github.com/pism/pism/issues/264>`__) and
   (`#262 <https://github.com/pism/pism/issues/262>`__).

Some changes that may break run scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Replace ``-boot_file foo.nc`` with ``-bootstrap -i foo.nc``
   (`#308 <https://github.com/pism/pism/issues/308>`__).
-  Rename ``-force_to_thk`` to ``-force_to_thickness``
   (`#268 <https://github.com/pism/pism/issues/268>`__).
-  ``-atmosphere searise_greenland`` uses
   ``-atmosphere_searise_greenland_file``
   (`#263 <https://github.com/pism/pism/issues/263>`__).
-  ``-surface ...,forcing`` requires ``ftt_mask``; floating ice is not
   affected.
-  Remove automatic vertical grid extension.

Minor changes
^^^^^^^^^^^^^

-  Remove the "preliminary" time step
   `#210 <https://github.com/pism/pism/issues/210>`__.
-  ``-topg_to_phi`` arguments are configuration parameters
   (`#289 <https://github.com/pism/pism/issues/289>`__).
-  Add ``energy_advection_ice_thickness_threshold``
   (`#292 <https://github.com/pism/pism/issues/292>`__).
-  Allow using a different ``ftt_alpha`` in "ice-free" areas. (`commit
   93c09d3d <https://github.com/pism/pism/commit/93c09d3d>`__)
-  Close `#254 <https://github.com/pism/pism/issues/254>`__. Adds
   ``-timestep_hit_multiples X``.
-  Close `#255 <https://github.com/pism/pism/issues/255>`__. (add
   ability to start ``-hydrology distributed`` model if ``bwp`` is
   missing)
-  Close `#259 <https://github.com/pism/pism/issues/259>`__ (new
   adaptive time-stepping decision output).
-  Fix discharge reporting (now it's the same for all calving
   mechanisms, cumulative only in 2D, cumulative and a "rate" as a
   scalar time-series).
-  Include ``pism_config`` in output files
   (`#270 <https://github.com/pism/pism/issues/270>`__).
-  Pressure is set to overburden in grounded-ice areas with empty
   subglacial aquifer, a change of "``-hydrology distributed``\ "
   behavior relative to v0.6 which avoids "sucking" water out at the
   margin. (The runs which went into Bueler & van Pelt (2015) used bmelt
   = -1 m/a to empty ice-free, but this is not necessary now.)
-  Remove Storglaciaren example.
-  Separate Glen exponents for SIA and SSA flow laws
   (`#285 <https://github.com/pism/pism/issues/285>`__).
-  Use latitude and longitude bounds names compatible with CDO.
-  Use the global attribute "``proj4``\ " instead of
   "``mapping:proj4``\ ".

Some bug fixes
~~~~~~~~~~~~~~

-  `#267 <https://github.com/pism/pism/issues/267>`__ (ensure that
   threshold thickness is non-negative).
-  `#277 <https://github.com/pism/pism/issues/277>`__ (the time axis
   length of ``ts_times``).
-  `#278 <https://github.com/pism/pism/issues/278>`__ (``-energy none``
   should start with either ``temp`` or ``enthalpy``).
-  `#281 <https://github.com/pism/pism/issues/281>`__ (a bug related to
   (now-removed) vertical grid extension).
-  `#283 <https://github.com/pism/pism/issues/283>`__ (unreasonable
   SSAFD velocities near the "cliff" at periodic boundary).
-  `#297 <https://github.com/pism/pism/issues/297>`__ (record bed
   deformation choice in a configuration parameter).
-  `#299 <https://github.com/pism/pism/issues/299>`__ (fix verbPrintf).
-  `#315 <https://github.com/pism/pism/issues/315>`__ (segfault in the
   PDD model).
-  unreported: Fix snow depth reset code in the PDD model `commit
   cae55774 <https://github.com/pism/pism/commit/cae55774>`__.

Code changes that may be of interest to developers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the changes in this category come from trying to make PISM more
library-like, which in the long run will make it easier to extend, test,
maintain, and couple to other models. This transformation is not nearly
done, but here's the current state of the code. I am not claiming that
this is great design, so please do let me know if something looks broken
to you.

Use exceptions to handle errors in PISM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM v0.7 uses exceptions to propagate information about runtime errors.

Use the ``RuntimeError`` class to signal an error:

::

        throw RuntimeError("message");
        // or
        throw RuntimeError::formatted("format string %s", "data");

Sometimes it is helpful to add context to an error message so that a
user can get more information about a failure. Here's a way to do that:

::

        try {
          // code that may fail
          foo();
        } catch (RuntimeError &e) {
          e.add_context("doing foo to %s", "data");
          throw;                        // don't forget to re-throw!
        }

Some benefits of using exceptions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  We can now use function return values instead of output arguments in
   most places.
-  Better resource management. (No half-allocated objects because we can
   allocate in constructors and report allocation errors.)
-  PISM code is easier to wrap with SWIG; the SWIG interface file is
   much easier to maintain.
-  Error propagation from PISM (C++) to Python and back **works**.
-  PISM's C++ code can be tested using Python scripts.

Errors in parallel code sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a computation fails on some (but not all) ranks in a communicator,
the next blocking MPI call will, well, block. This can make PISM hang
instead of stopping with an error message.

We try to prevent this by wrapping all ``i,j`` loops in

::

        ParallelSection loop(communicator);
        try {
          // for loop goes here
        } catch (...) {
          loop.failed();
        }
        loop.check();

``loop.failed()`` prints an error message and sets a flag indicating a
failure. Then ``loop.check()`` calls ``MPI_Allreduce()`` to tell **all**
ranks in the communicator that something failed and then throws an
exception on **all** ranks.

All loops containing code that might throw should be wrapped this way.

**Note:** This problem exists regardless of the error handling method in
use.

IceModelVec
^^^^^^^^^^^

First of all, ``const IceModelVec`` is usable now, so it's OK to return
a const reference to an internal field in a component's interface to
provide read-only access.

IceModelVec::AccessList
^^^^^^^^^^^^^^^^^^^^^^^

Accessing PETSc ``Vec`` arrays requires wrapping code in
``DMDAVecGetArray`` and ``DMDAVecRestoreArray`` calls. This is a problem
if we assume that all code can throw an exception.

::

        DMDAVecGetArray(...);
        // if do_work(...) throws, DMDAVecRestoreArray will not be called.
        do_work(i, j);
        DMDAVecRestoreArray(...);

This issue affects ``IceModelVec`` fields, too.

To get around this I created ``IceModelVec::AccessList`` which calls
``DMDAVecGetArray`` in the constructor and ``DMDAVecRestoreArray`` in
the destructor. This guarantees that ``DMDAVecRestoreArray`` is called
when we exit the scope. Here's how to use it:

::

        IceModelVec::AccessList list(f);
        list.add(g);

        f(i, j) = g(i, j);

Accessing "raw" PETSc Vecs
^^^^^^^^^^^^^^^^^^^^^^^^^^

To access "raw" PETSc Vecs and avoid the risk of not calling
``...RestoreArray...()``, use

-  ``pism::petsc::VecArray`` to access a ``Vec`` by calling
   ``VecGetArray`` and ``VecRestoreArray``,
-  ``pism::petsc::VecArray2D`` to use ``VecGetArray2d`` and
   ``VecRestoreArray2d`` calls,
-  ``pism::petsc::DMDAVecArray`` to use ``DMDAVecGetArray`` and
   ``DMDAVecRestoreArray``,
-  ``pism::petsc::DMDAVecArrayDOF`` to use ``DMDAVecGetArrayDOF`` and
   ``DMDAVecRestoreArrayDOF``.

New ``IceModelVec2S`` methods for moving data to/from processor 0.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``IceModelVec2S`` takes care of allocating copies on processor 0 and
communication.

See this code from ``IcebergRemover`` for an example:

::

        // In the constructor: allocate a copy on processor 0
        m_mask_p0 = m_iceberg_mask.allocate_proc0_copy();

        // Later: identify icebergs using serial code on processor 0:
        {
          m_iceberg_mask.put_on_proc0(*m_mask_p0);

          ParallelSection rank0(m_grid->com);
          try {
            if (m_grid->rank() == 0) {
              petsc::VecArray mask(*m_mask_p0);
              cc(mask.get(), m_grid->Mx(), m_grid->My(), true, mask_grounded_ice);
            }
          } catch (...) {
            rank0.failed();
          }
          rank0.check();

          m_iceberg_mask.get_from_proc0(*m_mask_p0);
        }

``IceModelVec::has_nan()`` is gone
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use PETSc's option ``-fp_trap`` to detect ``NaNs`` (on Linux).

Wrappers for all PETSc objects used in PISM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM has wrappers for all PETSc objects it uses: ``DM``, ``IS``,
``KSP``, ``Mat``, ``SNES``, ``Tao``, ``Vec``, ``VecScatter``,
``Viewer``.

These wrappers ensure that the corresponding ``...Destroy()`` function
is called when a ``DM``, ``Mat``, etc needs to be destroyed.

To use, write code similar to

::

        pism::petsc::Vec my_vec;
        ierr = VecCreateSeq(PETSC_COMM_SELF, size, my_vec.rawptr());
        PISM_CHK(ierr, "VecCreateSeq");
        // my_vec will be destroyed automatically when we exit the scope

Making PISM more library-like
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execution context ``Context`` (`#150 <https://github.com/pism/pism/issues/150>`__).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An execution context ``pism::Context`` collects objects that are shared
by several (possibly all) components of a PISM instance:

-  an MPI communicator
-  a UDUNITS unit system
-  a ``Config`` instance
-  an ``EnthalpyConverter`` instance
-  a ``Time`` manager object
-  a ``Profiling`` object
-  a ``Logger``

Putting them in a ``Context`` instead of using global objects will allow
running more than one PISM instance at the same time while preserving
consistency of modeling choices.

(Using PETSc's command-line options is still a problem, but this is a
step in the right direction.)

I imagine that we may have one or more ``Context``, with one or more
``IceGrid`` for each ``Context``, with multiple components for each
``IceGrid`` instance.

The ``pism::Logger`` class.
^^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM v0.6 uses the ``verbPrintf()`` function to print text to standard
out.

This may be a problem: - If two or more PISM instances run at the same
time in the same MPI process their outputs will get intermixed, which
will make this output useless. - If PISM is used as a library in a
bigger system we may want to suppress or redirect PISM's output.

To address this issue I created ``pism::Logger``, a simple class
wrapping ``verbPrintf``. Its default implementation does not add
anything new, but makes it possible to replace writing to ``stdout``
with writing to a file, for example. (Just write a class derived from
``Logger`` and use it to create a ``Context`` instance.)

New code for processing command-line options.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM v0.7 has new classes

-  ``pism::options::String``: an option with a string argument
-  ``pism::options::StringList``: an option with the argument which is a
   comma-separated list of strings, used as a vector of strings
-  ``pism::options::StringSet``: an option with the argument which is a
   comma-separated list of strings, used as a set of strings
-  ``pism::options::Keyword``: an option which takes one of pre-defined
   keywords as an argument
-  ``pism::options::Integer``: an option with an integer argument
-  ``pism::options::IntegerList``: takes a comma-separated list of
   integers, returned as a vector
-  ``pism::options::Real``: an option with an numerical argument
-  ``pism::options::RealList``: takes a comma-separated list of numbers,
   returned as a vector of doubles
-  ``pism::options::Bool``: is a function that returns ``true`` if an
   option was set and ``false`` if it was not (or if ``-no_...`` was
   set).

Here's a way to use ``options::Integer``, for example.

::

        int default_value = 100;
        options::Integer N("-N", "description", default_value);

        if (N.is_set()) {
          // -N was set
          int M = N + 1;                // N is automatically converted to int
         }

``Real`` is automatically converted to ``double``, ``String`` and
``Keyword`` to ``std::string``, ``StringList``, ``IntegerList``,
``RealList`` to ``std::vector`` of strings, integers, and doubles.

Config improvements
^^^^^^^^^^^^^^^^^^^

The ``pism::Config`` is an interface class now. It was re-written so as
to reduce dependencies on the rest of PISM.

It should now be easy to create a ``Config`` derived classes that get
parameter values from a model PISM is coupled to, for example.

PISM automatically processes command-line options corresponding to
configuration parameters: ``pism_config.cdl`` provides all the
information needed to associate a configuration parameter with an option
and process this command-line option (if ``set_config_from_options()``
is called):

For example:

::

        pism_config:bed_deformation_model_type = "keyword";
        pism_config:bed_deformation_model_option = "bed_def";
        pism_config:bed_deformation_model_choices = "none,iso,lc";
        pism_config:bed_deformation_model = "none";
        pism_config:bed_deformation_model_doc = "Selects a bed deformation model to use...";

Each configuration parameter has a corresponding command-line option,
either the one specified using ``..._option`` or the one that matches
the parameter name.

Overhaul pism::Vars
^^^^^^^^^^^^^^^^^^^

No need for unsightly ``dynamic_casts`` when getting a field from this
dictionary.

::

        const IceModelVec2S *field = vars.get_2d_scalar("field_name");

This throws an exception if a field is not present; use
``is_available()`` to check first.

Overhaul IceGrid
^^^^^^^^^^^^^^^^

Now ``IceGrid`` contains grid information only and cannot be changed
once it is allocated. (Previously ``IceGrid`` was a kind of a
"context".)

Each ``IceGrid`` instance still has a ``pism::Vars`` instance: fields
(variables) are stored on a particular grid and so correspond to this
grid.

Previously ``IceModel`` was responsible for getting grid parameters from
a file or command-line options and allocating a grid; now we can
allocate an ``IceGrid`` using one of these methods:

-  Fill all members of ``GridParameters`` and use the constructor of
   ``IceGrid``.
-  Create a shallow grid using ``IceGrid::Shallow`` (a static method).
-  Create a grid by getting its parameters from a variable in a file (or
   using a variable from a list of candidates) with
   ``IceGrid::FromFile``.
-  Create a grid by processing command-line options ``-i``,
   ``-bootstrap``, ``-Mx``, ``-My``, ``-Mz``, ``-Lx``, ``-Ly``, ``-Lz``,
   ``-x_range``, ``-y_range``, and ``-periodicity`` by calling
   ``IceGrid::FromOptions``. (This is what ``pismr`` does.)

(This makes it easier to run PISM's sub-models independently from
``IceModel``.)

Point iterators ``Points`` and ``PointsWithGhosts``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To simplify iterating over the part of the computational domain of the
current processor PISM has an iterator ``Points``. This code is
equivalent to the double for loop but hides the grid traversal order.

::

    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      field(i,j) = value;
    }

To update ghost points locally, use ``PointsWithGhosts``:

::

    for (PointsWithGhosts p(grid, ghost_width); p; p.next()) {
      const int i = p.i(), j = p.j();
      field(i,j) = value;
    }

Other
~~~~~

-  Add ``make retest`` (re-run failed tests), ``make test-python`` (run
   Python tests only), ``make pylint`` (run ``pylint``).
-  Reduce the number of ``-I`` flags needed to build code that uses PISM
   as a library.
-  PISM tries not to terminate the run by calling ``MPI_Abort()`` and
   such.

   Note: ``PISMEnd`` in PISM v0.6 called ``exit()``, and according to
   the C++ standard ``exit()`` does not unwind the stack, i.e. cleanup
   done in destructors never happens. So, we should avoid ``exit()``.

Less interesting internal changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Minor improvements
^^^^^^^^^^^^^^^^^^

-  Effective viscosity under-relaxation as a first recovery strategy in
   SSAFD (`#282 <https://github.com/pism/pism/issues/282>`__).
-  Isolate the basal melt rate computation
   (`#99 <https://github.com/pism/pism/issues/99>`__).
-  Moving towards stand-alone bed deformation model runs
   (`#181 <https://github.com/pism/pism/issues/181>`__).
-  Skip bootstrapping heuristics whenever possible
   (`#291 <https://github.com/pism/pism/issues/291>`__).
-  Consistent metadata in NetCDF calls
   (`#151 <https://github.com/pism/pism/issues/151>`__).
-  Re-use PETSc DMs to reduce memory footprint
   (`#132 <https://github.com/pism/pism/issues/132>`__).
-  Use GSL binary search to speed up vertical index lookup.
-  ``PIO``: make it possible to overwrite existing files without
   creating a backup copy
   (`#224 <https://github.com/pism/pism/issues/224>`__).
-  Class ``MaxTimestep`` helps with comparing time-step restrictions.
-  Class ``Profiling`` helps use PETSc's profiling; add ``-log_summary``
   to see results.
-  Use automatically-generated docstrings in Python bindings.

Regression testing
^^^^^^^^^^^^^^^^^^

-  Use ``nose`` and ``coverage`` Python modules to test PISM's Python
   code and keep track of code coverage.
-  Better regression tests
   (`#305 <https://github.com/pism/pism/issues/305>`__).
-  Add unit (regression) tests for enthalpy converters
   (`#272 <https://github.com/pism/pism/issues/272>`__).
-  Isolate vertical interpolation in the column (linear and quadratic;
   see ``ColumnInterpolation``). and horizontal linear interpolation
   code (``LinearInterpolation``); add regression tests.
-  Improve regression test coverage
   (`#294 <https://github.com/pism/pism/issues/294>`__).

Code organization and cleanup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Give less confusing names to flow law classes.
-  Remove some uses of STL ``upper_bound`` and ``lower_bound``
   (`#170 <https://github.com/pism/pism/issues/170>`__).
-  Replace ``NCSpatialVariable`` and ``NCVariable`` with
   ``VariableMetadata`` and ``SpatialVariableMetadata``.
-  Don't pass classes by value
   (`#269 <https://github.com/pism/pism/issues/269>`__).
-  Use namespaces to organize PISM code (see
   `#188 <https://github.com/pism/pism/issues/188>`__).
-  Reduce inter-dependencies of PISM code.
-  Code formatting (see
   ```coding_style.md`` <https://github.com/pism/pism/blob/dev/doc/browser/coding_style.md>`__).
-  ``autopep8`` the Python code, ``pylint`` support.
-  Clean up FEM code.
