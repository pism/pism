.. include:: ../global.txt

.. default-role:: literal

How do I...?
============

.. contents::

Create and use configuration flags and parameters
-------------------------------------------------

- Edit |config-cdl|. Each flag or parameter is stored as a NetCDF attribute and should
  have a corresponding "`_doc`" attribute describing its meaning and the "`_type`"
  defining the type. All scalar parameters have to have "`_units`" set as well.

.. code-block:: none

   pism_config:constants.standard_gravity = 9.81;
   pism_config:constants.standard_gravity_doc = "acceleration due to gravity on Earth geoid";
   pism_config:constants.standard_gravity_type = "scalar";
   pism_config:constants.standard_gravity_units = "meter second-2";

- One can access these parameters using the `Config` class. `IceModel` and all classes
  derived from `Component` have an pointer to an instance of this class as a data member
  `m_config`, so no additional code is necessary to initialize the configuration database.

To use a scalar parameter, do

.. code-block:: c++

   double g = m_config->get_number("constants.standard_gravity");

To use a flag, do

.. code-block:: c++

   bool compute_age = config->get_flag("age.enabled");

.. note::

   - It is best to avoid calling `m_config->get_...()` from within loops: looking up a
     parameter by its name is slow.
   - Please see :ref:`sec-parameter-list` for a list of flags and parameters currently
     used in PISM.

Create and use additional variables
-----------------------------------

Creating IceModelVec instances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM uses the following classes to manage 2D and 3D fields, their I/O and metadata:

.. list-table::

   * - `IceModelVec2S`
     - scalar 2D fields
   * - `IceModelVec2V`
     - vector 2D fields such as horizontal velocities; corresponds to 2 NetCDF variables
   * - `IceModelVec2Int`
     - 2D masks, such as the grounded/floating mask
   * - `IceModelVec2T`
     - 2D time-dependent fields (used to read and store forcing data)
   * - `IceModelVec3`
     - scalar 3D fields (usually within the ice)

Please see the documentation of these classes for more info. The base class `IceModelVec` is
a virtual class, so code should use the above derived classes.

To create a scalar field, for example, one needs to create an instance of `IceModelVec2S`
and then set its metadata:

.. code-block:: c++

   // land ice thickness
   IceModelVec2S ice_thickness(grid, "thk", WITH_GHOSTS, 2);
   ice_thickness.set_attrs("model_state", "land ice thickness",
                           "m", "land_ice_thickness");
   ice_thickness.metadata().set_number("valid_min", 0.0);

Here `grid` is an `IceGrid` instance, `thk` is the name of the NetCDF variable,
`WITH_GHOSTS` means that storage for "ghost" ("halo") points will be allocated, and "2" is
the number of ghosts (in other words: required stencil width).

The `IceModelVec::set_attrs()` call sets commonly used NetCDF variable attributes seen in
PISM output files:

.. list-table::

   * - `pism_intent`
     - variables that are a part of the model state of a sub-model should have
       `pism_intent` set to "model_state"
   * - `long_name`
     - the (descriptive) long name used for plotting, etc (a free-form string)
   * - `units`
     - units used *in the code*; does not have to match units in a file.
   * - `standard_name`
     - CF standard name, if defined, or an empty string.

The `ice_thickness.metadata()` call above allows accessing variable metadata and
adding arbitrary attributes. See `VariableMetadata` for details.

The CF convention covers some attribute semantics, including `valid_min` in this example.

PISM will automatically convert units from ones present in an input file into internal
units defines by the `set_attrs()` call above.

If you want PISM to save data in units other than internal ones, first set these
"glaciological" units:

.. code-block:: c++

   ice_thickness.metadata().set_string("glaciological_units", "km");

Read data from a file
---------------------

There are at least three cases of "reading data from a file":

1. reading a field stored in an input file on a grid matching the one used by the current
   run (restarting a run),
2. reading a field stored in an input file on a different grid, interpolating onto the
   current grid (bootstrapping), and
3. reading a field stored in a file **other** than the input file using interpolation
   (assuming that grids are compatible but not identical)

.. code-block:: c++

   // case 1, using a file name
   unsigned int time = 0;
   std::string filename = "filename.nc";
   ice_thickness.read(filename, time);

   // case 1, using an existing File instance (file is already open)
   File file(communicator, "guess_mode", filename, PISM_READONLY);
   ice_thickness.read(file, time);

   RegriddingFlag flag = OPTIONAL;
   double default_value = 0.0;
   // cases 2 and 3 (interpolation)
   ice_thickness.regrid(filename, flag, default_value);

   // cases 2 and 3 (interpolation) using an existing File instance
   ice_thickness.regrid(file, flag, default_value);

When interpolating ("regridding") a field, the ``flag`` specifies whether a variable is
required (``flag`` is ``CRITICAL`` or ``CRITICAL_FILL_MISSING``) or optional (``flag`` is
``OPTIONAL`` or ``OPTIONAL_FILL_MISSING``). PISM will stop with an error message if a
required variable is not found in an input file.

Flag values of ``CRITICAL_FILL_MISSING`` and ``OPTIONAL_FILL_MISSING`` replaces "missing"
values matching the ``_FillValue`` attribute by the default value.

If ``flag`` is ``OPTIONAL`` or ``OPTIONAL_FILL_MISSING`` PISM will fill the variable with
``default_value`` if it was not found in the file.

Write data to a file
--------------------

`IceModelVec::define()` will define all spatial dimensions used by a variable. However, we
usually need to "prepare" a file by defining the time dimension and appending a value to
the time variable.

.. code-block:: c++

   File file(m_grid->com, m_config->get_string("output.format"),
            filename, PISM_READWRITE_CLOBBER, m_grid->ctx()->pio_iosys_id());

   io::define_time(file, *m_grid->ctx());
   io::append_time(file, *m_grid->ctx()->config(), current_time);

When a file is opened with the mode `PISM_READWRITE_CLOBBER`, PISM checks if this file is
present overwrites if it is; to append to an existing file, use `PISM_READWRITE`. To move
the file aside (appending "`~`" to the file name), use `PISM_READWRITE_MOVE`.

A newly-created file is "empty" and contains no records. The `io::append_time()` call
creates a record corresponding to the `current_time` (in seconds).

To write a field to an already "prepared" file, call

.. code-block:: c++

   precip.write("filename.nc");
   // or, if the file is already open and a File instance is available:

   precip.write.(file);

Read scalar forcing data
------------------------

PISM uses instances of the `ScalarForcing` class to read scalar forcing data; please see
`pism::surface::Delta_T` for an example.


Read 2D forcing fields
----------------------

PISM uses instances of the `IceModelVec2T` class to read 2D forcing fields that
vary in time; please see `pism::surface::Given` for an example.
