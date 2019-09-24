.. include:: ../global.txt

How do I...?
============

.. default-role:: literal

.. contents::

Create and use configuration flags and parameters
-------------------------------------------------

- Edit |config-cdl|. Each flag or parameter is stored as a NetCDF attribute and
  should have a corresponding "_doc" attribute describing its meaning.

.. code-block:: none

   pism_config:standard_gravity = 9.81;
   pism_config:standard_gravity_doc = "m s-2; acceleration due to gravity on Earth geoid";

- One can access these parameters using the `Config` class. `IceModel`\ [#]_ has an
  instance of this class as a data member `m_config`, so no additional code is necessary to
  initialize the configuration database.

To use a parameter, do

.. code-block:: c++

   double g = m_config->get_double("standard_gravity");

To use a flag, do

.. code-block:: c++

   bool compute_age = config->get_boolean("do_age");

.. note::

   - It is a good idea to avoid calling `m_config->get_double()` and
     `m_config->get_boolean()` from within loops: looking up a parameter by its name is
     slow.
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
     - scalar 3D fields (within the ice)

Please see the documentation of these classes for more info. The base class `IceModelVec` is
a virtual class, so code should use the above derived classes. Only the derived classes
have `create()` methods, in particular.

To create a scalar field, for example, one needs to create an instance of one of the
classes listed above and then call "create" to allocate it.

.. code-block:: c++

   // land ice thickness
   ice_thickness.create(grid, "thk", WITH_GHOSTS, 2);
   ice_thickness.set_attrs("model_state", "land ice thickness",
                           "m", "land_ice_thickness");
   ice_thickness.metadata().set_double("valid_min", 0.0);

Here `grid` is an `IceGrid` instance, `thk` is the name of the NetCDF variable,
`WITH_GHOSTS` means that storage for "ghost" ("halo") points will be allocated, and "2" is
the number of ghosts (in other words: needed stencil width).

The `IceModelVec` destructor takes care of undoing all that's done by the `create()` call.
Therefore you don't need to explicitly de-allocate variables unless you dynamically
created the `IceModelVec` (or derived) instance using the C++ "`new`" operator. (In which
case "`delete`" should be used.)

The `IceModelVec::set_attrs()` call sets commonly used NetCDF variable attributes seen in
PISM output files:

.. list-table::

   * - `pism_intent`
     - the only important case is "model_state", see below
   * - `long_name`
     - the (descriptive) long name used for plotting, etc (a free-form string)
   * - `units`
     - units used *in the code*. Does not have to match units in a file
   * - `standard_name`
     - CF standard name, if defined, or an empty string.

The third call above `ice_thickness.metadata()` allows accessing variable metadata and
adding arbitrary named attributes. See `VariableMetadata` for details.

The CF convention covers some attribute semantics, including `valid_min` in this example.

PISM will automatically convert units from ones present in an input file into internal
units defines by the `set_attrs()` call (see above).

If you want PISM to save data in units other than internal ones, first set these
"glaciological" units:

.. code-block:: c++

   ice_thickness.metadata().set_string("glaciological_units", "km");

Read data from a file
---------------------

There are at least three cases of "reading data from a file":

- reading a field stored in an input file on a grid matching the one used by the current
  run (restarting a run) and
- reading a field stored in an input file on a different grid, interpolating onto the
  current grid (bootstrapping).
- reading a field stored in a file **other** than the input file using interpolation
  (assuming that grids are compatible but not identical)

FIXME

Write data to a file
--------------------

To write a field stored in an `IceModelVec` to an already "prepared" file, just call

.. code-block:: c++

   precip.write(filename);

The file referred to by "filename" here has to have the time "time" dimension created,
that is, it must be prepared (other dimensions are created automatically, if needed). No
action is needed to be able to write to an output ("-o") file, a snapshot file or the
like; IceModel has already prepared it.

If you do need to "prepare" a file, do:

.. code-block:: c++

   PIO nc(grid.com, grid.config.get_string("output_format"));

   std::string time_name = config.get_string("time_dimension_name");
   nc.open(filename, PISM_WRITE); // append == false
   nc.def_time(time_name, grid.time->calendar(),
               grid.time->CF_units_string());
   nc.append_time(time_name, grid.time->current());
   nc.close();

When a file is opened with the `PISM_WRITE` mode, PISM checks if this file is present and
moves it aside if it is; to append to an existing file, use

.. code-block:: c++

   nc.open(filename, PISM_WRITE, true); // append == true

A newly-created file is "empty" and contains no records. The `nc.append_time()` call
creates a record corresponding to a particular model year.

Create IceModelVec instances that are data members of a class
-------------------------------------------------------------

To add a new variable to IceModel, allocate it in the createVecs() method.

If `pism_intent` is set to `model_state` and a variable is added to the "variables"
dictionary (see PISMVars), IceModel will automatically read this variable from a file it
is re-starting from and will always write it to an output, snapshot and backup files.

.. code-block:: c++

   variables.add(ice_thickness);

Read scalar forcing data
------------------------

PISM uses instances of the `Timeseries` class to read scalar forcing data; please see
`PA_delta_T` for an example.

The following snippet illustrates creating a `Timeseries` instance and reading data from a
file.

.. code-block:: c++

   std::string offset_name = "delta_T";

   offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
   offset->set_string("units", "Kelvin");
   offset->set_dimension_units(grid.time->units_string(), "");

   verbPrintf(2, g.com,
              "  reading %s data from forcing file %s...\n",
              offset->short_name.c_str(), filename.c_str());

   PIO nc(g.com, "netcdf3", grid.get_unit_system());
   nc.open(filename, PISM_NOWRITE);
   {
     offset->read(nc, grid.time);
   }
   nc.close();

   // use offset

   delete offset; // when done


To use `offset`, call

.. code-block:: c++

   double data = (*offset)(time);

to get the value corresponding to the time `time` in seconds. The value returned will be
computed using linear interpolation.

Read 2D forcing fields
----------------------

PISM uses instances of the `IceModelVec2T` class to read 2D forcing fields that
vary in time; please see PSDirectForcing for an example.

The following snippet from PSDirectForcing::init() illustrates creating an
IceModelVec2T object and reading data from a file.

.. code-block:: c++

   IceModelVec2T temperature;
   temperature.set_n_records((unsigned int) config.get_double("climate_forcing_buffer_size"));
   temperature.create(grid, "artm", WITHOUT_GHOSTS);
   temperature.set_attrs("climate_forcing",
                         "temperature of the ice at the ice surface but below firn processes",
                         "Kelvin", "");
   temperature.init(filename);


@section using_vars Using fields managed by IceModel in a surface model to implement a parameterization

Please see PA_SeaRISE_Greenland::init() and PA_SeaRISE_Greenland::update() for an example.

@section writing_components Managing I/O in a PISMComponent derived class

A PISM component needs to implement the following I/O methods:

- init(). It is not an I/O method per se, but most PISM components
  read their input fields there; see PA_SeaRISE_Greenland::init().
- add_vars_to_output(), which adds variable names to the list of
  fields that need to be written. See
  PSTemperatureIndex::add_vars_to_output_impl() for an example.
- define_variables(), which defines variables. (See
  PSTemperatureIndex::define_variables().)
- write_variables(), which writes data; see
  PSTemperatureIndex::write_variables().

Why are all these methods needed? In PISM we separate defining and writing NetCDF
variables because defining all the NetCDF variables before writing data is a lot faster
than defining a variable, writing it, defining the second variable, etc. (See `The NetCDF
Users' Guide <netcdf-classic-format_>`_ for a technical explanation.)

Within `IceModel` the following steps are done to write 2D and 3D fields to an
output file:

- Assemble the list of variables to be written (see
  IceModel::set_output_size()); calls add_vars_to_output()
- Create a NetCDF file
- Define all the variables in the file (see IceModel::write_variables());
  calls define_variables()
- Write all the variables to the file (same method); calls write_variables().

Add a new "diagnostic" quantity to an atmosphere model
------------------------------------------------------

To add a new "diagnostic" quantity (i.e. a 2D or 3D field that needs to be saved to an
output file but is not permanently stored in memory and is computed on demand), do the
following.

Create the class implementing the diagnostic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming that PA_foo is the atmosphere model we're working with and "bar" is
the name of the new quantity, we need to add this code

.. code-block:: c++

   class PA_foo_bar : public PISMDiag<PA_foo>
   {
   public:
     PA_foo_bar(PA_foo *m, IceGrid &g, PISMVars &my_vars);
     virtual IceModelVec::Ptr compute();
   };

to a header (`.hh`) file and implement it in a `.cc` file.

See `IceModel_diagnostics.cc` for examples. Generally speaking, in every class
implementing a "diagnostic" quantity

- the constructor sets metadata
- you have access to a data member "var" of an atmosphere model as
  "model->var"; you might need to add

.. code-block:: c++

   friend class PA_foo_bar;

to the definition of PA_foo if you need to access a private data member (see the
definition of IceModel for one example).

- when working with atmosphere models and other classes derived from PISMComponent, you
  can access the config database as "model->config".
- the `PISMDiagnostic::compute()` method allocates memory and performs the computation.
- to use a field managed by IceModel, use "variables":

.. code-block:: c++

   const IceModelVec2S *surface = variables.get_2d_scalar("surface_altitude");

- the **caller** of the `PISMDiagnostic::compute()` method has to de-allocate the field
  returned by `PISMDiagnostic::compute()`

Note that in almost every (current) implementation of `PISMDiagnostic::compute()` you see

.. code-block:: c++

   IceModelVec::Ptr ...::compute() {
     const PetscScalar fillval = -0.01;

     <...>

     // 1
     IceModelVec2S::Ptr result(new IceModelVec2S);
     result->create(grid, "hardav", WITHOUT_GHOSTS);
     result->metadata(0) = m_vars[0];

     <...>

     // 2
     return result;
   }


The block marked "1" allocates a 2D field and copies metadata stored in `m_vars[0]`, while
the block marked "2" returns a pointer to a 2D field, which gets cast to a pointer to a
"generic" field.

This allows us to have the same interface for both 2D and 3D diagnostics.

"Register" the new diagnostic.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To make the new diagnostic field available (i.e. to be able to use the new `PA_foo_bar`
class), implement `PA_foo::diagnostics_impl()` or `PA_foo::ts_diagnostics_impl()`.

.. code-block:: c++

   std::map<std::string, Diagnostic::Ptr> PA_foo_bar::diagnostics_impl() const {
     return {{"name", Diagnostic::Ptr(new PA_foo_bar(this))}};
   }

Note that if you are implementing this method in a "modifier" of a surface model, you need
to remember to call

.. code-block:: c++

   input_model->diagnostics();

.. rubric:: Footnotes

.. [#] And all classes derived from `Component`.

