# How do I...? {#howto}

[TOC]

@section adding_parameters Creating and using configuration flags and parameters

- Edit 'src/pism_config.cdl'. Each flag or parameter is stored as a NetCDF attribute and should have a corresponding "_doc" attribute describing its meaning.

~~~
    pism_config:standard_gravity = 9.81;
    pism_config:standard_gravity_doc = "m s-2; acceleration due to gravity on Earth geoid";
~~~

- One can access these parameters using the NCConfigVariable class. IceModel
has an instance of this class as a data member 'config', so no additional code is
necessary to initialize the configuration database.

To use a parameter, do

~~~
    double g = config.get("standard_gravity");
~~~

To use a flag, do

~~~
    bool compute_age = config.get_flag("do_age");
~~~

@note
- It is a good idea to avoid calling `config.get()` and `config.get_flag()`
from within loops: looking up a parameter by its name is slow.
- Please see @ref config "this page" for a list of flags and parameters
currently used in PISM.

@section adding_variables Creating and using additional variables.

@subsection standalone_code Creating IceModelVec instances.

PISM uses the following classes to manage 2D and 3D fields, their I/O and
metadata:

- IceModelVec2S -- scalar 2D fields
- IceModelVec2V -- vector 2D fields (such as horizontal velocities;
corresponds to 2 NetCDF variables)
- IceModelVec2Int -- 2D masks, such as the grounded/floating mask
- IceModelVec2T -- 2D time-dependent fields (used to read and store forcing data)
- IceModelVec3 -- scalar 3D fields (within the ice)

Please see the documentation of these classes for more info.  The base class
IceModelVec is a virtual class, so code should use the above derived classes.
Only the derived classes have create() methods, in particular.

To create a scalar field, for example, one needs to create an instance of one
of the classes listed above and then call "create" to allocate it.

~~~
    // land ice thickness
    ierr = vH.create(grid, "thk", true, 2); CHKERRQ(ierr);
    ierr = vH.set_attrs("model_state", "land ice thickness",
                        "m", "land_ice_thickness"); CHKERRQ(ierr);
    ierr = vH.set_attr("valid_min", 0.0); CHKERRQ(ierr);
~~~

Here `grid` is an IceGrid instance, "thk" is the name of the NetCDF variable,
"true" corresponds to "has ghosts" and "2" is the number of ghosts (in other
words: needed stencil width).

The IceModelVec destructor takes care of undoing all that's done by the
create() call.  Therefore you don't need to explicitly de-allocate variables
unless you dynamically created the IceModelVec (or derived) instance using the
C++ "new" operator.  (In which case "delete" should be used.)

The set_attrs() call sets variable metadata you will see in an output file
- **pism_intent**  -- the only important case is "model_state", see below
- **long_name** -- the (descriptive) long name used for plotting, etc
- **units** -- units used *in* *the* *code*. Does not have to match units
in a file
- **standard_name** -- CF standard name.

The third call above using set_attr() allows any named attribute to be added to
the variable.  The CF convention covers some attribute semantics, including
"valid_min" in this example.

If you want PISM to automatically convert data from the units used in an input
file to the units used internally, you need to add

~~~
    ierr = vH.set_glaciological_units("km"); CHKERRQ(ierr);
~~~

If, in addition to this, you want PISM to *write* data using "glaciological"
units, add
~~~
vH.write_in_glaciological_units = true;
~~~

@subsubsection reading_more Reading data from a file.

There are at least three cases of "reading data from a file":
- reading a field stored in an input file on a grid matching the one used by the current run
(restarting a run) and
- reading a field stored in an input file on a different grid, interpolating onto the current
grid (bootstrapping).
- reading a field stored in a file **other** than the input file using
interpolation (assuming that grids are compatible but not identical)

This snippet from PAYearlyCycle::init() covers the first two cases (at least for
surface, atmosphere and ocean models; %i.e. for all the classes derived from PISMComponent)

~~~
    ierr = find_pism_input(precip_filename, regrid, start); CHKERRQ(ierr);

    // read precipitation rate from file
    ierr = verbPrintf(2, grid.com, 
                      "    reading mean annual ice-equivalent precipitation rate 'precip'\n"
                      "      from %s ... \n",
                      precip_filename.c_str()); CHKERRQ(ierr); 
    if (regrid) {
      ierr = precip.regrid(precip_filename.c_str(), true); CHKERRQ(ierr); // fails if not found!
    } else {
      ierr = precip.read(precip_filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
    }
    string precip_history = "read from " + precip_filename + "\n";

    ierr = precip.set_attr("history", precip_history); CHKERRQ(ierr);

    delete lic;
~~~

Here
- "start" is an index of a record *within* *a* *file* that we need to
  read. PISM almost always reads the last record.

Please see PISMComponent::find_pism_input() to see how to compute the "start"
index "by hand".

The snippet below is an example of case 3 (for 2D fields, reading the last record).

~~~
    ierr = variable.regrid(filename, true); CHKERRQ(ierr); // fails if not found!
~~~

@subsubsection writing_more Writing data to a file.

To write a field stored in an IceModelVec to an already "prepared" file, just call

~~~
    ierr = precip.write(filename); CHKERRQ(ierr);
~~~

The file referred to by "filename" here has to have "t", "x", "y", "z", "zb"
dimensions created, that is, it must be prepared.  No action is needed to be
able to write to an output ("-o") file, a snapshot file or the like; IceModel
has already prepared it.

If you do need to "prepare" a file, do:

~~~
    PISMIO nc(&grid);

    ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
    ierr = nc.append_time(grid.year); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
~~~

Regarding the second and third arguments to nc.open_for_writing():
- append == false means that if a file exists already, it will be moved
aside
- check_dims == true means tells PISM to create dimensions. 

A newly-created file is "empty" and contains no records. The nc.append_time()
call creates a record corresponding to a particular model year.

@subsection inside_icemodel Creating IceModelVec instances that are data members of IceModel or a derived class.

To add a new variable to IceModel, allocate it in the createVecs() method.

If "pism_intent" is set to "model_state" and a variable is added to the
"variables" dictionary (see PISMVars), IceModel will automatically read this
variable from a file it is re-starting from and will always write it to an
output, snapshot and backup files.

~~~
    ierr = variables.add(vH); CHKERRQ(ierr);
~~~

@section reading_scalar_forcing Reading scalar forcing data

PISM uses instances of the Timeseries class to read scalar forcing data; please
see PAForcing or PA_SeaRISE_Greenland for an example.

The following snippet from PAForcing::init() illustrates creating a Timeseries
object and reading data from a file.

~~~
    delta_T = new Timeseries(grid.com, grid.rank, "delta_T", "t");
    ierr = delta_T->set_units("Kelvin", ""); CHKERRQ(ierr);
    ierr = delta_T->set_dimension_units("years", ""); CHKERRQ(ierr);
    ierr = delta_T->set_attr("long_name", "near-surface air temperature offsets");
    CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
		      "  reading delta T data from forcing file %s...\n", dT_file);
		      CHKERRQ(ierr);
	 
    ierr = delta_T->read(dT_file); CHKERRQ(ierr);
~~~
Call

~~~
    double offset = (*delta_T)(time);
~~~
to get the value corresponding to the time "time", in this case in years. The
value returned will be computed using linear interpolation.

@section reading_2d_forcing Reading 2D forcing fields

PISM uses instances of the IceModelVec2T class to read 2D forcing fields that
vary in time; please see PSDirectForcing for an example. See the section @ref
adding_variables if you need to read and store **one** record of a 2D field.

The following snippet from PSDirectForcing::init() illustrates creating an
IceModelVec2T object and reading data from a file.

~~~
    temperature.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
    ierr = temperature.create(grid, "artm", false); CHKERRQ(ierr);
    ierr = temperature.set_attrs("climate_forcing",
                                 "temperature of the ice at the ice surface but below firn processes",
                                 "Kelvin", ""); CHKERRQ(ierr);
    ierr = temperature.init(filename); CHKERRQ(ierr);
~~~

@section get_variables Using fields managed by IceModel in a surface model to implement a parameterization.

Please see PA_EISMINT_Greenland::init() and PA_EISMINT_Greenland::update() for an example.

@section component_output Managing I/O in a PISMComponent derived class

A PISM component needs to implement the following I/O methods:

- init(). It is not an I/O method per se, but most PISM components read their
input fields there; see PA_EISMINT_Greenland::init().
- add_vars_to_output(), which adds variable names to the list of fields that need
to be written. See PSTemperatureIndex::add_vars_to_output() for an example.
- define_variables(), which defines variables. (See PSTemperatureIndex::define_variables().)
- write_variables(), which writes data; see PSTemperatureIndex::write_variables().

Why are all these methods needed? In PISM we separate defining and writing NetCDF variables because defining all
the NetCDF variables before writing data is a lot faster than defining a
variable, writing it, defining the second variable, etc. (See <a
href="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Parts-of-a-NetCDF-Classic-File.html#Parts-of-a-NetCDF-Classic-Filel">The
NetCDF Users' Guide</a> for a technical explanation.)

Within IceModel the following steps are done to write 2D and 3D fields to an
output file:

- Assemble the list of variables to be written (see
IceModel::set_output_size()); calls add_vars_to_output()
- Create a NetCDF file
- Define all the variables in the file (see IceModel::write_variables());
calls define_variables()
- Write all the variables to the file (same method); calls write_variables().

@section diag Adding a new "diagnostic" quantity to an atmosphere model.

To add a new "diagnostic" quantity (i.e. a 2D or 3D field that needs to be
saved to an output file but is not permanently stored in memory and is computed
on demand), do the following.

@subsection create_diag_class Create the class implementing the diagnostic

Assuming that PA_foo is the atmosphere model we're working with and "bar" is
the name of the new quantity, we need to add this code

~~~
    class PA_foo_bar : public PISMDiag<PA_foo>
    {
    public:
      PA_foo_bar(PA_foo *m, IceGrid &g, PISMVars &my_vars);
      virtual PetscErrorCode compute(IceModelVec* &result);
    };
~~~
to a header (.hh) file and implement it in a .cc file.

See IceModel_diagnostics.cc for examples. Generally speaking, in every class
implementing a "diagnostic" quantity
- the constructor sets metadata
- you have access to a data member "var" of an atmosphere model as
"model->var"; you might need to add

~~~
    friend class PA_foo_bar;
~~~

to the definition of PA_foo if you need to access a private data member (see
the definition of IceModel for one example).
- when working with atmosphere models and other classes derived from
PISMComponent, you can access the config database as "model->config".
- the PISMDiagnostic::compute() method allocates memory and performs the computation.
- to use a field managed by IceModel, use "variables":

~~~
    IceModelVec2S *surface;
    surface = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
    if (surface == NULL) SETERRQ(1, "surface_altitude is not available");
~~~

- the **caller** of the PISMDiagnostic::compute() method has to de-allocate the field
allocated by PISMDiagnostic::compute()

Note that in almost every (current) implementation of PISMDiagnostic::compute() you see

~~~
    PetscErrorCode ...::compute(IceModelVec* &output) {
      PetscErrorCode ierr;
      const PetscScalar fillval = -0.01;

      <...>

      // 1
      IceModelVec2S *result = new IceModelVec2S;
      ierr = result->create(grid, "hardav", false); CHKERRQ(ierr);
      ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

      <...>

      // 2
      output = result;
      return 0;
    }
~~~

The block marked "1" allocates a 2D field and copies metadata stored in
vars[0], while the block marked "2" casts a pointer to a 2D field to a pointer
to a "generic" field.

This allows us to have the same interface for both 2D and 3D diagnostics.

@subsection add_diag_to_the_list "Register" the new diagnostic.

To make the new diagnostic field available (i.e. to be able to use the new
PA_foo_bar class), implement PA_foo::get_diagnostics().

~~~
    void PA_foo::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
      dict["bar"] = new PA_foo_bar(this, grid, *variables);
    }
~~~

Note that if you are implementing this method in a "modifier" of a surface
model, you need to remember to call
~~~
    input_model->get_diagnostics(dict);
~~~
@note PISM's handling of scalar diagnostic quantities is different and should
be improved.
