// Copyright (C) 2008-2016 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __Component_hh
#define __Component_hh

#include <string>
#include <set>
#include <map>

#include "base/util/io/IO_Flags.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMUnits.hh"
#include "base/util/Logger.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMDiagnostic.hh"

namespace pism {

class MaxTimestep;
class PIO;
class IceModelVec;

enum InitializationType {INIT_RESTART, INIT_BOOTSTRAP, INIT_OTHER};

struct InputOptions {
  //! initialization type
  InitializationType type;
  //! name of the input file (if applicable)
  std::string filename;
  //! index of the record to re-start from
  unsigned int record;
};

InputOptions process_input_options(MPI_Comm com);

//! \brief A class defining a common interface for most PISM sub-models.
/*!
  \section pism_components PISM's model components and their interface

  We've found that many sub-models in PISM share some tasks: they need to be
  "initialized", "updated", asked for diagnostic quantities, asked to write the
  model state...

  Component and its derived classes were created to have a common interface
  for PISM sub-models, such as surface, atmosphere, ocean and bed deformation
  models.

  \subsection pismcomponent_init Initialization

  Component::init() should contain all the initialization code,
  excluding memory-allocation. (We might need to "re-initialize" a
  component.)

  Many PISM sub-models read data from the same file the rest of PISM reads
  from. Component::find_pism_input() checks options `-i` and `-bootstrap`
  options to simplify finding this file.

  \subsection pismcomponent_output Writing to an output file

  A PISM component needs to implement the following I/O methods:

  - define_model_state_impl()
  - write_model_state_impl()

  Why are all these methods needed? In PISM we separate defining and writing
  NetCDF variables because defining all the NetCDF variables before writing
  data is a lot faster than defining a variable, writing it, defining the
  second variable, etc. (See <a
  href="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Parts-of-a-NetCDF-Classic-File.html#Parts-of-a-NetCDF-Classic-File">The
  NetCDF Users' Guide</a> for a technical explanation.)

  Within IceModel the following steps are done to write 2D and 3D fields to an
  output file:

  - Assemble the list of variables to be written (see
  IceModel::set_output_size()); calls add_vars_to_output()
  - Create a NetCDF file
  - Define all the variables in the file (see IceModel::write_variables());
  calls define_variables()
  - Write all the variables to the file (same method); calls write_variables().

  \subsection pismcomponent_timestep Restricting time-steps

  Implement Component_TS::max_timestep() to affect PISM's adaptive time-stepping mechanism.
*/
class Component {
public:

  /** Create a Component instance given a grid. */
  Component(IceGrid::ConstPtr g);
  virtual ~Component();

  //! Add pointers to available diagnostic quantities to a dictionary.
  std::map<std::string, Diagnostic::Ptr> diagnostics() const;
  std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics() const;

  IceGrid::ConstPtr grid() const;

  void define_model_state(const PIO &output) const;
  void write_model_state(const PIO &output) const;

protected:
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;
  virtual std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;

  /** @brief This flag determines whether a variable is read from the
      `-regrid_file` file even if it is not listed among variables in
      `-regrid_vars`.
  */
  enum RegriddingFlag { REGRID_WITHOUT_REGRID_VARS, NO_REGRID_WITHOUT_REGRID_VARS };
  virtual void regrid(const std::string &module_name, IceModelVec &variable,
                      RegriddingFlag flag = NO_REGRID_WITHOUT_REGRID_VARS);
protected:
  //! grid used by this component
  const IceGrid::ConstPtr m_grid;
  //! configuration database used by this component
  const Config::ConstPtr m_config;
  //! unit system used by this component
  const units::System::Ptr m_sys;
  //! logger (for easy access)
  const Logger::ConstPtr m_log;
};

//! \brief An abstract class for time-stepping PISM components. Created to
//! simplify creating basic surface, snow, atmosphere, ocean... models for
//! PISM.
class Component_TS : public Component {
public:
  /** Create an instance of Component_TS given a grid. */
  Component_TS(IceGrid::ConstPtr g);
  virtual ~Component_TS();

  //! @brief Reports the maximum time-step the model can take at t.
  MaxTimestep max_timestep(double t) const;

  //! Update the *state* of a component, if necessary.
  /**
   * Defines the common interface of time-stepping components.
   *
   * Derived classes should use this method to perform computations
   * needed to step from `t` to `t + dt`. This could be
   * a no-op.
   *
   * Time-step length `dt` should never be zero.
   *
   * This method will be called only once with a given `t`, `dt` pair;
   * the value of `t` in a particular call should be equal to `t + dt`
   * in the previous call (i.e. there should be no "gaps" or "repeats").
   *
   * One unfortunate exception is the initialization stage:
   * IceModel::bootstrapFromFile() and IceModel::model_state_setup()
   * may need to "know" the state of the model at the beginning of
   * the run, which might require a non-trivial computation and so
   * requires an update() call. Because of this *currently* update()
   * may get called twice at the beginning of the run.
   *
   * Other interface methods
   * (SurfaceModel::ice_surface_temperature() is an example)
   * should use cached values if the corresponding computation is
   * expensive. Methods like
   * SurfaceModel::ice_surface_temperature() might be called
   * multiple times per time-step.
   *
   * TemperatureIndex is an example of a component that does a
   * fairly expensive computation in TemperatureIndex::update() and
   * uses cached values in
   * TemperatureIndex::ice_surface_mass_flux_impl().
   *
   * *Who* calls this depends on the kind of the component in
   * question, but all calls originate from IceModel::step() and the
   * initialization methods mentioned above.
   *
   * @param[in] t time corresponding to the beginning of the time-step, in seconds
   * @param[in] dt length of the time-step, in seconds
   *
   */
  void update(double t, double dt);

protected:
  virtual MaxTimestep max_timestep_impl(double t) const = 0;
  virtual void update_impl(double t, double dt) = 0;
protected:
  //! Last time used as an argument for the update() method.
  double m_t;
  //! Last time-step used as an argument for the update() method.
  double m_dt;
};

void init_step(Component_TS &model, const Time& time);

} // end of namespace pism

#endif // __Component_hh
