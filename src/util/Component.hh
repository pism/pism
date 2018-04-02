// Copyright (C) 2008-2018 Ed Bueler and Constantine Khroulev
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

#include "pism/util/io/IO_Flags.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Units.hh"
#include "pism/util/Logger.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Diagnostic.hh"

namespace pism {

class MaxTimestep;
class PIO;
class IceModelVec;

enum InitializationType {INIT_RESTART, INIT_BOOTSTRAP, INIT_OTHER};

struct InputOptions {
  InputOptions(InitializationType t, const std::string &file, unsigned int index);
  //! initialization type
  InitializationType type;
  //! name of the input file (if applicable)
  std::string filename;
  //! index of the record to re-start from
  unsigned int record;
};

InputOptions process_input_options(MPI_Comm com, Config::ConstPtr config);

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
  IceModel::output_variables()); calls add_vars_to_output()
  - Create a NetCDF file
  - Define all the variables in the file (see IceModel::write_variables());
  calls define_variables()
  - Write all the variables to the file (same method); calls write_variables().

  \subsection pismcomponent_timestep Restricting time-steps

  Implement Component::max_timestep() to affect PISM's adaptive time-stepping mechanism.
*/
class Component {
public:

  /** Create a Component instance given a grid. */
  Component(IceGrid::ConstPtr g);
  virtual ~Component();

  DiagnosticList diagnostics() const;
  TSDiagnosticList ts_diagnostics() const;

  IceGrid::ConstPtr grid() const;

  void define_model_state(const PIO &output) const;
  void write_model_state(const PIO &output) const;

  //! Reports the maximum time-step the model can take at time t.
  MaxTimestep max_timestep(double t) const;

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

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

} // end of namespace pism

#endif // __Component_hh
