// Copyright (C) 2008-2014 Ed Bueler and Constantine Khroulev
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

#include <petscsys.h>
#include <gsl/gsl_math.h>
#include <string>
#include <set>
#include <map>

#include "PIO.hh"

namespace pism {

class IceGrid;
class Config;
class NCSpatialVariable;
class Diagnostic;
class TSDiagnostic;
class Vars;
class IceModelVec;

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
  from. Component::find_pism_input() checks -i and -boot_file command-line
  options and simplifies finding this file.

  \subsection pismcomponent_output Writing to an output file

  A PISM component needs to implement the following I/O methods:

  - add_vars_to_output(), which adds variable names to the list of fields that need
  to be written.
  - define_variables(), which defines variables to be written and writes variable metadata.
  - write_variables(), which writes data itself.
  
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
  Component(const IceGrid &g)
    : m_grid(g), m_config(g.config) {}
  virtual ~Component() {}

  //! \brief Adds more variable names to result (to let sub-models respect
  //! -o_size or -save_size).
  /*!
    Keyword can be one of "small", "medium" or "big".
  */
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result) = 0;

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                IO_Type nctype) = 0;

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void write_variables(const std::set<std::string> &vars, const PIO& nc) = 0;

  //! Add pointers to available diagnostic quantities to a dictionary.
  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict)
  {
    (void)dict;
    (void)ts_dict;
  }

  const IceGrid& get_grid() const;

protected:
  virtual void find_pism_input(std::string &filename, bool &regrid, int &start);
  const IceGrid &m_grid;
  const Config &m_config;

  /** @brief This flag determines whether a variable is read from the
      `-regrid_file` file even if it is not listed among variables in
      `-regrid_vars`.
  */
  enum RegriddingFlag { REGRID_WITHOUT_REGRID_VARS, NO_REGRID_WITHOUT_REGRID_VARS };
  virtual void regrid(const std::string &module_name, IceModelVec *variable,
                                RegriddingFlag flag = NO_REGRID_WITHOUT_REGRID_VARS);
};

//! \brief An abstract class for time-stepping PISM components. Created to
//! simplify creating basic surface, snow, atmosphere, ocean... models for
//! PISM.
class Component_TS : public Component
{
public:
  /** Create an instance of Component_TS given a grid. */
  Component_TS(const IceGrid &g)
    : Component(g)
  { m_t = m_dt = GSL_NAN; }
  virtual ~Component_TS() {}

  //! \brief Reports the maximum time-step the model can take at t. Sets
  //! dt to -1 if any time-step is OK.
  virtual void max_timestep(double t, double &dt, bool &restrict)
  {
    (void)t;
    dt = -1;
    restrict = false;
  }

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
   * PSTemperatureIndex is an example of a component that does a
   * fairly expensive computation in PSTemperatureIndex::update() and
   * uses cached values in
   * PSTemperatureIndex::ice_surface_mass_flux().
   *
   * *Who* calls this depends on the kind of the component in
   * question, but all calls originate from IceModel::step() and the
   * initialization methods mentioned above.
   *
   * @param[in] t time corresponding to the beginning of the time-step, in seconds
   * @param[in] dt length of the time-step, in seconds
   *
   * @return 0 on success
   */
  virtual void update(double t, double dt) = 0;

protected:
  double m_t,                   //!< Last time used as an argument for the update() method.
    m_dt;                               //!< Last time-step used as an argument for the update() method.
};

//! \brief This template allows creating Component_TS (AtmosphereModel,
//! SurfaceModel and OceanModel) modifiers with minimum effort.
/*!
 * A specialization of this template will implement all important methods
 * except init(). This means that to create a complete modifier, one needs to
 * re-implement interesting methods, without worrying about preserving
 * modifier's "transparency".
 */
template<class Model>
class Modifier : public Model
{
public:
  Modifier(IceGrid &g, Model* in)
    : Model(g), input_model(in) {}
  virtual ~Modifier()
  {
    if (input_model != NULL) {
      delete input_model;
    }
  }

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result)
  {
    if (input_model != NULL) {
      input_model->add_vars_to_output(keyword, result);
    }
  }

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype)
  {
    if (input_model != NULL) {
      input_model->define_variables(vars, nc, nctype);
    }
  }

  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc)
  {
    if (input_model != NULL) {
      input_model->write_variables(vars, nc);
    }
  }

  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict)
  {
    if (input_model != NULL) {
      input_model->get_diagnostics(dict, ts_dict);
    }
  }

  virtual void max_timestep(double my_t, double &my_dt, bool &restrict)
  {
    if (input_model != NULL) {
      input_model->max_timestep(my_t, my_dt, restrict);
    } else {
      my_dt    = -1;
      restrict = false;
    }
  }

  virtual void update(double my_t, double my_dt)
  {
    Model::m_t = my_t;
    Model::m_dt = my_dt;
    if (input_model != NULL) {
      input_model->update(my_t, my_dt);
    }
  }

protected:
  Model *input_model;
};

} // end of namespace pism

#endif // __Component_hh
