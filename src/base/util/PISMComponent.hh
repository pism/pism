// Copyright (C) 2008-2012 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef __PISMComponent_hh
#define __PISMComponent_hh

#include <petscsys.h>
#include <gsl/gsl_math.h>
#include <string>
#include <set>
#include <map>

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

#include "PIO.hh"

class IceGrid;
class NCConfigVariable;
class NCSpatialVariable;
class PISMDiagnostic;
class PISMVars;

//! \brief A class defining a common interface for most PISM sub-models.
/*!
  \section pism_components PISM's model components and their interface

  We've found that many sub-models in PISM share some tasks: they need to be
  "initialized", "updated", asked for diagnostic quantities, asked to write the
  model state...

  PISMComponent and its derived classes were created to have a common interface
  for PISM sub-models, such as surface, atmosphere, ocean and bed deformation
  models.

  There are two kinds of PISM's components:

  \li diagnostic components (PISMComponent_Diag) and
  \li time-stepping components (PISMComponent_TS).

  The main difference is that diagnostic components do not need to know the
  model time to perform an update, while time-stepping ones need to know the
  time-step to update for (usually given as the my_t, my_dt pair defining
  the (my_t, my_t + my_dt) interval) and may impose restrictions on a
  time-step that is possible at a given time during a run.

  \subsection pismcomponent_init Initialization

  PISMComponent::init() should contain all the initialization code, preferably
  including memory-allocation. This makes it possible to separate the PISM
  initialization stage at which the choice of a sub-model is made from its own
  initialization, so we can avoid initializing a sub-model just to throw it away.

  Many PISM sub-models read data from the same file the rest of PISM reads
  from. PISMComponent::find_pism_input() checks -i and -boot_file command-line
  options and simplifies finding this file.

  \subsection pismcomponent_output Writing to an output file

  A PISM component needs to implement the following I/O methods:

  \li add_vars_to_output(), which adds variable names to the list of fields that need
  to be written.
  \li define_variables(), which defines variables to be written and writes variable metadata.
  \li write_variables(), which writes data itself.
  
  Why are all these methods needed? In PISM we separate defining and writing
  NetCDF variables because defining all the NetCDF variables before writing
  data is a lot faster than defining a variable, writing it, defining the
  second variable, etc. (See <a
  href="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Parts-of-a-NetCDF-Classic-File.html#Parts-of-a-NetCDF-Classic-Filel">The
  NetCDF Users' Guide</a> for a technical explanation.)

  Within IceModel the following steps are done to write 2D and 3D fields to an
  output file:

  \li Assemble the list of variables to be written (see
  IceModel::set_output_size()); calls add_vars_to_output()
  \li Create a NetCDF file
  \li Define all the variables in the file (see IceModel::write_variables());
  calls define_variables()
  \li Write all the variables to the file (same method); calls write_variables().

  \subsection pismcomponent_timestep Restricting time-steps

  Implement PISMComponent_TS::max_timestep() to affect PISM's adaptive time-stepping mechanism.
 */
class PISMComponent {
public:
  PISMComponent(IceGrid &g, const NCConfigVariable &conf)
    : grid(g), config(conf) {}
  virtual ~PISMComponent() {}

  virtual PetscErrorCode init(PISMVars &vars) = 0;

  //! \brief Adds more variable names to result (to let sub-models respect
  //! -o_size or -save_size).
  /*!
    Keyword can be one of "small", "medium" or "big".
   */
  virtual void add_vars_to_output(string /*keyword*/,
                                  map<string,NCSpatialVariable> &/*result*/) = 0;

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode define_variables(set<string> /*vars*/, const PIO &/*nc*/,
                                          PISM_IO_Type /*nctype*/) = 0;

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_variables(set<string> /*vars*/, string /*filename*/) = 0;

  //! Add pointers to available diagnostic quantities to a dictionary.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &/*dict*/) {}

  // //! Add pointers to scalar diagnostic quantities to a dictionary.
  // virtual void get_scalar_diagnostics(map<string, PISMDiagnostic_Scalar*> &/*dict*/) {}
protected:
  virtual PetscErrorCode find_pism_input(string &filename, bool &regrid, int &start);
  IceGrid &grid;
  const NCConfigVariable &config;
};

//! \brief An abstract class for "diagnostic" components (such as stress
//! balance modules).
/*!
 * Here "diagnostic" means "one that performs a computation which does not
 * involve time-stepping".
 */
class PISMComponent_Diag : public PISMComponent
{
public:
  PISMComponent_Diag(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent(g, conf) {}
  virtual ~PISMComponent_Diag() {}

  virtual PetscErrorCode update(bool /*fast*/)
  { return 0; }
};

//! \brief An abstract class for time-stepping PISM components. Created to
//! simplify creating basic surface, snow, atmosphere, ocean... models for
//! PISM.
class PISMComponent_TS : public PISMComponent
{
public:
  PISMComponent_TS(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent(g, conf)
  { t = dt = GSL_NAN; }
  virtual ~PISMComponent_TS() {}

  //! \brief Reports the maximum time-step the model can take at my_t. Sets
  //! my_dt to -1 if any time-step is OK.
  virtual PetscErrorCode max_timestep(PetscReal /*my_t*/, PetscReal &my_dt, bool &restrict)
  { my_dt = -1; restrict = false; return 0; }

  //! Update a model, if necessary.
  virtual PetscErrorCode update(PetscReal /*my_t*/, PetscReal /*my_dt*/) = 0;

protected:
  PetscReal t,			//!< Last time used as an argument for the update() method.
    dt;				//!< Last time-step used as an argument for the update() method.
};

//! \brief This template allows creating PISMComponent_TS (PISMAtmosphereModel,
//! PISMSurfaceModel and PISMOceanModel) modifiers with minimum effort.
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
  Modifier(IceGrid &g, const NCConfigVariable &conf, Model* in)
    : Model(g, conf), input_model(in) {}
  virtual ~Modifier()
  {
    if (input_model != NULL) {
      delete input_model;
    }
  }

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result)
  {
    if (input_model != NULL) {
      input_model->add_vars_to_output(keyword, result);
    }
  }

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype)
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode write_variables(set<string> vars, string filename)
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->write_variables(vars, filename); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict)
  {
    if (input_model != NULL) {
      input_model->get_diagnostics(dict);
    }
  }

  virtual PetscErrorCode max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict)
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->max_timestep(my_t, my_dt, restrict); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  {
    Model::t = my_t;
    Model::dt = my_dt;
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->update(my_t, my_dt); CHKERRQ(ierr);
    }
    return 0;
  }

protected:
  Model *input_model;
};

#endif // __PISMComponent_hh
