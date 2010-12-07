// Copyright (C) 2008-2010 Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include <gsl/gsl_math.h>
#include "NCVariable.hh"
#include "PISMVars.hh"
#include "grid.hh"
#include "LocalInterpCtx.hh"
#include "PISMDiagnostic.hh"

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
  virtual void add_vars_to_output(string /*keyword*/, set<string> &/*result*/) {}

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode define_variables(set<string> /*vars*/, const NCTool &/*nc*/)
  { return 0; }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_variables(set<string> /*vars*/, string /*filename*/)
  { return 0; }

  //! Add pointers to available diagnostic quantities to a dictionary.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &/*dict*/) {}

  // //! Add pointers to scalar diagnostic quantities to a dictionary.
  // virtual void get_scalar_diagnostics(map<string, PISMDiagnostic_Scalar*> &/*dict*/) {}
protected:
  virtual PetscErrorCode find_pism_input(string &filename, LocalInterpCtx* &lic,
					 bool &regrid, int &start);
  IceGrid &grid;
  const NCConfigVariable &config;
};

//! \brief An abstract class for "diagnostic" components (such as stress
//! balance modules).
class PISMComponent_Diag : public PISMComponent
{
public:
  PISMComponent_Diag(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent(g, conf) {}
  virtual ~PISMComponent_Diag() {}

  virtual PetscErrorCode update(bool /*fast*/)
  { return 0; }

  virtual PetscErrorCode write_model_state(string /*filename*/)
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

  //! \brief Reports the maximum time-step the model can take at t_years. Sets
  //! dt_years to -1 if any time-step is OK.
  virtual PetscErrorCode max_timestep(PetscReal /*t_years*/, PetscReal &dt_years)
  { dt_years = -1; return 0; }

  //! Update a model, if necessary.
  virtual PetscErrorCode update(PetscReal /*t_years*/, PetscReal /*dt_years*/)
  { return 0; }

  //! \brief Writes fields that were read from an input file and are necessary
  //! for restarting.
  virtual PetscErrorCode write_model_state(PetscReal /*t_years*/, PetscReal /*dt_years*/,
					   string /*filename*/)
  { return 0; }

protected:
  PetscReal t,			//!< Last time used as an argument for the update() method.
    dt;				//!< Last time-step used as an argument for the update() method.
};

#endif // __PISMComponent_hh
