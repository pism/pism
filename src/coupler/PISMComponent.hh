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
#include "../base/NCVariable.hh"
#include "../base/PISMVars.hh"
#include "../base/grid.hh"
#include "../base/LocalInterpCtx.hh"
#include <gsl/gsl_math.h>

//! \brief An abstract class intended to simplify creating basic surface,
//! snow, atmosphere, ocean... models for PISM.
class PISMComponent {
public:
  PISMComponent(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PISMComponent() {};

  virtual PetscErrorCode init(PISMVars &vars) = 0;

  //! \brief Reports the maximum time-step the model can take at t_years. Sets
  //! dt_years to -1 if any time-step is OK.
  virtual PetscErrorCode max_timestep(PetscReal /*t_years*/, PetscReal &dt_years)
  { dt_years = -1; return 0; }

  //! \brief Writes fields that were read from an input file and are necessary
  //! for restarting.
  virtual PetscErrorCode write_input_fields(PetscReal /*t_years*/, PetscReal /*dt_years*/,
					    string /*filename*/)
  { return 0; }

  //! \brief Updates the model and writes all the internal fields (for testing
  //! and debugging).
  virtual PetscErrorCode write_diagnostic_fields(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						 string /*filename*/)
  { return 0; }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_fields(set<string> /*vars*/, PetscReal /*t_years*/,
				      PetscReal /*dt_years*/, string /*filename*/)
  { return 0; }

  //! Update a model, if necessary.
  virtual PetscErrorCode update(PetscReal /*t_years*/, PetscReal /*dt_years*/)
  { return 0; }

protected:
  virtual PetscErrorCode find_pism_input(string &filename, LocalInterpCtx* &lic,
					 bool &regrid, int &start);
  IceGrid &grid;
  const NCConfigVariable &config;
  PetscReal t,			//!< Last time used as an argument for the update() method.
    dt;				//!< Lasr time-step used as an argument for the update() method.
};

#endif // __PISMComponent_hh
