// Copyright (C) 2010, 2011 Constantine Khroulev
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

#include "PISMSurface.hh"

const int TAG_EBM_RUN     = 1;
const int TAG_EBM_STOP    = 2;
const int TAG_EBM_COMMAND = 3;
const int TAG_EBM_STATUS  = 4;

const int EBM_STATUS_READY  = 1;
const int EBM_STATUS_FAILED = 2;

//! \brief A class executing an external energy balance model.
class EBM_driver {
public:
  EBM_driver(MPI_Comm com);
  int run();
protected:
  MPI_Comm inter_comm;          // communicator used to send messages to PISM
  string command;               // EBM command
  int run_ebm();
};

//! \brief A derived class created to couple PISM to an energy balance model (through
//! files).
/*!
  Uses an atmospheric lapse rate to correct temperatures at the ice surface.
 */ 
class PSExternal : public PISMSurfaceModel {
public:
  PSExternal(IceGrid &g, const NCConfigVariable &conf, MPI_Comm my_inter_comm)
    : PISMSurfaceModel(g, conf)
  {
    gamma = 0;                  // essentially disables the lapse rate correction
    update_interval = 1;        // years
    last_update = GSL_NAN;

    inter_comm = my_inter_comm;
  }

  virtual ~PSExternal();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &/*dict*/) {}

  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2S &result);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years);
protected:
  double gamma, update_interval, last_update;
  IceModelVec2S acab, artm, artm_0, *usurf, *topg;
  MPI_Comm inter_comm;
  string ebm_command, ebm_input, ebm_output;

  virtual PetscErrorCode update_artm();
  virtual PetscErrorCode update_acab();
  virtual PetscErrorCode run();
  virtual PetscErrorCode wait();
  virtual PetscErrorCode write_coupling_fields();
};

