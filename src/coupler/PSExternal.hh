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

#ifndef _PSEXTERNAL_H_
#define _PSEXTERNAL_H_

#include "PISMSurface.hh"

const int TAG_EBM_RUN     = 1;
const int TAG_EBM_STOP    = 2;
const int TAG_EBM_COMMAND = 3;
const int TAG_EBM_STATUS  = 4;

const int EBM_STATUS_READY  = 1;
const int EBM_STATUS_FAILED = 2;

//! \brief A class running an external energy balance model.
class EBM_driver {
public:
  EBM_driver(MPI_Comm com);
  int run();
protected:
  MPI_Comm inter_comm;          // communicator used to send messages to PISM
  string command;               // EBM command
  int run_ebm(double year);
};

//! \brief A derived class created to couple PISM to an energy balance model (through
//! files).
class PSExternal : public PISMSurfaceModel {
public:
  PSExternal(IceGrid &g, const NCConfigVariable &conf, MPI_Comm my_inter_comm)
    : PISMComponent_TS(g, conf), PISMSurfaceModel(g, conf)
  {
    update_interval = 1;        // years
    ebm_update_interval = 0.5 * update_interval;
    last_ebm_update_year = GSL_NAN;
    last_bc_update_year = GSL_NAN;
    ebm_is_running = false;
    inter_comm = my_inter_comm;
  }

  virtual ~PSExternal();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc,
                                          nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &/*dict*/) {}

  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years);
protected:
  double gamma, update_interval, ebm_update_interval, last_ebm_update_year, last_bc_update_year;
  IceModelVec2S acab, artm;
  MPI_Comm inter_comm;
  string ebm_command, ebm_input, ebm_output;
  bool ebm_is_running;
  vector<IceModelVec*> ebm_vars;

  virtual PetscErrorCode update_artm();
  virtual PetscErrorCode update_acab();
  virtual PetscErrorCode run(double t_years);
  virtual PetscErrorCode wait();
  virtual PetscErrorCode write_coupling_fields();
};

//! \brief A less-generic coupling class implementing a lapse-rate correction
//! of the annual temperature at the top of the ice.
/*!
 * ALR stands for "Atmospheric Lapse Rate"; this class was written at the
 * request of Nick Golledge of the Antarctic Research Centre, VUW, New Zealand.
 */
class PSExternal_ALR : public PSExternal
{
public:
  PSExternal_ALR(IceGrid &g, const NCConfigVariable &conf, MPI_Comm my_inter_comm)
    : PISMComponent_TS(g, conf), PSExternal(g, conf, my_inter_comm)
  {
    gamma = 0;                  // essentially disables the lapse rate correction
  }

  virtual ~PSExternal_ALR() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual void add_vars_to_output(string keyword, set<string> &result);
protected:
  IceModelVec2S artm_0, *usurf;
  PetscReal gamma;

  virtual PetscErrorCode update_artm();
};

#endif /* _PSEXTERNAL_H_ */
