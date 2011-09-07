// Copyright (C) 2011 Constantine Khroulev
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

// Implementation of surface and atmosphere models reading spatially-variable
// time-dependent B.C. data from a file (-surface given and -atmosphere given).

#ifndef _PASDIRECTFORCING_H_
#define _PASDIRECTFORCING_H_

#include "PISMSurface.hh"
#include "PISMAtmosphere.hh"
#include "iceModelVec2T.hh"

template <class Model>
class PDirectForcing : public Model
{
public:
  PDirectForcing(IceGrid &g, const NCConfigVariable &conf)
    : Model(g, conf) {}

  virtual ~PDirectForcing() {}

  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years)
  {
    PetscReal max_dt = -1;

    // "Periodize" the climate:
    t_years = my_mod(t_years);

    dt_years = temp.max_timestep(t_years);

    max_dt = smb.max_timestep(t_years);

    if (dt_years > 0) {
      if (max_dt > 0)
        dt_years = PetscMin(max_dt, dt_years);
    }
    else dt_years = max_dt;

    // If the user asked for periodized climate, limit time-steps so that PISM
    // never tries to average data over an interval that begins in one period and
    // ends in the next one.
    if (bc_period > 0.01)
      dt_years = PetscMin(dt_years, bc_period - t_years);

    return 0;
  }

  virtual void add_vars_to_output(string keyword, set<string> &result)
  {
    if (keyword != "small") {
      result.insert(temp_name);
      result.insert(smb_name);
    }
  }

  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype)
  {
    PetscErrorCode ierr;

    if (set_contains(vars, temp_name)) {
      ierr = temp.define(nc, nctype); CHKERRQ(ierr); 
    }

    if (set_contains(vars, smb_name)) {
      ierr = smb.define(nc, nctype); CHKERRQ(ierr);
    }

    return 0;
  }

  virtual PetscErrorCode write_variables(set<string> vars, string fname)
  {
    PetscErrorCode ierr;

    if (set_contains(vars, temp_name)) {
      ierr = temp.write(fname.c_str()); CHKERRQ(ierr);
    }

    if (set_contains(vars, smb_name)) {
      ierr = smb.write(fname.c_str()); CHKERRQ(ierr);
    }

    return 0;
  }

protected:
  IceModelVec2T temp, smb;
  string filename, temp_name, smb_name;

  PetscReal bc_period, bc_reference_year;
  bool enable_time_averaging;

  PetscErrorCode process_options()
  {
    PetscErrorCode ierr;
    bool bc_file_set, bc_period_set, bc_ref_year_set;
  
    bc_period = 0;
    bc_reference_year = 0;
    enable_time_averaging = false;

    ierr = PetscOptionsBegin(Model::grid.com, "", "Direct forcing options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-bc_file", "Specifies a file with top-surface boundary conditions",
                               filename, bc_file_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-bc_period", "Specifies the length of the climate data period",
                             bc_period, bc_period_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-bc_reference_year", "Boundary condition reference year",
                             bc_reference_year, bc_ref_year_set); CHKERRQ(ierr);
      ierr = PISMOptionsIsSet("-bc_time_average", "Enable time-averaging of boundary condition data",
                              enable_time_averaging); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (bc_file_set == false) {
      PetscPrintf(Model::grid.com, "PISM ERROR: option -bc_file is required.\n");
      PISMEnd();
    }

    return 0;
  }

  PetscErrorCode set_vec_parameters(string temp_std_name, string smb_std_name)
  {
    PetscErrorCode ierr;
    unsigned int buffer_size = (unsigned int) Model::config.get("climate_forcing_buffer_size"),
      temp_n_records = 1, smb_n_records = 1;

    NCTool nc(Model::grid.com, Model::grid.rank);
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    ierr = nc.get_nrecords(temp_name, temp_std_name, temp_n_records); CHKERRQ(ierr);
    ierr = nc.get_nrecords(smb_name,  smb_std_name,  smb_n_records);  CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    temp_n_records = PetscMin(temp_n_records, buffer_size);
    smb_n_records  = PetscMin(smb_n_records, buffer_size);

    if (temp_n_records < 1) {
      PetscPrintf(Model::grid.com, "PISM ERROR: Can't find '%s' (%s) in %s.\n",
                  temp_name.c_str(), temp_std_name.c_str(), filename.c_str());
      PISMEnd();

    }

    if (smb_n_records < 1) {
      PetscPrintf(Model::grid.com, "PISM ERROR: Can't find '%s' (%s) in %s.\n",
                  smb_name.c_str(), smb_std_name.c_str(), filename.c_str());
      PISMEnd();

    }

    temp.set_n_records(temp_n_records);
    smb.set_n_records(smb_n_records);

    temp.strict_timestep_limit = ! enable_time_averaging;
    smb.strict_timestep_limit   = ! enable_time_averaging;

    return 0;
  }

  virtual PetscErrorCode update_internal(PetscReal t_years, PetscReal dt_years)
  {
    PetscErrorCode ierr;

    // "Periodize" the climate:
    t_years = my_mod(t_years);

    if ((fabs(t_years - Model::t) < 1e-12) &&
        (fabs(dt_years - Model::dt) < 1e-12))
      return 0;

    Model::t  = t_years;
    Model::dt = dt_years;

    ierr = temp.update(Model::t, Model::dt); CHKERRQ(ierr);
    ierr = smb.update(Model::t, Model::dt); CHKERRQ(ierr);

    return 0;
  }

  PetscReal my_mod(PetscReal input)
  {
    if (bc_period < 0.01) return input;

    // compute time since the reference year:
    PetscReal delta = input - bc_reference_year;

    // compute delta mod bc_period:
    return delta - floor(delta / bc_period) * bc_period;
  }
};

class PSDirectForcing : public PDirectForcing<PISMSurfaceModel>
{
public:
  PSDirectForcing(IceGrid &g, const NCConfigVariable &conf)
    : PDirectForcing<PISMSurfaceModel>(g, conf)
  {
    temp_name = "artm";
    smb_name = "acab";
  }
  virtual ~PSDirectForcing() {}

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);

  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
};

class PADirectForcing : public PDirectForcing<PISMAtmosphereModel>
{
public:
  PADirectForcing(IceGrid &g, const NCConfigVariable &conf)
    : PDirectForcing<PISMAtmosphereModel>(g, conf)
  {
    temp_name = "artm";
    smb_name  = "precip";
  }

  virtual ~PADirectForcing() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);

  virtual PetscErrorCode mean_precip(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result); 
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);

  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();  
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
protected:
  vector<PetscReal> ts_mod;
};

#endif /* _PASDIRECTFORCING_H_ */
