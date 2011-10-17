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
#include "PISMTime.hh"
#include "NCTool.hh"

template <class Model>
class PDirectForcing : public Model
{
public:
  PDirectForcing(IceGrid &g, const NCConfigVariable &conf)
    : Model(g, conf) {}

  virtual ~PDirectForcing() {}

  virtual PetscErrorCode max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict)
  {
    PetscReal max_dt = -1;

    // "Periodize" the climate:
    my_t = Model::grid.time->mod(my_t - bc_reference_time, bc_period);

    my_dt = temp.max_timestep(my_t);

    max_dt = mass_flux.max_timestep(my_t);

    if (my_dt > 0) {
      if (max_dt > 0)
        my_dt = PetscMin(max_dt, my_dt);
    }
    else my_dt = max_dt;

    // If the user asked for periodized climate, limit time-steps so that PISM
    // never tries to average data over an interval that begins in one period and
    // ends in the next one.
    if (bc_period > 0.01)
      my_dt = PetscMin(my_dt, bc_period - my_t);

    if (my_dt > 0)
      restrict = true;
    else
      restrict = false;

    return 0;
  }

  virtual void add_vars_to_output(string keyword, set<string> &result)
  {
    if (keyword != "small") {
      result.insert(temp_name);
      result.insert(mass_flux_name);
    }
  }

  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype)
  {
    PetscErrorCode ierr;

    if (set_contains(vars, temp_name)) {
      ierr = temp.define(nc, nctype); CHKERRQ(ierr); 
    }

    if (set_contains(vars, mass_flux_name)) {
      ierr = mass_flux.define(nc, nctype); CHKERRQ(ierr);
    }

    return 0;
  }

  virtual PetscErrorCode write_variables(set<string> vars, string fname)
  {
    PetscErrorCode ierr;

    if (set_contains(vars, temp_name)) {
      ierr = temp.write(fname.c_str()); CHKERRQ(ierr);
    }

    if (set_contains(vars, mass_flux_name)) {
      ierr = mass_flux.write(fname.c_str()); CHKERRQ(ierr);
    }

    return 0;
  }

protected:
  IceModelVec2T temp, mass_flux;
  string filename, temp_name, mass_flux_name, option_prefix;

  PetscReal bc_period,          // in seconds
    bc_reference_time;          // in seconds
  bool enable_time_averaging;

  PetscErrorCode process_options()
  {
    PetscErrorCode ierr;
    bool bc_file_set, bc_period_set, bc_ref_year_set;
  
    PetscReal bc_period_years = 0,
      bc_reference_year = 0;

    enable_time_averaging = false;
    bc_period = 0;
    bc_reference_time = 0;

    ierr = PetscOptionsBegin(Model::grid.com, "", "Direct forcing options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString(option_prefix + "_file",
                               "Specifies a file with boundary conditions",
                               filename, bc_file_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_period",
                             "Specifies the length of the climate data period",
                             bc_period_years, bc_period_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_reference_year",
                             "Boundary condition reference year",
                             bc_reference_year, bc_ref_year_set); CHKERRQ(ierr);
      ierr = PISMOptionsIsSet(option_prefix + "_time_average",
                              "Enable time-averaging of boundary condition data",
                              enable_time_averaging); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (bc_file_set == false) {
      PetscPrintf(Model::grid.com, "PISM ERROR: option %s_file is required.\n", option_prefix.c_str());
      PISMEnd();
    }

    if (bc_ref_year_set) {
      bc_reference_time = Model::grid.time->years_to_seconds(bc_reference_year);
    }

    if (bc_period_set) {
      bc_period = Model::grid.time->years_to_seconds(bc_period_years);
    }

    return 0;
  }

  PetscErrorCode set_vec_parameters(string temp_std_name, string mass_flux_std_name)
  {
    PetscErrorCode ierr;
    unsigned int buffer_size = (unsigned int) Model::config.get("climate_forcing_buffer_size"),
      temp_n_records = 1, mass_flux_n_records = 1;

    NCTool nc(Model::grid.com, Model::grid.rank);
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    ierr = nc.get_nrecords(temp_name, temp_std_name, temp_n_records); CHKERRQ(ierr);
    ierr = nc.get_nrecords(mass_flux_name,  mass_flux_std_name,  mass_flux_n_records);  CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    temp_n_records = PetscMin(temp_n_records, buffer_size);
    mass_flux_n_records  = PetscMin(mass_flux_n_records, buffer_size);

    if (temp_n_records < 1) {
      PetscPrintf(Model::grid.com, "PISM ERROR: Can't find '%s' (%s) in %s.\n",
                  temp_name.c_str(), temp_std_name.c_str(), filename.c_str());
      PISMEnd();

    }

    if (mass_flux_n_records < 1) {
      PetscPrintf(Model::grid.com, "PISM ERROR: Can't find '%s' (%s) in %s.\n",
                  mass_flux_name.c_str(), mass_flux_std_name.c_str(), filename.c_str());
      PISMEnd();

    }

    temp.set_n_records(temp_n_records);
    mass_flux.set_n_records(mass_flux_n_records);

    return 0;
  }

  virtual PetscErrorCode update_internal(PetscReal my_t, PetscReal my_dt)
  {
    PetscErrorCode ierr;

    // "Periodize" the climate:
    my_t = Model::grid.time->mod(my_t - bc_reference_time, bc_period);

    if ((fabs(my_t - Model::t) < 1e-12) &&
        (fabs(my_dt - Model::dt) < 1e-12))
      return 0;

    Model::t  = my_t;
    Model::dt = my_dt;

    ierr = temp.update(Model::t, Model::dt); CHKERRQ(ierr);
    ierr = mass_flux.update(Model::t, Model::dt); CHKERRQ(ierr);

    return 0;
  }
};

class PSDirectForcing : public PDirectForcing<PISMSurfaceModel>
{
public:
  PSDirectForcing(IceGrid &g, const NCConfigVariable &conf)
    : PDirectForcing<PISMSurfaceModel>(g, conf)
  {
    temp_name = "artm";
    mass_flux_name = "acab";
    option_prefix = "-surface";
  }
  virtual ~PSDirectForcing() {}

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

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
    mass_flux_name  = "precip";
    option_prefix = "-atmosphere";
  }

  virtual ~PADirectForcing() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

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
