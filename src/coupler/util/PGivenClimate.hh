// Copyright (C) 2011, 2012 PISM Authors
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

#ifndef _PGIVENCLIMATE_H_
#define _PGIVENCLIMATE_H_

#include "iceModelVec2T.hh"
#include "PISMTime.hh"
#include "PIO.hh"
#include "pism_options.hh"

template <class Model, class Input>
class PGivenClimate : public Model
{
public:
  PGivenClimate(IceGrid &g, const NCConfigVariable &conf, Input *in)
    : Model(g, conf, in) {}

  virtual ~PGivenClimate() {}

  virtual PetscErrorCode max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict)
  {
    PetscReal mass_flux_max_dt = -1;

    // "Periodize" the climate:
    my_t = Model::grid.time->mod(my_t - bc_reference_time, bc_period);

    my_dt = temp.max_timestep(my_t);

    mass_flux_max_dt = mass_flux.max_timestep(my_t);

    if (my_dt > 0) {
      if (mass_flux_max_dt > 0)
        my_dt = PetscMin(mass_flux_max_dt, my_dt);
    }
    else my_dt = mass_flux_max_dt;

    // If the user asked for periodized climate, limit time-steps so that PISM
    // never tries to average data over an interval that begins in one period and
    // ends in the next one.
    if (bc_period > 0.01)
      my_dt = PetscMin(my_dt, bc_period - my_t);

    // my_dt is fully determined (in the case input_model == NULL). Now get
    // max_dt from an input model:

    if (Model::input_model != NULL) {
      PetscReal input_dt;
      bool input_restrict;

      // Note: we use "periodized" t here:
      PetscErrorCode ierr = Model::input_model->max_timestep(my_t, input_dt, input_restrict); CHKERRQ(ierr);

      if (input_restrict) {
        if (my_dt > 0)
          my_dt = PetscMin(input_dt, my_dt);
        else
          my_dt = input_dt;
      }
      // else my_dt is not changed

    }

    if (my_dt > 0)
      restrict = true;
    else
      restrict = false;

    return 0;
  }

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result)
  {
    if (Model::input_model != NULL) {
      Model::input_model->add_vars_to_output(keyword, result);
    }

    result[temp_name] = temp.get_metadata();
    result[mass_flux_name] = mass_flux.get_metadata();
  }

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype)
  {
    PetscErrorCode ierr;

    if (set_contains(vars, temp_name)) {
      ierr = temp.define(nc, nctype); CHKERRQ(ierr);
      vars.erase(temp_name);
    }

    if (set_contains(vars, mass_flux_name)) {
      ierr = mass_flux.define(nc, nctype); CHKERRQ(ierr);
      vars.erase(mass_flux_name);
    }

    if (Model::input_model != NULL) {
      ierr = Model::input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    return 0;
  }

  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc)
  {
    PetscErrorCode ierr;

    if (set_contains(vars, temp_name)) {
      ierr = temp.write(nc); CHKERRQ(ierr);
      vars.erase(temp_name);
    }

    if (set_contains(vars, mass_flux_name)) {
      ierr = mass_flux.write(nc); CHKERRQ(ierr);
      vars.erase(mass_flux_name);
    }

    if (Model::input_model != NULL) {
      ierr = Model::input_model->write_variables(vars, nc); CHKERRQ(ierr);
    }

    return 0;
  }

protected:
  IceModelVec2T temp, mass_flux;
  string filename, temp_name, mass_flux_name, option_prefix;

  PetscReal bc_period,          // in seconds
    bc_reference_time;          // in seconds

  PetscErrorCode process_options()
  {
    PetscErrorCode ierr;
    bool bc_file_set, bc_period_set, bc_ref_year_set;

    PetscReal bc_period_years = 0,
      bc_reference_year = 0;

    bc_period = 0;
    bc_reference_time = 0;

    ierr = PetscOptionsBegin(Model::grid.com, "", "Climate forcing options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString(option_prefix + "_file",
                               "Specifies a file with boundary conditions",
                               filename, bc_file_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_period",
                             "Specifies the length of the climate data period (in years)",
                             bc_period_years, bc_period_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_reference_year",
                             "Boundary condition reference year",
                             bc_reference_year, bc_ref_year_set); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (bc_file_set == false) {
      // find PISM input file to read data from:
      bool regrid; int start;   // will be ignored
      ierr = Model::find_pism_input(filename, regrid, start); CHKERRQ(ierr);

      ierr = verbPrintf(2, Model::grid.com,
                        "  - Option %s_file is not set. Trying the input file '%s'...\n",
                        option_prefix.c_str(), filename.c_str());
      CHKERRQ(ierr);

    } else {
      ierr = verbPrintf(2, Model::grid.com,
                        "  - Reading boundary conditions from '%s'...\n",
                        filename.c_str()); CHKERRQ(ierr);
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

    PIO nc(Model::grid.com, Model::grid.rank, "netcdf3");
    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_nrecords(temp_name, temp_std_name, temp_n_records); CHKERRQ(ierr);
    ierr = nc.inq_nrecords(mass_flux_name,  mass_flux_std_name,  mass_flux_n_records);  CHKERRQ(ierr);
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

    if (Model::input_model != NULL) {
      ierr = Model::input_model->update(Model::t, Model::dt); CHKERRQ(ierr);
    }

    return 0;
  }
};

#endif /* _PGIVENCLIMATE_H_ */
