// Copyright (C) 2011, 2012, 2013 PISM Authors
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

  virtual ~PGivenClimate() {
    map<string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      delete k->second;
      ++k;
    }
  }

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result)
  {
    map<string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      result[k->first] = (k->second)->get_metadata();
      ++k;
    }

    if (Model::input_model != NULL) {
      Model::input_model->add_vars_to_output(keyword, result);
    }

  }

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype)
  {
    PetscErrorCode ierr;

    map<string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      if (set_contains(vars, k->first)) {
        ierr = (k->second)->define(nc, nctype);
        vars.erase(k->first);
      }
      ++k;
    }

    if (Model::input_model != NULL) {
      ierr = Model::input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    return 0;
  }

  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc)
  {
    PetscErrorCode ierr;

    map<string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {

      if (set_contains(vars, k->first)) {
        ierr = (k->second)->write(nc); CHKERRQ(ierr);
      }

      ++k;
    }

    if (Model::input_model != NULL) {
      ierr = Model::input_model->write_variables(vars, nc); CHKERRQ(ierr);
    }

    return 0;
  }

protected:
  map<string, IceModelVec2T*> m_fields;
  string filename, option_prefix;

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

  PetscErrorCode set_vec_parameters(map<string, string> standard_names)
  {
    PetscErrorCode ierr;
    unsigned int buffer_size = (unsigned int) Model::config.get("climate_forcing_buffer_size");

    PIO nc(Model::grid.com, Model::grid.rank, "netcdf3",
           Model::grid.get_unit_system());
    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);

    map<string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      unsigned int n_records = 0;
      string short_name = k->first,
        standard_name = standard_names[short_name];

      ierr = nc.inq_nrecords(short_name, standard_name, n_records); CHKERRQ(ierr);

      // If -..._period is not set, make ..._n_records the minimum of the
      // buffer size and the number of available records. Otherwise try
      // to keep all available records in memory.
      if (bc_period == 0.0)
        n_records = PetscMin(n_records, buffer_size);

      if (n_records < 1) {
        PetscPrintf(Model::grid.com, "PISM ERROR: Can't find '%s' (%s) in %s.\n",
                    short_name.c_str(), standard_name.c_str(), filename.c_str());
        PISMEnd();
      }

      (k->second)->set_n_records(n_records);

      (k->second)->set_n_evaluations_per_year((unsigned int)Model::config.get("climate_forcing_evaluations_per_year"));

      ++k;
    }

    ierr = nc.close(); CHKERRQ(ierr);

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

    if (Model::input_model != NULL) {
      ierr = Model::input_model->update(Model::t, Model::dt); CHKERRQ(ierr);
    }

    map<string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      ierr = (k->second)->update(Model::t, Model::dt); CHKERRQ(ierr);

      ++k;
    }

    return 0;
  }
};

#endif /* _PGIVENCLIMATE_H_ */
