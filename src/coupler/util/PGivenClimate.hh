// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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
#include "PISMConfig.hh"
#include "PIO.hh"
#include "pism_options.hh"

#include "error_handling.hh"

namespace pism {

template <class Model, class Input>
class PGivenClimate : public Model
{
public:
  PGivenClimate(IceGrid &g, Input *in)
    : Model(g, in) {}

  virtual ~PGivenClimate() {
    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      delete k->second;
      ++k;
    }
  }

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result)
  {
    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      result.insert(k->first);
      ++k;
    }

    if (Model::input_model != NULL) {
      Model::input_model->add_vars_to_output(keyword, result);
    }

  }

  virtual void define_variables(const std::set<std::string> &vars_input, const PIO &nc, IO_Type nctype)
  {
    std::set<std::string> vars = vars_input;
    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      if (set_contains(vars, k->first)) {
        (k->second)->define(nc, nctype);
        vars.erase(k->first);
      }
      ++k;
    }

    if (Model::input_model != NULL) {
      Model::input_model->define_variables(vars, nc, nctype);
    }
  }

  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc) {

    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {

      if (set_contains(vars, k->first)) {
        (k->second)->write(nc);
      }

      ++k;
    }

    if (Model::input_model != NULL) {
      Model::input_model->write_variables(vars, nc);
    }
  }

protected:
  std::map<std::string, IceModelVec2T*> m_fields;
  std::string filename, option_prefix;

  unsigned int bc_period;       // in (integer) years
  double bc_reference_time;  // in seconds

  void process_options()
  {
    options::String file(option_prefix + "_file",
                         "Specifies a file with boundary conditions");
    if (file.is_set()) {
      filename = file;
      verbPrintf(2, Model::m_grid.com,
                 "  - Reading boundary conditions from '%s'...\n",
                 filename.c_str());
    } else {
      // find PISM input file to read data from:
      bool do_regrid; int start;   // will be ignored
      Model::find_pism_input(filename, do_regrid, start);

      verbPrintf(2, Model::m_grid.com,
                 "  - Option %s_file is not set. Trying the input file '%s'...\n",
                 option_prefix.c_str(), filename.c_str());
    }

    options::Integer period(option_prefix + "_period",
                            "Specifies the length of the climate data period (in years)", 0);
    if (period.value() < 0.0) {
      throw RuntimeError::formatted("invalid %s_period %d (period length cannot be negative)",
                                    option_prefix.c_str(), period.value());
    }
    bc_period = (unsigned int)period;

    options::Integer ref_year(option_prefix + "_reference_year",
                              "Boundary condition reference year", 0);
    if (ref_year.is_set()) {
      bc_reference_time = Model::m_grid.convert(ref_year, "years", "seconds");
    } else {
      bc_reference_time = 0;
    }
  }

  void set_vec_parameters(std::map<std::string, std::string> standard_names)
  {
    unsigned int buffer_size = (unsigned int) Model::m_config.get("climate_forcing_buffer_size");

    PIO nc(Model::m_grid.com, "netcdf3", Model::m_grid.config.get_unit_system());
    nc.open(filename, PISM_READONLY);

    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      unsigned int n_records = 0;
      std::string short_name = k->first,
        standard_name = standard_names[short_name];

      n_records = nc.inq_nrecords(short_name, standard_name);

      // If -..._period is not set, make ..._n_records the minimum of the
      // buffer size and the number of available records. Otherwise try
      // to keep all available records in memory.
      if (bc_period == 0) {
        n_records = std::min(n_records, buffer_size);
      }

      if (n_records < 1) {
        // If the variable was not found we allocate storage for one
        // record. This is needed to be able to allocate and then
        // discard an "-atmosphere given" model (PAGivenClimate) when
        // "-surface given" (PSGivenClimate) is selected.
        n_records = 1;
      }

      (k->second)->set_n_records(n_records);

      (k->second)->set_n_evaluations_per_year((unsigned int)Model::m_config.get("climate_forcing_evaluations_per_year"));

      ++k;
    }

    nc.close();
  }

  virtual void update_internal(double my_t, double my_dt)
  {
    // "Periodize" the climate:
    my_t = Model::m_grid.time->mod(my_t - bc_reference_time, bc_period);

    if ((fabs(my_t - Model::m_t) < 1e-12) &&
        (fabs(my_dt - Model::m_dt) < 1e-12)) {
      return;
    }

    Model::m_t  = my_t;
    Model::m_dt = my_dt;

    if (Model::input_model != NULL) {
      Model::input_model->update(Model::m_t, Model::m_dt);
    }

    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      (k->second)->update(Model::m_t, Model::m_dt);

      ++k;
    }
  }
};

} // end of namespace pism

#endif /* _PGIVENCLIMATE_H_ */
