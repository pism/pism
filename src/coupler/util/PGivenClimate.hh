// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec2T.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/PISMComponent.hh"

namespace pism {

template <class Model, class Input>
class PGivenClimate : public Model
{
public:
  PGivenClimate(IceGrid::ConstPtr g, Input *in)
    : Model(g, in) {}

  virtual ~PGivenClimate() {
    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      delete k->second;
      ++k;
    }
  }

protected:
  virtual MaxTimestep max_timestep_impl(double t) const {
    (void) t;
    return MaxTimestep();
  }

  virtual void define_model_state_impl(const PIO &output) const {
    std::map<std::string, IceModelVec2T*>::const_iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      k->second->define(output);
      ++k;
    }
  }

  virtual void write_model_state_impl(const PIO &output) const {
    std::map<std::string, IceModelVec2T*>::const_iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      k->second->write(output);
      ++k;
    }
  }

  void process_options()
  {
    options::String file(m_option_prefix + "_file",
                         "Specifies a file with boundary conditions");
    if (file.is_set()) {
      m_filename = file;
      Model::m_log->message(2,
                 "  - Reading boundary conditions from '%s'...\n",
                 m_filename.c_str());
    } else {
      m_filename = process_input_options(Model::m_grid->com).filename;

      Model::m_log->message(2,
                            "  - Option %s_file is not set. Trying the input file '%s'...\n",
                            m_option_prefix.c_str(), m_filename.c_str());
    }

    options::Integer period(m_option_prefix + "_period",
                            "Specifies the length of the climate data period (in years)", 0);
    if (period.value() < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid %s_period %d (period length cannot be negative)",
                                    m_option_prefix.c_str(), period.value());
    }
    m_bc_period = (unsigned int)period;

    options::Integer ref_year(m_option_prefix + "_reference_year",
                              "Boundary condition reference year", 0);
    if (ref_year.is_set()) {
      m_bc_reference_time = units::convert(Model::m_sys, ref_year, "years", "seconds");
    } else {
      m_bc_reference_time = 0;
    }
  }

  void set_vec_parameters(const std::map<std::string, std::string> &standard_names)
  {
    unsigned int buffer_size = (unsigned int) Model::m_config->get_double("climate_forcing.buffer_size");

    PIO nc(Model::m_grid->com, "netcdf3", m_filename, PISM_READONLY);

    std::map<std::string, IceModelVec2T*>::const_iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      unsigned int n_records = 0;
      const std::string &short_name = k->first;
      std::string standard_name;
      if (standard_names.find(short_name) != standard_names.end()) {
        standard_name = standard_names.find(short_name)->second;
      }
      // else leave standard_name empty

      n_records = nc.inq_nrecords(short_name, standard_name,
                                  Model::m_grid->ctx()->unit_system());

      // If -..._period is not set, make ..._n_records the minimum of the
      // buffer size and the number of available records. Otherwise try
      // to keep all available records in memory.
      if (m_bc_period == 0) {
        n_records = std::min(n_records, buffer_size);
      }

      if (n_records < 1) {
        // If the variable was not found we allocate storage for one
        // record. This is needed to be able to allocate and then
        // discard an "-atmosphere given" model (atmosphere::Given) when
        // "-surface given" (Given) is selected.
        n_records = 1;
      }

      (k->second)->set_n_records(n_records);

      (k->second)->set_n_evaluations_per_year((unsigned int)Model::m_config->get_double("climate_forcing.evaluations_per_year"));

      ++k;
    }

    nc.close();
  }

  virtual void update_internal(double my_t, double my_dt)
  {
    // "Periodize" the climate:
    my_t = Model::m_grid->ctx()->time()->mod(my_t - m_bc_reference_time, m_bc_period);

    if ((fabs(my_t - Model::m_t) < 1e-12) &&
        (fabs(my_dt - Model::m_dt) < 1e-12)) {
      return;
    }

    Model::m_t  = my_t;
    Model::m_dt = my_dt;

    if (Model::m_input_model != NULL) {
      Model::m_input_model->update(Model::m_t, Model::m_dt);
    }

    std::map<std::string, IceModelVec2T*>::iterator k = m_fields.begin();
    while(k != m_fields.end()) {
      (k->second)->update(Model::m_t, Model::m_dt);

      ++k;
    }
  }
protected:
  std::map<std::string, IceModelVec2T*> m_fields;
  std::string m_filename, m_option_prefix;

  unsigned int m_bc_period;       // in (integer) years
  double m_bc_reference_time;  // in seconds
};

} // end of namespace pism

#endif /* _PGIVENCLIMATE_H_ */
