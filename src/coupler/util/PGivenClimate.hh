// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include <memory>

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Component.hh"

namespace pism {

class Geometry;

template <class Model>
class PGivenClimate : public Model
{
public:
  PGivenClimate(IceGrid::ConstPtr g, std::shared_ptr<Model> in)
    : Model(g, in) {}

  virtual ~PGivenClimate() {
    // empty
  }

protected:
  virtual MaxTimestep max_timestep_impl(double t) const {
    (void) t;
    return MaxTimestep();
  }

  virtual void define_model_state_impl(const PIO &output) const {
    for (auto f : m_fields) {
      f.second->define(output);
    }
  }

  virtual void write_model_state_impl(const PIO &output) const {
    for (auto f : m_fields) {
      f.second->write(output);
    }
  }

  std::string process_options(const std::string &option_prefix)
  {
    std::string filename;
    {
      options::String file(option_prefix + "_file",
                           "Specifies a file with boundary conditions");
      if (file.is_set()) {
        filename = file;
        Model::m_log->message(2,
                              "  - Reading boundary conditions from '%s'...\n",
                              filename.c_str());
      } else {
        filename = process_input_options(Model::m_grid->com).filename;

        Model::m_log->message(2,
                              "  - Option %s_file is not set. Trying the input file '%s'...\n",
                              option_prefix.c_str(), filename.c_str());
      }
    }

    {
      options::Integer period(option_prefix + "_period",
                              "Specifies the length of the climate data period (in years)", 0);
      if (period.value() < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid %s_period %d (period length cannot be negative)",
                                      option_prefix.c_str(), period.value());
      }
      m_bc_period = (unsigned int)period;
    }

    {
      options::Integer ref_year(option_prefix + "_reference_year",
                                "Boundary condition reference year", 0);
      if (ref_year.is_set()) {
        m_bc_reference_time = units::convert(Model::m_sys, ref_year, "years", "seconds");
      } else {
        m_bc_reference_time = 0;
      }
    }

    return filename;
  }

  IceModelVec2T::Ptr allocate(const PIO &file,
                              units::System::Ptr unit_system,
                              const std::string &short_name,
                              const std::string &standard_name,
                              int max_buffer_size,
                              bool periodic) {
    int evaluations_per_year = Model::m_config->get_double("climate_forcing.evaluations_per_year");

    int n_records = file.inq_nrecords(short_name, standard_name, unit_system);

    if (not periodic) {
      n_records = std::min(n_records, max_buffer_size);
    }
    // In the periodic case we try to keep all the records in RAM.

    // Allocate storage for one record if the variable was not found. This is needed to be
    // able to cheaply allocate and then discard an "-atmosphere given" model
    // (atmosphere::Given) when "-surface given" (Given) is selected.
    n_records = std::max(n_records, 1);

    return IceModelVec2T::Ptr(new IceModelVec2T(Model::m_grid, short_name, n_records,
                                                evaluations_per_year));
  }

  virtual void update_internal(const Geometry &geometry, double t, double dt)
  {
    if (Model::m_input_model) {
      Model::m_input_model->update(geometry, t, dt);
    }

    for (auto f : m_fields) {
      f.second->update(t, dt);
    }
  }
protected:
  std::map<std::string, std::shared_ptr<IceModelVec2T> > m_fields;
  std::string m_filename;

  unsigned int m_bc_period;       // in (integer) years
  double m_bc_reference_time;  // in seconds
};

} // end of namespace pism

#endif /* _PGIVENCLIMATE_H_ */
