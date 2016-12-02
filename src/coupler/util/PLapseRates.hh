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

#ifndef _PLAPSERATES_H_
#define _PLAPSERATES_H_

#include <cassert>

#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMTime.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec2T.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"

namespace pism {

template <class Model, class Mod>
class PLapseRates : public Mod {
public:
  PLapseRates(IceGrid::ConstPtr g, Model* in)
    : Mod(g, in) {
    m_temp_lapse_rate = 0.0;
  }

  virtual ~PLapseRates() {
    // empty
  }

protected:
  virtual MaxTimestep max_timestep_impl(double t) const {
    // "Periodize" the climate:
    t = Mod::m_grid->ctx()->time()->mod(t - m_bc_reference_time, m_bc_period);

    MaxTimestep input_max_dt = Mod::m_input_model->max_timestep(t);
    MaxTimestep surface_max_dt = m_reference_surface.max_timestep(t);

    if (input_max_dt.finite()) {
      return std::min(surface_max_dt, input_max_dt);
    } else {
      return surface_max_dt;
    }
  }

  virtual void update_impl(double my_t, double my_dt) {
    // a convenience
    double &t = Mod::m_t;
    double &dt = Mod::m_dt;

    // "Periodize" the climate:
    my_t = Mod::m_grid->ctx()->time()->mod(my_t - m_bc_reference_time,  m_bc_period);

    if ((fabs(my_t - t) < 1e-12) &&
        (fabs(my_dt - dt) < 1e-12)) {
      return;
    }

    t  = my_t;
    dt = my_dt;

    // NB! Input model uses original t and dt
    Mod::m_input_model->update(my_t, my_dt);

    m_reference_surface.update(t, dt);

    m_reference_surface.interp(t + 0.5*dt);
  }

  virtual void init_internal() {
    IceGrid::ConstPtr g = Mod::m_grid;

    options::String file(m_option_prefix + "_file",
                         "Specifies a file with top-surface boundary conditions");

    options::Integer period(m_option_prefix + "_period",
                            "Specifies the length of the climate data period", 0.0);

    options::Real reference_year(m_option_prefix + "_reference_year",
                                 "Boundary condition reference year", 0.0);

    options::Real T_lapse_rate("-temp_lapse_rate",
                               "Elevation lapse rate for the temperature, in K per km",
                               m_temp_lapse_rate);
    m_temp_lapse_rate = T_lapse_rate;

    if (not file.is_set()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "command-line option %s_file is required.",
                                    m_option_prefix.c_str());
    }

    if (reference_year.is_set()) {
      m_bc_reference_time = units::convert(Model::m_sys, reference_year, "years", "seconds");
    } else {
      m_bc_reference_time = 0;
    }

    if (period.value() < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid %s_period %d (period length cannot be negative)",
                                    m_option_prefix.c_str(), period.value());
    }
    m_bc_period = (unsigned int)period;

    if (not m_reference_surface.was_created()) {
      unsigned int buffer_size = (unsigned int) Mod::m_config->get_double("climate_forcing.buffer_size"),
        ref_surface_n_records = 1;

      PIO nc(g->com, "netcdf3", file, PISM_READONLY);
      ref_surface_n_records = nc.inq_nrecords("usurf", "surface_altitude",
                                              Mod::m_sys);
      nc.close();

      // if -..._period is not set, make n_records the minimum of the
      // buffer size and the number of available records. Otherwise try
      // to keep all available records in memory.
      if (not period.is_set()) {
        ref_surface_n_records = std::min(ref_surface_n_records, buffer_size);
      }

      if (ref_surface_n_records == 0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "can't find reference surface elevation (usurf) in %s.\n",
                                      file->c_str());
      }

      m_reference_surface.set_n_records(ref_surface_n_records);
      m_reference_surface.create(g, "usurf");
      m_reference_surface.set_attrs("climate_forcing",
                                  "reference surface for lapse rate corrections",
                                  "m", "surface_altitude");
      m_reference_surface.set_n_evaluations_per_year((unsigned int)Mod::m_config->get_double("climate_forcing.evaluations_per_year"));
    }

    Mod::m_log->message(2,
               "    reading reference surface elevation from %s ...\n",
               file->c_str());

    m_reference_surface.init(file, m_bc_period, m_bc_reference_time);
  }

  void lapse_rate_correction(IceModelVec2S &result, double lapse_rate) const {
    if (fabs(lapse_rate) < 1e-12) {
      return;
    }

    const IceModelVec2S
      &surface = *Mod::m_grid->variables().get_2d_scalar("surface_altitude");

    IceModelVec::AccessList list;
    list.add(surface);
    list.add(m_reference_surface);
    list.add(result);

    for (Points p(*Mod::m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      result(i, j) -= lapse_rate * (surface(i,j) - m_reference_surface(i, j));
    }
  }
protected:
  // "mutable" is needed here because some methods (average, interp) change the state of an
  // "IceModelVec2T"
  mutable IceModelVec2T m_reference_surface;
  unsigned int m_bc_period;
  double m_bc_reference_time,          // in seconds
    m_temp_lapse_rate;
  std::string m_option_prefix;
};


} // end of namespace pism

#endif /* _PLAPSERATES_H_ */
