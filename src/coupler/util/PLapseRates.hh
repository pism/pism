// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#include "IceGrid.hh"
#include "iceModelVec2T.hh"
#include "pism_options.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "PISMConfig.hh"
#include <assert.h>

#include "error_handling.hh"

namespace pism {

template <class Model, class Mod>
class PLapseRates : public Mod
{
public:
  PLapseRates(const IceGrid &g, Model* in)
    : Mod(g, in)
  {
    surface = thk = NULL;
    temp_lapse_rate = 0.0;
  }

  virtual ~PLapseRates() {}

  virtual void update(double my_t, double my_dt)
  {
    // a convenience
    double &t = Mod::m_t;
    double &dt = Mod::m_dt;

    // "Periodize" the climate:
    my_t = Mod::m_grid.time->mod(my_t - bc_reference_time,  bc_period);

    if ((fabs(my_t - t) < 1e-12) &&
        (fabs(my_dt - dt) < 1e-12)) {
      return;
    }

    t  = my_t;
    dt = my_dt;

    // NB! Input model uses original t and dt
    Mod::input_model->update(my_t, my_dt);

    reference_surface.update(t, dt);

    reference_surface.interp(t + 0.5*dt);
  }

  virtual void max_timestep(double t, double &dt, bool &restrict) {
    double max_dt = -1;

    // "Periodize" the climate:
    t = Mod::m_grid.time->mod(t - bc_reference_time, bc_period);

    Mod::input_model->max_timestep(t, dt, restrict);

    max_dt = reference_surface.max_timestep(t);

    if (restrict == true) {
      if (max_dt > 0) {
        dt = std::min(max_dt, dt);
      }
    } else {
      dt = max_dt;
    }

    if (dt > 0) {
      restrict = true;
    } else {
      restrict = false;
    }
  }

protected:
  IceModelVec2T reference_surface;
  const IceModelVec2S *surface, *thk;
  unsigned int bc_period;
  double bc_reference_time,          // in seconds
    temp_lapse_rate;
  std::string option_prefix;

  virtual void init_internal()
  {
    const IceGrid &g = Mod::m_grid;

    options::String file(option_prefix + "_file",
                         "Specifies a file with top-surface boundary conditions");

    options::Integer period(option_prefix + "_period",
                            "Specifies the length of the climate data period", 0.0);

    options::Real reference_year(option_prefix + "_reference_year",
                                 "Boundary condition reference year", 0.0);

    options::Real T_lapse_rate("-temp_lapse_rate",
                               "Elevation lapse rate for the temperature, in K per km",
                               temp_lapse_rate);
    temp_lapse_rate = T_lapse_rate;

    if (not file.is_set()) {
      throw RuntimeError::formatted("command-line option %s_file is required.",
                                    option_prefix.c_str());
    }

    if (reference_year.is_set()) {
      bc_reference_time = Model::m_grid.convert(reference_year, "years", "seconds");
    } else {
      bc_reference_time = 0;
    }

    if (period.value() < 0.0) {
      throw RuntimeError::formatted("invalid %s_period %d (period length cannot be negative)",
                                    option_prefix.c_str(), period.value());
    }
    bc_period = (unsigned int)period;

    if (reference_surface.was_created() == false) {
      unsigned int buffer_size = (unsigned int) Mod::m_config.get("climate_forcing_buffer_size"),
        ref_surface_n_records = 1;

      PIO nc(g.com, "netcdf3", g.config.get_unit_system());
      nc.open(file, PISM_READONLY);
      ref_surface_n_records = nc.inq_nrecords("usurf", "surface_altitude");
      nc.close();

      // if -..._period is not set, make n_records the minimum of the
      // buffer size and the number of available records. Otherwise try
      // to keep all available records in memory.
      if (not period.is_set()) {
        ref_surface_n_records = std::min(ref_surface_n_records, buffer_size);
      }

      if (ref_surface_n_records == 0) {
        throw RuntimeError::formatted("can't find reference surface elevation (usurf) in %s.\n",
                                      file->c_str());
      }

      reference_surface.set_n_records(ref_surface_n_records);
      reference_surface.create(g, "usurf", false);
      reference_surface.set_attrs("climate_forcing",
                                  "reference surface for lapse rate corrections",
                                  "m", "surface_altitude");
      reference_surface.set_n_evaluations_per_year((unsigned int)Mod::m_config.get("climate_forcing_evaluations_per_year"));
    }

    verbPrintf(2, g.com,
               "    reading reference surface elevation from %s ...\n",
               file->c_str());

    reference_surface.init(file, bc_period, bc_reference_time);

    surface = g.variables().get_2d_scalar("surface_altitude");
    thk     = g.variables().get_2d_scalar("land_ice_thickness");
  }

  void lapse_rate_correction(IceModelVec2S &result, double lapse_rate)
  {
    if (fabs(lapse_rate) < 1e-12) {
      return;
    }

    IceModelVec::AccessList list;
    list.add(*thk);
    list.add(*surface);
    list.add(reference_surface);
    list.add(result);

    for (Points p(Mod::m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((*thk)(i,j) > 0) {
        const double correction = lapse_rate * ((*surface)(i,j) - reference_surface(i,j));
        result(i,j) -= correction;
      }
    }
  }
};


} // end of namespace pism

#endif /* _PLAPSERATES_H_ */
