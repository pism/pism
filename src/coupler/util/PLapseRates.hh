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

#ifndef _PLAPSERATES_H_
#define _PLAPSERATES_H_

#include "IceGrid.hh"
#include "iceModelVec2T.hh"
#include "pism_options.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"

template <class Model, class Mod>
class PLapseRates : public Mod
{
public:
  PLapseRates(IceGrid &g, const NCConfigVariable &conf, Model* in)
    : Mod(g, conf, in)
  {
    surface = thk = NULL;
    temp_lapse_rate = 0;
  }

  virtual ~PLapseRates() {}

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  {
    PetscErrorCode ierr;

    // a convenience
    PetscReal &m_t = Mod::t;
    PetscReal &m_dt = Mod::dt;

    // "Periodize" the climate:
    my_t = Mod::grid.time->mod(my_t - bc_reference_time,  bc_period);

    if ((fabs(my_t - m_t) < 1e-12) &&
        (fabs(my_dt - m_dt) < 1e-12))
      return 0;

    m_t = my_t;
    m_dt = my_dt;

    ierr = Mod::input_model->update(my_t, my_dt); CHKERRQ(ierr);

    ierr = reference_surface.update(m_t, m_dt); CHKERRQ(ierr);

    ierr = reference_surface.at_time(m_t); CHKERRQ(ierr);

    return 0;
  }

  virtual PetscErrorCode max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict) {
    PetscErrorCode ierr;
    PetscReal max_dt = -1;

    // "Periodize" the climate:
    my_t = Mod::grid.time->mod(my_t - bc_reference_time, bc_period);

    ierr = Mod::input_model->max_timestep(my_t, my_dt, restrict); CHKERRQ(ierr);

    max_dt = reference_surface.max_timestep(my_t);

    if (my_dt > 0) {
      if (max_dt > 0)
        my_dt = PetscMin(max_dt, my_dt);
    }
    else my_dt = max_dt;

    // If the user asked for periodized climate, limit time-steps so that PISM
    // never tries to average data over an interval that begins in one period and
    // ends in the next one.
    if (bc_period > 1e-6)
      my_dt = PetscMin(my_dt, bc_period - my_t);

    if (my_dt > 0)
      restrict = true;
    else
      restrict = false;

    return 0;
  }

protected:
  IceModelVec2T reference_surface;
  IceModelVec2S *surface, *thk;
  PetscReal bc_period,          // in seconds
    bc_reference_time,          // in seconds
    temp_lapse_rate;
  string option_prefix;

  virtual PetscErrorCode init_internal(PISMVars &vars)
  {
    PetscErrorCode ierr;
    string filename;
    bool bc_file_set, bc_period_set, bc_ref_year_set, temp_lapse_rate_set;

    IceGrid &g = Mod::grid;

    PetscReal bc_period_years = 0,
      bc_reference_year = 0;

    ierr = PetscOptionsBegin(g.com, "", "Lapse rate options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString(option_prefix + "_file",
                               "Specifies a file with top-surface boundary conditions",
                               filename, bc_file_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_period",
                             "Specifies the length of the climate data period",
                             bc_period_years, bc_period_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_reference_year",
                             "Boundary condition reference year",
                             bc_reference_year, bc_ref_year_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-temp_lapse_rate",
                             "Elevation lapse rate for the temperature, in K per km",
                             temp_lapse_rate, temp_lapse_rate_set); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (bc_file_set == false) {
      PetscPrintf(Model::grid.com, "PISM ERROR: option %s_file is required.\n", option_prefix.c_str());
      PISMEnd();
    }

    if (bc_ref_year_set) {
      bc_reference_time = Model::grid.time->years_to_seconds(bc_reference_year);
    } else {
      bc_reference_time = 0;
    }

    if (bc_period_set) {
      bc_period = Model::grid.time->years_to_seconds(bc_period_years);
    } else {
      bc_period = 0;
    }

    unsigned int buffer_size = (unsigned int) Mod::config.get("climate_forcing_buffer_size"),
      ref_surface_n_records = 1;

    PIO nc(g.com, g.rank, "netcdf3");
    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_nrecords("usurf", "surface_altitude", ref_surface_n_records); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ref_surface_n_records = PetscMin(ref_surface_n_records, buffer_size);

    if (ref_surface_n_records == 0) {
      PetscPrintf(g.com, "PISM ERROR: can't find reference surface elevation (usurf) in %s.\n",
                  filename.c_str());
      PISMEnd();
    }

    ierr = verbPrintf(2,g.com,
                      "    reading reference surface elevation from %s ...\n",
                      filename.c_str()); CHKERRQ(ierr);

    reference_surface.set_n_records(ref_surface_n_records);
    ierr = reference_surface.create(g, "usurf", false); CHKERRQ(ierr);
    ierr = reference_surface.set_attrs("climate_forcing",
                                       "reference surface for lapse rate corrections",
                                       "m", "surface_altitude"); CHKERRQ(ierr);
    ierr = reference_surface.init(filename); CHKERRQ(ierr);

    surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
    if (!surface) SETERRQ(g.com, 1, "ERROR: 'usurf' is not available or is wrong type in dictionary");

    thk = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
    if (!thk) SETERRQ(g.com, 1, "ERROR: 'ice thickness' is not available or is wrong type in dictionary");

    return 0;
  }

  PetscErrorCode lapse_rate_correction(IceModelVec2S &result, PetscReal lapse_rate)
  {
    PetscErrorCode ierr;

    if (PetscAbs(lapse_rate) < 1e-12)
      return 0;

    ierr = thk->begin_access(); CHKERRQ(ierr);
    ierr = surface->begin_access(); CHKERRQ(ierr);
    ierr = reference_surface.begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access(); CHKERRQ(ierr);

    IceGrid &g = Mod::grid;

    for (PetscInt   i = g.xs; i < g.xs + g.xm; ++i) {
      for (PetscInt j = g.ys; j < g.ys + g.ym; ++j) {
        if ((*thk)(i,j) > 0)
          result(i,j) -= lapse_rate * ((*surface)(i,j) - reference_surface(i,j));
      }
    }

    ierr = result.end_access(); CHKERRQ(ierr);
    ierr = reference_surface.end_access(); CHKERRQ(ierr);
    ierr = surface->end_access(); CHKERRQ(ierr);
    ierr = thk->end_access(); CHKERRQ(ierr);

    return 0;
  }
};


#endif /* _PLAPSERATES_H_ */
