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

#ifndef _PASLAPSERATES_H_
#define _PASLAPSERATES_H_

#include "PISMAtmosphere.hh"
#include "PISMSurface.hh"

template <class Model, class Modifier>
class PLapseRates : public Modifier
{
public:
  PLapseRates(IceGrid &g, const NCConfigVariable &conf, Model* in)
    : PISMComponent_TS(g, conf), Modifier(g, conf, in), input(in)
  {
    surface = thk = NULL;
    temp_lapse_rate = 0;
  }

  virtual ~PLapseRates() {}

  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years)
  {
    PetscErrorCode ierr;

    PetscReal &m_t = Modifier::t;
    PetscReal &m_dt = Modifier::dt;

    // "Periodize" the climate:
    t_years = my_mod(t_years);

    if ((fabs(t_years - m_t) < 1e-12) &&
        (fabs(dt_years - m_dt) < 1e-12))
      return 0;

    m_t = t_years;
    m_dt = dt_years;

    ierr = input->update(t_years, dt_years); CHKERRQ(ierr);

    ierr = reference_surface.update(m_t, m_dt); CHKERRQ(ierr);

    if (enable_time_averaging) {
      ierr = reference_surface.interp(m_t + 0.5*m_dt); CHKERRQ(ierr);
    } else {
      ierr = reference_surface.get_record_years(m_t); CHKERRQ(ierr);
    }

    return 0;
  }

  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years) {
    PetscErrorCode ierr;
    PetscReal max_dt = -1;

    // "Periodize" the climate:
    t_years = my_mod(t_years);

    ierr = input->max_timestep(t_years, dt_years); CHKERRQ(ierr);

    max_dt = reference_surface.max_timestep(t_years);

    if (dt_years > 0) {
      if (max_dt > 0)
        dt_years = PetscMin(max_dt, dt_years);
    }
    else dt_years = max_dt;

    // If the user asked for periodized climate, limit time-steps so that PISM
    // never tries to average data over an interval that begins in one period and
    // ends in the next one.
    if (bc_period > 1e-6)
      dt_years = PetscMin(dt_years, bc_period - t_years);

    return 0;
  }

  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
    PetscErrorCode ierr = input->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode write_variables(set<string> vars, string filename) {
    PetscErrorCode ierr = input->write_variables(vars, filename); CHKERRQ(ierr);
    return 0;
  }

protected:
  IceModelVec2T reference_surface;
  IceModelVec2S *surface, *thk;
  PetscReal bc_period, bc_reference_year, temp_lapse_rate;
  bool enable_time_averaging;
  Model *input;

  virtual PetscErrorCode common_init(PISMVars &vars)
  {
    PetscErrorCode ierr;
    string filename;
    bool bc_file_set, bc_period_set, bc_ref_year_set, temp_lapse_rate_set;

    IceGrid &g = Modifier::grid;

    enable_time_averaging = false;

    ierr = PetscOptionsBegin(g.com, "", "Lapse rate options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-bc_file", "Specifies a file with top-surface boundary conditions",
                               filename, bc_file_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-bc_period", "Specifies the length of the climate data period",
                             bc_period, bc_period_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-bc_reference_year", "Boundary condition reference year",
                             bc_reference_year, bc_ref_year_set); CHKERRQ(ierr);
      ierr = PISMOptionsIsSet("-bc_time_average", "Enable time-averaging of boundary condition data",
                              enable_time_averaging); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-temp_lapse_rate",
                             "Elevation lapse rate for the temperature, in K per km",
                             temp_lapse_rate, temp_lapse_rate_set); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (bc_file_set == false) {
      PetscPrintf(g.com, "PISM ERROR: option -bc_file is required.\n");
      PISMEnd();
    }

    if (bc_period_set == false) {
      bc_period = 0;
    }

    if (bc_ref_year_set == false) {
      bc_reference_year = 0;
    }

    unsigned int buffer_size = (unsigned int) Modifier::config.get("climate_forcing_buffer_size"),
      ref_surface_n_records = 1;

    NCTool nc(g.com, g.rank);
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    ierr = nc.get_nrecords("usurf", "surface_altitude", ref_surface_n_records); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ref_surface_n_records = PetscMin(ref_surface_n_records, buffer_size);

    ierr = verbPrintf(2,g.com,
                      "    reading reference surface elevation from %s ...\n",
                      filename.c_str()); CHKERRQ(ierr);

    reference_surface.set_n_records(ref_surface_n_records);
    ierr = reference_surface.create(g, "usurf", false); CHKERRQ(ierr);
    ierr = reference_surface.set_attrs("climate_forcing",
                                       "reference surface for lapse rate corrections",
                                       "m", "surface_altitude"); CHKERRQ(ierr);
    ierr = reference_surface.init(filename); CHKERRQ(ierr);

    reference_surface.strict_timestep_limit = ! enable_time_averaging;

    surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
    if (!surface) SETERRQ(1, "ERROR: 'usurf' is not available or is wrong type in dictionary");

    thk = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
    if (!thk) SETERRQ(1, "ERROR: 'ice thickness' is not available or is wrong type in dictionary");

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
    
    IceGrid &g = Modifier::grid;

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

  //! \brief Computes year modulo bc_period if bc_period is active.
  PetscReal my_mod(PetscReal in) {
    if (bc_period < 1e-6) return in;

    // compute time since the reference year:
    PetscReal delta = in - bc_reference_year;

    // compute delta mod bc_period:
    return delta - floor(delta / bc_period) * bc_period;}
};

class PSLapseRates : public PLapseRates<PISMSurfaceModel,PSModifier>
{
public:
  PSLapseRates(IceGrid &g, const NCConfigVariable &conf, PISMSurfaceModel* in)
    : PISMComponent_TS(g, conf), PLapseRates<PISMSurfaceModel,PSModifier>(g, conf, in)
  {
    smb_lapse_rate = 0;
  }

  virtual ~PSLapseRates() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_liquid_water_fraction(IceModelVec2S &result);
  virtual PetscErrorCode mass_held_in_surface_layer(IceModelVec2S &result);
  virtual PetscErrorCode surface_layer_thickness(IceModelVec2S &result);
protected:
  PetscReal smb_lapse_rate;
};

class PALapseRates : public PLapseRates<PISMAtmosphereModel,PAModifier>
{
public:
  PALapseRates(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in)
    : PISMComponent_TS(g, conf), PLapseRates<PISMAtmosphereModel,PAModifier>(g, conf, in)
  {
    precip_lapse_rate = 0;
  }

  virtual ~PALapseRates() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode mean_precip(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result);

  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();

  virtual PetscErrorCode temp_time_series(int i, int j, int N,
                                          PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);
protected:
  PetscReal precip_lapse_rate;
};

#endif /* _PASLAPSERATES_H_ */
