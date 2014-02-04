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

#include "PSTemperatureIndex.hh"
#include "localMassBalance.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "PISMAtmosphere.hh"
#include "Mask.hh"

#include <assert.h>

///// PISM surface model implementing a PDD scheme.

PSTemperatureIndex::PSTemperatureIndex(IceGrid &g, const PISMConfig &conf)
  : PISMSurfaceModel(g, conf),
    ice_surface_temp(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_PSTemperatureIndex(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PSTemperatureIndex::~PSTemperatureIndex() {
  delete mbscheme;
  delete faustogreve;
}

PetscErrorCode PSTemperatureIndex::allocate_PSTemperatureIndex() {
  PetscErrorCode ierr;
  bool flag;

  mbscheme		= NULL;
  faustogreve		= NULL;
  base_ddf.snow		= config.get("pdd_factor_snow");
  base_ddf.ice		= config.get("pdd_factor_ice");
  base_ddf.refreezeFrac = config.get("pdd_refreeze");
  base_pddStdDev	= config.get("pdd_std_dev");
  base_pddThresholdTemp = config.get("pdd_positive_threshold_temp");

  ierr = PetscOptionsBegin(grid.com, "",
                           "Temperature-index (PDD) scheme for surface (snow) processes", "");
                           CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-pdd_rand",
                            "Use a PDD implementation based on simulating a random process",
			    randomized); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-pdd_rand_repeatable",
                            "Use a PDD implementation based on simulating a repeatable random process",
			    randomized_repeatable); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-pdd_fausto",
                            "Set PDD parameters using formulas (6) and (7) in [Faustoetal2009]",
			    fausto_params); CHKERRQ(ierr);

    ierr = PISMOptionsReal("-pdd_factor_snow", "PDD snow factor",
                           base_ddf.snow, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_factor_ice", "PDD ice factor",
                           base_ddf.ice, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_refreeze", "PDD refreeze fraction",
                           base_ddf.refreezeFrac, flag); CHKERRQ(ierr);

    ierr = PISMOptionsReal("-pdd_std_dev", "PDD standard deviation",
                           base_pddStdDev, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_positive_threshold_temp",
                           "PDD uses this temp in K to determine 'positive' temperatures",
                           base_pddThresholdTemp, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (randomized_repeatable) {
    mbscheme = new PDDrandMassBalance(config, true);
  } else if (randomized) {
    mbscheme = new PDDrandMassBalance(config, false);
  } else {
    mbscheme = new PDDMassBalance(config);
  }

  if (fausto_params) {
    faustogreve = new FaustoGrevePDDObject(grid, config);
  }

  ierr = climatic_mass_balance.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_attrs("diagnostic",
					 "instantaneous ice-equivalent surface mass balance (accumulation/ablation) rate",
					 "m s-1",  // m *ice-equivalent* per second
					 "land_ice_surface_specific_mass_balance");  // CF standard_name
  CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  climatic_mass_balance.write_in_glaciological_units = true;
  climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  ierr = accumulation_rate.create(grid, "saccum", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = accumulation_rate.set_attrs("diagnostic",
                                     "instantaneous ice-equivalent surface accumulation rate"
                                     " (precipitation minus rain)",
                                     "m s-1",
                                     ""); CHKERRQ(ierr);
  ierr = accumulation_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  accumulation_rate.write_in_glaciological_units = true;

  ierr = melt_rate.create(grid, "smelt", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = melt_rate.set_attrs("diagnostic",
                             "instantaneous ice-equivalent surface melt rate",
                             "m s-1",
                             ""); CHKERRQ(ierr);
  ierr = melt_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  melt_rate.write_in_glaciological_units = true;

  ierr = runoff_rate.create(grid, "srunoff", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = runoff_rate.set_attrs("diagnostic",
                               "instantaneous ice-equivalent surface meltwater runoff rate",
                               "m s-1",
                               ""); CHKERRQ(ierr);
  ierr = runoff_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  runoff_rate.write_in_glaciological_units = true;

  ierr = snow_depth.create(grid, "snow_depth", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = snow_depth.set_attrs("diagnostic",
			      "snow cover depth (set to zero once a year)",
			      "m", ""); CHKERRQ(ierr);
  ierr = snow_depth.set(0.0); CHKERRQ(ierr);

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                  "ice temperature at the ice surface");
  ierr = ice_surface_temp.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSTemperatureIndex::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = PISMSurfaceModel::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
    "* Initializing the default temperature-index, PDD-based surface processes scheme.\n"
    "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
    "  Surface mass balance and ice upper surface temperature are outputs.\n"
    "  See PISM User's Manual for control of degree-day factors.\n");
    CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
    "  Computing number of positive degree-days by: "); CHKERRQ(ierr);
  if (randomized_repeatable) {
    ierr = verbPrintf(2, grid.com, "simulation of a random process.\n"); CHKERRQ(ierr);
  } else if (randomized) {
    ierr = verbPrintf(2, grid.com, "repeatable simulation of a random process.\n"); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "an expectation integral.\n"); CHKERRQ(ierr);
  }

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  assert(mask != NULL);

  if ((config.get("pdd_std_dev_lapse_lat_rate") != 0.0) || fausto_params) {
    lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
    assert(lat != NULL);
  } else {
    lat = NULL;
  }

  if (fausto_params) {
    ierr = verbPrintf(2, grid.com,
       "  Setting PDD parameters from [Faustoetal2009] ...\n");
       CHKERRQ(ierr);

    base_pddStdDev = 2.53;

    lon = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
    assert(lon != NULL);

    usurf = dynamic_cast<IceModelVec2S*>(vars.get("usurf"));
    assert(usurf != NULL);
  } else {
    // generally, this is the case in which degree day factors do not depend
    //   on location; we use base_ddf
    lon = NULL;
    usurf = NULL;
  }

  std::string input_file;
  bool do_regrid = false;
  int start = -1;
  
  // find PISM input file to read data from:
  ierr = find_pism_input(input_file, do_regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate from file
  ierr = verbPrintf(2, grid.com,
		    "    reading snow depth (ice equivalent meters) from %s ... \n",
		    input_file.c_str()); CHKERRQ(ierr);
  ierr = snow_depth.regrid(input_file, OPTIONAL, 0.0); CHKERRQ(ierr);

  m_next_balance_year_start = compute_next_balance_year_start(grid.time->current());

  return 0;
}

PetscErrorCode PSTemperatureIndex::max_timestep(double my_t, double &my_dt, bool &restrict) {
  PetscErrorCode ierr;

  ierr = atmosphere->max_timestep(my_t, my_dt, restrict); CHKERRQ(ierr);

  return 0;
}

double PSTemperatureIndex::compute_next_balance_year_start(double time) {
    // compute the time corresponding to the beginning of the next balance year
    double
      balance_year_start_day = config.get("pdd_balance_year_start_day"),
      one_day                = grid.convert(1.0, "days", "seconds"),
      year_start             = grid.time->calendar_year_start(time),
      balance_year_start     = year_start + (balance_year_start_day - 1.0) * one_day;

    if (balance_year_start > time) {
      return balance_year_start;
    }
    return grid.time->increment_date(balance_year_start, 1);
}


PetscErrorCode PSTemperatureIndex::update(double my_t, double my_dt) {
  PetscErrorCode ierr;

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12))
    return 0;

  m_t  = my_t;
  m_dt = my_dt;

  // upate to ensure that temperature and precipitation time series
  // are correct:
  ierr = atmosphere->update(my_t, my_dt); CHKERRQ(ierr);

  // set up air temperature and precipitation time series
  int Nseries = mbscheme->get_timeseries_length(my_dt);

  const double dtseries = my_dt / Nseries;
  std::vector<double> ts(Nseries), T(Nseries), P(Nseries), PDDs(Nseries);
  for (int k = 0; k < Nseries; ++k)
    ts[k] = my_t + k * dtseries;

  MaskQuery m(*mask);
  ierr = mask->begin_access(); CHKERRQ(ierr);

  if (lat != NULL) {
    ierr = lat->begin_access(); CHKERRQ(ierr);
  }

  if (faustogreve != NULL) {
    assert(lat != NULL && lon != NULL && usurf != NULL);
    ierr = lon->begin_access(); CHKERRQ(ierr);
    ierr = usurf->begin_access(); CHKERRQ(ierr);
    ierr = faustogreve->update_temp_mj(usurf, lat, lon); CHKERRQ(ierr);
  }

  const double sigmalapserate = config.get("pdd_std_dev_lapse_lat_rate"),
    sigmabaselat   = config.get("pdd_std_dev_lapse_lat_base");
  if (sigmalapserate != 0.0) {
    assert(lat != NULL);
  }

  DegreeDayFactors  ddf = base_ddf;

  ierr = atmosphere->begin_pointwise_access(); CHKERRQ(ierr);
  ierr = climatic_mass_balance.begin_access(); CHKERRQ(ierr);

  ierr = accumulation_rate.begin_access(); CHKERRQ(ierr);
  ierr = melt_rate.begin_access(); CHKERRQ(ierr);
  ierr = runoff_rate.begin_access(); CHKERRQ(ierr);
  ierr = snow_depth.begin_access(); CHKERRQ(ierr);

  ierr = atmosphere->init_timeseries(&ts[0], Nseries); CHKERRQ(ierr);

  for (int i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j<grid.ys+grid.ym; ++j) {

      // the temperature time series from the PISMAtmosphereModel and its modifiers
      ierr = atmosphere->temp_time_series(i, j, &T[0]); CHKERRQ(ierr);

      // the precipitation time series from PISMAtmosphereModel and its modifiers
      ierr = atmosphere->precip_time_series(i, j, &P[0]); CHKERRQ(ierr);

      if (faustogreve != NULL) {
	// we have been asked to set mass balance parameters according to
	//   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
	ierr = faustogreve->setDegreeDayFactors(i, j, (*usurf)(i, j),
                                                (*lat)(i, j), (*lon)(i, j), ddf);
        CHKERRQ(ierr);
      }

      // Use temperature time series, the "positive" threshhold, and
      // the standard deviation of the daily variability to get the
      // number of positive degree days (PDDs)
      double sigma = base_pddStdDev;
      if (sigmalapserate != 0.0) {
        sigma += sigmalapserate * ((*lat)(i,j) - sigmabaselat);
      }
      mbscheme->get_PDDs(sigma, dtseries, &T[0], Nseries, &PDDs[0]);

      // Use temperature time series to remove rainfall from precipitation
      mbscheme->get_snow_accumulation(&P[0], // precipitation rate (input-output)
                                      &T[0], // air temperature (input)
                                      Nseries);

      // Use degree-day factors, and number of PDDs, and the snow
      // precipitation, to get surface mass balance (and diagnostics:
      // accumulation, melt, runoff)
      {
        double next_snow_depth_reset = m_next_balance_year_start;
        accumulation_rate(i,j)     = 0.0;
        melt_rate(i,j)             = 0.0;
        runoff_rate(i,j)           = 0.0;
        climatic_mass_balance(i,j) = 0.0;
        for (int k = 0; k < Nseries; ++k) {
          if (ts[k] >= next_snow_depth_reset) {
            snow_depth(i,j)       = 0.0;
            next_snow_depth_reset = grid.time->increment_date(next_snow_depth_reset, 1);
          }

          double accumulation     = P[k] * dtseries;
          accumulation_rate(i,j) += accumulation;

          mbscheme->step(ddf, PDDs[k], accumulation,
                         snow_depth(i,j), melt_rate(i,j), runoff_rate(i,j), climatic_mass_balance(i,j));
        }
        accumulation_rate(i,j)     /= m_dt;
        melt_rate(i,j)             /= m_dt;
        runoff_rate(i,j)           /= m_dt;
        climatic_mass_balance(i,j) /= m_dt;
      }

      if (m.ocean(i,j)) {
        snow_depth(i,j) = 0.0;  // snow over the ocean does not stick
      }
    } // j-loop
  } // i-loop

  ierr = accumulation_rate.end_access(); CHKERRQ(ierr);
  ierr = melt_rate.end_access(); CHKERRQ(ierr);
  ierr = runoff_rate.end_access(); CHKERRQ(ierr);
  ierr = snow_depth.end_access(); CHKERRQ(ierr);

  ierr = climatic_mass_balance.end_access(); CHKERRQ(ierr);
  ierr = atmosphere->end_pointwise_access(); CHKERRQ(ierr);

  ierr = mask->end_access(); CHKERRQ(ierr);

  if (lat != NULL) {
    ierr = lat->end_access(); CHKERRQ(ierr);
  }

  if (faustogreve != NULL) {
    ierr = lon->end_access(); CHKERRQ(ierr);
    ierr = usurf->end_access(); CHKERRQ(ierr);
  }

  m_next_balance_year_start = compute_next_balance_year_start(grid.time->current());
  
  return 0;
}


PetscErrorCode PSTemperatureIndex::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = climatic_mass_balance.copy_to(result); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PSTemperatureIndex::ice_surface_temperature(IceModelVec2S &result) {

  PetscErrorCode ierr = atmosphere->mean_annual_temp(result); CHKERRQ(ierr);

  return 0;
}

void PSTemperatureIndex::add_vars_to_output(std::string keyword, std::set<std::string> &result) {

  PISMSurfaceModel::add_vars_to_output(keyword, result);

  result.insert("snow_depth");
  
  if (keyword == "medium" || keyword == "big") {
    result.insert("climatic_mass_balance");
    result.insert("ice_surface_temp");
  }

  if (keyword == "big") {
    result.insert("saccum");
    result.insert("smelt");
    result.insert("srunoff");
  }
}

PetscErrorCode PSTemperatureIndex::define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = ice_surface_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "saccum")) {
    ierr = accumulation_rate.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "smelt")) {
    ierr = melt_rate.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "srunoff")) {
    ierr = runoff_rate.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "snow_depth")) {
    ierr = snow_depth.define(nc, nctype); CHKERRQ(ierr);
  }

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode PSTemperatureIndex::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = ice_surface_temp;

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);
    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.write(nc); CHKERRQ(ierr);
    vars.erase("climatic_mass_balance");
  }

  if (set_contains(vars, "saccum")) {
    ierr = accumulation_rate.write(nc); CHKERRQ(ierr);
    vars.erase("saccum");
  }

  if (set_contains(vars, "smelt")) {
    ierr = melt_rate.write(nc); CHKERRQ(ierr);
    vars.erase("smelt");
  }

  if (set_contains(vars, "srunoff")) {
    ierr = runoff_rate.write(nc); CHKERRQ(ierr);
    vars.erase("srunoff");
  }

  if (set_contains(vars, "snow_depth")) {
    ierr = snow_depth.write(nc); CHKERRQ(ierr);
    vars.erase("snow_depth");
  }

  ierr = PISMSurfaceModel::write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}
