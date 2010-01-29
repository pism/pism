// Copyright (C) 2008-2010 Ed Bueler, Constantine Khroulev, Gudfinna Adalgeirsdottir, and Andy Aschwanden
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


#include <petscda.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModelVec.hh"
#include "../base/LocalInterpCtx.hh"
#include "../base/nc_util.hh"
#include "../base/NCVariable.hh"
#include "../base/Timeseries.hh"
#include "localMassBalance.hh"
#include "iceModelVec2T.hh"
#include "pccoupler.hh"
#include "pGreenlandAtmosCoupler.hh"


PISMGreenlandAtmosCoupler::PISMGreenlandAtmosCoupler() : PISMAtmosphereCoupler() {
  snowtempmaps = NULL;
  snowprecipmaps = NULL;
  mbscheme = NULL;
}


PISMGreenlandAtmosCoupler::~PISMGreenlandAtmosCoupler() {
  delete snowtempmaps;
  delete snowprecipmaps;
  delete mbscheme;
}


PetscErrorCode PISMGreenlandAtmosCoupler::setLMBScheme(LocalMassBalance *usethisscheme) {
  delete mbscheme;
  mbscheme = usethisscheme;
  return 0;
}


//! Initializes the snow process model.
/*!
Performs the many actions which are needed to initialize the snow process model:

-# finds the input file,
-# initializes pointers to IceModel members which have information needed for parameterizing climate inputs (namely, latitude, longitude, surface elevation),
-# reads vsnowprecip from file,
-# initializes annual mean snow temperature anomalies if -anomaly_temp_ma is given
-# initializes snow precipitation anomalies if -anomaly_precip is given
-# chooses among mass balance schemes (from options \c -pdd, \c -pdd_rand, \c -pdd_rand_repeatable),
-# and reads vsurftemp from file.

Regarding the major variables, vsurfmassflux is a computed result of the actions
of this derived class, vsnowtemp_ma and vsnowtemp_mj are intermediate quantities to determine annual
cycle of temperature in PDD scheme (i.e. LocalMassBalance), and vsnowprecip a time-invariant map
which is read from the input file.  The map vsurftemp is also time-invariant, and is read from file.
 */
PetscErrorCode PISMGreenlandAtmosCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  PetscTruth   optSet;

  ierr = PISMAtmosphereCoupler::initFromOptions(g, variables); CHKERRQ(ierr);

  char  filename[PETSC_MAX_PATH_LEN];
  LocalInterpCtx* lic;
  bool regrid;
  int start;
  // allocates and initializes lic (if necessary)
  ierr = findPISMInputFile((char*)filename, lic, regrid, start); CHKERRQ(ierr);

  ierr = verbPrintf(2, g->com, 
       "  initializing atmospheric climate coupler with a snow process model ...\n"); CHKERRQ(ierr); 

  surfelev = dynamic_cast<IceModelVec2*>(variables.get("surface_altitude"));
  if (!surfelev) SETERRQ(1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2*>(variables.get("latitude"));
  if (!lat) SETERRQ(1, "ERROR: latitude is not available");

  lon = dynamic_cast<IceModelVec2*>(variables.get("longitude"));
  if (!lon) SETERRQ(1, "ERROR: longitude is not available");

  // clear out; will be overwritten by mass balance model
  ierr = vsurfmassflux.set_attr("pism_intent","climate_diagnostic"); CHKERRQ(ierr);
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr);

  // check on whether we should read snow temperature anomalies
  char anomalies_file[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-anomaly_temp_ma", 
             anomalies_file, PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
  if (optSet) {
    // stop if -dTforcing is set:
    PetscTruth dTforcing_set;
    ierr = check_option("-dTforcing", dTforcing_set); CHKERRQ(ierr);
    if (dTforcing_set) {
      ierr = PetscPrintf(g->com, "PISM ERROR: option -anomaly_temp_ma is incompatible with -dTforcing.\n");
      PetscEnd();
    }

    ierr = verbPrintf(2,grid->com,
       "    reading snow temperature anomalies from %s ...\n",
       anomalies_file); CHKERRQ(ierr);

    snowtempmaps = new IceModelVec2T;
    snowtempmaps->set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
    ierr = snowtempmaps->create(*grid, "delta_snowtemp", false); CHKERRQ(ierr);
    ierr = snowtempmaps->set_attrs("climate_forcing", "snow surface temperature anomalies",
				   "Kelvin", ""); CHKERRQ(ierr);
    ierr = snowtempmaps->init(anomalies_file); CHKERRQ(ierr);
  }

  // check on whether we should read snow precipitation anomalies
  ierr = PetscOptionsGetString(PETSC_NULL, "-anomaly_precip", 
             anomalies_file, PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
  if (optSet) {
    ierr = verbPrintf(2,grid->com,
       "    reading ice-equivalent snow precipitation rate anomalies from %s ...\n",
       anomalies_file); CHKERRQ(ierr);

    snowprecipmaps = new IceModelVec2T;
    snowprecipmaps->set_n_records((unsigned int) config.get("climate_forcing_buffer_size")); 
    ierr = snowprecipmaps->create(*grid, "delta_snowprecip", false); CHKERRQ(ierr);
    ierr = snowprecipmaps->set_attrs("climate_forcing",
				     "anomalies of ice-equivalent snow precipitation rate",
				     "m s-1", ""); CHKERRQ(ierr);
    ierr = snowprecipmaps->init(anomalies_file); CHKERRQ(ierr);
    ierr = snowprecipmaps->set_glaciological_units("m year-1");
    snowprecipmaps->write_in_glaciological_units = true;
  }
  
  // create mean annual ice equivalent snow precipitation rate (before melt, and not including rain)
  ierr = vsnowprecip.create(*g, "snowprecip", false); CHKERRQ(ierr);
  ierr = vsnowprecip.set_attrs("climate_state", 
			       "mean annual ice-equivalent snow precipitation rate",
			       "m s-1", 
			       "");  // no CF standard_name ??
  CHKERRQ(ierr);
  ierr = vsnowprecip.set_glaciological_units("m year-1");
  vsnowprecip.write_in_glaciological_units = true;

  // read snow precipitation rate from file
  ierr = verbPrintf(2, g->com, 
		    "    reading mean annual ice-equivalent snow precipitation rate 'snowprecip'\n"
		    "      from %s ... \n",
		    filename); CHKERRQ(ierr); 
  if (regrid) {
    ierr = vsnowprecip.regrid(filename, *lic, true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = vsnowprecip.read(filename, start); CHKERRQ(ierr); // fails if not found!
  }


  ierr = verbPrintf(2,grid->com,
		    "    using default snow-surface temperature parameterization\n"); CHKERRQ(ierr);
  ierr = vsnowtemp_mj.create(*g, "snowtemp_mj", false); CHKERRQ(ierr);
  ierr = vsnowtemp_mj.set_attrs("climate_state",
				"mean July snow-surface temperature used in mass balance scheme",
				"Kelvin",
				""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = vsnowtemp_mj.set(273.15); CHKERRQ(ierr);  // merely a default value

  // check if user wants default or random PDD; initialize mbscheme to one of these PDDs
  if (mbscheme == NULL) { // only read user options if scheme is not chosen yet;
                          //   derived class could choose
    PetscTruth  pddRandSet, pddRepeatableSet;
    ierr = check_option("-pdd_rand", pddRandSet); CHKERRQ(ierr);
    ierr = check_option("-pdd_rand_repeatable", pddRepeatableSet); CHKERRQ(ierr);
    if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
      ierr = verbPrintf(2,grid->com,
         "    mass balance scheme chosen: PDD with simulated random daily variability ...\n");
         CHKERRQ(ierr);
      mbscheme = new PDDrandMassBalance(config,(pddRepeatableSet == PETSC_TRUE));
    } else { // default case
      ierr = verbPrintf(2,grid->com,
         "    mass balance scheme chosen: PDD with expected value for daily variability ...\n");
         CHKERRQ(ierr);
      mbscheme = new PDDMassBalance(config);
    }
  }
  ierr = mbscheme->init(); CHKERRQ(ierr);

  // check if artm=vsurftemp is available in input file, and if so warn user that it
  //   will be ignored and overwritten by parameterization
  NCTool nc(grid);
  bool surftemp_exists;
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  ierr = nc.find_variable("artm", "", NULL, surftemp_exists); CHKERRQ(ierr);
  ierr = nc.close();
  if (surftemp_exists) {
    ierr = verbPrintf(1, g->com, 
        "  WARNING: 'artm' found in file %s; IGNORING IT and using parameterization instead ...\n",
        filename); CHKERRQ(ierr);
  }

  delete lic;

  return 0;
}

PetscErrorCode PISMGreenlandAtmosCoupler::writeCouplingFieldsToFile(
                              PetscScalar t_years, const char *filename) {
  PetscErrorCode ierr;
  IceModelVec2 tmp;  // no work vectors, so we create a new one
  ierr = tmp.create(*grid, "temporary_vector", false); CHKERRQ(ierr);
  
  ierr = PISMAtmosphereCoupler::writeCouplingFieldsToFile(t_years,filename); CHKERRQ(ierr);

  // paleo-precipitation:
  double dT_offset = 0.0;
  if (dTforcing != NULL) {
    dT_offset = (*dTforcing)(t_years);

    // *if* dTforcing then modify precipitation, which is assumed to be present-day, to be 
    //   modeled paleo-precipitation; see formula for P_G(t,x,y) on SeaRISE
    //   "Model Initialization" page http://websrv.cs.umt.edu/isis/index.php/Model_Initialization
    //   the version here assumes Delta T_SC = 0.0; write out with new name

    ierr = tmp.set_name("snowprecip_paleo"); CHKERRQ(ierr);
    ierr = tmp.set_attrs("climate_diagnostic", 
			 "mean annual, time-dependent, paleo, ice-equivalent snow precipitation rate",
			 "m s-1", 
			 ""); CHKERRQ(ierr);
    ierr = tmp.set_glaciological_units("m year-1");
    tmp.write_in_glaciological_units = true;

    ierr = tmp.copy_from(vsnowprecip); CHKERRQ(ierr);
    const PetscScalar precipexpfactor
       = exp( config.get("precip_exponential_factor_for_temperature") * dT_offset );
    ierr = tmp.scale(precipexpfactor); CHKERRQ(ierr);

    ierr = tmp.write(filename, NC_FLOAT); CHKERRQ(ierr);
  } // end of if (dTforcing != NULL)

  // temperature snapshot:
  ierr = tmp.set_name("snowtemp"); CHKERRQ(ierr);
  ierr = tmp.set_attrs("climate_diagnostic",
		       "time-dependent snow-surface temperature used in mass balance scheme",
		       "K",
		       ""); CHKERRQ(ierr);  // use no CF standard_name

  PetscScalar **T, **T_ma, **T_mj;
  const PetscScalar radpersec = 2.0 * pi / secpera, // radians per second frequency for annual cycle
    sperd = 8.64e4, // exact number of seconds per day
    julydaysec = sperd * config.get("snow_temp_july_day"),
    t = (t_years - floor(t_years)) * secpera,
    cosfactor = cos(radpersec * (t - julydaysec));

  // note that T_ma already has dT_offset or anomalies applied by
  // updateSurfTempAndProvide()
  ierr = tmp.get_array(T);  CHKERRQ(ierr);
  ierr = vsurftemp.get_array(T_ma); CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj); CHKERRQ(ierr);
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      T[i][j] = T_ma[i][j] + (T_mj[i][j] - T_ma[i][j]) * cosfactor;
    }
  }
  ierr = vsurftemp.end_access(); CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  ierr = tmp.end_access();  CHKERRQ(ierr);

  ierr = tmp.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // write the current anomaly to the file (if applicable):
  if (snowtempmaps != NULL) {
    ierr = snowtempmaps->update(t_years, 0); CHKERRQ(ierr);
    ierr = snowtempmaps->interp(t_years); CHKERRQ(ierr);
    ierr = snowtempmaps->write(filename, NC_FLOAT); CHKERRQ(ierr);
  }

  // precipitation:
  if (snowprecipmaps != NULL) {
    ierr = snowprecipmaps->update(t_years, 0); CHKERRQ(ierr);
    ierr = snowprecipmaps->interp(t_years); CHKERRQ(ierr);
    ierr = snowprecipmaps->write(filename, NC_FLOAT); CHKERRQ(ierr);

    ierr = tmp.copy_from(vsnowprecip); CHKERRQ(ierr);
    ierr = tmp.add(1.0, *snowprecipmaps); CHKERRQ(ierr);

    ierr = tmp.set_name("snowprecip_inst"); CHKERRQ(ierr);
    ierr = tmp.set_attrs("climate_diagnostic", 
			 "instantaneous ice-equivalent snow precipitation rate",
			 "m s-1", ""); CHKERRQ(ierr);
    ierr = tmp.set_glaciological_units("m year-1");
    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(filename, NC_FLOAT); CHKERRQ(ierr);
  }

  // duplicate original precipitation map:
  ierr = vsnowprecip.write(filename, NC_FLOAT); CHKERRQ(ierr);

  return 0;
}


/*!
The default method here is the Fausto et al parameterization scheme
appropriate to the Greenland ice sheet.  The parameterization depends linearly
on surface elevation, latitude, and longitude.
  
See formulas (1) and (2) and Table 3 in [\ref Faustoetal2009].
  
The default values for lapse rates, which are stored in src/pism_config.cdl, use values 
from the lines 'Best annual fit: This study with \f$\kappa_{\text{ma}}\f$' 
and 'Best July fit: This study with \f$\kappa_{\text{ma}}\f$' from Table 3.
But any scheme that has the same kind of dependence on the same quantities
can be used just by changing the "snow_temp_fausto_..." configuration parameters 
in src/pism_config.cdl.  Other temperature parameterization schemes will need to
re-implement this procedure.
 */
PetscErrorCode PISMGreenlandAtmosCoupler::parameterizedUpdateSnowSurfaceTemp(
		PetscScalar /*t_years*/, PetscScalar /*dt_years*/) {
  PetscErrorCode ierr;
  const PetscScalar 
    d_ma     = config.get("snow_temp_fausto_d_ma"),      // K
    gamma_ma = config.get("snow_temp_fausto_gamma_ma"),  // K m-1
    c_ma     = config.get("snow_temp_fausto_c_ma"),      // K (degN)-1
    kappa_ma = config.get("snow_temp_fausto_kappa_ma"),  // K (degW)-1
    d_mj     = config.get("snow_temp_fausto_d_mj"),      // SAME UNITS as for _ma ...
    gamma_mj = config.get("snow_temp_fausto_gamma_mj"),
    c_mj     = config.get("snow_temp_fausto_c_mj"),
    kappa_mj = config.get("snow_temp_fausto_kappa_mj");
  
  PetscScalar **lat_degN, **lon_degE, **h, **T_ma, **T_mj;
  ierr = surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = lon->get_array(lon_degE); CHKERRQ(ierr);
  ierr = vsurftemp.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);

  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      T_ma[i][j] = d_ma + gamma_ma * h[i][j] + c_ma * lat_degN[i][j] + kappa_ma * (-lon_degE[i][j]);
      T_mj[i][j] = d_mj + gamma_mj * h[i][j] + c_mj * lat_degN[i][j] + kappa_mj * (-lon_degE[i][j]);
    }
  }
  
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = lon->end_access(); CHKERRQ(ierr);
  ierr = vsurftemp.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  return 0;
}


/*! Note that in this case surface temperature always comes from the
  parameterization, so shifting it around is safe (it will be re-created next
  time it's needed).
 */
PetscErrorCode PISMGreenlandAtmosCoupler::updateSurfTempAndProvide(
                  PetscScalar t_years, PetscScalar dt_years, 
                  IceModelVec2* &pvst) {
  PetscErrorCode ierr;

  ierr = parameterizedUpdateSnowSurfaceTemp(t_years, dt_years); CHKERRQ(ierr);

  if (dTforcing != NULL) {
    TsOffset = (*dTforcing)(t_years);
    ierr = vsurftemp.shift(TsOffset); CHKERRQ(ierr);  // apply the new offset
  }

  if (snowtempmaps != NULL) {
    ierr = snowtempmaps->update(t_years, 0); CHKERRQ(ierr);
    ierr = snowtempmaps->interp(t_years); CHKERRQ(ierr);

    ierr = vsurftemp.add(1.0, *snowtempmaps); CHKERRQ(ierr);
  }

  pvst = &vsurftemp;
  
  return 0;
}


//! Compute the ice surface mass flux from the snow temperature, a stored map of snow precipication rate, and a choice of mass balance scheme (typically a PDD).
/*! First this method recomputes the snow temperature, either from a
parameterization or from interpolating stored snow temperature maps. At each
point on the surface a temperature time series with short (weekly or less) time
steps is generated. There is also a parameterized yearly cycle of temperature and
additional weather-related variability according to a normally distributed
random temperature change for each week and grid point.

At each point on the ice surface a temperature time series, for the duration of the 
time step specified in calling this routine, is used by a LocalMassBalance object
to compute the number of positive degree days.  There are two such schemes, derived 
classes of LocalMassBalance.  There is a deterministic
(expectation) default method and random (monte carlo) method.  The standard deviation of
the temperature variation change can be controlled by option <tt>-pdd_std_dev</tt>.  
The deterministic method computes only the expected number of positive degree days
for that amount of variability \ref CalovGreve05; it is chosen by option <tt>-pdd</tt>.
The random method uses pseudo-random numbers to simulate the variability and then directly 
sums the number of positive degree days.  It is chosen by either <tt>-pdd_rand</tt> or 
<tt>-pdd_rand_repeatable</tt>.

The surface mass balance is computed from the stored map of snow precipitation rate
by a call to the getMassFluxFromTemperatureTimeSeries() method of the LocalMassBalance object.
 */
PetscErrorCode PISMGreenlandAtmosCoupler::updateSurfMassFluxAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsmf) {
  PetscErrorCode ierr;

  // constants related to standard yearly cycle and precipitation offset
  const PetscScalar
    radpersec = 2.0 * pi / secpera, // radians per second frequency for annual cycle
    sperd = 8.64e4, // exact number of seconds per day
    julydaysec = sperd * config.get("snow_temp_july_day"),
    precipexpfactor = config.get("precip_exponential_factor_for_temperature");
 
  // set up snow temperature time series
  PetscInt     Nseries;
  ierr = mbscheme->getNForTemperatureSeries(t_years * secpera,
					    dt_years * secpera, Nseries); CHKERRQ(ierr);
  PetscScalar  *Tseries = new PetscScalar[Nseries];

  // times at which we compute snow temps using the parameterization are
  //    tseries, tseries + dtseries, ..., tseries + (Nseries-1) * dtseries;
  const PetscScalar tseries = (t_years - floor(t_years)) * secpera,
    dtseries = (dt_years * secpera) / ((PetscScalar) Nseries);

  // times for the snow temperature time-series (in years, used to get mean
  // annual temperature anomalies):
  vector<PetscScalar> ts(Nseries);
  for (PetscInt k = 0; k < Nseries; ++k)
    ts[k] = t_years + k * dt_years / Nseries;

  // use the parameterization to update T_ma and T_mj:
  ierr = parameterizedUpdateSnowSurfaceTemp(t_years,dt_years); CHKERRQ(ierr);

  if (snowprecipmaps != NULL) {
    ierr = snowprecipmaps->update(t_years, dt_years); CHKERRQ(ierr);
    ierr = snowprecipmaps->begin_access(); CHKERRQ(ierr);
  }

  if (snowtempmaps != NULL) {
    ierr = snowtempmaps->update(t_years, dt_years); CHKERRQ(ierr);
    ierr = snowtempmaps->begin_access(); CHKERRQ(ierr);
  }
  
  PetscScalar **T_ma, **T_mj, **snowrate, **lat_degN, **smflux;
  ierr = vsurfmassflux.get_array(smflux);  CHKERRQ(ierr);
  ierr = vsurftemp.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);
  ierr = vsnowprecip.get_array(snowrate);  CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      double precip = 0;

      // get the precipitation during the time interval (t_years, t_years +
      // dt_years); note that addition and averaging commute:
      if (snowprecipmaps != NULL) {
	ierr = snowprecipmaps->average(i, j, t_years, dt_years, precip); CHKERRQ(ierr);
      }
      
      precip += snowrate[i][j];

      // *if* dTforcing then modify precipitation, which is assumed to be
      //   present-day, to be modeled paleo-precipitation; use dTforcing at
      //   midpoint of interval [t_years,t_years+dt_years]; see formula for
      //   P_G(t,x,y) on SeaRISE "Model Initialization" page
      //   http://websrv.cs.umt.edu/isis/index.php/Model_Initialization the
      //   version here assumes Delta T_SC = 0.0
      if (dTforcing != NULL) {
	precip *= exp( precipexpfactor * (*dTforcing)(t_years + 0.5 * dt_years) );
      }

      // done with precipitation

      // build temperature time-series at the point i,j by filling Tseries with
      // anomalies and then adding the yearly cycle:
      if (snowtempmaps != NULL) {
	ierr = snowtempmaps->interp(i, j, Nseries, &ts[0], Tseries); CHKERRQ(ierr);
      } else {
	// apply the temperature offset if needed
	if (dTforcing != NULL) {

	  for (PetscInt k = 0; k < Nseries; ++k)
	    Tseries[k] = (*dTforcing)(ts[k]);

	} else {

	  for (PetscInt k = 0; k < Nseries; ++k)
	    Tseries[k] = 0.0;

	}
      }

      // add the standard yearly cycle using corrected formula (4) in [\ref Faustoetal2009]
      for (PetscInt k = 0; k < Nseries; ++k) {
	double tk = tseries + k * dtseries;
	Tseries[k] += T_ma[i][j] + (T_mj[i][j] - T_ma[i][j]) * cos(radpersec * (tk - julydaysec));
      }

      // done with temperature time-series

      // if we can, set mass balance parameters according to formula (6) in
      // [\ref Faustoetal2009]
      PDDMassBalance* pddscheme = dynamic_cast<PDDMassBalance*>(mbscheme);
      if (pddscheme != NULL) {
	ierr = pddscheme->setDegreeDayFactorsFromSpecialInfo(lat_degN[i][j],T_mj[i][j]); CHKERRQ(ierr);
      }

      // get surface mass balance at point i,j (units of day^-1 K-1)
      smflux[i][j] = mbscheme->getMassFluxFromTemperatureTimeSeries(tseries, dtseries, Tseries,
								    Nseries, precip);
    }
  }
  ierr = vsurfmassflux.end_access(); CHKERRQ(ierr);
  ierr = vsurftemp.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  ierr = vsnowprecip.end_access();  CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);

  if (snowtempmaps != NULL) {
    ierr = snowtempmaps->end_access(); CHKERRQ(ierr);

    // add the snow temperature anomaly to artm so that its value is the same
    // after updateSurfTempAndProvide() and updateSurfMassFluxAndProvide():
    ierr = snowtempmaps->interp(t_years); CHKERRQ(ierr);
    ierr = vsurftemp.add(1.0, *snowtempmaps); CHKERRQ(ierr);
  }

  if (snowprecipmaps != NULL) {
    ierr = snowprecipmaps->end_access(); CHKERRQ(ierr);
  }

  delete [] Tseries;

  // post-update in base class to add on -force_to_thk, etc.
  ierr = PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
     t_years, dt_years, pvsmf); CHKERRQ(ierr);

  return 0;
}


//! \brief Computes the maximum time-step (in seconds) this climate coupler can
//! take. Sets dt_years to a negative number if there is no restriction.
PetscErrorCode PISMGreenlandAtmosCoupler::max_timestep(PetscScalar t_years,
						       PetscScalar &dt_years) {
  PetscErrorCode ierr;
  
  dt_years = -1.0;
  
  // time-step restriction if reading climate from stored maps
  if ((snowtempmaps != NULL) && (snowprecipmaps != NULL)) {
    double snowtemp_dt, snowprecip_dt;
    snowtemp_dt   = snowtempmaps->max_timestep(t_years);
    snowprecip_dt = snowprecipmaps->max_timestep(t_years);

    if (snowtemp_dt < 0)   dt_years = snowprecip_dt;
    if (snowprecip_dt < 0) dt_years = snowtemp_dt;
    else                   dt_years = PetscMin(snowtemp_dt, snowprecip_dt);
  }
  
  // ask base class for time-step restriction (e.g. -force_to_thk)
  PetscScalar pac_dt_years;
  ierr = PISMAtmosphereCoupler::max_timestep(t_years, pac_dt_years); CHKERRQ(ierr);

  // resolve
  if (pac_dt_years >= 0.0) {
    if (dt_years < 0.0) {
      dt_years = pac_dt_years;
    } else {
      // case in which both PISMGreenlandAtmosCoupler and PISMAtmosphereCoupler
      //    have time-step restrictions; return minimum
      if (pac_dt_years < dt_years)  dt_years = pac_dt_years;
    }
  }
  return 0;
}

