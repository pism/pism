// Copyright (C) 2008-2010 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

#ifndef __PISMAtmosphere_hh
#define __PISMAtmosphere_hh

#include "PISMComponent.hh"
#include "../base/iceModelVec.hh"
#include "../base/Timeseries.hh"
#include "iceModelVec2T.hh"

//! A purely virtual class defining the interface of a PISM Atmosphere Model.
class PISMAtmosphereModel : public PISMComponent {
public:
  PISMAtmosphereModel(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars) : PISMComponent(g, conf, vars) {};

  //! \brief Sets result to the mean precipitation over the time interval
  //! (t_years, t_years + dt_years), in m/s ice equivalent.
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2 &result) = 0;

  //! \brief Sets result to the mean annual near-surface air temperature corresponding
  //! to the time interval (t_years, t_years + dt_years), in degrees Kelvin.
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2 &result) = 0;

  virtual PetscErrorCode begin_pointwise_access() = 0;
  virtual PetscErrorCode end_pointwise_access() = 0;
  
  //! \brief Sets a pre-allocated N-element array "values" to the time-series
  //! of near-surface air temperature (degrees Kelvin) at the point i,j on the
  //! grid. Times (in years) are specified in ts. NB! Has to be surrounded by
  //! begin_pointwise_access() and end_pointwise_access()
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values) = 0;
  //! \brief Sets result to a snapshot of temperature for the interval
  //! (t_years, dt_years).
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2 &result) = 0;
};

//! \brief A class implementing a constant-in-time atmosphere model. Reads data
//! from a PISM input file.
class PAConstant : public PISMAtmosphereModel {
public:
  PAConstant(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars)
    : PISMAtmosphereModel(g, conf, vars) {};
  virtual PetscErrorCode init();
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2 &result);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2 &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode write_input_fields(PetscReal t_years, PetscReal dt_years,
					    string filename);
  virtual PetscErrorCode write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						 string filename);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2 &result);
protected:
  string input_file;
  IceModelVec2 snowprecip, temperature;
};

//! A class implementing an atmosphere model cobmining the Fausto [\ref
//! Faustoetal2009] present-day temperature parameterization for Greenland and
//! stored precipitation data.
class PAFausto : public PISMAtmosphereModel {
public:
  PAFausto(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars)
    : PISMAtmosphereModel(g, conf, vars) {};
  virtual PetscErrorCode init();
  virtual PetscErrorCode write_input_fields(PetscReal /*t_years*/,
					    PetscReal /*dt_years*/,
					    string filename);
  virtual PetscErrorCode write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						 string filename);
  virtual PetscErrorCode write_fields(set<string> vars, PetscReal t_years,
				      PetscReal dt_years, string filename);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2 &result);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2 &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2 &result);
protected:
  string reference, snowprecip_filename;
  IceModelVec2 temp_ma, temp_mj, snowprecip;
  IceModelVec2 *lat, *lon, *surfelev;
};

class PAModifier : public PISMAtmosphereModel {
public:
  PAModifier(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars)
    : PISMAtmosphereModel(g, conf, vars)
  { input_model = NULL; }

  virtual ~PAModifier()
  { delete input_model; }

  virtual void attach_input(PISMAtmosphereModel *input);
protected:
  PISMAtmosphereModel *input_model;
};

//! \brief A class implementing an "atmosphere model" applying forcing data
//! (anomalies, temperature offsets...) to results of another PISM atmosphere
//! model.
/*! Processes command-line options -dTforcing, -temp_ma_anomaly,
  -snowprecip_anomaly, -paleo_precip.
 */
class PAForcing : public PAModifier {
public:
  PAForcing(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars);
  virtual ~PAForcing();
  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years);
  virtual PetscErrorCode init();
  virtual PetscErrorCode write_input_fields(PetscReal t_years, PetscReal dt_years,
					    string filename);
  virtual PetscErrorCode write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						 string filename);
  virtual PetscErrorCode write_fields(set<string> vars, PetscReal t_years,
				      PetscReal dt_years, string filename); // needs work
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2 &result);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2 &result); 
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();  
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2 &result);
protected:
  Timeseries *dTforcing;
  DiagnosticTimeseries *delta_T; //!< for debugging
  bool paleo_precipitation_correction;
  IceModelVec2T *temp_ma_anomaly, *snowprecip_anomaly;

  PISMAtmosphereModel *input_model;
};

/*
//! A class implementing a simple atmospheric lapse rate model.
class PALapseRates : public PAModifier {
public:
  PALapseRates(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars);
  virtual ~PALapseRates();
  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years); 

  virtual PetscErrorCode init(); 
  virtual PetscErrorCode write_fields(set<string> vars, PetscReal t_years,
				      PetscReal dt_years, string filename); 
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years); 
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2 &result); 
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2 &result); 
  virtual PetscErrorCode begin_pointwise_access(); 
  virtual PetscErrorCode end_pointwise_access();   
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values); 

protected:
  PetscReal temperature_lapse_rate;
  IceModelVec2 surfelev_initial;
};
*/

#endif	// __PISMAtmosphere_hh
