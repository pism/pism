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
#include "iceModelVec.hh"
#include "Timeseries.hh"
#include "iceModelVec2T.hh"


///// PISMAtmosphereModel: models which provide precipitation and temperature
/////                      to the PISMSurfaceModel below

//! A purely virtual class defining the interface of a PISM Atmosphere Model.
class PISMAtmosphereModel : public PISMComponent_TS {
public:
  PISMAtmosphereModel(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf) {};

  //! \brief Sets result to the mean precipitation over the time interval
  //! (t_years, t_years + dt_years), in m/s ice equivalent.
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2S &result) = 0;

  //! \brief Sets result to the mean annual near-surface air temperature corresponding
  //! to the time interval (t_years, t_years + dt_years), in degrees Kelvin.
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2S &result) = 0;

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
				       IceModelVec2S &result) = 0;
};


//! \brief A class implementing a constant-in-time atmosphere model. Reads data
//! from a PISM input file.
class PAConstant : public PISMAtmosphereModel {
public:
  PAConstant(IceGrid &g, const NCConfigVariable &conf)
    : PISMAtmosphereModel(g, conf) {};
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2S &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					    string filename);
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2S &result);
protected:
  string input_file;
  IceModelVec2S precip, temperature;
  NCSpatialVariable airtemp_var;
};


//! A class containing an incomplete implementation of an atmosphere model
//! based on a temperature parameterization using mean annual and mean July
//! (mean summer) temperatures and a cosine yearly cycle. Uses a stored
//! (constant in time) precipitation field.
class PAYearlyCycle : public PISMAtmosphereModel {
public:
  PAYearlyCycle(IceGrid &g, const NCConfigVariable &conf)
    : PISMAtmosphereModel(g, conf) {}
  virtual PetscErrorCode init(PISMVars &vars);	      // nb
  virtual PetscErrorCode write_model_state(PetscReal /*t_years*/,
					    PetscReal /*dt_years*/,
					    string filename);
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  //! This method implements the parameterization.
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years) = 0;
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2S &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2S &result);
protected:
  PetscScalar snow_temp_july_day;
  string reference, precip_filename;
  IceModelVec2S temp_ma, temp_mj, precip;
};


//! \brief A modification of PAYearlyCycle tailored for the
//! SeaRISE-Greenland assessment. Uses the Fausto [\ref Faustoetal2009]
//! present-day temperature parameterization and stored precipitation data.
//! Adds the precipitation correction for spin-ups.
class PA_SeaRISE_Greenland : public PAYearlyCycle {
public:
  PA_SeaRISE_Greenland(IceGrid &g, const NCConfigVariable &conf)
    : PAYearlyCycle(g, conf)
  {
    paleo_precipitation_correction = false;
    dTforcing = NULL;
  }

  virtual ~PA_SeaRISE_Greenland()
  {
    delete dTforcing;
  }
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2S &result);
protected:
  bool paleo_precipitation_correction;
  Timeseries *dTforcing;
  IceModelVec2S *lat, *lon, *surfelev;
};


//! A class implementing a simple atmospheric lapse rate model.
/*!
  Let \f$T\f$ be the temperature and \f$h\f$ the surface elevation. Then the lapse rate \f$\gamma\f$ is
\f[
\gamma = -\frac{dT}{dh}.
\f]
This equation can be solved exactly, to obtain
\f[
T(h) = -\gamma \cdot h + f,
\f]
or
\f[
T(x,y) = -\gamma \cdot h(x,y) + f(x,y),
\f]
where \f$f(x,y)\f$ is the initial condition. We have
\f[
f(x,y) = T_0(x,y) + \gamma\cdot h_0(x,y).
\f]

Class PALapseRates implements this lapse-rate correction mechanism.
 */
class PALapseRates : public PAConstant {
public:
  PALapseRates(IceGrid &g, const NCConfigVariable &conf)
    : PAConstant(g, conf)
  {
    gamma = 0;
    usurf = NULL;
  }
  virtual ~PALapseRates() {}
  virtual PetscErrorCode init(PISMVars &vars); 
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2S &result); 
  virtual PetscErrorCode begin_pointwise_access(); 
  virtual PetscErrorCode end_pointwise_access();   
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values); 
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					    string filename);
protected:
  PetscReal gamma;
  IceModelVec2S f, *usurf;
};


///// PAModifier: classes which modify outputs of PISMAtmosphereModel and its derived classes

class PAModifier : public PISMAtmosphereModel {
public:
  PAModifier(IceGrid &g, const NCConfigVariable &conf)
    : PISMAtmosphereModel(g, conf)
  { input_model = NULL; }

  virtual ~PAModifier()
  { delete input_model; }

  virtual void attach_input(PISMAtmosphereModel *input) {
    if (input_model != NULL)  delete input_model;
    input_model = input;
  }

  virtual void add_vars_to_output(string key, set<string> &result) {
    if (input_model != NULL)  input_model->add_vars_to_output(key, result);
  }

protected:
  PISMAtmosphereModel *input_model;
};

//! \brief A class implementing an "atmosphere modifier" model applying forcing data
//! (anomalies, temperature offsets...) to results of another PISM atmosphere model.
/*!
Processes command-line options -dTforcing, -anomaly_temp, -anomaly_precip.

The temperature anomaly should be interpreted as a change to the \e air temperature,
for example at 2m.  An underlying PISMSurfaceModel is in charge of producing an
actual ice fluid upper surface temperature.  The precipitation anomaly could
include rain, and again the underlying PISMSurfaceModel is in charge of removing
rain (if that's the model ...) using (for example) the 2m air temperature.
 */
class PAForcing : public PAModifier {
public:
  PAForcing(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PAForcing();
  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years);
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					    string filename);
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2S &result); 
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();  
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2S &result);
protected:
  Timeseries *dTforcing;
  DiagnosticTimeseries *delta_T; //!< for diagnostic time series output (-ts_vars?)
  IceModelVec2T *temp_anomaly, *precip_anomaly;
  NCSpatialVariable airtemp_var;
};


#endif	// __PISMAtmosphere_hh
