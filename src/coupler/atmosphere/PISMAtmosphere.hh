// Copyright (C) 2008-2011 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

  //! \brief Sets result to the mean precipitation, in m/s ice equivalent.
  virtual PetscErrorCode mean_precip(IceModelVec2S &result) = 0;

  //! \brief Sets result to the mean annual near-surface air temperature, in degrees Kelvin.
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result) = 0;

  virtual PetscErrorCode begin_pointwise_access() = 0;
  virtual PetscErrorCode end_pointwise_access() = 0;
  
  //! \brief Sets a pre-allocated N-element array "values" to the time-series
  //! of near-surface air temperature (degrees Kelvin) at the point i,j on the
  //! grid. Times (in years) are specified in ts. NB! Has to be surrounded by
  //! begin_pointwise_access() and end_pointwise_access()
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values) = 0;
  //! \brief Sets result to a snapshot of temperature for the current time.
  //! (For diagnostic purposes.)
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result) = 0;
};


//! \brief A class implementing a constant-in-time atmosphere model. Reads data
//! from a PISM input file.
class PAConstant : public PISMAtmosphereModel {
public:
  PAConstant(IceGrid &g, const NCConfigVariable &conf)
    : PISMAtmosphereModel(g, conf) {};
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  { t = my_t; dt = my_dt; return 0; } // do nothing
  virtual PetscErrorCode mean_precip(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);
protected:
  string input_file;
  IceModelVec2S precip, temperature;
  NCSpatialVariable airtemp_var;
};

class PAConstantPIK : public PAConstant
{
public:
  PAConstantPIK(IceGrid &g, const NCConfigVariable &conf)
    : PAConstant(g, conf) {};
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
protected:
  IceModelVec2S *usurf, *lat;
};

//! A class containing an incomplete implementation of an atmosphere model
//! based on a temperature parameterization using mean annual and mean July
//! (mean summer) temperatures and a cosine yearly cycle. Uses a stored
//! (constant in time) precipitation field.
class PAYearlyCycle : public PISMAtmosphereModel {
public:
  PAYearlyCycle(IceGrid &g, const NCConfigVariable &conf)
    : PISMAtmosphereModel(g, conf) {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  //! This method implements the parameterization.
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt) = 0;
  virtual PetscErrorCode mean_precip(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);
protected:
  PISMVars *variables;
  PetscScalar snow_temp_july_day;
  string reference, precip_filename;
  IceModelVec2S temp_ma, temp_mj, precip;
  NCSpatialVariable airtemp_var;
};

class PAModifier : public Modifier<PISMAtmosphereModel>
{
public:
  PAModifier(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in)
    : Modifier<PISMAtmosphereModel>(g, conf, in) {}
  virtual ~PAModifier() {}

  virtual PetscErrorCode mean_precip(IceModelVec2S &result)
  {
    PetscErrorCode ierr = input_model->mean_precip(result); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result)
  {
    PetscErrorCode ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode begin_pointwise_access()
  {
    PetscErrorCode ierr = input_model->begin_pointwise_access(); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode end_pointwise_access()
  {
    PetscErrorCode ierr = input_model->end_pointwise_access(); CHKERRQ(ierr);
    return 0;
  }
  
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values)
  {
    PetscErrorCode ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result)
  {
    PetscErrorCode ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);
    return 0;
  }
};

//! \brief A class implementing an "atmosphere modifier" model applying forcing data
//! (anomalies, temperature offsets...) to results of another PISM atmosphere model.
/*!
Processes command-line options -anomaly_temp, -anomaly_precip.

The temperature anomaly should be interpreted as a change to the \e air temperature,
for example at 2m.  An underlying PISMSurfaceModel is in charge of producing an
actual ice fluid upper surface temperature.  The precipitation anomaly could
include rain, and again the underlying PISMSurfaceModel is in charge of removing
rain (if that's the model ...) using (for example) the 2m air temperature.
 */
class PAForcing : public PAModifier {
public:
  PAForcing(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel *input);
  virtual ~PAForcing();
  virtual PetscErrorCode max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict);
  virtual PetscErrorCode init(PISMVars &vars);
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
  virtual PetscErrorCode mean_precip(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result); 
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();  
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);
protected:
  IceModelVec2T *temp_anomaly, *precip_anomaly;
  NCSpatialVariable airtemp_var;
};


#endif	// __PISMAtmosphere_hh
